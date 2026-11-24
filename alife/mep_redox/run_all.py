"""Full sweep for mep_redox.  Usage:  python3 -I run_all.py [--workers 3] [--quick]
Writes results/*.tsv, synthetic/<run_id>/ (observer files) and truth/<run_id>/ (keys)."""
import argparse
import json
import os
import sys
import time
import numpy as np
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from model import Model, FixedEnzyme, default_params, TYPES, RXN, T  # noqa: E402
import analysis as A  # noqa: E402
import mep  # noqa: E402
import observe as O  # noqa: E402
import infer as I  # noqa: E402

RES = os.path.join(HERE, "results")
SYN = os.path.join(HERE, "synthetic")
TRU = os.path.join(HERE, "truth")
T_END = 1500.0
SEEDS = [0, 1, 2]
KO_HORIZON = 600.0

MECH = {  # switchable kinetic symmetry breakers
    "neutral": dict(),
    "crowd": dict(alpha=1.3),
    "stackadv": dict(alpha=0.8),
    "constitutive": dict(beta=0.0),
    "channel": dict(phi=0.8),
    "toxic_channel": dict(K_I_NO2=2.0, phi=0.8),
}
OTHER = {  # non-MEP / forcing controls (full pool)
    "electron_yield": dict(yield_mode="electron"),
    "dispersal": dict(cell_mix=30.0),
    "forced": dict(force_amp=0.8, force_period=30.0),
}
SPLIT_SPECIALISTS = ["NAR", "NIR", "AOB", "NOB", "FER", "MET", "SRBH", "SNR"]
STACKED = ["DEN", "AERDEN", "CMX", "FERMET", "SRBX", "SOBX"]
POOLS = {
    "full": list(TYPES),
    "split_only": [t for t in TYPES if t not in STACKED],
    "stack_only": [t for t in TYPES if t not in ("NAR", "NIR", "AOB", "NOB", "FER", "MET",
                                                 "AERNAR")],
}


def params(scen, seed, pool):
    kw = dict(MECH.get(scen, {}))
    kw.update(OTHER.get(scen, {}))
    types = POOLS[pool]
    rng = np.random.default_rng(1000 + seed)
    kap_all = dict(zip(TYPES, rng.lognormal(0, 0.15, len(TYPES))))
    kap_all["ALG"] = 1.0
    kw.update(types=types, kappa=[kap_all[t] for t in types])
    return default_params(**kw)


def run_id(scen, seed, pool):
    return f"{scen}__{pool}__s{seed}"


def job_population(args):
    scen, seed, pool = args
    rid = run_id(scen, seed, pool)
    t0 = time.time()
    p = params(scen, seed, pool)
    m = Model(p)
    y0 = m.y0(np.random.default_rng(seed))
    te = np.unique(np.concatenate([np.linspace(0, 100, 51), np.arange(100, T_END + 1, 5),
                                   [T_END - 20, T_END - 10]]))
    sol = m.run(y0, T_END, t_eval=te)
    sig = np.array([m.sigma_parts(sol.y[:, k], sol.t[k])[:2] for k in range(len(sol.t))])
    np.savez_compressed(os.path.join(RES, "runs", rid + ".npz"), t=sol.t, y=sol.y[:, -1],
                        y_casts=sol.y[:, [list(sol.t).index(T_END - 20),
                                          list(sol.t).index(T_END - 10), len(sol.t) - 1]],
                        sigma=sig)
    return dict(run_id=rid, scenario=scen, pool=pool, seed=seed, runtime_s=round(time.time() - t0, 1))


def load_run(rid):
    z = np.load(os.path.join(RES, "runs", rid + ".npz"))
    return z["t"], z["y"], z["y_casts"], z["sigma"]


def job_analyse(args):
    """Ground truth + MEP yardsticks + synthetic observations for a full-pool run."""
    scen, seed, pool, do_mep, blind = args
    rid = run_id(scen, seed, pool)
    p = params(scen, seed, pool)
    m = Model(p)
    t, y, yc, sig = load_run(rid)
    out = dict(run_id=rid, scenario=scen, seed=seed, pool=pool)
    # time-averaged sigma over last 300 d (covers 10 forcing periods)
    sel = t >= T_END - 300
    out["sigma_chem"] = float(np.trapezoid(sig[sel, 0], t[sel]) / (t[sel][-1] - t[sel][0]))
    out["sigma_alg"] = float(np.trapezoid(sig[sel, 1], t[sel]) / (t[sel][-1] - t[sel][0]))
    out["alg_share_sigma"] = out["sigma_alg"] / (out["sigma_alg"] + out["sigma_chem"])
    out["npp"] = float(np.mean([m.npp(yc[:, k], T_END - 20 + 10 * k) for k in range(3)]))
    out["sigma_chem_per_npp"] = out["sigma_chem"] / max(out["npp"], 1e-12)
    C, B = m.unpack(y)
    bio = [dict(run_id=rid, type=tn, column_biomass=float(B[j].sum() * p["dz"]),
                depth_centroid=float((B[j] * m.z).sum() / max(B[j].sum(), 1e-12)))
           for j, tn in enumerate(m.types)]
    sp = m.sigma_parts(y)[3]
    for j, b in enumerate(bio):
        b["sigma_chem_share"] = float(sp[j] / max(sp.sum(), 1e-12))
    tr = A.attribution(m, y)
    for x in tr:
        x["run_id"] = rid
    prow, edges = A.partner_counts(m, y, tr)
    for x in prow:
        x["run_id"] = rid
    # algal share of each recipient's intake (all resources, excluding abiotic sources)
    for x in prow:
        if x["resource"] == "*":
            rn = x["recipient"]
            mine = [q for q in tr if q["recipient"] == rn and q["kind"] == "pool"
                    and q["donor"] not in ("ATM", "LAT")]
            tot = sum(q["flux"] for q in mine)
            x["alg_share_intake"] = (sum(q["flux"] for q in mine if q["donor"] == "ALG") / tot
                                     if tot > 0 else np.nan)
    fe_rows = []
    fe = A.free_energy_budget(m, y, T_END)
    fe.update(run_id=rid, config="emergent:" + scen, npp=m.npp(y, T_END))
    fe_rows.append(fe)
    dol = A.division_of_labour(m, y)
    for x in dol:
        x["run_id"] = rid
    mrow = None
    if do_mep:
        C0, E, Balg, necro = mep.fixed_inputs(m, y)
        g_mean, g_last, _ = mep.sigma_fixed(p, E, Balg, necro, C0, 1000)
        gr_mean, gr_last, gr = mep.sigma_fixed(p, E, Balg, necro, C0, 1000, greedy=20.0)
        Eopt, s_opt, s_start, info = mep.optimise(p, E, Balg, necro, C0, 1000, maxfev=250,
                                                  steady=True)
        rs = []
        rfields = mep.random_fields(E, np.random.default_rng(seed + 77), n=6)
        for Er in rfields:
            a, b_, s_ = mep.sigma_fixed(p, Er, Balg, necro, C0, 1000)
            rs.append(b_)
        # common feasible set: the reference cannot represent channelling and holds B12 and
        # algae fixed; emergent sigma is comparable only where the emergent state lies in
        # the reference family (gauge copy == emergent) -> flag
        gauge_rel = abs(g_last - out["sigma_chem"]) / out["sigma_chem"]
        mrow = dict(run_id=rid, sigma_emergent=out["sigma_chem"], sigma_gauge_fixed=g_last,
                    gauge_rel_diff=gauge_rel, fixed_field_sigma_within_0p5pct=bool(gauge_rel < 0.005),
                    sigma_greedy=gr_last, sigma_best_found=s_opt, sigma_ref_start=s_start,
                    sigma_random_mean=float(np.mean(rs)),
                    sigma_random_sd=float(np.std(rs)), sigma_random_max=float(np.max(rs)),
                    emergent_over_best_found=out["sigma_chem"] / s_opt,
                    gauge_over_best_found=g_last / s_opt, **info)
        np.savez_compressed(os.path.join(RES, "runs", rid + "__mep.npz"), Eopt=Eopt, E=E,
                            Balg=Balg, necro=necro, C0=C0)
    # synthetic observations
    cat = O.build_catalogue()
    ttimes = [T_END - 20, T_END - 10, T_END]
    casts = [O.cast_from_population(m, yc[:, k], ttimes[k]) for k in range(3)]
    O.emit(os.path.join(SYN, blind[rid]), os.path.join(TRU, rid), blind[rid], casts, cat,
           seed=seed + 5)
    truth_typed = {(a, b, r) for a, b, r in edges}
    _write_truth(rid, tr, m, y)
    key = {g["mag_id"]: tn for tn, g in cat.items()}
    sc = I.score(os.path.join(SYN, blind[rid]), truth_typed, key)
    sc.update(run_id=rid, dataset_id=blind[rid], config="emergent:" + scen,
              sigma_chem=out["sigma_chem"])
    sigs = [sc]
    if do_mep:
        # observations from fixed yardstick configurations (same medium)
        for lab, Ef, greedy in (("best_found_ref", Eopt, None), ("greedy", E, 20.0),
                                ("random1", rfields[0], None), ("random2", rfields[1], None)):
            f = FixedEnzyme(p, Ef, Balg, necro, greedy_beta=greedy)
            mean, last, res_ = mep.sigma_fixed(p, Ef, Balg, necro, C0, 1000, greedy=greedy)
            if res_ is None:
                continue
            yl = res_[2]
            fe2 = A.free_energy_budget(f, yl, T_END)
            fe2.update(run_id=rid + "__" + lab, config=lab, npp=f.npp(yl, T_END))
            fe_rows.append(fe2)
            cs = [O.cast_from_fixed(f, yl, T_END)] * 1
            rid2 = rid + "__" + lab
            O.emit(os.path.join(SYN, blind[rid2]), os.path.join(TRU, rid2), blind[rid2], cs, cat,
                   seed=seed + 9)
            tr2 = A.attribution(f, yl)
            # map reaction pseudo-types to carrier types for truth edges
            mp = dict(O.REACTION_CARRIER)
            mp["OMNI"] = "OMNI"
            mp["ALG"] = "ALG"
            e2 = set()
            intake = {}
            for q in tr2:
                if q["kind"] == "pool":
                    intake[(q["recipient"], q["resource"])] = intake.get((q["recipient"], q["resource"]), 0) + q["flux"]
            for q in tr2:
                if q["kind"] == "pool" and q["donor"] not in ("ATM", "LAT") and q["donor"] != q["recipient"]:
                    if q["flux"] >= 0.1 * intake[(q["recipient"], q["resource"])]:
                        e2.add((mp.get(q["donor"], q["donor"]), mp.get(q["recipient"], q["recipient"]), q["resource"]))
            sc2 = I.score(os.path.join(SYN, blind[rid2]), e2, key)
            sc2.update(run_id=rid2, dataset_id=blind[rid2], config=lab, sigma_chem=last)
            sigs.append(sc2)
    return dict(out=out, bio=bio, tr=tr, prow=prow, dol=dol, mrow=mrow, sigs=sigs, fe=fe_rows)


def _write_truth(rid, tr, m, y):
    d = os.path.join(TRU, rid)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "realised_transfers.tsv"), "w") as fh:
        fh.write("resource\tdonor\trecipient\tflux_mmol_m2_d\tkind\n")
        for x in tr:
            fh.write(f"{x['resource']}\t{x['donor']}\t{x['recipient']}\t{x['flux']:.6g}\t{x['kind']}\n")
    C, B = m.unpack(y)
    P = m.tendencies(0.0, y, parts=True)
    Jt = P["J"] + P["Jint"]
    with open(os.path.join(d, "flux_by_depth.tsv"), "w") as fh:
        fh.write("type\treaction\tdepth_m\tflux_mmol_m3_d\tdG_kJ_mol\n")
        for j, tn in enumerate(m.types):
            for r in RXN:
                ri = RXN.index(r)
                if m.M[j, ri]:
                    for k, z in enumerate(m.z):
                        if Jt[j, ri, k] > 1e-6:
                            fh.write(f"{tn}\t{r}\t{z}\t{Jt[j, ri, k]:.5g}\t{P['dG'][ri, k]:.4g}\n")


def job_knockout(args):
    scen, seed = args
    rid = run_id(scen, seed, "full")
    p = params(scen, seed, "full")
    m = Model(p)
    t, y, yc, sig = load_run(rid)
    tr = A.attribution(m, y)
    rows, runs = A.knockouts(m, y, tr, horizon=KO_HORIZON)
    for r in rows:
        r["run_id"] = rid
    return rows


def job_horizon(args):
    """Static enzyme fields optimised for mean sigma over [0, H] from a sterile start, then
    evaluated over every horizon; plus greedy and the emergent assembly from the same start."""
    scen, seed = args
    rid = run_id(scen, seed, "full")
    p = params(scen, seed, "full")
    m = Model(p)
    t, y, yc, sig = load_run(rid)
    C0, E, Balg, necro = mep.fixed_inputs(m, y)
    # sterile start: algae and necromass as emergent, no chemotroph enzyme, 300 d
    f0 = FixedEnzyme(p, E * 0.0, Balg, necro)
    s0 = f0.run(C0.ravel(), 300.0)
    Cs = s0.y[:, -1]
    H_list = [1.0, 10.0, 100.0, 1000.0]
    rows = []
    hinfo = []
    pols = {}
    for H in H_list:
        Eh, sh, s_start, info = mep.optimise(p, E, Balg, necro, Cs, H, maxfev=100,
                                             recheck_fev=40)
        pols[f"static_best_H{int(H)}"] = (Eh, None)
        hinfo.append(dict(run_id=rid, policy=f"static_best_H{int(H)}", opt_horizon_d=H,
                          sigma_best_found=sh, **info))
    pols["emergent_field_static"] = (E, None)
    pols["greedy_instantaneous"] = (E, 20.0)
    for lab, (Eh, gr) in pols.items():
        f = FixedEnzyme(p, Eh, Balg, necro, greedy_beta=gr)
        te = np.unique(np.concatenate([np.linspace(0, 1, 21), np.linspace(1, 10, 37),
                                       np.linspace(10, 100, 91), np.linspace(100, 1000, 91)]))
        sol = f.run(Cs, 1000.0, t_eval=te, rtol=1e-4, atol=1e-6)
        ss = np.array([f.sigma_parts(sol.y[:, k], sol.t[k])[0] for k in range(len(sol.t))])
        for He in H_list:
            sel = sol.t <= He
            rows.append(dict(run_id=rid, policy=lab, eval_horizon_d=He,
                             mean_sigma=float(np.trapezoid(ss[sel], sol.t[sel]) / He)))
    # emergent assembly from the sterile start (seed biomass 0.01 uM C per type and layer)
    mm = Model(p)
    yA = mm.y0(np.random.default_rng(seed))
    Cm, Bm = mm.unpack(yA)
    Cm[:] = Cs.reshape(Cm.shape)
    Bm[mm.ialg] = Balg
    te = np.unique(np.concatenate([np.linspace(0, 1, 21), np.linspace(1, 10, 37),
                                   np.linspace(10, 100, 91), np.linspace(100, 1000, 91)]))
    sol = mm.run(yA, 1000.0, t_eval=te)
    ss = np.array([mm.sigma_parts(sol.y[:, k], sol.t[k])[0] for k in range(len(sol.t))])
    for He in H_list:
        sel = sol.t <= He
        rows.append(dict(run_id=rid, policy="population_assembly", eval_horizon_d=He,
                         mean_sigma=float(np.trapezoid(ss[sel], sol.t[sel]) / He)))
    return rows, hinfo


def make_blind_key(pairs):
    """Anonymous dataset ids for synthetic/ (scenario names stripped).  The permutation uses
    OS entropy, so ids cannot be recomputed from scenario names; the key is written to
    truth/unblinding_key.tsv only."""
    import secrets
    ids = []
    for s, sd in pairs:
        rid = run_id(s, sd, "full")
        ids.append(rid)
        if s != "forced":
            ids += [rid + "__" + lab for lab in ("best_found_ref", "greedy", "random1", "random2")]
    rng = np.random.default_rng(secrets.randbits(64))
    nums = rng.permutation(len(ids)) + 1
    key = {rid: f"D{n:03d}" for rid, n in zip(ids, nums)}
    os.makedirs(TRU, exist_ok=True)
    with open(os.path.join(TRU, "unblinding_key.tsv"), "w") as fh:
        fh.write("dataset_id\trun_id\n")
        for rid, d in sorted(key.items(), key=lambda x: x[1]):
            fh.write(f"{d}\t{rid}\n")
    return key


def write_tsv(name, rows, cols=None):
    if not rows:
        return
    cols = cols or list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(RES, name), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            vals = []
            for c in cols:
                v = r.get(c, "")
                if isinstance(v, float):
                    v = f"{v:.6g}"
                vals.append(str(v))
            fh.write("\t".join(vals) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", type=int, default=3)
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    os.makedirs(os.path.join(RES, "runs"), exist_ok=True)
    seeds = [0] if a.quick else SEEDS
    t0 = time.time()
    pop_jobs = [(s, sd, pool) for s in MECH for sd in seeds for pool in POOLS]
    pop_jobs += [(s, sd, "full") for s in OTHER for sd in seeds]
    with Pool(a.workers) as pool:
        runs = pool.map(job_population, pop_jobs, chunksize=1)
        print(f"population runs done {time.time() - t0:.0f}s", flush=True)
        blind = make_blind_key([(s, sd) for s in list(MECH) + list(OTHER) for sd in seeds])
        an_jobs = [(s, sd, "full", s != "forced", blind) for s in list(MECH) + list(OTHER)
                   for sd in seeds]
        ko_jobs = [("neutral", sd) for sd in seeds] + \
            [(s, 0) for s in ("crowd", "constitutive", "channel", "electron_yield")]
        hz_jobs = [("neutral", 0)]
        r_an = pool.map_async(job_analyse, an_jobs, chunksize=1)
        r_ko = pool.map_async(job_knockout, ko_jobs, chunksize=1)
        r_hz = pool.map_async(job_horizon, hz_jobs, chunksize=1)
        an = r_an.get()
        print(f"analysis done {time.time() - t0:.0f}s", flush=True)
        ko = r_ko.get()
        print(f"knockouts done {time.time() - t0:.0f}s", flush=True)
        hz = r_hz.get()
        print(f"horizon done {time.time() - t0:.0f}s", flush=True)
    write_tsv("population_runs.tsv", runs)
    write_tsv("sigma.tsv", [x["out"] for x in an])
    write_tsv("biomass.tsv", [b for x in an for b in x["bio"]])
    write_tsv("realised_transfers.tsv", [b for x in an for b in x["tr"]],
              ["run_id", "resource", "donor", "recipient", "flux", "kind"])
    write_tsv("partners_realised.tsv", [b for x in an for b in x["prow"]],
              ["run_id", "recipient", "resource", "intake", "shares", "d_real", "donors",
               "mean_hill", "c_real", "scc", "alg_share_intake"])
    write_tsv("division_of_labour.tsv", [b for x in an for b in x["dol"]])
    write_tsv("mep_yardsticks.tsv", [x["mrow"] for x in an if x["mrow"]])
    write_tsv("signatures_inference.tsv", [b for x in an for b in x["sigs"]])
    write_tsv("free_energy_budget.tsv", [b for x in an for b in x["fe"]])
    write_tsv("necessity_knockouts.tsv", [r for rows in ko for r in rows])
    write_tsv("horizon.tsv", [r for rows, hi in hz for r in rows])
    write_tsv("horizon_convergence.tsv", [r for rows, hi in hz for r in hi])
    # pool contrast: sigma of full vs split-only vs stack-only (steady, last 300 d)
    pc = []
    for s in MECH:
        for sd in seeds:
            row = dict(scenario=s, seed=sd)
            for pl in POOLS:
                t, y, yc, sig = load_run(run_id(s, sd, pl))
                sel = t >= T_END - 300
                row["sigma_" + pl] = float(np.trapezoid(sig[sel, 0], t[sel]) / (t[sel][-1] - t[sel][0]))
                mm_ = Model(params(s, sd, pl))
                row["npp_" + pl] = mm_.npp(y, T_END)
                row["perC_" + pl] = row["sigma_" + pl] / max(row["npp_" + pl], 1e-12)
            pc.append(row)
    write_tsv("pool_contrast.tsv", pc)
    with open(os.path.join(RES, "sweep_meta.json"), "w") as fh:
        json.dump(dict(seeds=seeds, T_END=T_END, workers=a.workers,
                       wall_s=round(time.time() - t0, 1), mech=MECH,
                       other={k: v for k, v in OTHER.items()}, pools=POOLS), fh, indent=1)
    print(f"total {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
