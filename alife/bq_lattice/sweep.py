"""Parameter sweep for bq_lattice: run + infer each config x seed, then aggregate tables.

usage: python3 -I sweep.py OUTROOT [NWORKERS]     (aggregate only: python3 -I sweep.py OUTROOT agg)
"""
import sys, os, json, time, itertools
import multiprocessing as mp
from collections import defaultdict
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

BASE = dict(leak=0.5, cost=0.03, D=0.3, iota=1.0, phi=0.2, euk_frac=0.0, o_resp=0.2)


def configs():
    """Revised reduced sweep (after cross-family review): block A prokaryotes only at iota 0/1 for
    base, leak 0.2/0.8, cost 0.08, D 0.05, phi 0/0.4; block E with the eukaryote at base."""
    out = []
    for io in (0.0, 1.0):
        out.append(("A", dict(BASE, iota=io), 2))
        for k, v in [("leak", 0.2), ("leak", 0.8), ("cost", 0.08), ("D", 0.05), ("phi", 0.0), ("phi", 0.4)]:
            out.append(("A", dict(BASE, iota=io, **{k: v}), 2))
        out.append(("E", dict(BASE, iota=io, euk_frac=0.1), 2))
    return out


def cname(block, c):
    return block + "_" + "_".join(f"{k}{c[k]:g}" for k in ("iota", "leak", "cost", "D", "phi", "o_resp", "euk_frac"))


def job(args):
    root, name, seed, params = args
    rd = f"{root}/runs/{name}/s{seed}"
    if os.path.exists(f"{rd}/inference.json"):
        return rd, "cached"
    import run as R, infer as I
    t0 = time.time()
    try:
        R.main(rd, seed, params)
        I.main(rd)
        import gzip, shutil
        for p in ("obs/per_mag_cds.tsv", "truth/per_mag_cds_complete.tsv", "obs/expression.tsv"):
            with open(f"{rd}/{p}", "rb") as fi, gzip.open(f"{rd}/{p}.gz", "wb") as fo:
                shutil.copyfileobj(fi, fo)
            os.remove(f"{rd}/{p}")
        return rd, f"ok {time.time() - t0:.0f}s"
    except Exception as e:  # keep the sweep going; record failure
        return rd, f"FAIL {e!r}"


def read_tsv(p):
    with open(p) as f:
        h = f.readline().rstrip("\n").split("\t")
        return [dict(zip(h, l.rstrip("\n").split("\t"))) for l in f if l.strip()]


def summarise(rd):
    meta = json.load(open(f"{rd}/meta.json"))
    inf = json.load(open(f"{rd}/inference.json"))
    s = dict(final_cells=meta["final_cells"], runtime_s=meta["runtime_s"])
    ts = read_tsv(f"{rd}/truth/timeseries.tsv")[-1]
    s["n_genotypes"] = int(ts["n_genotypes"]); s["euk_cell_frac"] = int(ts["n_euk_cells"]) / max(int(ts["n_cells"]), 1)
    for n in ("L0", "OX", "ANX", "P0"):
        s[f"carry_{n}"] = float(ts[f"carry_{n}"])
    # genome-level: cells with no leaky need, max leaky-compatible functions in a common genotype,
    # co-carriage of OX and ANX (last snapshot)
    z = np.load(f"{rd}/truth/lattice_snapshots.npz")
    g = z["geno"][-1]; g = g[g >= 0].astype(np.int64)
    prok = g[(g >> 11) & 1 == 0]
    s["frac_genome_prototroph"] = float(np.mean((prok & 127) == 127)) if len(prok) else np.nan
    s["frac_cells_OX_and_ANX"] = float(np.mean(((prok >> 5) & 1) & ((prok >> 6) & 1))) if len(prok) else np.nan
    cnt = np.bincount(prok, minlength=2048)
    common = np.nonzero(cnt >= 20)[0]
    s["max_leaky_in_common_genotype"] = int(max([bin(x & 31).count("1") for x in common], default=0))
    s["mean_leaky_carried"] = float(np.mean([bin(x & 31).count("1") for x in prok])) if len(prok) else np.nan
    # necessity regime (P-all over tested kept sets), cell-weighted over necessity time points
    nec = [r for r in read_tsv(f"{rd}/truth/necessity_lineage.tsv") if r["scope"] == "all"]
    w = defaultdict(float); tot = 0
    for r in nec:
        n = int(r["n_cells"]); tot += n
        st = r["status"]
        if st == "plexi":
            w[r.get("plexi_type") or "fanin"] += n
        else:
            w[st] += n
    for k in ("autonomous", "not_plexi", "fanin", "cyclic", "mixed", "unresolved_no_feasible_set",
              "not_viable_in_situ"):
        s[f"nec_{k}"] = w[k] / tot if tot else np.nan
    s["nec_plexi_any"] = (w["fanin"] + w["cyclic"] + w["mixed"]) / tot if tot else np.nan
    cs = inf.get("capacity_vs_survival", {})
    ga = sum(v for k, v in cs.items() if k.startswith("genome_autonomous=True"))
    ga_s = cs.get("genome_autonomous=True|survives_alone=True", 0)
    gn_s = cs.get("genome_autonomous=False|survives_alone=True", 0)
    gn = sum(v for k, v in cs.items() if k.startswith("genome_autonomous=False"))
    s["genome_auto_rows"] = ga; s["genome_auto_survive_frac"] = ga_s / ga if ga else np.nan
    s["genome_dep_rows"] = gn; s["genome_dep_survive_alone_frac"] = gn_s / gn if gn else np.nan
    # realised fan-in (cell-weighted, scope all, last 3 sample times)
    real = [r for r in read_tsv(f"{rd}/truth/realised_lineage.tsv") if r["scope"] == "all"]
    tt = sorted({int(r["time"]) for r in real})[-3:]
    rr = [r for r in real if int(r["time"]) in tt]
    tn = sum(int(r["n_cells"]) for r in rr)
    s["real_d2"] = sum(int(r["n_cells"]) for r in rr if int(r["d_real"]) >= 2) / tn if tn else np.nan
    s["real_c3"] = sum(int(r["n_cells"]) for r in rr if int(r.get("c_real") or 1) >= 3) / tn if tn else np.nan
    # abundance tail
    ab = read_tsv(f"{rd}/truth/abundance_fit.tsv")
    vz = [float(r["vuong_z"]) for r in ab if r.get("vuong_z") not in (None, "", "nan")]
    s["abund_vuong_z_mean"] = float(np.mean(vz)) if vz else np.nan
    s["abund_frac_z_gt2"] = float(np.mean([x > 2 for x in vz])) if vz else np.nan
    s["abund_frac_z_lt-2"] = float(np.mean([x < -2 for x in vz])) if vz else np.nan
    s["abund_alpha"] = float(np.mean([float(r["alpha_pl"]) for r in ab if r.get("alpha_pl")])) if ab else np.nan
    s["abund_top10_share"] = float(np.mean([float(r["top10_share"]) for r in ab if r.get("top10_share")])) if ab else np.nan
    # inference
    for vn in ("complete_prok", "incomplete_prok", "incomplete_all", "complete_all", "incomplete_prok_core"):
        sc = inf["score_necessity"][vn]
        for m_ in ("res", "low", "high"):
            s[f"{vn}_prec_{m_}"] = sc.get(f"precision_{m_}"); s[f"{vn}_rec_{m_}"] = sc.get(f"recall_{m_}")
        ev = sum(sc.get(f"res_{x}", 0) for x in ("tp", "fp", "fn", "tn"))
        s[f"{vn}_frac_unres"] = sc.get("unresolved", 0) / max(ev + sc.get("unresolved", 0), 1)
        bt = sc.get("recall_by_truth_type", {})
        for k in ("fanin", "cyclic", "mixed"):
            v = bt.get(k)
            s[f"{vn}_rec_{k}"] = v["called"] / v["n"] if v and v["n"] else np.nan
    s["inc_d_bias"] = inf["incompleteness"]["mean_d_bias"]
    s["inc_frac_d_inflated"] = inf["incompleteness"]["frac_d_inflated"]
    s["false_absence_hq"] = inf["incompleteness"]["false_absence_rate_hq"]
    s["false_absence_all"] = inf["incompleteness"]["false_absence_rate_all"]
    for vn in ("incomplete_prok", "oracle_genotypes"):
        sg = inf["sigma"][vn]
        s[f"I_{vn}"] = sg["point"]["I"]
        ci = sg["ci"].get("I", [np.nan, np.nan])
        s[f"I_{vn}_lo"], s[f"I_{vn}_hi"] = ci
        s[f"Sigma_leaky_{vn}"] = sg["point"]["leaky"]; s[f"Sigma_private_{vn}"] = sg["point"]["private"]
        s[f"Sigma_inc_{vn}"] = sg["point"]["inc"]
        s[f"prev_leaky_{vn}"] = sg["prevalence"].get("leaky"); s[f"prev_private_{vn}"] = sg["prevalence"].get("private")
    for vn in ("complete_prok", "incomplete_prok"):
        oe = inf["fanin_oe"][vn]
        for mode_ in ("uncond", "cond"):
            n = sum(v[f"n_{mode_}"] for v in oe.values())
            for q in ("obs", "exp", "oe"):
                s[f"{q}_{mode_}_{vn}"] = (sum((v[f"{q}_{mode_}"] or 0) * v[f"n_{mode_}"] for v in oe.values()
                                              if v[f"n_{mode_}"]) / n) if n else np.nan
    p4 = inf["p4"]["incomplete_prok"]
    for k in ("n_units", "obs_P_minus_W", "obs_donor_standardised", "mean_donors_W", "mean_donors_P",
              "n_truth_units", "truth_P_minus_W", "truth_fanin_P_minus_W"):
        s[f"p4_{k}"] = p4.get(k)
    e = inf["eukaryote"]
    s["euk_share_uptake"] = e["cell_weighted_euk_share_of_uptake"]
    s["euk_frac_sets"] = e["frac_feasible_sets_with_euk"]
    for mode in ("emergent", "lognormal", "powerlaw"):
        r = inf["rare_donors"][mode]
        s[f"member_detected_{mode}"] = r["frac_detected_at_time"]
        s[f"member_in_core_{mode}"] = r["frac_in_core"]
        s[f"witness_outside_core_{mode}"] = r["frac_witness_member_outside_core"]
    s["member_in_catalogue"] = inf["rare_donors"]["emergent"]["frac_in_catalogue"]
    s["n_mags"] = inf["n_mags"]
    return s


def auroc(pos, neg):
    pos, neg = [x for x in pos if x == x], [x for x in neg if x == x]
    if not pos or not neg:
        return np.nan
    return float(np.mean([(p > n) + 0.5 * (p == n) for p in pos for n in neg]))


def aggregate(root):
    rows = []
    for block, c, nseed in configs():
        name = cname(block, c)
        for sd in range(nseed):
            seed = seed_of(name, sd)
            rd = f"{root}/runs/{name}/s{seed}"
            if not os.path.exists(f"{rd}/inference.json"):
                continue
            rows.append(dict(block=block, config=name, seed=seed, **c, **summarise(rd)))
    os.makedirs(f"{root}/results", exist_ok=True)
    keys = list(rows[0].keys())
    with open(f"{root}/results/runs.tsv", "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join(fmt(r.get(k)) for k in keys) + "\n")
    # config means
    byc = defaultdict(list)
    for r in rows:
        byc[r["config"]].append(r)
    num = [k for k in keys if isinstance(rows[0][k], (int, float, np.floating)) and k != "seed"]
    with open(f"{root}/results/config_means.tsv", "w") as f:
        f.write("config\tn_seeds\t" + "\t".join(num) + "\n")
        for cn, rs in byc.items():
            f.write(cn + f"\t{len(rs)}\t" + "\t".join(fmt(np.nanmean([np.nan if r[k] is None else r[k] for r in rs]))
                                                    for k in num) + "\n")
    # signature analysis: AUROC of each observable for iota=1 vs iota=0 (block A matched configs),
    # and its Spearman correlation with true fan-in among iota=0 runs only
    from scipy.stats import spearmanr
    A = [r for r in rows if r["block"] == "A" and r["iota"] in (0.0, 1.0)]
    obs_keys = ["Sigma_inc_incomplete_prok", "I_incomplete_prok", "I_oracle_genotypes",
                "oe_cond_incomplete_prok", "oe_uncond_incomplete_prok", "obs_cond_incomplete_prok",
                "max_leaky_in_common_genotype", "mean_leaky_carried", "frac_cells_OX_and_ANX",
                "p4_obs_P_minus_W"]
    sig = []
    for k in obs_keys:
        pos = [r[k] for r in A if r["iota"] == 1.0 and r[k] is not None]
        neg = [r[k] for r in A if r["iota"] == 0.0 and r[k] is not None]
        z = [(r[k], r["nec_plexi_any"]) for r in A if r["iota"] == 0.0 and r[k] is not None and r[k] == r[k]
             and r["nec_plexi_any"] == r["nec_plexi_any"]]
        rho = spearmanr(*zip(*z))[0] if len(z) > 4 else np.nan
        z1 = [(r[k], r["nec_plexi_any"]) for r in A if r[k] is not None and r[k] == r[k]
              and r["nec_plexi_any"] == r["nec_plexi_any"]]
        rho_all = spearmanr(*zip(*z1))[0] if len(z1) > 4 else np.nan
        sig.append((k, auroc(pos, neg), np.nanmean(pos) if pos else np.nan, np.nanmean(neg) if neg else np.nan,
                    rho, rho_all))
    with open(f"{root}/results/signatures.tsv", "w") as f:
        f.write("observable\tAUROC_iota1_vs_iota0\tmean_iota1\tmean_iota0\tspearman_with_true_plexi_iota0_only\tspearman_with_true_plexi_all\n")
        for s_ in sig:
            f.write("\t".join(fmt(x) for x in s_) + "\n")
    return rows


def seed_of(name, rep):
    import zlib
    return 20261118 + (zlib.crc32(name.encode()) % 100000) * 10 + rep


def fmt(x):
    if x is None:
        return "nan"
    if isinstance(x, (float, np.floating)):
        return "nan" if x != x else f"{x:.4g}"
    return str(x)


if __name__ == "__main__":
    root = sys.argv[1]
    os.makedirs(root, exist_ok=True)
    if len(sys.argv) > 2 and sys.argv[2] == "agg":
        aggregate(root); sys.exit()
    nw = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    jobs = []
    with open(f"{root}/seeds.tsv", "w") as f:
        f.write("config\tblock\tseed\tparams_json\n")
        for block, c, nseed in configs():
            name = cname(block, c)
            for rep in range(nseed):
                sd = seed_of(name, rep)
                jobs.append((root, name, sd, c))
                f.write(f"{name}\t{block}\t{sd}\t{json.dumps(c)}\n")
    # longest first (D=1.0 has more diffusion substeps)
    jobs.sort(key=lambda j: -j[3]["D"])
    t0 = time.time()
    with mp.Pool(nw) as pool:
        for i, (rd, st) in enumerate(pool.imap_unordered(job, jobs)):
            print(f"[{time.time() - t0:6.0f}s] {i + 1}/{len(jobs)} {rd} {st}", flush=True)
    aggregate(root)
    print(f"done in {time.time() - t0:.0f}s")
