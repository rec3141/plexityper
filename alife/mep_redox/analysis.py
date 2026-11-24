"""Ground-truth analysis of a population state: realised transfers (tagged-source
attribution with transport), division of labour, graded partner counts and cycles, and
counterfactual necessity by knock-out."""
import itertools
import numpy as np
from model import (Model, NS, NR, NU, S, SPECIES, RXN, RI, TYPES, T)

LIVING_EXCLUDE = {"ATM", "LAT"}          # abiotic source nodes


def _transport_matrix(m, s, lam):
    """Steady-state operator A for species s: A c = sources (per volume)."""
    p = m.p
    NZ, dz = p["NZ"], p["dz"]
    A = np.zeros((NZ, NZ))
    for i in range(NZ - 1):
        k = m.Kf[i] / dz / dz
        A[i, i] += k
        A[i + 1, i + 1] += k
        A[i, i + 1] -= k
        A[i + 1, i] -= k
    if m.top_mask[s]:
        A[0, 0] += m.Ktop / (dz * dz / 2)
    A[np.arange(NZ), np.arange(NZ)] += lam + p["flush"]
    if s == S["POC"]:
        w = p["w_sink"] / dz
        A[np.arange(NZ), np.arange(NZ)] += w + p["k_hyd"]
        A[np.arange(1, NZ), np.arange(NZ - 1)] -= w
    return A


def attribution(m, y):
    """Realised transfers.  Returns list of dict(resource, donor, recipient, flux) with flux
    in mmol m^-2 d^-1 (column integrated), plus self-use and channelled flows.
    Donor nodes: type names (incl. ALG), 'ATM' (top boundary), 'LAT' (lateral inflow).
    Organic C is traced to its ultimate producer through POC (necromass) and hydrolysis."""
    p = m.p
    C, B = m.unpack(y)
    P = m.tendencies(0.0, y, parts=True)
    Jtot = P["J"] + P["Jint"]
    NZ, dz = p["NZ"], p["dz"]
    names = m.types
    # per-type production and consumption of each species (pool exchange only)
    prod = np.einsum("jrz,rs->jsz", Jtot, np.maximum(NU, 0))
    cons = np.einsum("jrz,rs->jsz", Jtot, np.maximum(-NU, 0))
    chan = np.zeros_like(prod)
    for (sp, js, amt) in P["chan"]:
        chan[js, sp] += amt
    prod -= chan
    cons -= chan
    ia = m.ialg
    if ia is not None:
        prod[ia, S["DOC"]] += P["ex"]
        prod[ia, S["DMSP"]] += P["dm"]
        prod[ia, S["O2"]] += P["G"] - P["resp"]
        cons[ia, S["B12"]] += p["q_B12"] * p["eps_alg"] * P["G"]
    prod[:, S["B12"]] += P["b12rel"]
    if ia is not None:                      # (v3) only algal necromass enters POC
        prod[ia, S["POC"]] += p["f_necro"] * P["mortB"][ia]
        cons[ia, S["NH4"]] += P["u4"]
        cons[ia, S["NO3"]] += P["u3"]
        cons[ia, S["SO4"]] += 0.2 * P["dm"]
        prod[ia, S["O2"]] += 0.6 * P["dm"] + 2.0 * P["u3"]
    out = []
    tags = {}
    order = ["POC", "DOC"] + [s for s in SPECIES if s not in ("POC", "DOC")]
    for sp in order:
        s = S[sp]
        c = np.maximum(C[s], 1e-12)
        ctot = cons[:, s].sum(0)
        extra = np.zeros(NZ)
        if sp == "B12":
            extra = p["k_B12dec"] * np.ones(NZ)
        lam = ctot / c + extra
        A = _transport_matrix(m, s, lam)
        src = {}
        for j, nm in enumerate(names):
            if prod[j, s].sum() > 0:
                src[nm] = prod[j, s].copy()
        if m.top_mask[s] and m.top_val[s] > 0:
            v = np.zeros(NZ)
            v[0] = m.Ktop * m.top_val[s] / (dz * dz / 2)
            src["ATM"] = v
        if m.ref[s] > 0:
            src["LAT"] = p["flush"] * m.ref[s] * np.ones(NZ)
        if sp == "DOC":
            # hydrolysis and sediment return carry POC tags
            hyd = P["hyd"]
            sinkb = p["w_sink"] * np.maximum(C[S["POC"]][-1], 0) / dz
            for nm, frac in tags["POC"].items():
                v = hyd * frac
                v = v.copy()
                v[-1] += sinkb * frac[-1]
                src[nm] = src.get(nm, 0) + v
        sol = {nm: np.linalg.solve(A, v) for nm, v in src.items()}
        tot = sum(sol.values()) if sol else np.zeros(NZ)
        frac = {nm: v / np.maximum(tot, 1e-30) for nm, v in sol.items()}
        tags[sp] = frac
        for j, rn in enumerate(names):
            u = cons[j, s]
            if u.sum() <= 0:
                continue
            for nm, fr in frac.items():
                fl = float((u * fr).sum() * dz)
                if fl > 1e-9:
                    out.append(dict(resource=sp, donor=nm, recipient=rn, flux=fl, kind="pool"))
        # channelled (intra-genome) flows
        for j, rn in enumerate(names):
            if chan[j, s].sum() > 0:
                out.append(dict(resource=sp, donor=rn, recipient=rn,
                                flux=float(chan[j, s].sum() * dz), kind="channel"))
    return out


def partner_counts(m, y, transfers, thr=0.1, min_bio_frac=0.002):
    """Graded realised partner structure per recipient type: shares per input, d_real (number
    of distinct living non-self donors with share >= thr of some input), Hill-1 effective
    donor numbers, and c (size of the strongly connected component in the thresholded
    realised graph).  NOT necessity: see knockouts()."""
    C, B = m.unpack(y)
    bio = B.sum(1) * m.p["dz"]
    tot = bio.sum()
    rows = []
    edges = set()
    for j, rn in enumerate(m.types):
        if bio[j] < min_bio_frac * tot:
            continue
        mine = [t for t in transfers if t["recipient"] == rn and t["kind"] == "pool"]
        res = sorted({t["resource"] for t in mine})
        donors = set()
        hill = []
        for r in res:
            fl = {}
            for t in mine:
                if t["resource"] == r:
                    fl[t["donor"]] = fl.get(t["donor"], 0) + t["flux"]
            ssum = sum(fl.values())
            sh = {k: v / ssum for k, v in fl.items()}
            pr = np.array([v for v in sh.values() if v > 0])
            hill.append(float(np.exp(-(pr * np.log(pr)).sum())))
            for k, v in sh.items():
                if v >= thr and k not in LIVING_EXCLUDE and k != rn:
                    donors.add(k)
                    edges.add((k, rn, r))
            rows.append(dict(recipient=rn, resource=r, intake=ssum,
                             shares=";".join(f"{k}:{v:.3f}" for k, v in
                                             sorted(sh.items(), key=lambda x: -x[1]))))
        rows.append(dict(recipient=rn, resource="*", intake=float(bio[j]),
                         shares="", d_real=len(donors), donors=",".join(sorted(donors)),
                         mean_hill=float(np.mean(hill)) if hill else 0.0))
    # SCCs of realised graph
    nodes = sorted({e[0] for e in edges} | {e[1] for e in edges})
    adj = {n: set() for n in nodes}
    for a, b, r in edges:
        adj[a].add(b)
    scc = _scc(nodes, adj)
    cmap = {n: len(c) for c in scc for n in c}
    for rr in rows:
        if rr["resource"] == "*":
            rr["c_real"] = cmap.get(rr["recipient"], 1)
            rr["scc"] = ",".join(sorted(next((c for c in scc if rr["recipient"] in c), [rr["recipient"]])))
    return rows, edges


def _scc(nodes, adj):
    index, low, onst, st, out = {}, {}, set(), [], []
    counter = [0]

    def strong(v):
        index[v] = low[v] = counter[0]
        counter[0] += 1
        st.append(v)
        onst.add(v)
        for w in adj.get(v, ()):
            if w not in index:
                strong(w)
                low[v] = min(low[v], low[w])
            elif w in onst:
                low[v] = min(low[v], index[w])
        if low[v] == index[v]:
            comp = []
            while True:
                w = st.pop()
                onst.discard(w)
                comp.append(w)
                if w == v:
                    break
            out.append(comp)
    for v in nodes:
        if v not in index:
            strong(v)
    return out


PATHWAYS = [  # (label, producer step, consumer step, intermediate)
    ("denitrification NO3->NO2->N2", "NAR", "NIR", "NO2"),
    ("nitrification NH4->NO2->NO3", "AMO", "NXR", "NO2"),
    ("fermentation->methanogenesis (H2)", "FER", "MET", "H2"),
    ("fermentation->H2 sulfate reduction", "FER", "HSR", "H2"),
    ("sulfate reduction->sulfide oxidation", "SRB", "SOX", "H2S"),
    ("methanogenesis->methanotrophy", "MET", "MOX", "CH4"),
]


def division_of_labour(m, y):
    """Per pathway: share of consumer-step flux done by genomes that also carry the
    producer step (stacked), flux-weighted depth centroids of the two steps, and the
    fraction of intermediate uptake that is self-supplied (same genome) vs transferred."""
    C, B = m.unpack(y)
    P = m.tendencies(0.0, y, parts=True)
    Jt = P["J"] + P["Jint"]
    z = m.z
    rows = []
    for lab, a, b, sp in PATHWAYS:
        ia_, ib_ = RI[a], RI[b]
        Ja = Jt[:, ia_].sum(1)
        Jb = Jt[:, ib_].sum(1)
        if Jb.sum() <= 0 or Ja.sum() <= 0:
            rows.append(dict(pathway=lab, flux_prod=float(Ja.sum()), flux_cons=float(Jb.sum()),
                             stacked_share=np.nan, depth_prod=np.nan, depth_cons=np.nan,
                             depth_sep=np.nan, self_supply=np.nan, top_cons_types=""))
            continue
        stacked = m.M[:, ia_] * m.M[:, ib_]
        sh = float((Jb * stacked).sum() / Jb.sum())
        za = float((Jt[:, ia_].sum(0) * z).sum() / Jt[:, ia_].sum())
        zb = float((Jt[:, ib_].sum(0) * z).sum() / Jt[:, ib_].sum())
        # self-supply: per type and depth min(production, consumption) of the intermediate
        s = S[sp]
        pr = np.einsum("jz,->jz", Jt[:, ia_], max(NU[ia_, s], 0))
        cn = np.einsum("jz,->jz", Jt[:, ib_], max(-NU[ib_, s], 0))
        selfs = float(np.minimum(pr, cn).sum() / cn.sum())
        order = np.argsort(-Jb)[:3]
        rows.append(dict(pathway=lab, flux_prod=float(Ja.sum() * m.p["dz"]),
                         flux_cons=float(Jb.sum() * m.p["dz"]), stacked_share=sh,
                         depth_prod=za, depth_cons=zb, depth_sep=zb - za, self_supply=selfs,
                         top_cons_types=",".join(f"{m.types[k]}:{Jb[k] / Jb.sum():.2f}"
                                                 for k in order if Jb[k] > 0)))
    return rows


def knockouts(m, y, transfers, horizon=1000.0, fail_frac=0.1, min_bio_frac=0.002,
              min_share=0.02):
    """Counterfactual necessity (revised after Astra 248992826624753664).
    Same medium (forcing, boundaries) and horizon for every run, compared with an
    UNPERTURBED CONTROL integrated from the same state over the same horizon.  A recipient
    'fails' a removal if its column biomass at the horizon is < fail_frac of the control's.
    Singles: each living donor removed alone.  Groups: for each recipient x input resource,
    all living donors of that input removed together (OR-group).
    Reported per recipient: failures by direct donors (with a realised transfer to the
    recipient) and indirect ones (no realised transfer), solver failures separately, and
    d_lb_direct = minimum hitting set of the failing OR-groups of the recipient's own inputs
    plus failing direct donors; d_lb_upstream adds indirect (upstream) failures.  Both are
    lower bounds on the size of any viable support CONDITIONAL on monotone rescue (adding
    partners never hurts); neither is a demonstrated viable support, and only d_lb_direct
    speaks to fan-in (plexitrophy) rather than serial chains."""
    C, B = m.unpack(y)
    bio = B.sum(1) * m.p["dz"]
    tot = bio.sum()
    keep = [j for j in range(m.NT) if bio[j] >= min_bio_frac * tot]
    kname = [m.types[j] for j in keep]
    donors_of = {}
    for t in transfers:
        if t["kind"] != "pool" or t["donor"] in LIVING_EXCLUDE or t["donor"] == t["recipient"]:
            continue
        if t["recipient"] not in kname or t["donor"] not in kname:
            continue
        donors_of.setdefault((t["recipient"], t["resource"]), {}).setdefault(t["donor"], 0.0)
        donors_of[(t["recipient"], t["resource"])][t["donor"]] += t["flux"]
    intake = {}
    for t in transfers:
        if t["kind"] == "pool":
            intake[(t["recipient"], t["resource"])] = intake.get((t["recipient"], t["resource"]), 0) + t["flux"]
    groups = {}
    for (rn, res), d in donors_of.items():
        g = tuple(sorted(k for k, v in d.items() if v >= min_share * intake[(rn, res)]))
        if g:
            groups.setdefault(g, []).append((rn, res))

    def integrate(ko):
        mm = Model(m.p)
        mm.knockout = [m.types.index(k) for k in ko] if ko else None
        y0 = y.copy()
        if ko:
            Cb, Bb = mm.unpack(y0)
            Bb[mm.knockout] = 0.0
        try:
            sol = mm.run(y0, horizon, t_eval=[horizon])
            return sol.y[:, -1]
        except RuntimeError:
            return None
    control = integrate(())
    ctrl_bio = None if control is None else m.unpack(control)[1].sum(1) * m.p["dz"]
    runs = {}
    for g in sorted({(k,) for k in kname} | set(groups)):
        runs[g] = integrate(g)
    rows = []
    for j in keep:
        rn = m.types[j]
        direct = {d for (r2, res), dd in donors_of.items() if r2 == rn for d in dd}
        base = ctrl_bio[j] if ctrl_bio is not None else np.nan
        fail_single, solver_fail = [], []
        for k in kname:
            if k == rn:
                continue
            if runs[(k,)] is None or ctrl_bio is None:
                solver_fail.append(k)
                continue
            if m.unpack(runs[(k,)])[1][j].sum() * m.p["dz"] < fail_frac * base:
                fail_single.append(k)
        group_res = []
        for g, lst in groups.items():
            for (r2, res) in lst:
                if r2 != rn or rn in g:
                    continue
                if runs[g] is None or ctrl_bio is None:
                    group_res.append((res, g, "solver_fail"))
                    continue
                ok = m.unpack(runs[g])[1][j].sum() * m.p["dz"] >= fail_frac * base
                group_res.append((res, g, "persists" if ok else "fails"))
        fd = [k for k in fail_single if k in direct]
        fi = [k for k in fail_single if k not in direct]
        group_fail = [set(g) for (res, g, st) in group_res if st == "fails"]
        # (v3, Carbon 248994691424894976) direct-donor lower bound uses only failing OR-groups
        # of the recipient's own inputs and failing DIRECT donors; upstream closure (indirect
        # failures, e.g. a serial chain ROS -> ALG -> AOB) is reported separately.
        d_lb_direct = _min_hitting_set(group_fail + [{k} for k in fd])
        d_lb_upstream = _min_hitting_set(group_fail + [{k} for k in fail_single])
        rows.append(dict(recipient=rn, control_biomass=float(base),
                         baseline_biomass=float(bio[j]),
                         n_fail_single=len(fail_single), fail_single=",".join(fail_single),
                         n_fail_direct=len(fd), fail_direct=",".join(fd),
                         n_fail_indirect=len(fi), fail_indirect=",".join(fi),
                         n_solver_fail=len(solver_fail) + sum(1 for x in group_res if x[2] == "solver_fail"),
                         failed_groups=";".join(f"{res}:{'|'.join(g)}" for res, g, st in group_res if st == "fails"),
                         tested_groups=";".join(f"{res}:{'|'.join(g)}={st}" for res, g, st in group_res),
                         d_lb_direct=d_lb_direct, d_lb_upstream=d_lb_upstream))
    return rows, {",".join(g): (None if v is None else float(m.unpack(v)[1].sum() * m.p["dz"]))
                  for g, v in runs.items()}


def _min_hitting_set(sets):
    if not sets:
        return 0
    univ = sorted(set().union(*sets))
    for k in range(1, len(univ) + 1):
        for comb in itertools.combinations(univ, k):
            cs = set(comb)
            if all(cs & s for s in sets):
                return k
    return len(univ)


def free_energy_budget(m, y, t=0.0):
    """(v3, after Astra 248994794025959424 / Carbon 248994691424894976) Species-resolved
    chemical-potential budget of the column (kJ m^-2 d^-1 = T * J K^-1 m^-2 d^-1 / 1).
    System: all tracked solutes, POC and algal biomass; fixed species (CO2, N2, H2O, H+) are
    constant-potential reservoirs inside the system.  Organic pools carry Redfield organic N
    valued at the local NH4 potential, so remineralisation is non-energetic (consistent with
    the catabolic dG, which omits it).  B12 is excluded.  Chemotroph biomass is catalyst only.
    Terms: algal capture (chemical potential stored by photosynthesis: the photon free energy
    minus algal dissipation), boundary inputs/outputs per species, storage change, and the
    dissipation by chemotroph catalysis (T*sigma_chem), diffusion, sinking, hydrolysis,
    mortality conversion and sediment return.  residual = capture + in - out - storage -
    all dissipation terms."""
    from model import GF, FIXED, RN, ALL_SPECIES as ALLS, NU_ALL
    p = m.p
    NZ, dz = p["NZ"], p["dz"]
    C, B = m.unpack(y)
    P = m.tendencies(t, y, parts=True)
    dCv, dBv = m.tendencies(t, y)
    dC = dCv.reshape(NS, -1)
    RTl = 8.314e-3 * T

    def mu_of(sp, c):
        return GF[sp] + RTl * np.log(np.maximum(c, 1e-9) * 1e-6)
    mu = {sp: mu_of(sp, C[S[sp]]) for sp in SPECIES if sp != "B12"}
    for sp in ("DOC", "POC"):
        mu[sp] = mu[sp] + RN * mu["NH4"]
    ia = m.ialg
    Ba = np.maximum(B[ia], 0)
    mu_bio = mu_of("DOC", Ba) + RN * mu["NH4"]
    muf = {k: GF[k] + RTl * np.log(v) for k, v in FIXED.items()}
    sps = [sp for sp in SPECIES if sp != "B12"]
    # --- algal capture
    G, resp, dm, ex, u3, u4 = P["G"], P["resp"], P["dm"], P["ex"], P["u3"], P["u4"]
    grow = p["eps_alg"] * G
    o2net = (G - dm) - resp + 1.6 * dm + 2.0 * u3
    ch2o_net = G - dm - resp
    cap = (grow * mu_bio + ex * mu["DOC"] + dm * mu["DMSP"] + o2net * mu["O2"]
           - u4 * mu["NH4"] - u3 * mu["NO3"] - 0.2 * dm * mu["SO4"]
           - (G - resp) * muf["CO2"] - (ch2o_net + 0.8 * dm + u3) * muf["H2O"]
           - (0.4 * dm + 2 * u3) * muf["H+"]).sum() * dz
    # --- boundaries
    bin_, bout = {}, {}
    sig_diff = {}
    for sp in sps:
        s_ = S[sp]
        c = C[s_]
        mu_s = mu[sp]
        F = -m.Kf * (c[1:] - c[:-1]) / dz                  # downward diffusive flux
        sd = float((F * (mu_s[:-1] - mu_s[1:])).sum())
        if m.top_mask[s_]:
            val = m.top_val[s_]
            mval = mu_of(sp, np.array([val]))[0] + (RN * mu["NH4"][0] if sp in ("DOC", "POC") else 0)
            Ft = m.Ktop * (val - c[0]) / (dz / 2)
            sd += Ft * (mval - mu_s[0])
            if Ft >= 0:
                bin_[sp] = bin_.get(sp, 0) + Ft * mval
            else:
                bout[sp] = bout.get(sp, 0) - Ft * mval
        if m.ref[s_] > 0:
            mref = mu_of(sp, np.array([m.ref[s_]]))[0]
            bin_[sp] = bin_.get(sp, 0) + p["flush"] * m.ref[s_] * mref * NZ * dz
        bout[sp] = bout.get(sp, 0) + float((p["flush"] * np.maximum(c, 0) * mu_s).sum() * dz)
        sig_diff[sp] = sd
    # algal biomass transport / exports
    Fa = -m.Kf * (Ba[1:] - Ba[:-1]) / dz
    sig_diff["ALG"] = float((Fa * (mu_bio[:-1] - mu_bio[1:])).sum())
    wa = p["alg_sink"] * Ba
    sig_sink = float((wa[:-1] * (mu_bio[:-1] - mu_bio[1:])).sum())
    wp = p["w_sink"] * np.maximum(C[S["POC"]], 0)
    sig_sink += float((wp[:-1] * (mu["POC"][:-1] - mu["POC"][1:])).sum())
    sig_sed = float(wp[-1] * (mu["POC"][-1] - mu["DOC"][-1]))
    sig_hyd = float((P["hyd"] * (mu["POC"] - mu["DOC"])).sum() * dz)
    mort = P["mortB"][ia]
    sig_mort = float((p["f_necro"] * mort * (mu_bio - mu["POC"])).sum() * dz)
    bout["ALG"] = float(((1 - p["f_necro"]) * mort * mu_bio).sum() * dz + wa[-1] * mu_bio[-1]
                        + (p["flush"] * Ba * mu_bio).sum() * dz)
    bin_["ALG"] = float((p["immig"] * mu_bio).sum() * dz)
    # storage
    stor = sum(float((dC[S[sp]] * mu[sp]).sum() * dz) for sp in sps)
    if dBv is not None:                      # frozen algae (FixedEnzyme): no storage term
        stor += float((dBv.reshape(m.NT, -1)[ia] * mu_bio).sum() * dz)
    Jt = (P["J"] + P["Jint"]).sum(0)
    # fixed-species reservoirs are inside the system: their potential changes with the net
    # production by chemotroph reactions minus algal consumption
    fx = {}
    for k in FIXED:
        fx[k] = float((NU_ALL[:, ALLS.index(k)] @ Jt).sum() * dz)
    fx["CO2"] -= float((G - resp).sum() * dz)
    fx["H2O"] -= float((ch2o_net + 0.8 * dm + u3).sum() * dz)
    fx["H+"] -= float((0.4 * dm + 2 * u3).sum() * dz)
    reservoir = sum(fx[k] * muf[k] for k in FIXED)
    Tsig_chem = float((Jt * np.maximum(-P["dG"], 0)).sum() * dz)
    # reactions with dG > 0 carry no flux (F_T = 0); remin and channel are inside Jt
    tot_in = sum(bin_.values())
    tot_out = sum(bout.values())
    diss = Tsig_chem + sum(sig_diff.values()) + sig_sink + sig_sed + sig_hyd + sig_mort
    resid = cap + tot_in - tot_out - stor - reservoir - diss
    return dict(capture=cap, boundary_in=tot_in, boundary_out=tot_out, storage=stor,
                reservoir_net=reservoir,
                reservoir_by_species=";".join(f"{k}:{fx[k] * muf[k]:.4g}" for k in FIXED),
                Tsigma_chem=Tsig_chem, Tsigma_diffusion=sum(sig_diff.values()),
                Tsigma_sinking=sig_sink, Tsigma_sediment=sig_sed, Tsigma_hydrolysis=sig_hyd,
                Tsigma_mortality=sig_mort, residual=resid,
                residual_rel=resid / max(abs(Tsig_chem), 1e-12),
                in_by_species=";".join(f"{k}:{v:.4g}" for k, v in bin_.items() if abs(v) > 1e-9),
                out_by_species=";".join(f"{k}:{v:.4g}" for k, v in bout.items() if abs(v) > 1e-9),
                diff_by_species=";".join(f"{k}:{v:.4g}" for k, v in sig_diff.items() if abs(v) > 1e-6))
