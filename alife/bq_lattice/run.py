"""One bq_lattice run: evolve the world, log ground truth, emit synthetic observations.

usage: python3 -I run.py OUTDIR SEED [JSON-params]
Writes OUTDIR/truth/* (ground truth, never read by the inference except for scoring) and
OUTDIR/obs/* (synthetic observations in the shape of our real tables).
"""
import sys, os, json, time, itertools
from collections import defaultdict
import numpy as np
from scipy import ndimage
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from world import World, K, FL, OX, ANX, PRIV, NAMES, EUK, FULL, NG, B12, LEAKY

RUN_DEFAULTS = dict(t0=2000, dt=300, W=80, nec_last=2, H=100, Hm=50, min_cells=20, tau=0.1, k_max=4,
                    n_st=6, half=8, depth=20000, det_reads=10, bg_share=0.5, n_bg=3000,
                    n_cand=7)


def gname(g):
    return ("E" if g & EUK else "P") + format(g & FULL, "03x")


P_ACTIVE = {}   # world parameters, set in main(); used for active (not just encoded) capability


def needs(g):
    """Leaky products genotype g must obtain from outside, by ACTIVE capability (revision after
    Carbon's review): an encoded ANX counts as absent when the genome also carries OX and the
    self-blocked internal ANX supply Y*(1-leak)*(1-iota) falls below demand theta."""
    if g & EUK:
        return 1 << B12
    J = (~g) & ((1 << FL) - 1)
    if P_ACTIVE and (g >> OX) & 1 and (g >> ANX) & 1:
        p = P_ACTIVE
        if p["Y"] * (1 - p["leak"]) * (1 - p["iota"]) < p["theta"]:
            J |= 1 << ANX
    return J


def min_cover(J, masks, cap=4):
    """Smallest number of masks whose union covers J (0 if J empty); cap+1 if none <= cap."""
    if J == 0:
        return 0
    ms = {m & J for m in masks if m & J}
    masks = [m for m in ms if not any(m != o and m & o == m for o in ms)]   # drop dominated
    for k in range(1, cap + 1):
        for c in itertools.combinations(masks, k):
            u = 0
            for m in c:
                u |= m
            if u & J == J:
                return k
    return cap + 1


def scc_size(nodes, edges, r):
    """Size of the strongly connected component of r in a small digraph."""
    fwd = {n: set() for n in nodes}
    bwd = {n: set() for n in nodes}
    for a, b in edges:
        if a in fwd and b in fwd:
            fwd[a].add(b); bwd[b].add(a)

    def reach(adj):
        seen, st = {r}, [r]
        while st:
            x = st.pop()
            for y in adj[x]:
                if y not in seen:
                    seen.add(y); st.append(y)
        return seen
    return len(reach(fwd) & reach(bwd))


# ---------------------------------------------------------------- truth: realised transfers
def realised_truth(w, t, cnt, rp, out_edges):
    """Per recipient genotype x scope (all/water/particle): realised suppliers from tagged flux."""
    F = w.flux                                   # (FL, NT, NG*2)
    rows, graph = [], []
    gs = np.nonzero(cnt >= rp["min_cells"])[0]
    for r in gs:
        Jr = needs(int(r))
        for scope in ("all", "water", "particle"):
            cols = {"all": [2 * r, 2 * r + 1], "water": [2 * r], "particle": [2 * r + 1]}[scope]
            fr = F[:, :, cols].sum(2)            # (FL, NT)
            ncs = sum(w.nsum[c] for c in cols)
            Sj, shares, partner_J = {}, {}, 0
            tot_partner = tot_euk = tot_all = 0.0
            for j in range(FL):
                selft = [tg for tg in range(3, w.NT) if int(w.tag_geno[j, tg]) == int(r)]
                selff = sum(fr[j, tg] for tg in selft)
                if not (Jr >> j) & 1:
                    # v3 (Carbon 248999088913563648): encoded but short (O2-inhibited ANX,
                    # light-limited OX, high leak) -> counted if non-self uptake >= U_MIN/cell-step
                    if ncs <= 0 or (fr[j].sum() - fr[j, 0] - selff) / ncs < U_MIN:
                        continue
                tot = fr[j].sum() - selff
                if tot <= 0:
                    continue
                env, part_ = fr[j, 1], fr[j, 2:].sum() - selff
                tot_all += tot
                tot_partner += part_
                for tag in range(3, w.NT):
                    if w.tag_geno[j, tag] >= 0 and w.tag_geno[j, tag] & EUK:
                        tot_euk += fr[j, tag]
                if part_ / max(part_ + env, 1e-30) < 0.5:
                    continue                      # ambient-rescued for this product
                partner_J |= 1 << j
                srcs = set()
                for tag in range(2, w.NT):
                    if tag in selft:
                        continue
                    s = fr[j, tag] / part_
                    if s >= rp["tau"]:
                        src = "OTHER" if tag == 2 else int(w.tag_geno[j, tag])
                        srcs.add(src)
                        shares[(src, j)] = s
                        if scope == "all":
                            out_edges.append((t, gname(int(r)), NAMES[j],
                                              src if src == "OTHER" else gname(src),
                                              f"{fr[j, tag]:.4g}", f"{s:.4f}"))
                Sj[j] = srcs
            srcs_all = sorted({s for v in Sj.values() for s in v}, key=str)
            masks = [sum(1 << j for j in Sj if s in Sj[j]) for s in srcs_all]
            d = min_cover(partner_J, masks)
            best = max([min(shares.get((s, j), 0.0) for j in Sj) for s in srcs_all], default=np.nan) \
                if Sj else np.nan
            if scope == "all":
                for s in srcs_all:
                    if s != "OTHER":
                        graph.append((int(s), int(r)))
            rows.append(dict(time=t, genotype=gname(int(r)), scope=scope, n_cells=int(cnt[r]),
                             needs="".join(NAMES[j] + "," for j in range(FL) if (Jr >> j) & 1).rstrip(","),
                             partner_supplied=bin(partner_J).count("1"), d_real=d,
                             best_single_share=best, n_realised_sources=len(srcs_all),
                             euk_share=tot_euk / tot_all if tot_all > 0 else np.nan,
                             partner_share=tot_partner / tot_all if tot_all > 0 else np.nan))
    nodes = set(int(x) for x in gs)
    for row in rows:
        if row["scope"] == "all":
            gi = [g for g in gs if gname(int(g)) == row["genotype"]][0]
            row["c_real"] = scc_size(nodes, graph, int(gi))
    return rows


# ---------------------------------------------------------------- truth: necessity by removal
def counterfactual(w, st, keep, rp, rng):
    """Run the kept genotypes for H steps from state st; all other cells are inert ghosts.

    Ghosts (revision after review): occupy their sites and take up the products they lack (so
    competition for the endpoint metabolites is unchanged), but export nothing, emit/respire no
    O2, never divide and never die.  Kept cells behave normally, no evolution.  Provenance tags
    (sources restricted to kept genotypes) accumulate over steps Hm..H.
    Returns (gsum, nsum, flux, tag_geno)."""
    w.restore(st)
    w.ghost = (w.geno >= 0) & ~np.isin(w.geno, list(keep))
    saved = (w.tags, w.evolve, w.rng, w.t)
    w.evolve, w.rng = False, rng
    w.start_tags(ntop=min(len(keep), 12), only=keep)
    w.gsum = np.zeros(NG * 2); w.nsum = np.zeros(NG * 2)
    w.track_growth = False
    for s in range(rp["H"]):
        if s == rp["Hm"]:
            w.track_growth = True
            w.accum = True
        w.step()
    w.track_growth = False
    out = (w.gsum.copy(), w.nsum.copy(), w.flux.copy(), w.tag_geno.copy())
    w.stop_tags()
    w.tags, w.evolve, w.rng, w.t = saved
    w.ghost = None
    w.restore(st)
    return out


def viable(gn, g, scope, p):
    """Endpoint: mean per-capita growth probability over the horizon tail >= death rate."""
    gs, ns = gn[0], gn[1]
    mu = p["mu_e"] if g & EUK else p["mu"]
    if scope == "all":
        a, b = gs[2 * g] + gs[2 * g + 1], ns[2 * g] + ns[2 * g + 1]
    else:
        h = int(scope == "particle")
        a, b = gs[2 * g + h], ns[2 * g + h]
    return b > 0 and a / b >= mu


def supply_sets(res, K, Y, scope, tau, mode="deficit"):
    """Realised suppliers in a counterfactual: {product j: {Z in K: share >= tau}} for products
    Y obtains from partners (partner flux >= ENV flux); ENV-rescued products are omitted.
    mode="needs": only products absent by active genome capability (v2).
    mode="deficit" (v3, after Carbon 248999088913563648): additionally any ENCODED product whose
    realised non-self uptake in this scope is >= u_min per cell-step of Y, i.e. a carrier left short
    by local O2 (ANX) or particle light (OX) that is fed by partners."""
    F, tg = res[2], res[3]
    cols = {"all": [2 * Y, 2 * Y + 1], "water": [2 * Y], "particle": [2 * Y + 1]}[scope]
    ncs = sum(res[1][c] for c in cols)
    out = {}
    J = needs(Y)
    for j in range(FL):
        fr = F[j][:, cols].sum(1)
        nonself = fr[1] + sum(fr[t] for t in range(2, len(fr)) if int(tg[j, t]) != Y)
        if not (J >> j) & 1:
            if mode == "needs" or ncs <= 0 or nonself / ncs < U_MIN:
                continue
        part = sum(fr[t] for t in range(2, len(fr)) if int(tg[j, t]) != Y)
        if part <= 0 or part < fr[1]:
            continue
        sup = {int(tg[j, t]) for t in range(3, len(fr)) if tg[j, t] >= 0 and fr[t] / part >= tau
               and int(tg[j, t]) != Y}
        if sup:
            out[j] = sup
    return out


U_MIN = 0.05     # realised non-self uptake per cell-step (demand theta = 1) counted as a deficit


def min_covers(Sj):
    """All minimum-size supplier sets covering every product in Sj (dict j -> set)."""
    if not Sj:
        return 0, [()]
    allz = sorted({z for v in Sj.values() for z in v})
    for k in range(1, len(allz) + 1):
        cov = [c for c in itertools.combinations(allz, k) if all(set(c) & v for v in Sj.values())]
        if cov:
            return k, cov
    return None, []


def simple_cycle_through(X, arcs, minlen=3):
    """True if X lies on a simple directed cycle with >= minlen nodes."""
    adj = defaultdict(set)
    for a, b in arcs:
        adj[a].add(b)

    def dfs(u, path):
        for v in adj[u]:
            if v == X and len(path) >= minlen:
                return True
            if v not in path and v != X:
                if dfs(v, path | {v}):
                    return True
        return False
    return dfs(X, {X})


def bad_wiring(X, K, res, scope, rp, mode="deficit"):
    """Counterexample test on one feasible kept set K (P-all semantics, plexitrophy-computation.md
    §1.5 I1/I2): is there a wiring, within the realised suppliers, with d(X) <= 1 and X on no
    simple cycle of length >= 3?  Other members are wired by one greedy min cover each (a
    heuristic restriction of the wiring quantifier); X's covers are enumerated.
    Returns (bad, d, kappa_needed)."""
    SX = supply_sets(res, K, X, scope, rp["tau"], mode)
    d, covers = min_covers(SX)
    if d is None:
        return False, None, False
    if d >= 2:
        return False, d, False
    arcs_other = []
    for Y in K:
        if Y == X:
            continue
        SY = supply_sets(res, K, Y, "all", rp["tau"], mode)
        k, cv = min_covers(SY)
        if cv:
            arcs_other += [(z, Y) for z in cv[0]]
    for c in covers:
        arcs = arcs_other + [(z, X) for z in c]
        if not simple_cycle_through(X, arcs):
            return True, d, False
    return False, d, True        # d <= 1 but every cover puts X on a simple cycle >= 3


def necessity_truth(w, t, cnt, rp, realised_rows, rng):
    """P-all plexitrophy over tested kept sets (revision after Astra's review).

    Candidates: the n_cand most abundant genotypes (>= min_cells), NOT filtered by genome
    coverage.  Every kept set K with |K| <= k_max drawn from the candidates is run once as a
    counterfactual (cached; a run gives viability and realised supply for all members).
    K is FEASIBLE iff every member is viable (endpoint) over the horizon.
    For recipient X: autonomous if {X} is feasible; otherwise X is plexitrophic iff some K
    containing X is feasible and no feasible K containing X admits a counterexample wiring
    (d(X) <= 1 and X on no simple cycle >= 3).  The verdict is restricted to the tested sets,
    the stated medium, endpoint and horizon."""
    st = w.snapshot_state()
    p = w.p
    cache = {}

    def test(K):
        K = frozenset(K)
        if K not in cache:
            cache[K] = counterfactual(w, st, K, rp, rng)
        return cache[K]
    present = [int(g) for g in np.nonzero(cnt >= 5)[0]]
    full = test(present)
    cands = [int(g) for g in np.argsort(-cnt)[:rp["n_cand"]] if cnt[g] >= rp["min_cells"]]
    sets = [K for k in range(1, rp["k_max"] + 1) for K in itertools.combinations(cands, k)]
    for K in sets:
        test(K)
    rows, sup_rows = [], []
    for X in cands:
        base = dict(time=t, genotype=gname(X), n_cells=int(cnt[X]),
                    active_needs=bin(needs(X)).count("1"))
        for sc in ("all", "water", "particle"):
            rec = {**base, "scope": sc}
            if not viable(full, X, sc, p):
                rows.append({**rec, "status": "not_viable_in_situ"}); continue
            feas = []
            for K in sets:
                if X not in K:
                    continue
                res = cache[frozenset(K)]
                okX = viable(res, X, sc, p)
                okall = okX and all(viable(res, Y, "all", p) for Y in K if Y != X)
                if okall:
                    feas.append(K)
            if (X,) in feas:
                rows.append({**rec, "status": "autonomous", "plexi_pall": False, "d_min": 0,
                             "n_feasible": len(feas), "witness": gname(X)}); continue
            if not feas:
                rows.append({**rec, "status": "unresolved_no_feasible_set", "n_feasible": 0}); continue
            # v2 (needs-based) verdict and edge count, kept for comparison
            bads2 = [K for K in feas if bad_wiring(X, K, cache[frozenset(K)], sc, rp, "needs")[0]]
            e2 = sum(len(set().union(*supply_sets(cache[frozenset(K)], K, X, sc, rp["tau"], "needs").values()))
                     if supply_sets(cache[frozenset(K)], K, X, sc, rp["tau"], "needs") else 0 for K in feas)
            e3 = sum(len(set().union(*supply_sets(cache[frozenset(K)], K, X, sc, rp["tau"], "deficit").values()))
                     if supply_sets(cache[frozenset(K)], K, X, sc, rp["tau"], "deficit") else 0 for K in feas)
            bads, ds, kap = [], [], False
            for K in feas:
                b, d, kn = bad_wiring(X, K, cache[frozenset(K)], sc, rp)
                if d is not None:
                    ds.append(d)
                kap |= kn
                if b:
                    bads.append(K)
                if sc == "all":
                    sup_rows.append((t, gname(X), "+".join(gname(g) for g in K), d, int(b), int(kn),
                                     int(any(g & EUK for g in K if g != X))))
            plexi = not bads
            typ = ("" if not plexi else "fanin" if not kap else
                   "cyclic" if all(d is not None and d <= 1 for d in ds) else "mixed")
            rows.append({**rec, "status": "plexi" if plexi else "not_plexi", "plexi_pall": plexi,
                         "status_v2_needs": "plexi" if not bads2 else "not_plexi",
                         "in_edges_v2_needs": e2, "in_edges_v3_deficit": e3,
                         "plexi_type": typ, "d_min": min(ds) if ds else np.nan,
                         "n_feasible": len(feas), "n_bad": len(bads),
                         "witness": "+".join(gname(g) for g in bads[0]) if bads else ""})
    return rows, sup_rows, len(cache)


# ---------------------------------------------------------------- observations
def sample_events(w, t, rp, rng):
    """Paired water / particle fractions from the same station window (one event id)."""
    L, half = w.L, rp["half"]
    p = w.p
    out = []
    psites = np.nonzero(w.part)[0]
    for k in range(rp["n_st"]):
        c = rng.choice(psites) if len(psites) else rng.integers(w.N)
        cx, cy = divmod(int(c), L)
        xs = (np.arange(cx - half, cx + half + 1) % L)
        ys = (np.arange(cy - half, cy + half + 1) % L)
        win = (xs[:, None] * L + ys[None, :]).ravel()
        ev = f"t{t:05d}_st{k}"
        for frac, sel in (("W", win[~w.part[win]]), ("P", win[w.part[win]])):
            g = w.geno[sel]
            sel = sel[g >= 0]; g = g[g >= 0].astype(np.int64)
            if len(g) < 5:
                continue
            # RNA-like activity per module: ANX efficiency depends on O2 and OX co-carriage
            o2 = w.O22[sel]
            iota = p["iota"]
            Bm = ((g[:, None] >> np.arange(K)) & 1).astype(float)
            act = Bm.copy()
            act[:, ANX] *= (1 - iota * Bm[:, OX]) / (1 + iota * o2 / p["K_O"])
            out.append(dict(sample_id=f"{ev}_{frac}", event_id=ev, time=t, station=k, fraction=frac,
                            n_cells=len(g), n_sites=len(sel), o2=float(o2.mean()),
                            counts=np.bincount(g, minlength=NG), act={int(x): act[g == x].sum(0)
                                                                     for x in np.unique(g)}))
    return out


def gene_layout(rng, n_hk=200, n_blank=145, op=5):
    """Ancestral prokaryotic chromosome: list of (kind, label). Modules are contiguous operons."""
    genes = [("hk", f"K2{i:04d}") for i in range(n_hk)] + [("blank", "") for _ in range(n_blank)]
    rng.shuffle(genes)
    pos = sorted(rng.choice(len(genes), K, replace=False), reverse=True)
    order = rng.permutation(K)
    for p_, j in zip(pos, order):
        genes[p_:p_] = [("mod", j, f"K1{j:02d}{i:02d}") for i in range(op)]
    return genes


def emit_observations(odir, samples, rp, rng, p):
    os.makedirs(odir, exist_ok=True)
    tot = np.zeros(NG)
    for s in samples:
        tot += s["counts"]
    layout = gene_layout(rng)
    glen = np.full(NG, len(layout), float)
    # --- reads under three abundance regimes with identical depth
    bg = {"emergent": None,
          "lognormal": rng.lognormal(0, 2.0, rp["n_bg"]),
          "powerlaw": rng.pareto(1.0, rp["n_bg"]) + 1}
    reads = {}
    for mode, wbg in bg.items():
        R = np.zeros((len(samples), NG), int)
        bgr = np.zeros(len(samples), int)
        for i, s in enumerate(samples):
            wl = s["counts"] * glen
            if wbg is None:
                pr = wl / wl.sum()
                R[i] = rng.multinomial(rp["depth"], pr)
            else:
                pl = (1 - rp["bg_share"]) * wl / wl.sum()
                pb = rp["bg_share"] * wbg / wbg.sum()
                draw = rng.multinomial(rp["depth"], np.concatenate([pl, pb]))
                R[i] = draw[:NG]; bgr[i] = draw[NG:].sum()
                if i == 0:
                    reads[mode + "_bgdetected"] = []
                reads[mode + "_bgdetected"].append(int((draw[NG:] >= rp["det_reads"]).sum()))
        reads[mode] = (R, bgr)
    R0 = reads["emergent"][0]
    det = R0 >= rp["det_reads"]
    cat = [int(g) for g in np.nonzero(det.sum(0) >= 2)[0]]
    cat = sorted(cat, key=lambda g: -R0[:, g].sum())
    magid = {g: f"mag_{i:04d}" for i, g in enumerate(cat)}
    # --- samples & abundance tables
    with open(f"{odir}/samples.tsv", "w") as f:
        f.write("sample_id\tevent_id\ttime\tstation\tfraction\tassays\tn_cells\tn_sites\tmean_o2\n")
        for s in samples:
            f.write(f"{s['sample_id']}\t{s['event_id']}\t{s['time']}\t{s['station']}\t{s['fraction']}"
                    f"\tDNA,RNA\t{s['n_cells']}\t{s['n_sites']}\t{s['o2']:.4f}\n")
    for mode in bg:
        R, bgr = reads[mode]
        with open(f"{odir}/abundance_{mode}.tsv", "w") as f:
            f.write("mag_id\t" + "\t".join(s["sample_id"] for s in samples) + "\n")
            for g in cat:
                f.write(magid[g] + "\t" + "\t".join(str(x) for x in R[:, g]) + "\n")
            f.write("unbinned_lattice\t" + "\t".join(
                str(R[i].sum() - R[i, cat].sum()) for i in range(len(samples))) + "\n")
            f.write("background_taxa\t" + "\t".join(str(x) for x in bgr) + "\n")
    # --- unbinned KO pool per sample (non-catalogue genotypes with >=1 read)
    with open(f"{odir}/unbinned_ko.tsv", "w") as f:
        f.write("sample_id\tko\n")
        for i, s in enumerate(samples):
            kos = set()
            for g in np.nonzero(R0[i] > 0)[0]:
                if int(g) in magid:
                    continue
                for gene in layout:
                    if gene[0] == "mod" and (int(g) >> gene[1]) & 1 and not (g & EUK):
                        kos.add(gene[2])
            for ko in sorted(kos):
                f.write(f"{s['sample_id']}\t{ko}\n")
    # --- per-MAG CDS tables (incomplete) and complete reference
    euk_extra = [("blank", "")] * 250 + [("hk", f"K3{i:04d}") for i in range(150)]
    fq = open(f"{odir}/mag_quality.tsv", "w")
    fq.write("mag_id\tdomain\tcompleteness\tcontamination\tn_contigs\tn_cds\n")
    fc = open(f"{odir}/per_mag_cds.tsv", "w")
    HDR = "mag_id\tcontig_id\tprotein_id\tfeature_index\tKO\n"   # p2_scorer schema
    fc.write(HDR)
    fl = open(f"{odir}/contig_lengths.tsv", "w")
    fl.write("mag_id\tcontig_id\tlength_bp\n")
    flc = open(f"{os.path.dirname(odir)}/truth/contig_lengths_complete.tsv", "w")
    flc.write("mag_id\tcontig_id\tlength_bp\n")
    fr = open(f"{os.path.dirname(odir)}/truth/per_mag_cds_complete.tsv", "w")
    fr.write(HDR)
    fm = open(f"{os.path.dirname(odir)}/truth/mag_map.tsv", "w")
    fm.write("mag_id\tgenotype\tdomain\ttrue_completeness\ttarget_completeness\n")
    for g in cat:
        m = magid[g]
        if g & EUK:
            genes = [("mod", j, f"K1{j:02d}{i:02d}") for j in range(K) if (g >> j) & 1 for i in range(5)]
            genes = genes + euk_extra
            keep = rng.random(len(genes)) < 0.35          # MetaEuk-like fragmentary calls
            ncont = 0
            for gi, gene in enumerate(genes):
                ko = gene[2] if gene[0] == "mod" else gene[1]
                fr.write(f"{m}\t{m}_c0\t{m}_c0_{gi + 1:05d}\t{gi + 1}\t{ko}\n")
                if keep[gi]:
                    fc.write(f"{m}\t{m}_c{ncont}\t{m}_c{ncont}_00001\t1\t{ko}\n")
                    fl.write(f"{m}\t{m}_c{ncont}\t{rng.integers(800, 3000)}\n"); ncont += 1
            flc.write(f"{m}\t{m}_c0\t{1000 * len(genes)}\n")
            fq.write(f"{m}\tEukaryota\tNA\tNA\t{ncont}\t{ncont}\n")
            fm.write(f"{m}\t{gname(g)}\tEukaryota\t{keep.mean():.3f}\tNA\n")
            continue
        genes = [gene for gene in layout if gene[0] != "mod" or (g >> gene[1]) & 1]
        n = len(genes)
        brk = np.nonzero(rng.random(n) < 1 / 20)[0]
        bounds = np.concatenate([[0], brk[brk > 0], [n]])
        contigs = [(bounds[i], bounds[i + 1]) for i in range(len(bounds) - 1) if bounds[i + 1] > bounds[i]]
        target = rng.uniform(0.9, 1.0) if rng.random() < 0.5 else rng.uniform(0.55, 0.9)
        keepc = np.ones(len(contigs), bool)
        for ci in rng.permutation(len(contigs)):
            size = contigs[ci][1] - contigs[ci][0]
            if (sum(b - a for (a, b), k in zip(contigs, keepc) if k) - size) / n < target:
                break
            keepc[ci] = False
        mk_tot = mk_kept = 0
        for ci, (a, b) in enumerate(contigs):
            for fi, gi in enumerate(range(a, b)):
                gene = genes[gi]
                ko = gene[2] if gene[0] == "mod" else gene[1]
                line = f"{m}\t{m}_c{ci}\t{m}_c{ci}_{fi + 1:05d}\t{fi + 1}\t{ko}\n"
                fr.write(line)
                if gene[0] == "hk" and int(gene[1][2:]) < 100:
                    mk_tot += 1; mk_kept += keepc[ci]
                if keepc[ci]:
                    fc.write(line)
        for ci, (a, b) in enumerate(contigs):       # ~1 kb per CDS
            flc.write(f"{m}\t{m}_c{ci}\t{1000 * (b - a)}\n")
            if keepc[ci]:
                fl.write(f"{m}\t{m}_c{ci}\t{1000 * (b - a)}\n")
        comp = mk_kept / mk_tot
        kept_cds = sum(b - a for (a, b), k in zip(contigs, keepc) if k)
        fq.write(f"{m}\tBacteria\t{100 * comp:.1f}\t0.0\t{int(keepc.sum())}\t{kept_cds}\n")
        fm.write(f"{m}\t{gname(g)}\tBacteria\t{kept_cds / n:.3f}\t{target:.3f}\n")
    for fh in (fq, fc, fr, fm, fl, flc):
        fh.close()
    # --- module definitions (labels from "chemistry", pre-registered, not from co-occurrence)
    with open(f"{odir}/modules.tsv", "w") as f:
        f.write("module\tclass\tkos\n")
        for j in range(K):
            cls = "leaky" if j < 5 else ("incompatible_O2_producing" if j == OX else
                                         "incompatible_O2_sensitive" if j == ANX else "private")
            f.write(f"{NAMES[j]}\t{cls}\t" + ",".join(f"K1{j:02d}{i:02d}" for i in range(5)) + "\n")
    # --- RNA-like module activity per sample x MAG (same sample/event ids as DNA)
    with open(f"{odir}/expression.tsv", "w") as f:
        f.write("sample_id\tevent_id\tmag_id\tmodule\tactivity\n")
        for s in samples:
            for g, a in s["act"].items():
                if g in magid:
                    for j in range(K):
                        if a[j] > 0:
                            f.write(f"{s['sample_id']}\t{s['event_id']}\t{magid[g]}\t{NAMES[j]}"
                                    f"\t{a[j] * rng.lognormal(0, 0.3):.3f}\n")
    # --- detection summary per regime for truth-side analyses
    summ = {m: {"bg_detected_mean": float(np.mean(reads.get(m + "_bgdetected", [0])))} for m in bg}
    return magid, {m: reads[m][0] for m in bg}, summ


# ---------------------------------------------------------------- abundance distribution
def tail_fit(x):
    """Compare power-law vs lognormal for abundances x >= xmin (continuous MLE, Vuong LR)."""
    from scipy import stats
    x = np.sort(np.asarray(x, float))
    x = x[x > 0]
    if len(x) < 10:
        return dict(n=len(x))
    xmin = x[0]
    a = 1 + len(x) / np.log(x / xmin).sum()
    ll_pl = np.log((a - 1) / xmin) - a * np.log(x / xmin)
    mu_, sd = np.log(x).mean(), np.log(x).std() + 1e-9
    # lognormal truncated at xmin
    z0 = (np.log(xmin) - mu_) / sd
    ll_ln = stats.norm.logpdf((np.log(x) - mu_) / sd) - np.log(x * sd) - np.log(max(1 - stats.norm.cdf(z0), 1e-12))
    dd = ll_ln - ll_pl
    R = dd.sum()
    z = R / (np.sqrt(len(x)) * dd.std() + 1e-12)
    ranks = np.arange(1, len(x) + 1)
    slope = np.polyfit(np.log(ranks), np.log(x[::-1]), 1)[0]
    return dict(n=len(x), alpha_pl=a, LR_lognormal_minus_powerlaw=R, vuong_z=z, rank_slope=slope,
                top10_share=x[::-1][:10].sum() / x.sum())


# ---------------------------------------------------------------- main
def main(outdir, seed, params):
    t_start = time.time()
    rp = {**RUN_DEFAULTS, **{k: v for k, v in params.items() if k in RUN_DEFAULTS}}
    wp = {k: v for k, v in params.items() if k not in RUN_DEFAULTS}
    ss = np.random.SeedSequence(seed)
    r_world, r_cf, r_obs = [np.random.default_rng(s) for s in ss.spawn(3)]
    w = World(wp, r_world)
    P_ACTIVE.clear(); P_ACTIVE.update(w.p)
    T = w.p["T"]
    os.makedirs(f"{outdir}/truth", exist_ok=True)
    os.makedirs(f"{outdir}/obs", exist_ok=True)
    stimes = list(range(rp["t0"], T + 1, rp["dt"]))
    nec_times = stimes[-rp["nec_last"]:]
    snaps, real_rows, edges, nec_rows, sup_rows, samples, ab_rows, ts_rows = [], [], [], [], [], [], [], []
    ncf = 0
    for t in range(T + 1):
        if t + rp["W"] in stimes:
            w.start_tags()
        if t + rp["W"] // 2 in stimes and w.tags:
            w.accum = True
            w.gsum = np.zeros(NG * 2); w.nsum = np.zeros(NG * 2); w.track_growth = True
        if t % 100 == 0:
            al = w.geno[w.geno >= 0]
            cnt = np.bincount(al, minlength=NG)
            nz = cnt[cnt > 0]
            ts_rows.append((t, len(al), len(nz), int((cnt[EUK:] > 0).sum()),
                            int(cnt[EUK:].sum()), float(np.exp(-(nz / nz.sum() * np.log(nz / nz.sum())).sum())),
                            *[float(((al >> j) & 1).mean()) if len(al) else 0 for j in range(K)]))
        if t in stimes:
            al = w.geno[w.geno >= 0]
            cnt = np.bincount(al, minlength=NG)
            snaps.append(w.geno.copy())
            real_rows += realised_truth(w, t, cnt, rp, edges)
            w.track_growth = False
            w.stop_tags()
            fit = tail_fit(cnt[cnt > 0])
            ab_rows.append(dict(time=t, **fit))
            samples += sample_events(w, t, rp, r_obs)
            if t in nec_times:
                rows, sr, nc = necessity_truth(w, t, cnt, rp, real_rows, r_cf)
                nec_rows += rows; sup_rows += sr; ncf += nc
        if t < T:
            w.step()
        if len(w.geno[w.geno >= 0]) == 0:
            break
    # ---- write truth
    np.savez_compressed(f"{outdir}/truth/lattice_snapshots.npz", times=np.array(stimes[:len(snaps)]),
                        geno=np.array(snaps, np.int32) if snaps else np.zeros((0, w.N), np.int32),
                        particle=w.part2, priv_map=w.priv_map)
    write_tsv(f"{outdir}/truth/realised_lineage.tsv", real_rows)
    with open(f"{outdir}/truth/realised_edges.tsv", "w") as f:
        f.write("time\trecipient\tproduct\tsource\tflux\tshare_of_partner_flux\n")
        for e in edges:
            f.write("\t".join(map(str, e)) + "\n")
    write_tsv(f"{outdir}/truth/necessity_lineage.tsv", nec_rows)
    with open(f"{outdir}/truth/necessity_supports.tsv", "w") as f:
        f.write("time\trecipient\tfeasible_set\td_X\tcounterexample\tcycle_needed\tset_has_eukaryote\n")
        for e in sup_rows:
            f.write("\t".join(map(str, e)) + "\n")
    write_tsv(f"{outdir}/truth/abundance_fit.tsv", ab_rows)
    with open(f"{outdir}/truth/timeseries.tsv", "w") as f:
        f.write("time\tn_cells\tn_genotypes\tn_euk_genotypes\tn_euk_cells\teff_n_genotypes\t"
                + "\t".join(f"carry_{n}" for n in NAMES) + "\n")
        for r in ts_rows:
            f.write("\t".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in r) + "\n")
    summ = {}
    if samples:
        magid, reads, summ = emit_observations(f"{outdir}/obs", samples, rp, r_obs, w.p)
        np.savez_compressed(f"{outdir}/truth/sample_counts.npz",
                            counts=np.array([s["counts"] for s in samples]),
                            sample_ids=np.array([s["sample_id"] for s in samples]))
    meta = dict(seed=seed, params=w.p, run_params=rp, runtime_s=time.time() - t_start,
                n_counterfactual_runs=ncf, final_cells=int((w.geno >= 0).sum()),
                extinct=bool((w.geno >= 0).sum() == 0), obs_summary=summ)
    json.dump(meta, open(f"{outdir}/meta.json", "w"), indent=1, default=float)
    return meta


def write_tsv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w") as f:
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join(fmt(r.get(k, "")) for k in keys) + "\n")


def fmt(x):
    if isinstance(x, (float, np.floating)):
        return "nan" if np.isnan(x) else f"{x:.4g}"
    return str(x)


if __name__ == "__main__":
    params = json.loads(sys.argv[3]) if len(sys.argv) > 3 else {}
    m = main(sys.argv[1], int(sys.argv[2]), params)
    print(json.dumps({k: m[k] for k in ("seed", "runtime_s", "final_cells", "extinct",
                                        "n_counterfactual_runs")}))
