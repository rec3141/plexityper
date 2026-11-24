"""Observer-side inference on the synthetic files only (no truth used until scoring):
(a) P6-like 'allocation follows power' test, (b) KO-rule exchange graph (P5 redox chains and
B12 public good), (c) co-occurrence graph and the depth-ordering ('lagged succession')
artefact, (d) candidate MEP signatures."""
import csv
import os
import numpy as np
from scipy.stats import spearmanr
from model import Model, default_params, SPECIES, S, NS, RXN, RI, NE, T

# observer's functional classes: marker KO -> reactions it may indicate (ambiguity kept)
CLASSES = {
    "coxA": ("K02274", ["AER"]), "dmdA": ("K17486", ["DMP"]),
    "narG": ("K00370", ["NAR", "NXR", "SNR"]), "nirK": ("K00368", ["NIR"]),
    "dsrA": ("K11180", ["SRB", "HSR"]), "hydA": ("K00532", ["FER"]),
    "mcrA": ("K00399", ["MET"]), "hyaA": ("K06282", ["HSR"]),
    "amoA/pmoA": ("K10944", ["AMO", "MOX"]), "hzsA": ("K20932", ["AMX"]),
    "soxB": ("K17224", ["SOX"]), "sqr": ("K17218", ["SOX", "SNR"]), "mmoX": ("K16157", ["MOX"]),
}
# exchange rules: intermediate -> (producer markers, consumer markers)
RULES = {
    "NO2": ({"narG", "amoA/pmoA"}, {"nirK", "hzsA", "narG"}),
    "H2": ({"hydA"}, {"mcrA", "hyaA"}),
    "H2S": ({"dsrA", "dmdA"}, {"soxB", "sqr"}),
    "CH4": ({"mcrA"}, {"amoA/pmoA", "mmoX"}),
}
REDOX = set(RULES)


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load(obsdir):
    d = {}
    for k in ["events", "chemistry", "abundance", "per_mag_cds", "mag_completeness", "rna_counts"]:
        d[k] = read_tsv(os.path.join(obsdir, k + ".tsv"))
    return d


def predicted_power(chem_rows):
    """Per sample: per-reaction power per unit enzyme (k*F_K*F_T*inh*(-dG)/n_e) and Monod
    only, computed from observed chemistry with the published-style rate law (parameter-free
    for the observer: same constants as the model's literature-style defaults)."""
    m = Model(default_params(NZ=len(chem_rows)))
    C = np.zeros((NS, len(chem_rows)))
    for k, r in enumerate(chem_rows):
        for sp in SPECIES:
            if sp in r:
                C[S[sp], k] = float(r[sp])
    q, dG, FT = m.specific_rates(C)
    power = q / NE[:, None] * np.maximum(-dG, 0)
    Cp = np.maximum(C, 0)
    mon = Cp / (Cp + m.ks[:, None])
    FK = np.ones_like(power)
    from model import SUBSTR
    for i in range(len(RXN)):
        for s in SUBSTR[i]:
            FK[i] *= mon[s]
    return power, FK


def mag_markers(cds_rows):
    mk = {}
    ko2cls = {}
    for c, (ko, _) in CLASSES.items():
        ko2cls.setdefault(ko, []).append(c)
    for r in cds_rows:
        for c in ko2cls.get(r["KO"], []):
            mk.setdefault(r["mag_id"], set()).add(c)
    cob = {r["mag_id"] for r in cds_rows if r["KO"] in ("K02232", "K02230", "K02233")}
    metH = {r["mag_id"] for r in cds_rows if r["KO"] == "K00548"}
    metE = {r["mag_id"] for r in cds_rows if r["KO"] == "K00549"}
    return mk, cob, metH, metE


def power_test(d):
    """Per sample Spearman between RNA class shares and (power, Monod-only, DNA share)."""
    chem = d["chemistry"]
    power, FK = predicted_power(chem)
    sid_idx = {r["sample_id"]: k for k, r in enumerate(chem)}
    classes = list(CLASSES)
    rna = np.zeros((len(chem), len(classes)))
    ko_cls = {v[0]: c for c, v in CLASSES.items()}
    for r in d["rna_counts"]:
        c = ko_cls.get(r["KO"])
        if c is not None and r["sample_id"] in sid_idx:
            rna[sid_idx[r["sample_id"]], classes.index(c)] += float(r["reads"])
    # DNA share per class: sum of abundance of MAGs carrying the marker
    mk, cob, metH, metE = mag_markers(d["per_mag_cds"])
    dna = np.zeros_like(rna)
    for r in d["abundance"]:
        if r["mag_id"] == "unbinned" or r["sample_id"] not in sid_idx:
            continue
        for c in mk.get(r["mag_id"], ()):
            dna[sid_idx[r["sample_id"]], classes.index(c)] += float(r["reads"])
    pw = np.zeros_like(rna)
    mo = np.zeros_like(rna)
    for ci, c in enumerate(classes):
        for rx in CLASSES[c][1]:
            pw[:, ci] += power[RI[rx]]
            mo[:, ci] += FK[RI[rx]]
    res = {"power": [], "monod": [], "dna": [], "power_x_dna": []}
    for k in range(len(chem)):
        if rna[k].sum() <= 0:
            continue
        y = rna[k] / rna[k].sum()
        for lab, x in (("power", pw[k]), ("monod", mo[k]), ("dna", dna[k]),
                       ("power_x_dna", pw[k] * dna[k])):
            if np.ptp(x) > 0 and np.ptp(y) > 0:
                res[lab].append(spearmanr(x, y).statistic)
    # unexploited power: share of predicted power density in classes with ~no DNA
    tot = (pw / (pw.sum(1, keepdims=True) + 1e-30))
    dshare = dna / (dna.sum(1, keepdims=True) + 1e-30)
    unexpl = float(np.mean((tot * (dshare < 1e-3)).sum(1)))
    return {k: float(np.nanmean(v)) if v else np.nan for k, v in res.items()}, unexpl


def present_mags(d, thr=1e-3):
    tot = {}
    for r in d["abundance"]:
        tot[r["sample_id"]] = tot.get(r["sample_id"], 0) + float(r["reads"])
    pres = {}
    prof = {}
    for r in d["abundance"]:
        if r["mag_id"] == "unbinned":
            continue
        f = float(r["reads"]) / max(tot[r["sample_id"]], 1)
        prof.setdefault(r["mag_id"], {})[r["sample_id"]] = f
        if f >= thr:
            pres[r["mag_id"]] = True
    return set(pres), prof


def rule_graph(d, thr=1e-3):
    """(v3, Carbon 248994691424894976) KO-rule edges are drawn only between MAGs that are
    both present (relative abundance >= thr) in the SAME sample, then aggregated over
    samples (union; per-edge sample counts kept).  v2 pooled presence over all depths."""
    mk, cob, metH, metE = mag_markers(d["per_mag_cds"])
    pres, prof = present_mags(d, thr)
    sids = sorted({s for g in prof for s in prof[g]})
    count = {}
    for s in sids:
        here = {g for g in pres if prof[g].get(s, 0) >= thr}
        for sp, (pm, cm) in RULES.items():
            P = [g for g in here if mk.get(g, set()) & pm]
            Cn = [g for g in here if mk.get(g, set()) & cm]
            for a in P:
                for b in Cn:
                    if a != b:
                        count[(a, b, sp)] = count.get((a, b, sp), 0) + 1
        aux = [g for g in here if g in metH and g not in metE and g not in cob]
        for a in here & cob:
            for b in aux:
                if a != b:
                    count[(a, b, "B12")] = count.get((a, b, "B12"), 0) + 1
    return set(count), pres, prof


def cooccurrence(prof, pres, rho_thr=0.8):
    mags = sorted(pres)
    sids = sorted({s for g in mags for s in prof[g]})
    X = np.array([[prof[g].get(s, 0) for s in sids] for g in mags])
    edges = set()
    if len(mags) < 3:
        return edges, mags, X
    rho = spearmanr(X.T).statistic
    for i in range(len(mags)):
        for j in range(i + 1, len(mags)):
            if rho[i, j] >= rho_thr:
                edges.add((mags[i], mags[j]))
    return edges, mags, X


def mean_path_length(edges):
    """Mean geodesic length over reachable ordered pairs on the SCC-condensed digraph."""
    from analysis import _scc
    nodes = sorted({a for a, b in edges} | {b for a, b in edges})
    if not nodes:
        return np.nan
    adj = {n: set() for n in nodes}
    for a, b in edges:
        adj[a].add(b)
    comps = _scc(nodes, adj)
    cid = {n: i for i, c in enumerate(comps) for n in c}
    cadj = {i: set() for i in range(len(comps))}
    for a, b in edges:
        if cid[a] != cid[b]:
            cadj[cid[a]].add(cid[b])
    L = []
    for s in cadj:
        dist = {s: 0}
        fr = [s]
        while fr:
            nf = []
            for u in fr:
                for v in cadj[u]:
                    if v not in dist:
                        dist[v] = dist[u] + 1
                        nf.append(v)
            fr = nf
        L += [v for k, v in dist.items() if k != s]
    return float(np.mean(L)) if L else 0.0


def score(obsdir, truth_edges_typed, mag_key):
    """truth_edges_typed: set of (donor_type, recipient_type, resource) realised (share>=0.1).
    mag_key: mag_id -> type."""
    d = load(obsdir)
    pt, unexpl = power_test(d)
    edges, pres, prof = rule_graph(d)
    t2m = {v: k for k, v in mag_key.items()}
    truth = {(t2m[a], t2m[b], r) for a, b, r in truth_edges_typed if a in t2m and b in t2m}
    rel = REDOX | {"B12"}
    truth_r = {e for e in truth if e[2] in rel}
    inf = edges
    tp = len(inf & truth_r)
    prec = tp / len(inf) if inf else np.nan
    rec = tp / len(truth_r) if truth_r else np.nan
    redox_inf = {(a, b) for a, b, r in inf if r in REDOX}
    redox_true = {(a, b) for a, b, r in truth_r if r in REDOX}
    co, mags, X = cooccurrence(prof, pres)
    und_true = {frozenset((a, b)) for a, b, r in truth}
    co_tp = sum(1 for e in co if frozenset(e) in und_true)
    # depth-ordering artefact: co-occurrence edges between depth-adjacent MAGs
    sids = sorted({s for g in mags for s in prof[g]})
    depth = np.array([int(s.rsplit("_z", 1)[1]) for s in sids])
    cent = {g: float((X[i] * depth).sum() / max(X[i].sum(), 1e-12)) for i, g in enumerate(mags)}
    rank = {g: k for k, g in enumerate(sorted(mags, key=lambda g: cent[g]))}
    adj = [e for e in co if abs(rank[e[0]] - rank[e[1]]) == 1]
    adj_false = [e for e in adj if frozenset(e) not in und_true]
    return dict(power_rho=pt["power"], monod_rho=pt["monod"], dna_rho=pt["dna"],
                power_x_dna_rho=pt["power_x_dna"], unexploited_power=unexpl,
                n_present=len(pres), rule_edges=len(inf), true_edges=len(truth_r),
                rule_precision=prec, rule_recall=rec,
                redox_pathlen_inferred=mean_path_length(redox_inf),
                redox_pathlen_true=mean_path_length(redox_true),
                cooc_edges=len(co), cooc_precision=(co_tp / len(co)) if co else np.nan,
                cooc_depth_adjacent=len(adj), cooc_depth_adjacent_false=len(adj_false))
