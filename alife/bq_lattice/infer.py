"""Apply simple versions of the team's inference to one run's synthetic observations and score
them against that run's ground truth.

usage: python3 -I infer.py RUNDIR        -> RUNDIR/inference.json, RUNDIR/inference_calls.tsv

Inference reads only obs/* (plus truth/per_mag_cds_complete.tsv for the explicitly labelled
"complete genome" oracle variants).  truth/* is read only in the scoring section.
"""
import sys, os, json, itertools, math
from collections import defaultdict
import numpy as np

FL = 7
LEAKY_COMPAT = [0, 1, 2, 3, 4]
PRIVATE = [7, 8, 9, 10]
INC = (5, 6)
DET = 10            # reads for detection
CORE_REL = 0.01     # abundant core: max relative abundance >= 1% in some sample at that time
rng = np.random.default_rng(20261118)


def read_tsv(path):
    import gzip
    if not os.path.exists(path) and os.path.exists(path + ".gz"):
        f = gzip.open(path + ".gz", "rt")
    else:
        f = open(path)
    with f:
        hdr = f.readline().rstrip("\n").split("\t")
        return [dict(zip(hdr, l.rstrip("\n").split("\t"))) for l in f if l.strip()]


def min_cover(J, masks, cap=4):
    if J == 0:
        return 0
    ms = {m & J for m in masks if m & J}
    masks = [m for m in ms if not any(m != o and m & o == m for o in ms)]
    for k in range(1, cap + 1):
        for c in itertools.combinations(masks, k):
            u = 0
            for m in c:
                u |= m
            if u & J == J:
                return k
    return None


def module_calls(cds_path, mods, euk=frozenset()):
    """mag -> (present, absent, unresolved) masks.  Prokaryotes: 80% / 25% KO rules.
    Eukaryote bins (MetaEuk-like, uncalibrated detectability; revision after Carbon's review):
    present if >= 80% of KOs, otherwise UNRESOLVED (never called absent)."""
    kos = defaultdict(set)
    for r in read_tsv(cds_path):
        ko = r.get("KO", r.get("ko", ""))
        if ko:
            kos[r["mag_id"]].add(ko)
    out = {}
    for m, ks in kos.items():
        pres = absn = unr = 0
        for j, (name, cls, mk) in enumerate(mods):
            f = len(ks & mk) / len(mk)
            if f >= 0.8:
                pres |= 1 << j
            elif f <= 0.25 and m not in euk:
                absn |= 1 << j
            else:
                unr |= 1 << j
        out[m] = (pres, absn, unr)
    return out


def logit(p):
    p = min(max(p, 1e-4), 1 - 1e-4)
    return math.log(p / (1 - p))


def p_fanin_null(J, q, N):
    """Prop. 1 (fan-in-null-theory.md) with independent functions, heterogeneous q_j, N donors."""
    js = [j for j in range(FL) if (J >> j) & 1]
    a = np.prod([q[j] for j in js])
    s = 0.0
    for k in range(1, len(js) + 1):
        for S in itertools.combinations(js, k):
            s += (-1) ** (k + 1) * np.prod([1 - q[j] for j in S]) ** N
    return (1 - a) ** N - s, 1 - (1 - a) ** N


def main(rd):
    O = f"{rd}/obs"
    mods = [(r["module"], r["class"], set(r["kos"].split(","))) for r in read_tsv(f"{O}/modules.tsv")]
    qual = {r["mag_id"]: r for r in read_tsv(f"{O}/mag_quality.tsv")}
    samples = read_tsv(f"{O}/samples.tsv")
    sid = [s["sample_id"] for s in samples]
    smeta = {s["sample_id"]: s for s in samples}
    mags = sorted(qual)
    euk = {m for m in mags if qual[m]["domain"] == "Eukaryota"}
    calls_inc = module_calls(f"{O}/per_mag_cds.tsv", mods, frozenset(euk))
    calls_cmp = module_calls(f"{rd}/truth/per_mag_cds_complete.tsv", mods)   # oracle: complete
    for m in mags:
        calls_inc.setdefault(m, (0, (1 << len(mods)) - 1, 0))
    comp = {m: float(qual[m]["completeness"]) if qual[m]["completeness"] != "NA" else np.nan for m in mags}

    def load_ab(mode):
        rows = read_tsv(f"{O}/abundance_{mode}.tsv")
        A = {r["mag_id"]: np.array([int(r[s]) for s in sid]) for r in rows}
        tot = sum(A.values())
        return {m: A[m] for m in mags}, tot
    ab = {mode: load_ab(mode) for mode in ("emergent", "lognormal", "powerlaw")}
    unb = defaultdict(set)
    for r in read_tsv(f"{O}/unbinned_ko.tsv"):
        unb[r["sample_id"]].add(r["ko"])
    unb_mask = {}
    for s in sid:
        mk = 0
        for j, (_, _, kos) in enumerate(mods):
            if len(unb[s] & kos) / len(kos) >= 0.8:
                mk |= 1 << j
        unb_mask[s] = mk
    times = sorted({int(s["time"]) for s in samples})
    # core per time (emergent abundance): max relative abundance >= CORE_REL at that time
    def core_at(mode, t):
        A, tot = ab[mode]
        cols = [i for i, s in enumerate(sid) if int(smeta[s]["time"]) == t]
        return {m for m in mags if max(A[m][i] / max(tot[i], 1) for i in cols) >= CORE_REL}

    LM = (1 << FL) - 1
    VARIANTS = {
        # name: (calls, include eukaryotes as donors, HQ recipient filter, core-only donors, mode)
        "complete_all":        (calls_cmp, True, False, False, "emergent"),
        "complete_prok":       (calls_cmp, False, False, False, "emergent"),
        "incomplete_all":      (calls_inc, True, True, False, "emergent"),
        "incomplete_prok":     (calls_inc, False, True, False, "emergent"),   # P2-like primary
        "incomplete_prok_core": (calls_inc, False, True, True, "emergent"),
        "incomplete_prok_lognormal": (calls_inc, False, True, False, "lognormal"),
        "incomplete_prok_powerlaw":  (calls_inc, False, True, False, "powerlaw"),
    }
    per_variant, call_rows = {}, []
    sample_calls = {}
    d_def = {}
    for vn, (calls, use_euk, hq, core_only, mode) in VARIANTS.items():
        A, tot = ab[mode]
        recs = []
        for i, s in enumerate(sid):
            t = int(smeta[s]["time"])
            det = [m for m in mags if A[m][i] >= DET and (use_euk or m not in euk)]
            donors_pool = det
            if core_only:
                core = core_at(mode, t)
                donors_pool = [m for m in det if m in core]
            for X in det:
                if X in euk and not vn.startswith("complete"):
                    continue          # eukaryote completeness cannot be scored (P2 scope)
                if hq and not (comp[X] >= 90):
                    continue
                pres, absn, unr = calls[X]
                J = absn & LM
                Ju = unr & LM
                if J == 0 and Ju == 0:
                    recs.append((X, s, t, "autonomous", 0, 0, 0)); continue
                if Ju:
                    recs.append((X, s, t, "unresolved_recipient", None, J, 0)); continue
                dm = [calls[D][0] for D in donors_pool if D != X]
                amb = [calls[D][0] | calls[D][2] for D in det if D != X] + [unb_mask[s]]
                d = min_cover(J, dm)
                d_def[(vn, X, s)] = d             # definite-evidence d (for the null comparison)
                if d == 1:
                    recs.append((X, s, t, "resolved", 1, J, len(dm))); continue
                # permissive single-supplier check BEFORE accepting d >= 2 (P2 upper bound):
                # one ambiguous MAG (partial / eukaryote-unknown modules) or the unbinned pool
                if any(a_ & J == J for a_ in amb):
                    recs.append((X, s, t, "ambig_single_supplier", None, J, len(dm))); continue
                if d is not None:
                    recs.append((X, s, t, "resolved", d, J, len(dm))); continue
                st = "ambig_donor" if min_cover(J, amb) is not None else "no_donor"
                recs.append((X, s, t, st, None, J, len(dm)))
        sample_calls[vn] = recs
        # aggregate per MAG x time
        agg = defaultdict(list)
        for X, s, t, st, d, J, N in recs:
            agg[(X, t)].append((st, d))
        mt = {}
        for (X, t), L in agg.items():
            res = [d for st, d in L if st in ("resolved", "autonomous")]
            if res:
                frac2 = np.mean([d >= 2 for d in res])
                mt[(X, t)] = dict(call="plexi" if frac2 >= 0.5 else "not", frac_d2=frac2,
                                  d_med=float(np.median(res)), n_res=len(res), n=len(L))
            else:
                mt[(X, t)] = dict(call="unresolved", frac_d2=np.nan, d_med=np.nan, n_res=0, n=len(L),
                                  why=max(set(st for st, _ in L), key=[st for st, _ in L].count))
        per_variant[vn] = mt
        for (X, t), v in mt.items():
            call_rows.append((vn, X, t, v["call"], v["frac_d2"], v["d_med"], v["n_res"], v["n"]))

    # ------------------------------------------------ fan-in null (Prop. 1), SAME denominators
    def fanin_oe(vn):
        """Two matched comparisons (revision after Astra's review), definite evidence only:
        unconditional: units = eligible recipient-samples with J != 0; obs = [d_def >= 2]
                       vs Prop. 1 P(d >= 2) (no-donor units count as not fan-in on both sides);
        conditional:   units with a definite cover; obs = [d_def >= 2] vs P(d>=2)/P(covered)."""
        calls, use_euk, mode = VARIANTS[vn][0], VARIANTS[vn][1], VARIANTS[vn][4]
        A, _ = ab[mode]
        out = defaultdict(lambda: [0, 0.0, 0.0, 0, 0.0, 0.0])
        for X, s, t, st, d, J, N in sample_calls[vn]:
            if J == 0 or st == "unresolved_recipient":
                continue
            i = sid.index(s)
            det = [m for m in mags if A[m][i] >= DET and m != X and (use_euk or m not in euk)]
            Nd = len(det)
            q = [np.mean([(calls[m][0] >> j) & 1 for m in det]) if det else 0 for j in range(FL)]
            p2, p1 = p_fanin_null(J, q, Nd)
            dd = d_def.get((vn, X, s))
            b = "<5" if Nd < 5 else "5-9" if Nd < 10 else "10-19" if Nd < 20 else "20-39" if Nd < 40 else "40+"
            o = out[b]
            o[0] += 1; o[1] += dd is not None and dd >= 2; o[2] += p2
            if dd is not None and p1 + p2 > 0:
                o[3] += 1; o[4] += dd >= 2; o[5] += p2 / (p1 + p2)
        return {b: dict(n_uncond=v[0], obs_uncond=v[1] / v[0], exp_uncond=v[2] / v[0],
                        oe_uncond=(v[1] - v[2]) / v[0],
                        n_cond=v[3], obs_cond=v[4] / v[3] if v[3] else np.nan,
                        exp_cond=v[5] / v[3] if v[3] else np.nan,
                        oe_cond=(v[4] - v[5]) / v[3] if v[3] else np.nan)
                for b, v in out.items() if v[0]}

    # ------------------------------------------------ leaky vs private co-carriage excess I
    def sigma_stats(vn, oracle_counts=None):
        """Per-donor-genome co-carriage (revision after Carbon's review): for each pair (j,k) in a
        class, O = number of donor genomes carrying both; E = n_j n_k / N, its exact expectation
        under fixed-count column permutation.  Sigma_C = log(sum O + .5) - log(sum E + .5) over
        pairs and samples of the events; I = Sigma_leaky - Sigma_private; event bootstrap."""
        calls, use_euk = VARIANTS[vn][0], VARIANTS[vn][1]
        A, _ = ab["emergent"]
        per_event = defaultdict(lambda: defaultdict(lambda: [0.0, 0.0]))
        prev = defaultdict(list)
        for i, s in enumerate(sid):
            if oracle_counts is not None:
                rows = oracle_counts[i]
            else:
                det = [m for m in mags if A[m][i] >= DET and (use_euk or m not in euk)]
                rows = [calls[m][0] for m in det]
            N = len(rows)
            if N < 3:
                continue
            ev = smeta[s]["event_id"]
            for cls, pairs in (("leaky", list(itertools.combinations(LEAKY_COMPAT, 2))),
                               ("private", list(itertools.combinations(PRIVATE, 2))),
                               ("inc", [INC])):
                for j, k in pairs:
                    nj = sum((r >> j) & 1 for r in rows)
                    nk = sum((r >> k) & 1 for r in rows)
                    if cls != "inc":
                        prev[cls] += [nj / N, nk / N]
                    O_ = sum(((r >> j) & 1) and ((r >> k) & 1) for r in rows)
                    e = per_event[ev][cls]
                    e[0] += O_; e[1] += nj * nk / N
        def stat(evs):
            res = {}
            for cls in ("leaky", "private", "inc"):
                o = sum(per_event[ev][cls][0] for ev in evs); e = sum(per_event[ev][cls][1] for ev in evs)
                res[cls] = math.log(o + 0.5) - math.log(e + 0.5) if e > 0 else np.nan
            res["I"] = res["leaky"] - res["private"]
            return res
        evs = sorted(per_event)
        point = stat(evs)
        boots = [stat(list(rng.choice(evs, len(evs)))) for _ in range(200)] if evs else []
        ci = {k: [float(np.nanpercentile([b[k] for b in boots], 2.5)),
                  float(np.nanpercentile([b[k] for b in boots], 97.5))] for k in point} if boots else {}
        return dict(point=point, ci=ci, prevalence={c: float(np.mean(v)) for c, v in prev.items()})

    tc = np.load(f"{rd}/truth/sample_counts.npz")
    oracle_rows = []
    for row in tc["counts"]:
        gs = np.nonzero(row >= 3)[0]
        oracle_rows.append([int(g) & ((1 << 11) - 1) for g in gs if not (int(g) >> 11) & 1])
    results = dict(fanin_oe={vn: fanin_oe(vn) for vn in ("complete_prok", "incomplete_prok", "incomplete_all")},
                   sigma={vn: sigma_stats(vn) for vn in ("complete_prok", "incomplete_prok")})
    results["sigma"]["oracle_genotypes"] = sigma_stats("complete_prok", oracle_counts=oracle_rows)

    # ================================================= SCORING (reads truth)
    mm = {r["mag_id"]: r["genotype"] for r in read_tsv(f"{rd}/truth/mag_map.tsv")}
    g2m = {g: m for m, g in mm.items()}
    nec = read_tsv(f"{rd}/truth/necessity_lineage.tsv")
    real = [r for r in read_tsv(f"{rd}/truth/realised_lineage.tsv") if r["scope"] == "all"]
    RES = ("autonomous", "not_plexi", "plexi")
    truth_nec = {(g2m.get(r["genotype"]), int(r["time"])): r for r in nec
                 if r["scope"] == "all" and r["status"] in RES}
    truth_real = {(g2m.get(r["genotype"]), int(r["time"])):
                  (int(r["d_real"]) >= 2 or int(r.get("c_real") or 1) >= 3) for r in real}

    def score(vn, truth, kind):
        """Resolved-only precision/recall plus bounds: unresolved calls set to not-plexi (low)
        or to plexi (high), so ambiguous donors are reported, not dropped."""
        c = defaultdict(int)
        sub = defaultdict(lambda: [0, 0])
        for key, tv in truth.items():
            if key[0] is None:
                c["truth_not_in_catalogue"] += 1; continue
            typ = None
            if kind == "nec":
                typ = tv.get("plexi_type") or None
                tv = tv["status"] == "plexi"
            v = per_variant[vn].get(key)
            if v is None:
                c["not_evaluated_plexi" if tv else "not_evaluated_not"] += 1; continue
            cl = v["call"]
            for mode_, call in (("res", None if cl == "unresolved" else cl == "plexi"),
                                ("low", cl == "plexi"), ("high", cl != "not")):
                if call is None:
                    continue
                c[f"{mode_}_tp"] += call and tv; c[f"{mode_}_fp"] += call and not tv
                c[f"{mode_}_fn"] += (not call) and tv; c[f"{mode_}_tn"] += (not call) and not tv
            c["unresolved"] += cl == "unresolved"
            if typ and cl != "unresolved":
                sub[typ][0] += cl == "plexi"; sub[typ][1] += 1
        out = dict(c)
        for mode_ in ("res", "low", "high"):
            tp, fp, fn = c[f"{mode_}_tp"], c[f"{mode_}_fp"], c[f"{mode_}_fn"]
            out[f"precision_{mode_}"] = tp / (tp + fp) if tp + fp else np.nan
            out[f"recall_{mode_}"] = tp / (tp + fn) if tp + fn else np.nan
        out["precision"], out["recall"] = out["precision_res"], out["recall_res"]
        out["recall_by_truth_type"] = {k: dict(called=v[0], n=v[1]) for k, v in sub.items()}
        return out
    results["score_necessity"] = {vn: score(vn, truth_nec, "nec") for vn in VARIANTS}
    results["score_realised"] = {vn: score(vn, truth_real, "real") for vn in VARIANTS}

    # genome capacity vs survival (Astra 248989851437154304): among truth rows, genome-autonomous
    # by ACTIVE capability vs. viable alone in the isolated-recipient counterfactual
    capsurv = defaultdict(int)
    for r in nec:
        if r["scope"] != "all" or r["status"] == "not_viable_in_situ":
            continue
        ga = r.get("active_needs", "") == "0"
        capsurv[f"genome_autonomous={ga}|survives_alone={r['status'] == 'autonomous'}"] += 1
    results["capacity_vs_survival"] = dict(capsurv)

    # incompleteness bias: same MAG x sample, incomplete vs complete (prokaryotes, no HQ filter)
    cmpd = {(X, s): (st, d) for X, s, t, st, d, J, N in sample_calls["complete_prok"]}
    diffs = []
    for X, s, t, st, d, J, N in sample_calls["incomplete_prok"]:
        c = cmpd.get((X, s))
        if c and st == "resolved" and c[0] == "resolved" and d is not None:
            diffs.append(d - c[1])
    nmod = nfa = 0
    hq_fa = []
    for m in mags:
        if m in euk:
            continue
        pi, ai, ui = calls_inc[m]; pc, ac, uc = calls_cmp[m]
        nmod += bin(pc).count("1"); nfa += bin(pc & ai).count("1")
        if comp[m] >= 90:
            hq_fa.append(bin(pc & ai).count("1") / max(bin(pc).count("1"), 1))
    results["incompleteness"] = dict(mean_d_bias=float(np.mean(diffs)) if diffs else np.nan,
                                     frac_d_inflated=float(np.mean([x > 0 for x in diffs])) if diffs else np.nan,
                                     n_pairs=len(diffs), false_absence_rate_all=nfa / max(nmod, 1),
                                     false_absence_rate_hq=float(np.mean(hq_fa)) if hq_fa else np.nan)

    # ------------------------------------------------ P4: same recipient x event, equal weights
    nec_scope = {(r["genotype"], int(r["time"]), r["scope"]): r for r in nec}

    def p4(vn):
        """Units = (recipient MAG X, event) with X resolved in BOTH fractions of the event (same
        recipient, same endpoint set J, same eligibility), each unit weight 1.
        observer: [d>=2 in P] - [d>=2 in W]; donor-standardised: the larger fraction's definite
        donor set is randomly subsampled to the smaller size (50 draws) before recomputing d.
        truth (on the subset of units at necessity time points): habitat-specific P-all verdicts
        for X's genotype (any plexi, and fan-in-only), same units, weight 1."""
        A, _ = ab[VARIANTS[vn][4]]
        calls, use_euk = VARIANTS[vn][0], VARIANTS[vn][1]
        by = defaultdict(dict)
        for X, s, t, st, d, J, N in sample_calls[vn]:
            if st in ("resolved", "autonomous"):
                by[(X, smeta[s]["event_id"])][smeta[s]["fraction"]] = (d, J, s, t)
        units = [(k, v) for k, v in by.items() if "W" in v and "P" in v]
        obs, std, tru, truf, nW, nP = [], [], [], [], [], []
        for (X, ev), v in units:
            dW, J, sW, t = v["W"]; dP, _, sP, _ = v["P"]
            obs.append((dP >= 2) - (dW >= 2))
            don = {}
            for fr, s_ in (("W", sW), ("P", sP)):
                i = sid.index(s_)
                don[fr] = [calls[m][0] for m in mags if A[m][i] >= DET and m != X and (use_euk or m not in euk)]
            nW.append(len(don["W"])); nP.append(len(don["P"]))
            if J:
                k = min(len(don["W"]), len(don["P"]))
                vals = {}
                for fr in ("W", "P"):
                    L = don[fr]
                    if len(L) == k:
                        dd = min_cover(J, L); vals[fr] = float(dd is not None and dd >= 2)
                    else:
                        r_ = []
                        for _ in range(50):
                            sub_ = [L[i_] for i_ in rng.choice(len(L), k, replace=False)]
                            dd = min_cover(J, sub_); r_.append(dd is not None and dd >= 2)
                        vals[fr] = float(np.mean(r_))
                std.append(vals["P"] - vals["W"])
            else:
                std.append(0.0)
            g = mm.get(X)
            rp_ = nec_scope.get((g, t, "particle")); rw_ = nec_scope.get((g, t, "water"))
            if rp_ and rw_ and rp_["status"] in RES and rw_["status"] in RES:
                tru.append((rp_["status"] == "plexi") - (rw_["status"] == "plexi"))
                truf.append((rp_.get("plexi_type") == "fanin") - (rw_.get("plexi_type") == "fanin"))
        f_ = lambda L: float(np.mean(L)) if L else np.nan
        return dict(n_units=len(units), obs_P_minus_W=f_(obs), obs_donor_standardised=f_(std),
                    mean_donors_W=f_(nW), mean_donors_P=f_(nP), n_truth_units=len(tru),
                    truth_P_minus_W=f_(tru), truth_fanin_P_minus_W=f_(truf))
    results["p4"] = {vn: p4(vn) for vn in ("complete_prok", "incomplete_prok", "incomplete_all")}

    # ------------------------------------------------ eukaryote contribution and rare donors
    sup = read_tsv(f"{rd}/truth/necessity_supports.tsv")
    plexi_keys = {(r["genotype"], r["time"]) for r in nec if r["scope"] == "all" and r["status"] == "plexi"}
    notp_keys = {(r["genotype"], r["time"]) for r in nec if r["scope"] == "all" and r["status"] == "not_plexi"}
    rw = [(float(r["euk_share"]) if r["euk_share"] != "nan" else np.nan, int(r["n_cells"])) for r in real]
    rw = [(a, n) for a, n in rw if not np.isnan(a)]
    results["eukaryote"] = dict(
        cell_weighted_euk_share_of_uptake=float(sum(a * n for a, n in rw) / max(sum(n for a, n in rw), 1)),
        frac_feasible_sets_with_euk=float(np.mean([s["set_has_eukaryote"] == "1" for s in sup])) if sup else np.nan,
        frac_plexi_recipients_all_sets_euk=float(np.mean([
            all(s["set_has_eukaryote"] == "1" for s in sup if (s["recipient"], s["time"]) == k)
            for k in plexi_keys])) if plexi_keys else np.nan)
    rare = {}
    for mode in ("emergent", "lognormal", "powerlaw"):
        A, tot = ab[mode]
        nd = nin = ncore = ndet = nd1 = nd1_out = 0
        for s_ in sup:
            t = int(s_["time"]); key = (s_["recipient"], s_["time"])
            cols = [i for i, s in enumerate(sid) if int(smeta[s]["time"]) == t]
            core = core_at(mode, t)
            for dg in s_["feasible_set"].split("+"):
                if dg == s_["recipient"]:
                    continue
                m = g2m.get(dg)
                inc = m is not None
                detd = inc and any(A[m][i] >= DET for i in cols)
                if key in plexi_keys:
                    nd += 1; nin += inc; ndet += detd; ncore += inc and m in core
                if key in notp_keys and s_["counterexample"] == "1":
                    nd1 += 1; nd1_out += not (inc and m in core)
        q = lambda a, b: a / b if b else np.nan
        rare[mode] = dict(n_member_slots_plexi=nd, frac_in_catalogue=q(nin, nd),
                          frac_detected_at_time=q(ndet, nd), frac_in_core=q(ncore, nd),
                          n_witness_slots=nd1, frac_witness_member_outside_core=q(nd1_out, nd1))
    results["rare_donors"] = rare
    results["n_mags"] = len(mags); results["n_euk_mags"] = len(euk)
    results["n_samples"] = len(sid)
    json.dump(results, open(f"{rd}/inference.json", "w"), indent=1, default=lambda x: None if x != x else float(x))
    with open(f"{rd}/inference_calls.tsv", "w") as f:
        f.write("variant\tmag_id\ttime\tcall\tfrac_d2\td_median\tn_resolved_samples\tn_samples\n")
        for r in call_rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    return results


if __name__ == "__main__":
    r = main(sys.argv[1])
    print(json.dumps({k: r[k] for k in ("score_necessity", "incompleteness")}, default=str)[:3000])
