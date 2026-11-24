"""P6 scorer (protocol protocols/P6-thermodynamic-allocation.md, frozen 0a9b44648826).

Per-event predictors are computed exactly as in Check 7 of scratch/build_and_verify_p6_table.py
(score_mep_redox_run), but without its per-run asserts, because eligibility and the coverage
gate are applied once at cohort level (table notes §5.1, §5.3). Regression test: with a fixed
phi = 0.05 and each run as its own block, D001 and D002 reproduce Check 7's equal-block scores.

Usage:
  python3 -I scratch/p6_scorer.py regress            # D001/D002 vs Check 7
  python3 -I scratch/p6_scorer.py benchmark          # pooled synthetic benchmark (27 blocks)
  python3 -I scratch/p6_scorer.py calibrate [reps]   # count-generating calibration
Outputs go to results/p6/.
"""
import csv
import hashlib
import json
import math
import os
import sys
import time

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.special import gammaln

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_and_verify_p6_table as bv  # noqa: E402

ROOT = os.path.dirname(HERE)
TABLE = os.path.join(ROOT, "protocols", "P6-pathway-table.tsv")
TABLE_SHA12 = "6c26362681b8"
SYN = os.path.join(ROOT, "alife", "mep_redox", "synthetic")
BLOCK_MAP = os.path.join(ROOT, "results", "p6", "synthetic_block_map.tsv")
OUT_DIR = os.path.join(ROOT, "results", "p6")

MODELS = ["M_cap", "M_capK", "M_capKT", "M_power"]
DELTA = 0.05          # frozen margin, nats per event
N_PANEL_MIN = 100
OXIC_UM = 10.0
SEED = 20261120
EPS_OBS = 1e-4
LOGPHI_BOUNDS = (math.log(1e-4), math.log(10.0))


def sha12(path):
    with open(path, "rb") as fh:
        return hashlib.sha256(fh.read()).hexdigest()[:12]


def load_table():
    assert sha12(TABLE) == TABLE_SHA12, "P6 pathway table hash mismatch"
    with open(TABLE) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


# ---------------------------------------------------------------- per-event extraction
def extract_run(run_dir, rows, sp):
    """Return list of event dicts: sample_id, cast_id, O2_uM, N_panel, counts, probs{model}."""
    active_tiers = {"core", "denit_tierC_neither", "h2_gated"}
    panel_rows = [r for r in rows if r["panel_tier"] in active_tiers]
    pid_to_idx = {r["pathway_id"]: i for i, r in enumerate(panel_rows)}
    R_len = len(panel_rows)
    assert R_len == 13

    prok_mags, mag_bp = set(), {}
    with open(os.path.join(run_dir, "mag_completeness.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row["domain"] != "Eukaryota":
                prok_mags.add(row["mag_id"])
                mag_bp[row["mag_id"]] = float(row["genome_size_bp"])

    mag_kos = {m: set() for m in prok_mags}
    cds_records = {}
    with open(os.path.join(run_dir, "per_mag_cds.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            m = row["mag_id"]
            if m in prok_mags:
                ko = row["KO"]
                if ko:
                    mag_kos[m].add(ko)
                cds_records[row["protein_id"]] = (m, ko, int(row["end"]) - int(row["start"]) + 1)

    classify_ko = bv.make_binned_ko_classifier(panel_rows)
    cds_to_pathway = {}
    hk_by_ko = {k: [] for k in bv.HOUSEKEEPING_PROK_KOS}
    for cds_id, (m, ko, length) in cds_records.items():
        if ko in hk_by_ko:
            hk_by_ko[ko].append((cds_id, m, length))
        pid, _ = classify_ko(ko, mag_kos[m])
        if pid is not None:
            cds_to_pathway[cds_id] = pid

    dna_mag_reads = {}
    with open(os.path.join(run_dir, "abundance.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            dna_mag_reads.setdefault(row["sample_id"], {})[row["mag_id"]] = float(row["reads"])

    rna_panel_counts = {}
    with open(os.path.join(run_dir, "rna_counts.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            pid = cds_to_pathway.get(row["protein_id"])
            if pid is not None:
                vec = rna_panel_counts.setdefault(row["sample_id"], [0] * R_len)
                vec[pid_to_idx[pid]] += int(row["reads"])

    def dG0_T(st, T):
        dG298 = sum(nu * sp[s][2] * bv.CAL_TO_KJ for s, nu in st.items())
        dH298 = sum(nu * sp[s][3] * bv.CAL_TO_KJ for s, nu in st.items())
        return (T / bv.T_REF) * dG298 + dH298 * (1.0 - T / bv.T_REF)

    events = []
    with open(os.path.join(run_dir, "chemistry.tsv")) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            sid = row["sample_id"]
            T = float(row["T_K"])
            pH = float(row["pH"])
            RT = bv.R_KJ * T
            aH = 10.0 ** (-pH)

            def get_uM(col, dl=0.01):
                v = float(row[col])
                return v if v > 0.0 else 0.5 * dl

            DOC_uM = get_uM("DOC")
            O2_uM = get_uM("O2")
            NO3_uM = get_uM("NO3")
            NO2_uM = get_uM("NO2")
            TAN_uM = get_uM("NH4")
            SO4_uM = get_uM("SO4")
            S_tot_uM = get_uM("H2S")
            H2_uM = get_uM("H2")
            CH4_uM = get_uM("CH4")
            Ka1_C = math.exp(-dG0_T({"CO2": -1, "H2O": -1, "HCO3-": 1, "H+": 1}, T) / RT)
            Ka2_C = math.exp(-dG0_T({"HCO3-": -1, "CO3-2": 1, "H+": 1}, T) / RT)
            Ka1_S = math.exp(-dG0_T({"H2S": -1, "HS-": 1, "H+": 1}, T) / RT)
            Ka_N = math.exp(-dG0_T({"NH4+": -1, "NH3": 1, "H+": 1}, T) / RT)
            a_CO2, a_HCO3 = bv.compute_carbonate_activities(row, aH, Ka1_C, Ka2_C)
            N2_M = float(row["N2"]) * 1e-6 if "N2" in row else 5e-4
            a_H2S = (S_tot_uM * 1e-6) / (1.0 + Ka1_S / aH)
            a_HS = a_H2S * Ka1_S / aH
            a_NH4 = (TAN_uM * 1e-6) / (1.0 + Ka_N / aH)
            act = {
                "H2O": 1.0, "H+": aH,
                "O2": O2_uM * 1e-6, "NO3-": NO3_uM * 1e-6, "NO2-": NO2_uM * 1e-6,
                "N2": N2_M, "NH4+": a_NH4, "SO4-2": SO4_uM * 1e-6, "HS-": a_HS,
                "H2": H2_uM * 1e-6, "CH4": CH4_uM * 1e-6, "HCO3-": a_HCO3,
                "CH3COO-": 0.05 * DOC_uM * 1e-6, "C6H12O6": (0.10 / 6.0) * DOC_uM * 1e-6,
            }
            conc_uM = {
                "CH3COO-": 0.05 * DOC_uM, "C6H12O6": (0.10 / 6.0) * DOC_uM,
                "O2": O2_uM, "NO3-": NO3_uM, "NO2-": NO2_uM, "TAN": TAN_uM,
                "SO4-2": SO4_uM, "SigmaH2S": S_tot_uM, "H2": H2_uM, "CH4": CH4_uM,
            }

            m_reads = dna_mag_reads[sid]
            mag_depth = {m: m_reads.get(m, 0.0) / mag_bp[m] for m in prok_mags}
            hk_depths = []
            for _, items in hk_by_ko.items():
                if items:
                    num = sum(mag_depth[m] * L for _, m, L in items)
                    den = sum(L for _, m, L in items)
                    hk_depths.append(num / den)
            hk_depths.sort()
            d_prok = hk_depths[len(hk_depths) // 2]
            assert d_prok > 0.0

            caps = [0.0] * R_len
            for cid, pid in cds_to_pathway.items():
                m, _, L = cds_records[cid]
                caps[pid_to_idx[pid]] += (mag_depth[m] * L) / d_prok

            w = {k: [] for k in MODELS}
            for i, r in enumerate(panel_rows):
                fk = 1.0
                for s_mon, ks_val in bv.parse_monod_spec(r["monod_substrates_uM"]):
                    c_val = conc_uM[s_mon]
                    fk *= c_val / (ks_val + c_val)
                if r["o2_inhibited"] == "true":
                    ki = float(r["Ki_O2_uM"])
                    fk *= ki / (ki + O2_uM)
                st = bv.parse_stoich_spec(r["stoichiometry_spec"])
                ne = float(r["n_e"])
                lnQ = sum(nu * math.log(act[s_q]) for s_q, nu in st.items())
                dGr_per_e = (dG0_T(st, T) + RT * lnQ) / ne
                arg = (dGr_per_e + float(r["dG_min_kJ_per_e"])) / (float(r["chi_per_e"]) * RT)
                ft = 0.0 if arg >= 0.0 else (1.0 - math.exp(arg))
                yield_e = max(0.0, -dGr_per_e)
                c_i = caps[i]
                w["M_cap"].append(c_i)
                w["M_capK"].append(c_i * fk)
                w["M_capKT"].append(c_i * fk * ft)
                w["M_power"].append(c_i * fk * ft * yield_e)

            counts = rna_panel_counts.get(sid, [0] * R_len)
            events.append({
                "sample_id": sid,
                "cast_id": row["event_id"],
                "O2_uM": O2_uM,
                "N_panel": int(sum(counts)),
                "counts": counts,
                "probs": {k: bv.safe_probs(w[k], caps, EPS_OBS) for k in MODELS},
            })
    return events


# ---------------------------------------------------------------- cohort arrays
class Cohort:
    """Eligible events as arrays: n (E,R), P (M,E,R), block index (E,), oxic (E,)."""

    def __init__(self, events, block_of_event):
        keep = [i for i, e in enumerate(events) if e["N_panel"] >= N_PANEL_MIN]
        self.n_events_total = len(events)
        ev = [events[i] for i in keep]
        self.n = np.array([e["counts"] for e in ev], dtype=float)
        self.P = np.array([[e["probs"][m] for e in ev] for m in MODELS], dtype=float)
        labels = [block_of_event[i] for i in keep]
        self.block_labels = sorted(set(labels))
        idx = {b: j for j, b in enumerate(self.block_labels)}
        self.block = np.array([idx[b] for b in labels], dtype=int)
        self.O2 = np.array([e["O2_uM"] for e in ev])
        self.sample_ids = [e["sample_id"] for e in ev]
        self.C = len(self.block_labels)

    def gate(self):
        oxic_blocks = len({b for b, o in zip(self.block, self.O2) if o >= OXIC_UM})
        anox_blocks = len({b for b, o in zip(self.block, self.O2) if o < OXIC_UM})
        E = len(self.block)
        det = (self.n > 0).sum(axis=0)
        n_cov = int((det >= 0.5 * E).sum())
        ok = self.C >= 6 and oxic_blocks >= 3 and anox_blocks >= 3 and E >= 20 and n_cov >= 4
        return {"events_total": self.n_events_total, "events_eligible": E, "blocks": self.C,
                "oxic_blocks": oxic_blocks, "anoxic_blocks": anox_blocks,
                "rna_cov50_pathways": n_cov, "gate_met": bool(ok)}


def dm_loglik_vec(n, P, phi):
    """Exact DM log-likelihood per event; n (E,R), P (...,E,R)."""
    a0 = 1.0 / phi
    N = n.sum(axis=-1)
    A = a0 * P
    return (gammaln(N + 1.0) + gammaln(a0) - gammaln(N + a0)
            + (gammaln(n + A) - gammaln(n + 1.0) - gammaln(A)).sum(axis=-1))


def fit_phi(n, P_capK):
    f = lambda lp: -dm_loglik_vec(n, P_capK, math.exp(lp)).sum()
    res = minimize_scalar(f, bounds=LOGPHI_BOUNDS, method="bounded", options={"xatol": 1e-6})
    return math.exp(res.x)


def heldout_scores(n, P, block, C):
    """Leave-one-block-out: returns (M,E) held-out log densities and per-fold phi."""
    S = np.empty((P.shape[0], n.shape[0]))
    phis = np.empty(C)
    iK = MODELS.index("M_capK")
    for c in range(C):
        te = block == c
        tr = ~te
        phi = fit_phi(n[tr], P[iK, tr])
        phis[c] = phi
        S[:, te] = dm_loglik_vec(n[te], P[:, te], phi)
    return S, phis


def block_means(x, block, C):
    return np.bincount(block, weights=x, minlength=C) / np.bincount(block, minlength=C)


def boot_t_ci(d, B, rng, alpha=0.05):
    """Studentised cluster bootstrap CI for the mean of per-block values d."""
    C = len(d)
    m = d.mean()
    se = d.std(ddof=1) / math.sqrt(C)
    idx = rng.integers(0, C, size=(B, C))
    db = d[idx]
    mb = db.mean(axis=1)
    sb = db.std(axis=1, ddof=1) / math.sqrt(C)
    ok = sb > 0
    t = (mb[ok] - m) / sb[ok]
    lo_t, hi_t = np.quantile(t, [1 - alpha / 2, alpha / 2])
    return float(m), float(m - lo_t * se), float(m - hi_t * se), float(se)


def verdict(lo, hi, lo_weight):
    if lo > DELTA:
        return "supported" + ("+power_allocation" if lo_weight > 0 else "(incremental_energy_weighting_unresolved)")
    if lo >= -DELTA and hi <= DELTA:
        return "equivalent"
    if hi < -DELTA:
        return "against"
    return "inconclusive"


def analyse(n, P, block, C, B, rng):
    S, phis = heldout_scores(n, P, block, C)
    bm = {m: block_means(S[i], block, C) for i, m in enumerate(MODELS)}
    contrasts = {
        "Delta": bm["M_power"] - bm["M_capK"],
        "Delta_gate": bm["M_capKT"] - bm["M_capK"],
        "Delta_weight": bm["M_power"] - bm["M_capKT"],
        "capK_minus_cap": bm["M_capK"] - bm["M_cap"],
    }
    out = {"S": {m: float(bm[m].mean()) for m in MODELS},
           "phi_folds": [float(phis.min()), float(np.median(phis)), float(phis.max())]}
    for k, d in contrasts.items():
        est, lo, hi, se = boot_t_ci(d, B, rng)
        out[k] = {"est": est, "lo95": lo, "hi95": hi, "se": se}
    out["verdict"] = verdict(out["Delta"]["lo95"], out["Delta"]["hi95"], out["Delta_weight"]["lo95"])
    return out


# ---------------------------------------------------------------- drivers
def load_synthetic(rows, sp, datasets=None):
    with open(BLOCK_MAP) as fh:
        bmap = {r["dataset_id"]: r["block"] for r in csv.DictReader(fh, delimiter="\t")}
    ds = sorted(datasets or bmap)
    events, blocks = [], []
    for d in ds:
        evs = extract_run(os.path.join(SYN, d), rows, sp)
        events += evs
        blocks += [bmap[d]] * len(evs)
    return events, blocks


def cmd_regress(rows, sp):
    expect = {"D001": {"M_cap": -64.41, "M_capK": -65.95, "M_capKT": -66.31, "M_power": -75.77},
              "D002": {"M_cap": -82.87, "M_capK": -83.03, "M_capKT": -87.26, "M_power": -109.86}}
    for d, exp in expect.items():
        evs = extract_run(os.path.join(SYN, d), rows, sp)
        co = Cohort(evs, [e["cast_id"] for e in evs])
        S = dm_loglik_vec(co.n, co.P, 0.05)
        got = {m: round(float(block_means(S[i], co.block, co.C).mean()), 2) for i, m in enumerate(MODELS)}
        assert got == exp, (d, got, exp)
        print(f"REGRESS {d} OK {got}")


def cmd_benchmark(rows, sp):
    t0 = time.time()
    events, blocks = load_synthetic(rows, sp)
    co = Cohort(events, blocks)
    g = co.gate()
    rng = np.random.default_rng(SEED)
    res = analyse(co.n, co.P, co.block, co.C, 9999, rng)
    out = {"protocol_sha12": sha12(os.path.join(ROOT, "protocols", "P6-thermodynamic-allocation.md")),
           "table_sha12": TABLE_SHA12, "block_map_sha12": sha12(BLOCK_MAP),
           "scorer_sha12": sha12(os.path.abspath(__file__)), "delta": DELTA, "seed": SEED,
           "gate": g, "result": res, "runtime_s": round(time.time() - t0, 1)}
    if not g["gate_met"]:
        out["result"]["verdict"] = "inconclusive (gate not met)"
    os.makedirs(OUT_DIR, exist_ok=True)
    with open(os.path.join(OUT_DIR, "synthetic_benchmark.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    # per-block contrasts for audit (blocks are anonymous labels)
    S, _ = heldout_scores(co.n, co.P, co.block, co.C)
    with open(os.path.join(OUT_DIR, "synthetic_benchmark_blocks.tsv"), "w") as fh:
        fh.write("block\tn_events\t" + "\t".join(MODELS) + "\n")
        cnt = np.bincount(co.block, minlength=co.C)
        bms = [block_means(S[i], co.block, co.C) for i in range(len(MODELS))]
        for c, lab in enumerate(co.block_labels):
            fh.write(f"{lab}\t{cnt[c]}\t" + "\t".join(f"{bm[c]:.4f}" for bm in bms) + "\n")
    print(json.dumps(out, indent=1))


def cmd_calibrate(rows, sp, reps):
    t0 = time.time()
    events, blocks = load_synthetic(rows, sp)
    co = Cohort(events, blocks)
    Nn = co.n.sum(axis=1).astype(int)  # panel depths only; RNA composition not used
    iK, iP = MODELS.index("M_capK"), MODELS.index("M_power")
    gens = {"M_capK": co.P[iK], "M_power": co.P[iP],
            "placebo": 0.98 * co.P[iK] + 0.02 * co.P[iP]}
    rng = np.random.default_rng(SEED + 1)
    table = []
    for phi in (0.01, 0.05, 0.2):
        for gname, Pi in gens.items():
            verdicts, deltas = {}, []
            for _ in range(reps):
                p = rng.gamma(Pi / phi)
                p /= p.sum(axis=1, keepdims=True)
                n = np.array([rng.multinomial(N, pe) for N, pe in zip(Nn, p)], dtype=float)
                # cohort-level coverage gate rechecked on each generated cohort (Astra 249690027885772800)
                det = (n > 0).sum(axis=0)
                if int((det >= 0.5 * n.shape[0]).sum()) < 4:
                    verdicts["inconclusive_gate"] = verdicts.get("inconclusive_gate", 0) + 1
                    continue
                r = analyse(n, co.P, co.block, co.C, 1999, rng)
                verdicts[r["verdict"]] = verdicts.get(r["verdict"], 0) + 1
                deltas.append(r["Delta"]["est"])
            # Monte-Carlo "true" Delta: expected held-out score difference under the generator,
            # at the generating phi, equal-block weighted, from 20 large draws.
            tv = []
            for _ in range(20):
                p = rng.gamma(Pi / phi)
                p /= p.sum(axis=1, keepdims=True)
                n = np.array([rng.multinomial(N, pe) for N, pe in zip(Nn, p)], dtype=float)
                s = dm_loglik_vec(n, co.P[[iK, iP]], phi)
                tv.append(block_means(s[1] - s[0], co.block, co.C).mean())
            # Procedure-level reference: under leave-one-block-out the held-out block's counts are
            # generated independently of the training blocks, so E[Delta_hat] over replicates is the
            # procedure-level target (training-only phi, independent held-out counts). Reported with
            # its Monte-Carlo SE; the oracle-phi contrast is kept separately labelled.
            dh = np.array(deltas)
            row = {"phi_gen": phi, "generator": gname, "reps": reps, "n_scored": int(len(dh)),
                   "Delta_ref_procedure": float(dh.mean()) if len(dh) else None,
                   "Delta_ref_procedure_mcse": float(dh.std(ddof=1) / np.sqrt(len(dh))) if len(dh) > 1 else None,
                   "Delta_oracle_phi_gen": float(np.mean(tv)),
                   "verdicts": verdicts}
            table.append(row)
            print(json.dumps(row), flush=True)
    out = {"protocol_sha12": sha12(os.path.join(ROOT, "protocols", "P6-thermodynamic-allocation.md")),
           "scorer_sha12": sha12(os.path.abspath(__file__)), "block_map_sha12": sha12(BLOCK_MAP),
           "delta": DELTA, "seed": SEED + 1, "cells": table, "runtime_s": round(time.time() - t0, 1)}
    with open(os.path.join(OUT_DIR, "calibration.json"), "w") as fh:
        json.dump(out, fh, indent=1)


def main():
    rows = load_table()
    sp = bv.load_obigt_species()
    cmd = sys.argv[1] if len(sys.argv) > 1 else "regress"
    if cmd == "regress":
        cmd_regress(rows, sp)
    elif cmd == "benchmark":
        cmd_benchmark(rows, sp)
    elif cmd == "calibrate":
        cmd_calibrate(rows, sp, int(sys.argv[2]) if len(sys.argv) > 2 else 200)
    else:
        sys.exit(__doc__)


if __name__ == "__main__":
    main()
