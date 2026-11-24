"""P6 synthetic benchmark, pre-registered secondary (run only after pooled scores were posted,
relay 249690023779549184): unblind and report Delta, Delta_gate, Delta_weight per scenario family and
emergent vs yardstick, using the same leave-one-block-out held-out event scores.
Cells are summarised as equal-dataset means of per-dataset mean differences, with the
dataset-level SD; no verdict bands are applied (secondary, descriptive)."""
import csv
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import p6_scorer as ps  # noqa: E402

KEY = os.path.join(ps.ROOT, "alife", "mep_redox", "truth", "unblinding_key.tsv")


def main():
    rows = ps.load_table()
    sp = ps.bv.load_obigt_species()
    with open(ps.BLOCK_MAP) as fh:
        bmap = {r["dataset_id"]: r["block"] for r in csv.DictReader(fh, delimiter="\t")}
    with open(KEY) as fh:
        key = {r["dataset_id"]: r["run_id"] for r in csv.DictReader(fh, delimiter="\t")}
    events, blocks, ds_of = [], [], []
    for d in sorted(bmap):
        evs = ps.extract_run(os.path.join(ps.SYN, d), rows, sp)
        events += evs
        blocks += [bmap[d]] * len(evs)
        ds_of += [d] * len(evs)
    co = ps.Cohort(events, blocks)
    keep_ds = [ds_of[i] for i, e in enumerate(events) if e["N_panel"] >= ps.N_PANEL_MIN]
    S, _ = ps.heldout_scores(co.n, co.P, co.block, co.C)
    iC, iK, iKT, iP = (ps.MODELS.index(m) for m in ps.MODELS)
    per_ds = {}
    for j, d in enumerate(keep_ds):
        per_ds.setdefault(d, []).append(j)

    def cls(run_id):
        parts = run_id.split("__")
        scen = parts[0]
        kind = parts[-1] if parts[-1] in ps_yard else "emergent"
        return scen, kind

    out = {}
    for d, idx in per_ds.items():
        scen, kind = cls(key[d])
        idx = np.array(idx)
        rec = {
            "Delta": float((S[iP, idx] - S[iK, idx]).mean()),
            "Delta_gate": float((S[iKT, idx] - S[iK, idx]).mean()),
            "Delta_weight": float((S[iP, idx] - S[iKT, idx]).mean()),
            "capK_minus_cap": float((S[iK, idx] - S[iC, idx]).mean()),
        }
        for k in [(scen, "all"), ("all", kind), (scen, kind)]:
            out.setdefault(k, []).append(rec)
    table = []
    for (scen, kind), recs in sorted(out.items()):
        row = {"scenario": scen, "kind": kind, "n_datasets": len(recs)}
        for m in ["Delta", "Delta_gate", "Delta_weight", "capK_minus_cap"]:
            v = np.array([r[m] for r in recs])
            row[m] = round(float(v.mean()), 3)
            row[m + "_sd"] = round(float(v.std(ddof=1)), 3) if len(v) > 1 else None
        table.append(row)
    path = os.path.join(ps.OUT_DIR, "synthetic_secondary.tsv")
    cols = list(table[0].keys())
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in table:
            fh.write("\t".join("" if r[c] is None else str(r[c]) for c in cols) + "\n")
    print(open(path).read())


ps_yard = ("best_found_ref", "mep_opt", "greedy", "random1", "random2")

if __name__ == "__main__":
    main()
