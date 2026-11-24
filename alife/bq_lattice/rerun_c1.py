"""v3 targeted rerun after Carbon 248999088913563648 / Astra 248999382242213888.

(a) C1 wiring: rerun 8 sweep_rev configurations (first seed each, identical seeds, unchanged world
    for these, so trajectories reproduce) and compare the v2 needs-based wiring with the v3
    deficit-based wiring on the same counterfactual runs: edge counts and verdict changes.
(b) Planted positive control for I: base, iota 0, stack_discount = 1 (2 seeds) vs stack_discount 0.
usage: python3 -I rerun_c1.py OUTROOT [NWORKERS]
"""
import sys, os, json, csv, time
import multiprocessing as mp
from collections import defaultdict
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import sweep as S

BASE = S.BASE


def jobs(root):
    out = []
    picks = [("A", dict(BASE, iota=io, **kv)) for io in (0.0, 1.0)
             for kv in ({}, {"cost": 0.08}, {"leak": 0.8})]
    picks += [("A", dict(BASE, iota=1.0, D=0.05)), ("E", dict(BASE, iota=1.0, euk_frac=0.1))]
    for b, c in picks:
        name = S.cname(b, c)
        out.append((root, name, S.seed_of(name, 0), c))
    for rep in range(2):
        c = dict(BASE, iota=0.0, stack_discount=1.0)
        name = "P_planted_stack1_" + S.cname("A", dict(BASE, iota=0.0))
        out.append((root, name, S.seed_of(name, rep), c))
    return out


def summarise(root):
    rows = []
    edges = defaultdict(int); flips = defaultdict(int)
    for rd in sorted(d for d in (os.path.join(dp) for dp, dn, fn in os.walk(f"{root}/runs")) if
                     os.path.exists(f"{d}/inference.json")):
        nec = list(csv.DictReader(open(f"{rd}/truth/necessity_lineage.tsv"), delimiter="\t"))
        inf = json.load(open(f"{rd}/inference.json"))
        name = rd.split("/runs/")[1]
        for r in nec:
            if r.get("status") not in ("plexi", "not_plexi"):
                continue
            sc = r["scope"]
            edges[(sc, "v2")] += int(r["in_edges_v2_needs"]); edges[(sc, "v3")] += int(r["in_edges_v3_deficit"])
            flips[(sc, r["status_v2_needs"], r["status"])] += 1
            if sc == "all":
                flips[(name.split("/")[0], r["status_v2_needs"], r["status"])] += 1
        sg = inf["sigma"]
        rows.append(dict(run=name, I_obs=sg["incomplete_prok"]["point"]["I"],
                         I_obs_ci=sg["incomplete_prok"]["ci"].get("I"),
                         I_oracle=sg["oracle_genotypes"]["point"]["I"],
                         I_oracle_ci=sg["oracle_genotypes"]["ci"].get("I"),
                         Sigma_leaky=sg["incomplete_prok"]["point"]["leaky"],
                         Sigma_private=sg["incomplete_prok"]["point"]["private"],
                         prev_leaky=sg["incomplete_prok"]["prevalence"].get("leaky"),
                         prev_private=sg["incomplete_prok"]["prevalence"].get("private")))
    os.makedirs(f"{root}/results", exist_ok=True)
    with open(f"{root}/results/c1_wiring_changes.tsv", "w") as f:
        f.write("key\tv2_status\tv3_status\tn\n")
        for k, v in sorted(flips.items(), key=str):
            f.write(f"{k[0]}\t{k[1]}\t{k[2]}\t{v}\n")
        for k, v in sorted(edges.items()):
            f.write(f"edges_{k[0]}\t{k[1]}\t\t{v}\n")
    with open(f"{root}/results/i_statistic.tsv", "w") as f:
        keys = list(rows[0].keys())
        f.write("\t".join(keys) + "\n")
        for r in rows:
            f.write("\t".join(str(r[k]) for k in keys) + "\n")
    print(open(f"{root}/results/c1_wiring_changes.tsv").read())
    print(open(f"{root}/results/i_statistic.tsv").read())


if __name__ == "__main__":
    root = sys.argv[1]
    nw = int(sys.argv[2]) if len(sys.argv) > 2 else 8
    if nw == 0:
        summarise(root); sys.exit()
    J = jobs(root)
    os.makedirs(root, exist_ok=True)
    with open(f"{root}/seeds.tsv", "w") as f:
        f.write("config\tseed\tparams_json\n")
        for _, n, sd, c in J:
            f.write(f"{n}\t{sd}\t{json.dumps(c)}\n")
    t0 = time.time()
    with mp.Pool(nw) as pool:
        for i, (rd, st) in enumerate(pool.imap_unordered(S.job, J)):
            print(f"[{time.time() - t0:5.0f}s] {i + 1}/{len(J)} {rd} {st}", flush=True)
    summarise(root)
    print(f"done in {time.time() - t0:.0f}s")
