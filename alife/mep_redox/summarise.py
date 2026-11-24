"""Aggregate results/*.tsv into results/summary.md tables (markdown).  python3 -I summarise.py"""
import csv
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
RES = os.path.join(HERE, "results")


def rd(name):
    p = os.path.join(RES, name)
    if not os.path.exists(p):
        return []
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


def ms(v, d=1):
    v = np.array([x for x in v if np.isfinite(x)])
    if len(v) == 0:
        return "na"
    return f"{v.mean():.{d}f} ± {v.std(ddof=1):.{d}f}" if len(v) > 1 else f"{v.mean():.{d}f}"


def table(head, rows):
    out = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def main():
    L = []
    sig = rd("sigma.tsv")
    mepr = {r["run_id"]: r for r in rd("mep_yardsticks.tsv")}
    scen = list(dict.fromkeys(r["scenario"] for r in sig))
    rows = []
    for s in scen:
        rr = [r for r in sig if r["scenario"] == s]
        mm = [mepr[r["run_id"]] for r in rr if r["run_id"] in mepr]
        rows.append([s, len(rr), ms([f(r["sigma_chem"]) for r in rr]),
                     ms([f(r["npp"]) for r in rr]), ms([f(r["sigma_chem_per_npp"]) for r in rr], 3),
                     ms([f(r["sigma_alg"]) for r in rr]),
                     ms([f(r["alg_share_sigma"]) for r in rr], 3),
                     ms([f(m["sigma_gauge_fixed"]) for m in mm]),
                     ms([f(m["sigma_greedy"]) for m in mm]),
                     ms([f(m["sigma_random_mean"]) for m in mm]),
                     ms([f(m["sigma_best_found"]) for m in mm]),
                     ms([f(m["emergent_over_best_found"]) for m in mm], 3),
                     ms([f(m["gauge_rel_diff"]) for m in mm], 4),
                     f"{sum(m['fixed_field_sigma_within_0p5pct'] == 'True' for m in mm)}/{len(mm)}",
                     f"{sum(m['powell_success'] == 'True' for m in mm)}/{len(mm)}",
                     ms([f(m["second_pass_rel_gain"]) for m in mm], 4),
                     ms([f(m["params_at_bound"]) for m in mm], 1)])
    L.append("## T1. Entropy production (J K^-1 m^-2 d^-1, mean ± sd over seeds; last 300 d)\n")
    L.append(table(["scenario", "n", "sigma_chemotroph", "NPP (mmol C m^-2 d^-1)",
                    "sigma_chem per NPP (J K^-1 per mmol C)", "sigma_algae", "algal share",
                    "fixed-enzyme copy (gauge)", "greedy (H->0)", "random alloc (mean)",
                    "best-found ref (steady)", "emergent / best-found",
                    "gauge rel. diff", "fixed-field sigma within 0.5%", "Powell converged",
                    "2nd-pass rel. gain", "params at bound"], rows))
    pc = rd("pool_contrast.tsv")
    rows = []
    for s in dict.fromkeys(r["scenario"] for r in pc):
        rr = [r for r in pc if r["scenario"] == s]
        full = np.array([f(r["sigma_full"]) for r in rr])
        sp = np.array([f(r["sigma_split_only"]) for r in rr])
        st = np.array([f(r["sigma_stack_only"]) for r in rr])
        pf = np.array([f(r["perC_full"]) for r in rr])
        psp = np.array([f(r["perC_split_only"]) for r in rr])
        pst = np.array([f(r["perC_stack_only"]) for r in rr])
        rows.append([s, ms(full), ms(sp), ms(st), ms(100 * (sp - st) / st, 2),
                     ms(pf, 3), ms(psp, 3), ms(pst, 3), ms(100 * (psp - pst) / pst, 2)])
    L.append("\n## T2. Packaging contrast: same system, pool restricted to split or stacked genomes\n")
    L.append(table(["scenario", "sigma full pool", "sigma split-only", "sigma stack-only",
                    "split vs stack (%)", "per-NPP full", "per-NPP split-only",
                    "per-NPP stack-only", "per-NPP split vs stack (%)"], rows))
    dol = rd("division_of_labour.tsv")
    paths = list(dict.fromkeys(r["pathway"] for r in dol))
    rows = []
    for s in scen:
        row = [s]
        for pth in paths[:3]:
            v = [f(r["stacked_share"]) for r in dol if r["pathway"] == pth and r["run_id"].startswith(s + "__")]
            row.append(ms(v, 2))
        for pth in paths:
            v = [f(r["depth_sep"]) for r in dol if r["pathway"] == pth and r["run_id"].startswith(s + "__")]
            row.append(ms(v, 1))
        rows.append(row)
    L.append("\n## T3. Division of labour (full pool). stacked share = fraction of consumer-step flux "
             "done by genomes that also carry the producer step; depth_sep = consumer minus producer "
             "flux-weighted depth (m)\n")
    L.append(table(["scenario"] + ["stacked: " + p for p in paths[:3]] +
                   ["sep: " + p for p in paths], rows))
    nk = rd("necessity_knockouts.tsv")
    pr = [r for r in rd("partners_realised.tsv") if r["resource"] == "*"]
    prd = {(r["run_id"], r["recipient"]): r for r in pr}
    rows = []
    for r in nk:
        q = prd.get((r["run_id"], r["recipient"]), {})
        rows.append([r["run_id"].replace("__full", ""), r["recipient"], q.get("d_real", ""),
                     q.get("c_real", ""), r["n_fail_single"], r["fail_direct"],
                     r["fail_indirect"], r["n_solver_fail"], r["d_lb_direct"], r["d_lb_upstream"],
                     r["failed_groups"],
                     f"{f(q.get('alg_share_intake')):.2f}"])
    L.append("\n## T4. Realised partners vs counterfactual necessity (knock-outs, same medium, "
             "endpoint = recipient column biomass >= 10% of an UNPERTURBED CONTROL at the same horizon; "
             "d_lb = min hitting set of failing removal sets (direct: own-input groups + direct donors; "
             "upstream: + indirect), lower bounds conditional on monotone rescue)\n")
    L.append(table(["run", "recipient", "d_real", "c_real", "n single removals failing",
                    "failing: direct donors", "failing: indirect", "solver failures",
                    "d_lb_direct", "d_lb_upstream", "failing OR-groups", "algal share of intake"], rows))
    # summary counts
    rows = []
    for rid in dict.fromkeys(r["run_id"] for r in nk):
        rr = [r for r in nk if r["run_id"] == rid]
        dr = [int(prd[(rid, r["recipient"])]["d_real"]) for r in rr if (rid, r["recipient"]) in prd]
        dq = [int(r["d_lb_direct"]) for r in rr]
        du = [int(r["d_lb_upstream"]) for r in rr]
        rows.append([rid, len(rr), sum(x >= 2 for x in dr), sum(x >= 2 for x in dq),
                     sum(x >= 2 for x in du), sum(1 for a, b in zip(dr, dq) if a >= 2 and b < 2)])
    L.append("\n### T4b. Realised d>=2 vs required d>=2 (recipients with >=0.2% biomass)\n")
    L.append(table(["run", "recipients", "realised d_real>=2", "d_lb_direct>=2", "d_lb_upstream>=2",
                    "realised>=2 but d_lb_direct<2"], rows))
    si = rd("signatures_inference.tsv")
    cfg = list(dict.fromkeys(r["config"] for r in si))
    rows = []
    for c in cfg:
        rr = [r for r in si if r["config"] == c]
        rows.append([c, len(rr), ms([f(r["sigma_chem"]) for r in rr]),
                     ms([f(r["power_rho"]) for r in rr], 2), ms([f(r["monod_rho"]) for r in rr], 2),
                     ms([f(r["dna_rho"]) for r in rr], 2),
                     ms([f(r["unexploited_power"]) for r in rr], 3),
                     ms([f(r["rule_precision"]) for r in rr], 2), ms([f(r["rule_recall"]) for r in rr], 2),
                     ms([f(r["redox_pathlen_inferred"]) for r in rr], 2),
                     ms([f(r["redox_pathlen_true"]) for r in rr], 2),
                     ms([f(r["cooc_precision"]) for r in rr], 2),
                     ms([f(r["cooc_depth_adjacent_false"]) for r in rr], 1)])
    L.append("\n## T5. Observer-side signatures and inference on synthetic data (by configuration)\n")
    L.append(table(["config", "n", "sigma_chem", "rho RNA~power", "rho RNA~Monod", "rho RNA~DNA",
                    "unexploited power", "KO-rule precision", "KO-rule recall",
                    "redox path (inferred)", "redox path (true)", "co-occ precision",
                    "co-occ depth-adjacent false edges"], rows))
    hz = rd("horizon.tsv")
    pols = list(dict.fromkeys(r["policy"] for r in hz))
    Hs = sorted({f(r["eval_horizon_d"]) for r in hz})
    rows = []
    for p_ in pols:
        row = [p_]
        for H in Hs:
            v = [f(r["mean_sigma"]) for r in hz if r["policy"] == p_ and f(r["eval_horizon_d"]) == H]
            row.append(f"{v[0]:.1f}" if v else "")
        rows.append(row)
    L.append("\n## T6. Horizon test (neutral, seed 0; sterile start). Mean chemotrophic sigma over "
             "[0, H]\n")
    L.append(table(["policy"] + [f"H={int(H)} d" for H in Hs], rows))
    bio = rd("biomass.tsv")
    types = list(dict.fromkeys(r["type"] for r in bio))
    rows = []
    for s in scen:
        row = [s]
        for t in types:
            v = [f(r["column_biomass"]) for r in bio if r["type"] == t and r["run_id"].startswith(s + "__")]
            row.append(f"{np.mean(v):.0f}")
        rows.append(row)
    L.append("\n## T7. Column biomass by type (mmol C m^-2, mean over seeds)\n")
    L.append(table(["scenario"] + types, rows))
    fe = rd("free_energy_budget.tsv")
    rows = []
    for c in dict.fromkeys(r["config"] for r in fe):
        rr = [r for r in fe if r["config"] == c]
        rows.append([c, len(rr), ms([f(r["Tsigma_chem"]) / f(r["npp"]) for r in rr], 1),
                     ms([f(r["capture"]) / f(r["npp"]) for r in rr], 1),
                     ms([(f(r["boundary_in"]) - f(r["boundary_out"])) / f(r["npp"]) for r in rr], 1),
                     ms([f(r["reservoir_net"]) / f(r["npp"]) for r in rr], 1),
                     ms([(f(r["Tsigma_diffusion"]) + f(r["Tsigma_sinking"]) + f(r["Tsigma_sediment"])
                          + f(r["Tsigma_hydrolysis"]) + f(r["Tsigma_mortality"])) / f(r["npp"]) for r in rr], 2),
                     ms([f(r["storage"]) / f(r["npp"]) for r in rr], 2),
                     ms([f(r["residual_rel"]) for r in rr], 5)])
    L.append("\n## T8. Species-resolved chemical-potential budget (kJ per mmol C of NPP; "
             "residual relative to T*sigma_chem). T*sigma_chem = capture + (in - out) - reservoir "
             "- storage - other dissipation - residual\n")
    L.append(table(["config", "n", "T*sigma_chem / NPP", "algal capture / NPP",
                    "boundary in - out / NPP", "net to CO2/N2/H2O/H+ reservoirs / NPP",
                    "transport+hydrolysis+mortality+sediment / NPP", "storage / NPP",
                    "residual / T*sigma_chem"], rows))
    open(os.path.join(RES, "summary.md"), "w").write("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
