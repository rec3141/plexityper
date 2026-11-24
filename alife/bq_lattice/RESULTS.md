# bq_lattice results: revised after cross-family review

**Status.** This file reports the **revised** model, truth and inference, run on a reduced sweep
(`sweep_rev/`).
- The first sweep (`sweep_out/`, 87 runs) used code with the defects listed in §6.
- **Do not cite `sweep_out/`.** Its tables are kept only as a record, and they cannot be
  regenerated with the current `sweep.py`.
- Every claim of the first version that is not restated below is withdrawn.
- Exploratory sweep, not pre-registered. Cross-family reviews so far: Astra (Codex,
  248988907492261888, 248989851437154304) and Carbon (Gemini, 248989754267713536), both on v1.
  v2 was re-reviewed by Carbon (248999088913563648, 249001216432324608) and Astra (248999382242213888, 249001252071325696); v3 changes are in §6. **v3 has not been re-reviewed.**

**Revised sweep.**
- 32 runs: 16 configurations × 2 seeds, 5,000 steps on a 64 × 64 lattice.
- Block A is prokaryotes only, at ι = 0 and 1, for:
  - base (leak 0.5, cost 0.03, D 0.3, φ 0.2);
  - leak 0.2 and 0.8;
  - cost 0.08;
  - D 0.05;
  - φ 0 and 0.4.
- Block E is base with the eukaryote, at ι = 0 and 1.
- 0 failures. Wall time 1,454 s on 8 workers. 290–400 s per run, including 197–200
  counterfactual runs (T7 in `sweep_rev/results/tables.md`).
- Seeds: `sweep_rev/seeds.tsv`.
- **Two seeds per cell and only 2 necessity time points per run.** Truth units per cell are few,
  so treat every number below as descriptive.

**How truth is defined now.**
- X is plexitrophic iff:
  - X is viable in situ;
  - at least one tested kept set containing X is **feasible**, meaning every member, not only X,
    is viable over the 100-step horizon (growth ≥ death over steps 50–100);
  - and **no feasible tested set** admits a counterexample wiring: d(X) ≤ 1 with X on **no simple
    directed cycle of ≥ 3**.
- This is P-all semantics with κ, from `discussion/deep/plexitrophy-computation.md` §1.5 I1 and
  I2.
- Tested sets are all subsets of 1–4 genotypes drawn from the 7 most abundant. They are not
  filtered by genome coverage, so indirect chain members are included.
- d and wiring come from **tagged realised flux inside each counterfactual**.
- Absent genotypes are inert ghosts: they keep their sites and take up only the products they
  lack, with no O2 flux, no export, no birth and no death.
- Needs use **active** capability: an OX-carrying genome's ANX counts as absent at ι = 1.
- The verdict is restricted to the tested sets, the stated medium (ambient 0.02 per site-step),
  the stated endpoint and the stated horizon.

## 1. Regimes (T1, cell-weighted, scope all, last two sample times)

| config | ι | autonomous | not plexitrophic (a d ≤ 1, κ = 0 witness exists) | plexitrophic (fan-in) | plexi range over seeds | unresolved (no feasible tested set) | realised d≥2 |
|---|---|---|---|---|---|---|---|
| base | 0 | 0.06 | 0.59 | 0.15 | 0.06–0.24 | 0.19 | 0.00 |
| base | 1 | 0.00 | 0.66 | 0.28 | 0.09–0.48 | 0.06 | 0.06 |
| leak 0.2 | 0 | **0.65** | 0.35 | 0.00 | – | 0.00 | 0.00 |
| leak 0.2 | 1 | 0.00 | 0.71 | 0.00 | – | 0.29 | 0.02 |
| leak 0.8 | 0 | 0.00 | 0.42 | 0.30 | 0.16–0.44 | 0.28 | 0.05 |
| leak 0.8 | 1 | 0.00 | 0.43 | 0.33 | 0.31–0.34 | 0.25 | 0.05 |
| cost 0.08 | 0 | 0.00 | 0.46 | **0.54** | 0.25–0.82 | 0.00 | 0.35 |
| cost 0.08 | 1 | 0.00 | 0.53 | 0.45 | 0.44–0.46 | 0.02 | 0.48 |
| D 0.05 | 0 | **0.60** | 0.40 | 0.00 | – | 0.00 | 0.00 |
| D 0.05 | 1 | 0.00 | 0.92 | 0.02 | 0.00–0.04 | 0.06 | 0.02 |
| φ 0 | 0 | 0.00 | 0.46 | 0.32 | 0.17–0.48 | 0.21 | 0.04 |
| φ 0 | 1 | 0.00 | 0.28 | 0.25 | 0.19–0.30 | 0.48 | 0.27 |
| φ 0.4 | 0 | 0.35 | 0.65 | 0.00 | – | 0.00 | 0.00 |
| φ 0.4 | 1 | 0.00 | 0.62 | 0.10 | 0.08–0.12 | 0.28 | 0.02 |
| E base | 0 | 0.03 | 0.59 | 0.22 (+0.04 cyclic) | 0.07–0.45 | 0.11 | 0.00 |
| E base | 1 | 0.00 | 0.69 | 0.31 | 0.00–0.63 | 0.00 | 0.43 |

What the simulation shows, under its assumptions:
- **Fan-in plexitrophy emerges with no coded dependencies.** It reaches 0.15–0.54 of cells where
  genes are leaky and costly. It is 0 where leakage, cost or diffusion is low (leak 0.2, D 0.05)
  and in high-particle worlds (φ 0.4) at ι = 0. Those worlds are 35–65% autonomous.
- **Single-supplier dependence (d ≤ 1 with no cycle witness) is the commonest state:** 0.28–0.92
  of cells.
- **Cyclic plexitrophy is essentially absent** under the simple-cycle criterion: 0 in block A and
  0.04 in one E cell.
- **What the O2 block does.**
  - It removes autonomy in every configuration that had any (D 0.05: −0.60; leak 0.2: −0.65;
    φ 0.4: −0.35).
  - That autonomy turns mostly into **single-supplier dependence** or unresolved states, not
    fan-in. Plexitrophy changes by +0.02, 0.00 and +0.10 in those three.
  - Across the 7 matched pairs, plexitrophy rises in 4 (base +0.13, φ 0.4 +0.10, leak 0.8 +0.03,
    D 0.05 +0.02), is unchanged in 1 (leak 0.2) and falls in 2 (cost 0.08 and φ 0, each −0.08).
    These are 2-seed means with wide ranges.
  - **The model therefore does not show that blocked stacking reliably increases plexitrophy.**
- **Genome capacity vs survival** (isolated-recipient counterfactual: same initial field,
  medium, horizon and endpoint).
  - Where any genome-autonomous genotypes existed (ι = 0), 100% of them survived alone.
  - Of the genome-dependent genotypes, 0–6% survived alone on the ambient medium.
  - At ι = 1, no genotype was autonomous by *active* capability, including OX + ANX carriers.
  - So in this world, genome capacity predicted isolated survival almost perfectly, and ambient
    rescue was rare (6% at D 0.05, ι = 0).
- **Use is not need, in both directions.**
  - Realised d ≥ 2 (two sources each ≥ 10% of partner flux, in situ) was 0.00–0.06 where
    plexitrophy was 0.15–0.33.
  - It was 0.35–0.48 where plexitrophy was 0.31–0.54 (cost 0.08, and E at ι = 1).
- **Unresolved fractions reach 0.29–0.48** in some ι = 1 cells. Those recipients had no feasible
  tested set of ≤ 4 genotypes, so they need more than 3 partners or partners outside the top 7.

## 2. Inference on the synthetic observations

### 2a. Genome-complementarity d (P2-like)
Setup:
- HQ recipients, prokaryotes only, 80/25% rules.
- The permissive single-supplier check is now applied before accepting d ≥ 2. One ambiguous MAG,
  an unknown-detectability eukaryote route, or the unbinned pool covering J makes the call
  unresolved.

Pooled over all 32 runs (MAG × time units with resolved truth):

| observer | TP | FP | FN | TN | unresolved | precision (resolved) | recall (resolved) | truth-plexi not evaluated |
|---|---|---|---|---|---|---|---|---|
| P2-like (incomplete, HQ, prok-only) | 27 | 9 | 22 | 113 | 8 | 0.75 | 0.55 | 46 |
| complete genomes, prok-only | 46 | 1 | 49 | 278 | 0 | 0.98 | 0.48 | 1 |
| incomplete, with fragmentary eukaryote | 24 | 9 | 25 | 113 | 8 | 0.73 | 0.49 | 46 |

Bounds:
- Unresolved calls counted as not plexitrophic: recall 27/50.
- Unresolved calls counted as plexitrophic: precision 28/44 = 0.64.

Readings:
- **Complete-genome complementarity is precise but finds about half of the plexitrophs.**
  Precision is 0.98 and recall 0.48. A genome that encodes every missing function often is not a
  feasible single supplier under removal, because of distance, quantity or its own unmet needs.
- **Incompleteness trades precision for apparent recall.** False absences (module false-absence
  rate 2.1–3.4% in HQ MAGs, 9.4–12% overall) create spurious fan-in: FP rose from 1 to 9.
  d_est rose by +0.00 to +0.11 on the same MAG × sample.
- **The HQ filter drops 46 of 95 true-plexitrophic units** from evaluation.
- Per-cell values are in T2. Most cells have fewer than 10 evaluable units.

### 2b. Fan-in null (Prop. 1) with matched denominators
Two comparisons were made, on definite evidence:
- **unconditional:** all eligible recipient-samples with J ≠ 0, where no-donor counts as not
  fan-in on both sides;
- **conditional on coverage:** observed vs P(d ≥ 2)/P(covered).

Results:
- O − E was negative in all 16 cells under both versions:
  - unconditional: −0.01 to −0.21;
  - conditional: −0.01 to −0.25.
- It did not track true plexitrophy at ι = 0 (Spearman −0.43 conditional, −0.36 unconditional).
- **The v1 conclusion stands with matched denominators:** observed minus expected against the
  independent closed form is not an indicator of plexitrophy here.

### 2c. Leaky-vs-private co-carriage I (revised estimator, prevalence-matched)
Changes:
- Private functions now have no background supply (`priv_bg` = 0) and patch coverage 0.95.
- Private prevalence was 0.30–0.64 vs leaky 0.36–0.64.
- Σ is per-donor co-carriage, log(observed/expected), against the fixed-count column-permutation
  expectation n_j n_k / N.

Results:
- **I was small everywhere:** −0.14 to +0.14 on observations, −0.17 to +0.11 on true genotypes.
- Its 95% intervals excluded 0 in both directions in different cells.
- It was unrelated to true plexitrophy (Spearman 0.00 at ι = 0, −0.28 overall).
- The v1 explanation (prevalence mismatch and saturation) is confirmed as an artefact.
- With both fixed, the near-zero pooled leaky-minus-private contrast is **consistent with the
  null in this sweep**. Within-donor co-carriage of leaky functions stayed close to the
  permutation expectation (Σ_leaky 0.00–0.12).
- The world has no phylogenetic linkage and no epistatic co-selection. Even so, independent
  mutations can still give correlated carriage after selection, so near-independence is an
  outcome here, not a design guarantee.
- **Specificity is not estimated.** That would need a pre-defined positive-call threshold and the
  fraction of known-null replicates that cross it.
- **Sensitivity is untested.** The planted positive control (§6, v3) did not provide a usable
  test.

### 2d. Rare donors and abundance (T6)
- Every member of feasible sets for plexitrophic recipients was in the catalogue, detected at
  that time point and in the abundant core (100% where defined).
- 0–1% of d ≤ 1 witness members lay outside the core.
- An imposed 50% background (lognormal or power-law, identical depth) did not change detection.
- Genotype richness was 52–122. The lognormal vs power-law Vuong z averaged −0.1 to +1.2, with
  z > 2 (favouring lognormal) in ≤ 14% of time points. **Tail shape is not identifiable at this richness.**
- The rare-donor result therefore does not carry over to 20k-ASV communities.

## 3. Signatures of blocked stacking

Seven matched block-A pairs (ι = 1 minus ι = 0, 2-seed means):

| observable | direction | Δ range | Spearman with true plexitrophy, ι = 0 only |
|---|---|---|---|
| Σ_inc, OX × ANX per-donor co-carriage (log O/E) | lower in 6/7 | −0.56 to +0.01 | +0.42 |
| fraction of cells co-carrying OX and ANX | lower in 6/7 | −0.38 to +0.03 | −0.72 |
| observed P(d_est ≥ 2 \| covered) | higher in 4/7 | −0.02 to +0.16 | **+0.84** |
| true plexitrophy | higher in 4/7 | −0.08 to +0.13 | (truth) |

Readings:
- The co-carriage deficit of the incompatible pair tracks blocking (AUROC 0.86 for ι = 1 lower).
  But it also correlates with plexitrophy when there is no block (ρ = +0.42 at ι = 0). The v1
  statement that it was "nearly unrelated to fan-in" (ρ = 0.24) **does not hold in v2**.
- Observed fan-in follows plexitrophy from any cause (ρ = 0.84).
- **No observable here uniquely indicates blocked stacking.**
- Blocking itself did not consistently raise true plexitrophy, so there is no model-backed
  "blocked stacking ⇒ plexitrophy" prediction to test from this world.
- What the world does predict: blocking removes autonomy and replaces it with *single-supplier*
  dependence. Testing that would need the necessity or autonomy readout, not a co-carriage
  statistic.

## 4. Particle enrichment (P4)
- **Particles are now anoxic.** Light attenuation (OX at 0.3 on particles) and 3× particulate
  respiration were added. In the φ 0.2 runs, mean O2 was 0.065 in particle samples vs 2.11 in
  water samples, and 80% of particle samples had O2 < 0.1.
- **Units:** the same recipient MAG × event, resolved in both fractions, weight 1. That fixes
  the recipient, the endpoint set J and eligibility. There are 288–464 units per cell.
- **Observer P − W:** −0.03 to +0.04, mostly within ±0.02. **No particle enrichment was
  inferred.**
- **Donor-count standardisation** subsampled the larger fraction's definite donor set to the
  smaller size. It lowered P − W in 14 of 14 cells, by 0.001–0.06 (result −0.08…+0.01).
  Particle samples had fewer donors: 10–18 vs 16–25.
  - This is consistent with richness raising d_est in the poorer fraction. It is a separate
    test, and it does not by itself attribute any sign.
- **Truth on the same units** (habitat-specific P-all verdicts, 2 necessity time points): 2–28
  units per cell, and P − W was 0 in 10 of 14 cells. The other four were −0.33, −0.21, −0.06 and
  +0.09. Particle and water verdicts for a genotype rarely differed.
- **Particle enrichment did not emerge in truth or in inference.** The v1 attribution of the
  observer's sign to donor count is **withdrawn**. v2 shows only that equalising donor count
  shifts the observer by −0.001 to −0.06.

## 5. Eukaryote (T5, block E, 2 seeds each)
Contribution:
- **ι = 0:** 5% of cells, 6% of true uptake, present in 49% of feasible sets.
- **ι = 1:** 15% of cells, 28% of uptake, present in 49% of feasible sets.
- The eukaryote now takes up only B12, so the v1 shares (23–49% at ι = 1) were partly scavenging.
  They are superseded.

Prokaryote-only inference:
- At ι = 1 it reached recall 0.43 on few units. **The v1 "100% unresolved" result is withdrawn.**
  It came from the v1 eukaryote scavenging and the OX block.
- Adding the fragmentary eukaryote bin (35% gene recovery, modules never called absent) did not
  help (recall 0.00 on those units).
- With complete eukaryote genes, precision was 1.00 and recall 0.55.
- **With 2 seeds and few units, these numbers are anecdotal.**

## 6. Revisions after cross-family review (Astra 248988907492261888)
These revisions also cover Astra 248989851437154304 and Carbon 248989754267713536.

| # | Finding | Fix (code) | Effect on claims |
|---|---|---|---|
| A1 | Ghosts produced and respired O2. The eukaryote took up every deficient product. | Ghosts are redox-inert: no O2 production, no respiration, no export, birth or death. They still take up the products they lack. Eukaryote uptake is restricted to B12 (`world.py`). | All truth regenerated. v1 eukaryote shares withdrawn. |
| A2 | Necessity tested only X's viability, filtered candidates by genome coverage, took member-minimal sets, and used SCC on possible transfers. | P-all over tested kept sets. Candidates are not coverage-filtered. Feasibility requires every member to be viable. Wiring comes from tagged realised flux within each counterfactual. κ is a simple cycle through X (`run.py` `necessity_truth`, `bad_wiring`). | "Necessary cycles rise with ι in 10/10" **withdrawn**: cyclic plexitrophy is about 0 under κ. "Blocking adds fan-in in 8/10" replaced by 4 up, 1 unchanged, 2 down of 7. |
| A3 | d ≥ 2 was accepted before the permissive single-supplier check. The null used a different denominator. | Single ambiguous MAG / eukaryote / unbinned-pool cover is checked first and makes the call unresolved. The null is compared unconditionally and conditionally on matched units (`infer.py`). | The null conclusion is unchanged. Precision and recall are now reported with [low, high] bounds. |
| A4 | P4 truth was cell-weighted and included cycles. The observer averaged different recipient sets. | Same recipient × event, equal weights, both fractions. Truth on the same units, any and fan-in-only. Separate donor-standardised observer. | v1 attribution of the sign to donor count **withdrawn**. No enrichment in truth or inference. |
| C1 | Needs used encoded bits, so silent OX + ANX genomes looked autonomous and acted as donors. | Active-capability needs. Every recipient gets the isolated test. Supplier identity comes from realised flux (`run.py` `needs`). | Genome capacity and survival are reported separately (§1). |
| C2 | The eukaryote at 35% recovery was called absent 43% of the time, and ambiguous donors were dropped. | Eukaryote modules are present at ≥ 80%, otherwise unresolved, never absent. Bounds are reported. | §5 |
| C3 | Particles were O2 hotspots. | Light attenuation (`light_p` = 0.3) and particulate respiration (`resp_p_mult` = 3). | Particles are anoxic (§4). |
| C4 | I failed through private prevalence (`priv_bg` = 0.85) and saturation of the any-co-carrier statistic. | `priv_bg` = 0 and `priv_cov` = 0.95 (prevalence matched). Per-donor co-carriage O/E against the permutation expectation. | I ≈ 0. The v1 "I flags unfinished loss" reading is withdrawn. I is not shown to detect anything (§2c). |
| C5 | CDS schema mismatch with `scratch/p2_scorer.py`. | `per_mag_cds.tsv` columns are `mag_id, contig_id, protein_id ({m}_c{ci}_{fi+1:05d}), feature_index, KO`, plus `contig_lengths.tsv`. | Format only. |

Claims from v1 that survive in v2 (re-measured):
- Fan-in emerges without coded dependencies.
- Single-supplier dependence is commonest.
- Use ≠ need.
- Genome complementarity is precise but low-recall, and incompleteness inflates d.
- The independent fan-in null is not indicative.
- Rare donors don't matter at this richness.
- Tail shape is not identifiable.

sha256, first 12 characters, v2:

| File | sha256 |
|---|---|
| world.py | 8ab7115300db |
| run.py | 075a53704772 |
| infer.py | b1b2ad8aa6c0 |
| sweep.py | eaa90a8d05de |
| make_tables.py | 080d36137b4c |
| sweep_rev/results/runs.tsv | 2f300072c3ea |
| sweep_rev/results/tables.md | 33cb3edaf56c |
| sweep_rev/results/signatures.tsv | 2ede2c442124 |
| sweep_rev/seeds.tsv | c0d828246c1a |

v1 hashes, for the record: world 7dabe801bea4, run 9176b9badaf9, infer 6f62909b0aba,
sweep_out/results/runs.tsv 25fe317736b4.


### v3 addendum (Carbon 248999088913563648 / 249001216432324608, Astra 248999382242213888 / 249001252071325696)

**C1 wiring (reopened, now fixed).**
- **The defect.** `supply_sets()` and `realised_truth()` considered only products absent by
  active genome capability. A carrier left short by realised conditions was skipped: ANX
  inhibited by local O2, OX dimmed on particles, or any gene at leak 0.8 (Y(1−leak) < θ). Its
  partner supply was therefore never examined.
- **Carbon's exposure counts** are rows *exposed* to the omission, not changed rows:
  - 1,035 of 3,966 v2 support rows for case (1);
  - 106 rows for case (2);
  - 179 rows for case (3).
- **v3 rule.** An encoded product also enters the wiring when the carrier's realised **non-self**
  uptake of it, in that scope, is ≥ 0.05 per cell-step (demand θ = 1). The partner ≥ ENV rule and
  the 10% supplier-share rule still apply. The same rule is applied to in-situ realised truth,
  with per-scope cell-steps now tracked in the tagging window. Self-tag uptake is excluded from
  supplier shares.
  - The 0.05 floor is a choice, and it is a **scope-averaged detection rule**. A small or
    intermittent contribution below it can still tip growth across the viability threshold.
    Sensitivity below 0.05 is unmeasured. Carbon suggested no floor.
- **Measured change** (Astra's correction applied: old and new rules compared **on the same
  rows**).
  - 10 targeted reruns (`rerun_c1.py`, `sweep_c1/`): 8 sweep_rev configurations, first seed each,
    plus 2 planted-control runs.
  - Both rules were evaluated on the same counterfactual runs. The table counts **paired
    nonautonomous resolved** recipient × time rows. A further 19 / 19 / 20 autonomous rows
    (all / water / particle) are excluded:

| scope | rows | v2 not → v3 plexi | v2 plexi → v3 not | in-edges into X, summed over feasible sets (v2 → v3) |
|---|---|---|---|---|
| all | 105 | 2 | 0 | 3,095 → 3,131 (+1.2%) |
| water | 106 | 3 | 0 | 3,101 → 3,135 (+1.1%) |
| particle | 90 | 1 | 0 | 1,939 → 2,003 (+3.3%) |

- The two scope-all changes were at leak 0.8, ι = 1 and E base, ι = 1.
- The omission therefore moved few verdicts, always toward plexitrophy. It added edges mainly on
  particles (light-limited OX).
- **Status of the sweep_rev tables.** They remain **v2 estimates**, and their v3 correction is
  unknown. The direction and size above were measured only on the targeted v3 trajectories,
  including the planted controls, with v3 self-tag handling in both arms. (An earlier draft of
  this addendum called them "slight underestimates"; that claim is withdrawn, per Astra
  249004336193388544.)
- v3 world runs are not bit-identical to sweep_rev for the same seed. Adding the planted-control
  term changed the float precision of the growth computation, and trajectories diverge
  chaotically. Comparisons above are within run.

**Planted positive control for I** (`stack_discount` = 1, base, ι = 0, 2 seeds).
- **Design.** A genome carrying all five leaky compatible functions pays for them as one
  function.
- **Outcome.**
  - 92–95% of cells carried all five, against 14% in the matched unplanted run. Leaky prevalence
    among detected MAGs was 0.79–0.81.
  - I was −0.015 [−0.026, −0.003] and +0.047 [+0.030, +0.062]. The unplanted runs gave
    −0.13 to +0.22.
- **Interpretation.** The planted effect turned nearly the whole community into stackers. That
  removed the helper vs recipient contrast that I is meant to see, and co-carriage relative to
  the permutation expectation stays near 1 at high prevalence.
- **This control is inconclusive.** It does not show that I lacks sensitivity, and it does not
  show sensitivity either. A usable control needs a planted effect that creates a minority of
  stacked helpers among recipients, at matched leaky and private prevalence. Not done.

**C2–C5:** closed by Carbon (249001216432324608). The E_iota0 CDS export parsed with
`scratch/p2_scorer.py`: 51,084 CDSs, 248 MAGs, 0 unknown.

sha256, first 12 characters, v3:

| File | sha256 |
|---|---|
| world.py | 184f53f37424 |
| run.py | fa4a7f83000c |
| infer.py | b1b2ad8aa6c0 (unchanged) |
| sweep.py | eaa90a8d05de (unchanged) |
| make_tables.py | 080d36137b4c (unchanged) |
| rerun_c1.py | 7e31927e52d2 |
| sweep_c1/results/c1_wiring_changes.tsv | f8b851f4b641 |
| sweep_c1/results/i_statistic.tsv | 1cd98275951f |
| sweep_c1/seeds.tsv | 11d71bd552c1 |

The v2 hashes above remain those of the code that produced `sweep_rev/`.

## 7. Limits
- Tested sets are ≤ 4 genotypes from the top 7.
- Wiring inside a set is realised flux with a 10% share threshold. The other members get one
  greedy minimum cover each, a heuristic restriction of P-all's "every wiring". X's covers are
  enumerated.
- Verdicts depend on the 100-step horizon and the growth ≥ death endpoint.
- Habitat-specific truth is per genotype × habitat, not per event.
- There is no phylogeny, no catabolic exchange (d_anab only) and no regain except by HGT.
- The world holds about 10² genotypes.
- There are 2 seeds per cell.
- **This shows what can happen under these assumptions, not what happens in lakes.**
