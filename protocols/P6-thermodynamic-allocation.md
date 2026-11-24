# P6: Does available reaction power predict expressed catabolic allocation? (FROZEN 2026-11-20)

## Status
- Requested by grid on 2026-11-19 in claude-plexi-2's session.
- Drafted by claude-plexi-2 (Opus), from discussion/deep/manifold-timescales.md §9 (being
  revised) with Astra's review applied in advance (relay 249247324022685696, point 2).
- Needs cross-family review (Astra, Carbon), then grid's go.
- No data have been inspected for this draft.
- **Frozen 2026-11-20 by claude-plexi-2** on grid's go to run (relay 249688876280889344).
  Freeze edits are listed under "Freeze record" at the end. No real RNA, and no synthetic RNA
  beyond the D001/D002 smoke tests already reported, has been scored for this protocol.

## Review status (2026-11-19)
- **Reaction / KO / constant table:** protocols/P6-pathway-table.tsv (6c26362681b8), with
  sources in notes/p6-pathway-table-sources.md (0ee6ff9add03). Built by Carbon and checked by
  scratch/build_and_verify_p6_table.py (b5e83a3742d6). **Cross-family review CLOSED** (Astra
  249264994596077568; Opus re-ran the verifier).
  - The table notes §5 freeze the likelihood edge cases (ε_obs floor, ambiguous reads) and the
    primary equal-block weighting.
  - Donor assignment is capacity-based: the respiratory donor actually used in a sample remains
    inferred, not measured.
- **Smoke tests:** the D001 and D002 synthetic runs score all four models finitely, but each
  is below the independent-block gate, so neither gives a verdict.
- **Before freeze (all resolved 2026-11-20):** grid's go (249688876280889344); δ frozen
  (Estimand section); the calibration and the synthetic benchmark are specified below. The data-gate
  check and the calibration *run* after freeze, outcome-blind, as specified.

## What P6 can and cannot show
- **It tests** whether a parameter-free, thermodynamically constrained rate law adds held-out
  predictive information about the community's *expressed* allocation among catabolic
  pathways. That information has to go beyond genomic capacity and ordinary (concentration-only)
  kinetics.
- **It does not** measure fluxes J, and it does not show that dissipation is maximised (MEP).
  A positive result means "thermodynamic driving force carries information about expression";
  it does not mean "the community maximises entropy production".
- It is the most direct test available of the thermodynamic reading (T) of grid's free-energy
  perspective (manifold note §1). The variational (V) and Bayesian (L) readings make no
  prediction here.

## Question
In samples with paired metatranscriptomes and chemistry, does the share of transcripts
allocated to each catabolic pathway r follow:

  capacity_r × F_K,r(c) × F_T,r(Δ_rG_r/n_e) × (−Δ_rG_r/n_e)

better than capacity_r × F_K,r(c) alone?

## Data gate (outcome-blind; do first)
1. **Metadata check.** Using sample metadata only, with no RNA counts opened, list every local
   sample with all of the following:
   - a paired DNA + RNA library from the same bottle, cast and depth;
   - temperature, pH and conductivity;
   - activities or concentrations, measured in the same water, for the electron donors and
     acceptors of every pathway in the panel. At minimum: O₂, NO₃⁻, NO₂⁻, NH₄⁺, SO₄²⁻, H₂S, CH₄
     and DIC, plus DOC or hydrocarbons for heterotrophy.
2. **Candidates to check:**
   - CMO mesocosms (diesel: is hydrocarbon or O₂ chemistry logged?);
   - Amundsen and Lake Winnipeg casts with a rosette chemistry file (e.g. the LWFC 2025 rosette
     chemistry already in grid's Downloads; check whether RNA is paired);
   - any redox-stratified series.
3. **Minimum to run** (as frozen in notes/p6-pathway-table-sources.md §5.3):
   - ≥ 20 sampling events (one cast × depth = one event; replicate libraries collapsed),
     each with ≥ 100 uniquely assigned panel RNA reads;
   - ≥ 6 independent blocks (block = cast_id, or mesocosm_id for mesocosm series), with
     ≥ 3 blocks containing an oxic event (O₂ ≥ 10 µM) and ≥ 3 blocks containing a suboxic or
     anoxic event (O₂ < 10 µM);
   - ≥ 4 panel pathways detected in RNA in ≥ 50% of events (coverage gate, applied once,
     without selecting pathways).
4. **If the gate fails,** P6 is reported as "not identifiable on local data", and only the
   precondition test (end of this file) runs. External data (e.g. redox-stratified series with
   paired chemistry and RNA) would need grid's choice.

## Pathway panel (frozen before any RNA is opened)
- **Catabolic pathways** (KO sets written into a frozen table, coordinate-joined per trap 3):
  - aerobic organic-C respiration (terminal oxidases);
  - denitrification steps (nar/nap, nir, nor, nos) as separate rows;
  - DNRA (nrf);
  - nitrification (amo, hao, nxr);
  - anammox (hzs/hdh);
  - dissimilatory sulfate reduction (dsr/apr);
  - sulfide oxidation (sqr/sox/fcc);
  - methanogenesis (mcr);
  - aerobic methanotrophy (pmo/mmo);
  - fermentative H₂ production and H₂ uptake (hyd groups).
- Each row gives the stoichiometric reaction (electron donor and acceptor), ΔG°′ at in-situ T,
  and the literature half-saturation constants K for its substrates. Sources are listed with
  page numbers, and **no constant is fitted to our data**.
- **Unbinned contigs are included,** through community annotation (trap 2).
- **Eukaryotic transcripts are excluded from the catabolic panel.** Their share of reads is
  reported separately (multiscale-sim closure artefact).

## Predictors (all parameter-free and fixed before unblinding)
For each event e and pathway r:
- **cap_r,e**: length-weighted DNA pathway capacity,
  Σ_genes∈r (depth_gene,e × effective_length_gene), normalised by prokaryotic single-copy
  marker depth so that eukaryotic and viral read shares don't move it.
  - The marker assignment is the same in DNA and RNA.
  - Example: two equal-length markers at equal depth give twice one marker's expected reads,
    not four times (Astra 195177c015c5).
- **F_K,r(c)**: the product of Monod terms over the substrates, using the frozen K values,
  **times frozen O₂ inhibition K_I,O2/(K_I,O2 + [O₂]) on strictly anaerobic rows**
  (nar/nap, nir, nor, nos, nrf, hzs, dsr-reductive, mcr, hyd-uptake). The same F_K is used in
  every model, so ordinary O₂ regulation (FNR/ArcA) can't masquerade as a thermodynamic win
  (Carbon 249248065344946176).
- **Gene length enters once, inside cap_r,e.** No separate L_r multiplier is applied, so
  subunits and length are counted exactly once for raw reads (Carbon, Astra).
- **n_e,r**: electrons transferred per reaction as written. **Every ΔG enters per electron
  (Δ_rG_r/n_e)**, so the arbitrary choice of how a reaction is written (e.g. per 1 vs 3 DOC)
  cannot rescale its power (Carbon 249247746036776960).
- **F_T,r**: the thermodynamic potential factor 1 − exp((Δ_rG_r/n_e + ΔG_min)/(χRT)), clipped at 0,
  with χ and ΔG_min fixed from the literature (LaRowe-type). This is Vallino-style.
- **Δ_rG_r** = ΔG°′(T) + RT ln Q, from in-situ activities.

**Models.** The predicted share is π_r,e ∝ the product shown, normalised over r.

| Model | π_r,e ∝ | Role |
|---|---|---|
| M_cap | cap | capacity explains expression (null 1) |
| M_capK | cap · F_K | capacity × ordinary kinetics (null 2; **the baseline P6 must beat**) |
| M_capKT | cap · F_K · F_T | adds the thermodynamic gate only (near-equilibrium / endergonic throttling) |
| M_power | cap · F_K · F_T · (−Δ_rG/n_e) | adds energy-weighted allocation per electron (primary) |

M_capKT separates "the reaction is near equilibrium, so it is throttled" from "more power,
so more allocation".

## Outcome and likelihood (frozen)
- **Outcome:** raw read counts n_r,e over the panel per event, plus an "other" category. No CLR.
- **Likelihood:** Dirichlet-multinomial with mean π_r,e and **one common dispersion φ**. φ is
  estimated once, on training folds only, under M_capK, and then held fixed for every model in
  that fold.
- **Score:** held-out log predictive density per event; the mean across events is S_M.

## Estimand and decision bands
- **Primary estimand:** Δ = S_M_power − S_M_capK, in nats per event.
- **Mechanism split (pre-registered secondary):** Δ_gate = S_M_capKT − S_M_capK and
  Δ_weight = S_M_power − S_M_capKT.
  - In alife/mep_redox v3 (T5), random catalyst fields and electron-yield selection scored
    nearly as high on RNA ~ power as energy-yield selection, because F_T and O₂ inhibition
    already gate which pathways can carry flux.
  - So a supported Δ that comes entirely from Δ_gate shows thermodynamic *feasibility*
    information, not energy-*proportional* allocation. Report which it is.
  - Neither component measures closeness to a global MEP optimum.
- **Folds:** blocked cross-validation by cast or date (all depths of a cast are held out
  together), with a leave-one-campaign-out check.
- **Interval:** a 95% cluster bootstrap over events (clusters = cast or campaign), studentised.
- **Margin: δ = 0.05 nats per event (FROZEN 2026-11-20, before any calibration output).**
  Rationale: on the equal-block mean score, 0.05 nats per event is a predictive likelihood
  ratio of e^0.05 ≈ 1.05 per event, i.e. about e^1 ≈ 2.7 across the 20-event minimum cohort
  (the cohort figure holds exactly only under equal block sizes, since the score is equal-block
  weighted).
  A smaller gain would not change which pathways one expects to be expressed. Simulations
  estimate operating characteristics only and do not move δ.
- **Folds, concretely:** leave-one-block-out. φ̂ is fitted by maximum DM likelihood under
  M_capK on the training blocks (bounded 1-D search over log φ ∈ [log 10⁻⁴, log 10]); all four
  models are then scored on the held-out block at that φ̂.
- **Interval, concretely:** equal-block weighting makes Δ the mean of per-block mean
  differences d_c; the 95% interval is a studentised (bootstrap-t) cluster bootstrap over
  blocks, 9,999 resamples, fixed seed 20261120. Held-out scores are not refitted inside the
  bootstrap.

| Verdict | Rule |
|---|---|
| **Supported** (thermodynamic information beyond kinetics) | lower bound of Δ > +δ. To claim that *allocation follows reaction power* (not just feasibility gating), also require the lower bound of Δ_weight = S_power − S_capKT to be > 0 (Carbon) |
| **Equivalent** (predictive difference within ±δ) | interval inside [−δ, +δ] |
| **Against** (the thermodynamic form predicts worse) | upper bound < −δ |
| **Inconclusive** | otherwise, or the data gate is not met |

**Secondary results** (reported, but no verdict depends on them):
- Δ_KT = S_M_capKT − S_M_capK;
- S_M_capK − S_M_cap;
- per-pathway calibration;
- the same analysis with DNA-only capacity replaced by MAG-resolved capacity.

## Design-only calibration before unblinding
- Use alife/mep_redox/synthetic (blinded D001–D123, which mixes 27 emergent runs with 96 fixed-enzyme yardstick runs, with the key withheld until the run is
  scored). Its world has known kinetics and thermodynamics, and its simulator found that RNA
  tracking of predicted power separates energy-based from electron-based selection
  (v3 T5: Spearman 0.76–0.80 for energy-yield selection vs 0.69; random catalyst fields 0.73, which shows that F_T gating alone gives much of the signal).
- **Operating characteristics come from count-generating controls** (Astra 249247782548193280,
  249247817151201280):
  - Simulate counts with known M_capK means, and separately with known M_power means, at matched
    depths and blocks, plus a known near-zero score contrast to calibrate "Equivalent".
  - If M_capK generates the counts, the expected Δ is minus a predictive KL divergence, so an
    "Against" verdict can be correct.
  - The F_T-off mep_redox world is only a mechanistic sensitivity check, because F_T off still
    leaves the −ΔG multiplier active.
- **δ is set as the smallest scientifically meaningful score difference, before any
  calibration outcome is seen.** The simulations estimate power and error rates only; they
  don't set δ.
- Note: that simulation's chemistry is synthetic and its per-sample genome presence is
  saturated (rule-graph caveat). Calibration there tells us whether the pipeline can detect
  the effect, not how large the effect is in nature.

## Synthetic benchmark (runs after freeze; scorer blind to scenario labels)
- **Data:** all 123 blinded datasets in alife/mep_redox/synthetic/ (each = one cast of 11
  depth events). Observer-side files only; scored with the frozen table and rules.
- **Independent block = provenance group** (Astra 249687826194939904): an emergent run and
  the fixed-enzyme yardstick variants derived from it form one block, so 27 blocks (8
  scenarios × 3 seeds, plus 3 forced runs). The dataset→block map is generated from
  truth/unblinding_key.tsv by a short audited script that writes **only** dataset_id →
  random block label; scenario names, seeds and emergent/yardstick status are not written or
  printed. The block map's hash is posted before scoring.
- **Scope (addendum before scoring, Astra 249689626339885056):** run_all.py:53,70 reuses the
  same three seed-specific catalyst draws and starting-state RNG seeds across scenarios, and
  lines 139,161 reuse the random-field and observation seeds. So a held-out scenario×seed block
  can have matched-seed relatives in training. The 27-block result is therefore a **27-cell
  conditional benchmark** (conditional on these three seed families), **not** 27 independent
  environmental replicates. Inference over fresh seeds would need whole seed families held
  out; that gives 3 families, below the 6-block gate, so it is not claimed. The analysis itself
  is unchanged.
- **Primary benchmark output:** Δ, Δ_gate, Δ_weight with the frozen bands, pooled over all
  27 blocks. **Secondary, after scores are posted:** unblind and report Δ per scenario family
  and emergent vs yardstick. This answers "can P6 tell energy-yield selection from the
  others?"; it is a pipeline check, not evidence about nature.
- **Count-generating calibration (Astra 249247782548193280):** using the benchmark's own
  per-event predictors and panel read depths (no RNA counts), simulate DM counts at
  φ ∈ {0.01, 0.05, 0.2} from (a) M_capK means, (b) M_power means, and (c) a placebo
  π = 0.98·π_capK + 0.02·π_power. Run the full frozen pipeline (same blocks, leave-one-block-out
  φ̂, bootstrap with 1,999 resamples for speed) for 200 replicate cohorts per cell. Report
  the frequency of each verdict and the Monte-Carlo true Δ per cell. This estimates operating
  characteristics only. A placebo or other verdict is called an error only against that cell's
  simulated true Δ, not against an assumed Δ ≈ 0 (Astra 249689626339885056).

## KO and direction collisions to resolve in the table (Carbon 249248065344946176)
- **amoABC and pmoABC share KOs** (K10944–K10946). Separate them by hao (K10535) or mmoX
  (K16157) co-carriage, or by clade-resolved MAG assignment.
- **narGH and nxrAB share KOs** (K00370, K00371), and the two run in opposite redox
  directions. Assign by MAG taxonomy (e.g. Nitrospira), and exclude these KOs from unbinned
  contigs; napAB (K02567, K02568) is an unambiguous reduction marker.
- **dsrAB and aprAB act reversibly.** Sulfur oxidisers (SUP05, Thiobacillus) run reverse
  rDsr. Split reductive from oxidative by MAG phylogeny or synteny (dsrEFH or sox
  co-carriage), and restrict unbinned reads to directional markers.
- **Terminal oxidases measure total O₂ respiration,** by lithotrophs (SOB, AOB, NOB, MOB) and
  organoheterotrophs alike, not organic-C oxidation alone. Either model O₂ respiration as one
  row with a donor mixture, or assign oxidases by MAG.

## Must be settled before freeze (Astra 249247782548193280) — RESOLVED
Status 2026-11-20: item 1 is settled in notes/p6-pathway-table-sources.md §§1–4 and the frozen
table; item 2 in §5.1–5.2 (ε_obs = 10⁻⁴ floor, panel-conditioned DM, N_panel ≥ 100, read
assignment); item 3 in §5.3 (block = cast_id/mesocosm_id, ≥ 6 blocks, ≥ 3 per regime,
equal-block weighting) and the Data gate above; item 4 in the decision bands. The text below is
kept as the record of what was asked.
1. **The actual reaction / marker / constant table,** with each of the following frozen:
   - reaction direction;
   - the electron donor for each pathway (terminal oxidases alone don't identify an organic
     donor, and bulk DOC has no defined activity: choose a model donor and state it);
   - the required product activities, adding H₂, NO and N₂O to the minimum chemistry for the
     steps that need them;
   - how to handle chemistry below detection;
   - the pH convention for ΔG°′, units, and the turnover basis.

   Per-electron weighting fixes the arbitrary scale of how a reaction is written only when the
   brake and the yield share that basis. Reading it as power per enzyme also needs turnover on
   that basis.
2. **Likelihood edge cases.**
   - "other" has counts but no model probability. Condition on uniquely allocated panel reads,
     or specify a common outside-panel mass.
   - cap = 0 or F_T = 0 with positive RNA gives a log score of −∞, and all-zero weights can't
     be normalised. Freeze a common observation-error floor or an explicit eligibility rule,
     and add tests for these cases.
   - Freeze how ambiguous reads are assigned.
3. **Independence and weighting.**
   - Twenty depths from one cast can't support a cast bootstrap. Freeze the independent block,
     the minimum number of blocks per redox regime, and the weighting.
   - Apply the RNA-detectability gate once, as a coverage gate, after freezing and without
     selecting pathways by outcome.
4. **The paired energy-magnitude contrast is Δ_weight = S_power − S_capKT.** It is
   pre-registered as above.

## Known threats (all reported, none fixable post hoc)
1. **Measurement location.** Activities measured in bulk water, not in particle or porewater
   microenvironments (P4). Report the share of events with O₂ near the detection limit.
2. **Expression is not flux.** There is post-transcriptional control, and constitutive
   expression (mep_redox switch).
3. **Overlapping KO sets** between pathways, e.g. nitrate reductase in DNRA and
   denitrification. Freeze the assignment rule and run a sensitivity analysis dropping shared
   KOs.
4. **Literature constants are uncertain.** Run a sensitivity analysis with all K × 2 and × 0.5,
   frozen in advance. The verdict uses the central values only.
5. **Eukaryotic read closure:** report the eukaryotic share per event and stratify by it.

## Precondition test (runs whatever the gate says; not a test of MEP)
- **Pooled layer ordering on CMO.** For each HQ MAG with enough coverage:
  D = log(τ̂_DNA / τ̂_RNA), from OU + noise fits after removing seasonal phase.
- **Prediction:** D > 0 (RNA faster).
- **Test:** NOT a sign test against 0.5. At CMO's n ≈ 20, unequal DNA and RNA measurement noise
  biases separate fits: in the manifold simulation the null gives P[D > 0] = 0.67 with median
  D = +0.38, so a sign test against 0.5 with 30 MAGs has about a 60% false-positive rate
  (Carbon 249247746036776960).
- **Instead,** fit a joint bivariate (DNA, RNA) state-space model per MAG with separate noise
  terms, test τ_DNA > τ_RNA against τ_DNA = τ_RNA by parametric bootstrap from the H0 joint fit,
  and generate the **full MAG panel jointly** at the observed dates in each bootstrap draw.
  Independent per-MAG draws on the same calendar would still erase cross-MAG dependence
  (Astra 5556e4aef76d).
- **Caveat (Astra 249247324022685696):** the simulated pooled power assumes MAG signs are
  independent. Power for this test comes from the joint parametric bootstrap above, which
  supersedes the earlier event-block permutation idea.
- **Interpretation:** a positive result shows only that layers are separable at CMO cadence.
  It is expected under almost any biological model. A negative result is informative against
  layer separability at that cadence.

## Compute
- **Table-scale.** It needs RNA counts per KO, coordinate-joined (existing illumina_rna
  outputs), the marker-normalised DNA coverage, and the chemistry table.
- **Cost:** under 1 CPU-hour, plus the design calibration on synthetic data (about 1
  CPU-hour). No cluster job is needed unless RNA must be remapped; that needs grid's go.

## Freeze record (2026-11-20, claude-plexi-2)
- Title and status set to FROZEN; grid's go 249688876280889344.
- δ fixed at 0.05 nats/event with rationale; fold and bootstrap details made concrete.
- Data-gate minimum aligned with the already-reviewed table notes §5.3 (it previously said
  "≥ 2 redox regimes").
- Synthetic benchmark and count-generating calibration specified (provenance blocks).
- "Must be settled" list marked resolved with pointers. No substantive model, predictor,
  likelihood or band was changed.
- Addendum before any scoring (same day): synthetic benchmark scope as a 27-cell conditional
  benchmark (shared seed families), placebo-error rule, and the equal-block caveat on the δ
  rationale (Astra 249689626339885056). Analysis unchanged.
