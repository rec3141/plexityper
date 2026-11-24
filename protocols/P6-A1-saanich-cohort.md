# P6-A1: Saanich Inlet external cohort and chemistry supply (PROPOSED, not approved)

- **Status.** Drafted 2026-11-24 by claude-plexi-2 at grid's request ("go ahead and draft the
  amendment"). It needs a cross-family review (codex-plexi) and grid's approval before any
  Saanich RNA count is opened. It does **not** change the frozen protocol
  `protocols/P6-thermodynamic-allocation.md` (9a47534487ba) or the pathway table (6c26362681b8).
  It only adds rules for one external cohort.
- **What has been inspected so far.** Sample inventory and chemistry metadata only: Hawley 2017
  Table S2 (MOESM151 732686f9cf36) and the Dryad doi:10.5061/dryad.nh035 zip (899e821a2bb4). No
  sequence or expression data has been downloaded or opened. Background:
  `notes/p6-external-data-scan.md`.

## 1. Cohort
- **Pairing is a cruise × depth proxy.**
  - The ISA DNA and RNA assay tables share **no** BioSample accessions (codex-plexi
    251092895373115392). Same-bottle or same-cast pairing is therefore *not established*. It
    needs an outcome-blind crosswalk of the sampling records.
  - The **primary estimand is explicitly a cruise × depth proxy pairing**: DNA and RNA taken
    on the same cruise at the same nominal depth. **This differs from frozen P6's
    "same bottle, cast and depth" requirement and needs grid's separate approval.**
  - The **strict set**, meaning pairs whose same-bottle or same-cast status and compatible
    filter fraction are documented by the crosswalk, is reported separately with all gates
    recomputed. It may be empty.
- **Events.** Saanich Inlet station S3 (SI03). An event is a cruise × depth that has both a
  metagenome and a metatranscriptome in Table S2.
  - That gives 44 sample IDs and 43 unique cruise × depth keys.
  - **Multiple IDs for one key** (cruise 75 at 135 m appears twice):
    - First keep only provenance-eligible pairs (fraction compatible, not flagged).
    - If more than one remains and they share a BioSample or bottle, collapse them as
      replicate libraries.
    - Otherwise take the first eligible ID in Table S2 order.
    - All of this is decided from metadata before any RNA is opened.
- **Blocks.** The independent block is the **cruise**, standing in for the cast. Cast IDs are not
  recorded, and the chemistry comes from one cast per cruise. The substitution is reported as
  such. There are 8 cruises: 42, 47, 48, 54, 72, 73, 74 and 75 (2010–2012).
- **Redox regime.** An event is oxic if CTD_O2 ≥ 10 µM and suboxic/anoxic otherwise, using
  Dryad Chemical `CTD_O2`.
- **Gate check on metadata.**
  - Every cruise has at least one oxic and at least one suboxic/anoxic event, giving 8 blocks of
    each.
  - Complete-case chemistry is 42 events. Cruise 73 at 200 m lacks H₂S, CH₄ and N₂O, and is
    excluded.
  - The frozen requirements are met: ≥ 20 events, ≥ 6 blocks, and ≥ 3 blocks in each regime.
  - The RNA coverage gate (≥ 4 pathways detected in ≥ 50% of events, each with ≥ 100 panel reads)
    is applied once, after approval, exactly as frozen.
- **Provenance audit, done before any RNA is opened.**
  - For each event, the metaG and metaT BioSamples (PRJNA247822) must share cruise, depth and
    collection date.
  - Any event whose same-bottle or same-cast status cannot be established is flagged, and
    sensitivity S-prov drops all flagged events.
  - Filter fractions are 0.22 µm Sterivex with a 2.7 µm GF/D prefilter (Hawley Methods). The
    metaG fraction must be confirmed to match, or the event is flagged.

## 2. Chemistry supply (fixed before RNA)

| Quantity | Source / rule |
|---|---|
| O₂ | `CTD_O2` (µM) |
| NO₃⁻, NO₂⁻ | `NO3`, `Mean_NO2` (µM) |
| TAN (NH₄⁺ + NH₃) | `Mean_NH4` (µM) |
| ΣH₂S | `Mean_H2S` (µM) |
| CH₄, N₂O | `Mean_CH4`, `Mean_N2O`. The README gives these in **nM**, so they are converted to µM. |
| N₂ | Air-saturation at in-situ T and S (table notes §1.5). `Mean_N2` is finite for 0 of the 42 retained events, so this applies to all of them. |
| Nondetects (`0`) | ½ DL, with a 0.01 µM fallback (table notes §3.3), applied identically to F_K and ln Q |
| T, S, ρ | Dryad CTD (Temperature, Salinity, Density) at the nearest pressure within ±2 dbar of the sample depth. |
| T, S, ρ: cruise 54 at 200 m | **Documented exception:** the CTD stops at 197 dbar, so the 197 dbar values are used (flagged). |
| T, S, ρ: cruise 42 | No CTD. Use the median at that depth (±2 dbar) over all S3 casts in the same calendar month of other years. If none exist, take the nearest calendar month that has a match. If still none, the event is excluded. Flagged. |
| SO₄²⁻ | Conservative with salinity: 28.24 mmol kg⁻¹ × S/35 (seawater SO₄/salinity ratio). Converted to µM as mmol kg⁻¹ × ρ (kg L⁻¹) × 1000. Sensitivity S-SO4 subtracts ½ × ΣH₂S, after both are in µM. |
| H₂ | Not measured, so the `h2_gated` rows (MCR_MET, HYD_FERM, HYD_UP) are **inactive**. |
| Denitrification tier | N₂O is measured and NO is not, so **Tier B_N2O** applies (DEN_NIRNOR + DEN_NOS), on the common-panel rule (table notes §1.1). |
| **pH, DIC, DOC** | **Not measured anywhere. They are supplied as a pre-registered assumption grid (§3).** |

## 3. Assumption grid for the unmeasured carbonate system and organic donor
None of these values comes from a measured Saanich profile; they are scenario values only.
- **pH, 4 scenarios.**
  - Three are depth-invariant: 7.4, 7.7 and 8.0.
  - One is depth-linear: 8.0 at 10 m falling to 7.4 at 200 m, interpolated on depth.
- **pH scale.** The scenario pH values are read as **total-scale** pH. a_H+ = 10^−pH is used as
  in the frozen code. This mixes scales for the frozen thermodynamic convention, so it is
  stated as an approximation.
- **DIC, 3 scenarios.** 1900, 2100 and 2300 µmol kg⁻¹. Converted to µM as × ρ (kg L⁻¹), using
  CTD density. HCO₃⁻ is speciated from DIC, pH and T with the table's carbonate constants.
- **DOC, 3 scenarios.** 50, 75 and 100 µM C. The model donor is acetate = 0.05 × DOC (table
  notes §3.3).
  - DOC is a substantive assumption. A fixed DOC changes heterotroph shares relative to the other
    pathways after normalisation (Astra 251077299830964224). It also changes ratios *within* the
    heterotroph rows, because K_s(acetate) is 5.0 µM for SRB_ORG but 15.6 µM for the others
    (Carbon 251077656854315008).
- **Central scenario.** Depth-linear pH, DIC 2100, DOC 75.
- **Robustness rule.** The full grid has 4 × 3 × 3 = 36 scenarios, and all of them are scored.
  The report gives band counts plus Δ, Δ_gate and Δ_weight with intervals for every scenario.
  - **All 36 give the same decisive band** (Supported, Equivalent or Against): that band is the
    verdict.
  - **All 36 are Inconclusive:** the verdict is "Inconclusive (insufficient precision across
    the tested grid)".
  - **The bands differ across scenarios:** the verdict is "Inconclusive (assumption-sensitive)".
  - The Δ_weight condition for "allocation follows power" must also hold in all 36 scenarios.

## 4. Analysis, as frozen except where noted
- **Models and scoring.** The four frozen models and the panel-conditioned Dirichlet-multinomial
  are unchanged. So are φ̂ fitted on training folds under M_capK, equal-block weighting,
  δ = 0.05 nats/event, and the studentised cluster bootstrap over cruises.
- **Folds.** Leave-one-cruise-out is primary.
  - Secondary and descriptive only: leave-one-calendar-year-out, with 3 folds (2010, 2011, 2012).
    Calendar years are not verified deep-water renewal cycles (Astra bb339d5d1ba1).
  - This secondary check carries no verdict.
- **Inputs.**
  - Capacity comes from each event's own metagenome: read mapping to a Saanich co-assembly or
    to the published JGI assemblies.
  - RNA counts come from the paired metatranscriptome.
  - Both use the frozen KO table and read-assignment rules (table notes §2, §5.2).
  - Assembly and mapping choices are fixed in a short run spec, posted with its hash before any
    RNA is opened.
- **Sensitivities**, reported without a verdict:
  - S-prov: drop events whose provenance is flagged.
  - S-SO4: correct SO₄ for sulfide.
  - S-SUP05: the frozen MAG-resolved `SOX_NO3` split. SUP05 couples HS⁻ oxidation to NO₃⁻
    reduction, so its nap/nar/nir/nor transcripts otherwise land on acetate-donor DEN rows.
    This donor mismatch is a known caveat (Carbon 251077656854315008).
  - S-K: all K values × 2 and × 0.5 (frozen threat 4).
  - S-eps: ε_obs of 10⁻⁵ and 10⁻³.
  - S-TSmeas: measured T and S only, which drops cruise 42. That leaves 7 cruises and 36
    chemistry-complete keys before the other gates.
  - S-strict: the strict same-bottle/cast set (§1), with gates recomputed.

## 5. Interpretation limits, stated in advance
- Saanich is one station, so a verdict generalises to this seasonally anoxic fjord, not to
  freshwater or the open ocean.
- The synthetic benchmark (results/p6/synthetic_benchmark.json, AGAINST) suggests M_power's
  energy weighting may be the wrong observable for transcripts that track turnover. A Saanich
  AGAINST would be consistent with that and would not by itself refute thermodynamic gating.
  Δ_gate is reported alongside.
- Expression is not flux, and bulk-water chemistry misses microenvironments (frozen threats 1–2).

## 6. Before running
1. Cross-family review of this addendum (codex-plexi), then grid's approval. The approval must
   explicitly cover the cruise × depth **proxy** pairing as the primary estimand.
2. Provenance audit (§1) from BioSample metadata only.
3. Run spec for assembly and mapping, with its hash posted.
4. Data download (SRA/JGI, PRJNA247822). This is large and needs grid's go for compute
   location and storage.
