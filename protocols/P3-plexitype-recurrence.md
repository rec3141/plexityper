# P3: Do designated members predict an unseen member in a new event? (DRAFT, not frozen)

> **Status 2026-11-18 (Opus), after the design-power review**
> (discussion/deep/p3-design-power.md a8f2417a7aed, closed by Astra 248975201261764608):
> - As drafted, P3 is safe against false support (≤ 1% simulated at the gated event counts;
>   5% at E = 5) but underpowered: support 0.23 at ρ = .7 and 0.05 at ρ = .5, with E = 20,
>   rotated targets, λ = 5 and δ = .02.
> - Under the null, clique discovery returns pseudo-cliques 85–93% of the time; true modules
>   are only partially recovered.
> - P3 is therefore reframed as a **precision pilot**: its primary output is the estimated
>   realised-2025 transfer gain with its interval, and verdict bands are secondary.
> - **Redesign to do before freeze:**
>   - compare the current discovery rule against event-blocked, cross-fitted environmental
>     residuals, with an environment-only negative control;
>   - re-run discovery inside the inner training folds when comparing selectors, since current
>     ranking scores carry selection optimism;
>   - adopt one common margin δ = 0.02 and the "N" (no clique) verdict;
>   - rotate masks;
>   - compute Mn as the mean of log scores;
>   - pre-register discovery-inclusive operating characteristics at the reconciled event counts.
>   - apply strong shrinkage **only to the added taxon slopes**, never to Me's environmental
>     covariates. In the v3 sim, shrinking the environmental covariates raised null outsider
>     inclusion from 0.74 to 0.86 (Carbon 248978424362090496).
>   - abundant core: a frozen 2024-defined core of about 200 lineages, with a depth-offset
>     NB/Poisson-lognormal count likelihood plus an occupancy (zero-hurdle) part for seasonal
>     bloomers. Avoid CLR, whose denominator leaks the masked target.

Status: **draft for cross-family review.**
- Drafted by claude-plexi from Astra's design in notes/plexitype-recurrence.md (c25bcc4ca7f7).
- Astra corrects; Carbon reviews. Frozen only after review and grid's go.
- No data have been inspected for this draft beyond Astra's metadata screen.

## Question
Do a candidate plexitype's designated members carry information about a masked member in an
independent event? The comparison is against season, measured water properties, a continuous
background-community gradient, and matched nonmembers. If they do, the candidate earns the
(descriptive) plexitype label at the stated resolution and scale. If not, we describe a seasonal
association or gradient instead.

Definition in force (discussion/plexi-concept-review.md): a plexitype is a recurring, possibly
overlapping **subset of taxa** within samples, not a whole-sample state. It carries no mechanism
claim.

## Data
- **Cohort A (primary): Lake Winnipeg / Namao DNA.** genice-sample-metadata
  genohub_metagenome_metadata.tsv (ca36581041af).
  - Astra's screen: 327 eligible library rows. 2024/25 has the largest literal station–month
    overlap (17 stations, 21 cells); 2023/24 has 12 cells. Counts are before depth, alias and
    protocol reconciliation (corrected per Astra 248921820019212288).
  - Because the same stations recur, the claim tested is **transfer to a new year at sampled
    stations**: not transfer to new stations, and not annual recurrence.
- **Gate before any fitting:**
  - Reconcile station aliases, depth units, field-log dates and library duplicates.
  - Report the number of independent events (cruise × station × depth, one per
    collection ID), with every exclusion and its reason.
  - 10 independent held-out events per year is a **feasibility floor only**. Same-cast depths
    count as one event. Account for station and cruise clustering, and check the detection and
    non-detection counts for each target.
  - Before outcomes are seen, assess expected precision by **design-only simulation** (or with
    external data). If the bands can't be resolved, report the result as *inconclusive for lack
    of design* without fitting.
- **Features:** frozen before inspection. Either a MAG coverage table (dereplicated catalogue,
  competitive mapping; existing outputs only) or an amplicon ASV table at one primer region. Not
  both in the same model.
  - Resolution: species-level MAG cluster (skani 95% ANI) or ASV, stated once.
  - Freshwater and marine cohorts are never pooled.

## Method
1. **Split before discovery.** The **forward** split (train 2024, test 2025) is primary. Keep
   technical replicates, DNA/RNA pairs and same-event samples together.
   - The reverse split may discover different taxa. It is reported as a **sensitivity
     analysis** only and cannot falsify the primary group.
   - Nothing is filtered, normalised or ordinated across the split.
2. **Discovery, in training only.**
   - Find candidate groups by conditional co-occurrence after environment, using a joint species
     distribution model with environmental covariates.
   - For each candidate, freeze the predictor members and **one or more outcome members**,
     chosen before looking at the test data.
   - **One primary group and control-selection rule** is frozen in advance. Any additional
     candidates get multiplicity-adjusted intervals before labels are assigned. All candidates
     found are listed.
2b. **Selection rule, deterministic from 2024 data alone** (so a second analyst picks the same
   test; Astra 248922002031034368):
   - **Eligible taxa:** detected in 20–80% of 2024 events, with at least 5 detections and at
     least 5 non-detections.
   - **Candidate groups:** maximal cliques of 3–8 eligible taxa in which every pair has a
     posterior residual correlation > 0.3, with the 95% interval excluding 0. These come from
     the JSDM with Me covariates fitted on 2024 data.
   - **Ranking:** candidates are ranked by leave-one-event-out cross-validated (Mm − Me) within
     2024 only.
     - Scored with the regularised JSDM conditional (no per-candidate refit), so larger cliques
       are not favoured by overfitting.
     - Tie-breaks, in order: *smaller* group; higher mean pairwise residual correlation;
       lexicographically smallest sorted taxon-ID list.
   - **Primary group:** rank 1.
   - **Outcome member:** the member whose 2024 prevalence is closest to 50%. Tie-break: smallest
     taxon ID. The remaining members are the predictors.
   - **Implementation:** a JSDM available on grid without installs, chosen at registration.
     Software installs wait for grid.
3. **Matched nonmembers.** For each outcome member, choose nonmember predictor sets of the same
   size, matched on prevalence, seasonal niche (month of peak detection) and detectability
   (mean depth-normalised abundance).
   - **Rule:** for each predictor member, take the eligible nonmember nearest in Mahalanobis
     distance on (logit prevalence, circular peak-month distance, log mean abundance) in 2024
     data. Matching is without replacement.
   - **Primary control pool:** all eligible taxa outside the frozen clique, *including* strong
     predictors of the outcome that fail a clique edge to another member (Astra
     248922173955555328). This is the fair test of membership specificity.
   - **Sensitivity check:** a pool screened to residual correlation ≤ 0.3 with the outcome.
   - If there are fewer than 2 × (number of predictors) eligible controls, specificity is
     **unresolved**.
   - Draw 20 matched sets by randomised order of predictor assignment (seed frozen) and average
     Mn over them.
4. **Models,** fit on training data and scored on the test events.
   - **Primary likelihood: Bernoulli detection** with a library-depth offset (Carbon
     248922965156478976). δ is calibrated on this scale. A count-model analysis is a secondary,
     separately margined sensitivity check.
   - **Scoring:** predictions use the JSDM's own regularised conditional predictive
     distribution p(y_outcome | y_predictors, x). No free weights are refit per candidate, in
     either the 2024 ranking or the 2025 scoring.
   - **Me:** month + water properties (temperature, conductivity, depth, nutrients where
     present) + fraction + library depth + batch + a latent gradient inferred from background
     taxa.
     - For each comparison of clique C1 against control set N_k, the background is
       V \ (C1 ∪ N_k), so neither predictor set helped build the gradient.
     - Loadings are fit on 2024 only, and 2025 is projected without refitting.
   - **Mn:** Me + matched nonmember predictors.
   - **Mm:** Me + designated member predictors.
5. **Score** the held-out log predictive density per eligible biological collection. Average
   per target, then per event, then across events. Uncertainty comes from a cluster bootstrap
   over events, with station and cruise as clusters.

## Decision bands
**δ is frozen at registration, before either split uses outcomes.** Proposal: δ = 0.05 nats per
collection, i.e. a geometric-mean predictive probability of the observed outcome about 5.1% higher
(exp(0.05) ≈ 1.051).

| Contrast | Rule |
|---|---|
| **Supported (conditional predictability):** both Mm−Me and Mm−Mn | Lower 95% bound > δ in the primary (forward) split |
| **Falsified (of a gain ≥ δ):** Mm−Me or Mm−Mn | Upper 95% bound < δ in the primary split. A gain of δ is ruled out (a smaller benefit is not); for Mm−Mn, the membership boundary is unsupported |
| **Inconclusive** | Otherwise; or fewer than 10 held-out events; or the gain is confounded with batch or cruise |

A supported result shows **conditional predictability**. It does not show interaction:
unmeasured shared resource pulses remain an explanation.

## Pitfalls to report
- Prevalence-dependent residual correlation (Zurell et al.).
- Effort-dependent detection (Caporaso et al.: one deep sample recovered 95.4% of the OTUs
  pooled over 72 shallow dates).
- Discreteness: this protocol makes **no** whole-sample class claim (Koren et al.).
- Prior art to beat: Wilkinson et al. 2021 (JSDM conditional prediction); Fuhrman et al. 2006;
  Venrick 1971 (abstract only).

## Compute
- Table-scale only, on one CPU node, using existing coverage or ASV tables.
- If remapping is needed: benchmark one representative library and multiply by the eligible
  libraries before asking grid. No cluster job is part of this draft.
