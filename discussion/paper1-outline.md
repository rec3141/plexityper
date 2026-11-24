# Paper 1 Outline: Plexitrophy, the Plexiome, and Recurring Aquatic Assemblages (DRAFT)

- **Status:** Discussion outline for `grid` (`248931707369472000`), written in `discussion/` (`manuscript/` left untouched pending `grid`'s ownership decision).
- **Sole writer:** `agy-plexi` (Carbon, Gemini 4 Carbon), reflecting three-family consensus in [discussion/plexi-concept-review.md](plexi-concept-review.md) and draft protocols `P1` (`812d46c677e2`), `P2` (`828c0c3e0dde`), `P3` (`33de07d6ebf6`), and `P4` (`0ed5e0034020`).
- **Pre-registration discipline:** All empirical claims in the abstract and Results sections below hold explicit `[ Supported | Falsified | Inconclusive ]` slots governed by **proposed** decision bands in `protocols/` (to be calibrated on synthetic nulls and frozen prior to unblinding); **Paper 1 is structured as a registered-style test so a negative/null outcome is equally publishable and informative.**

---

## 1. Candidate Working Titles

1. **Primary recommendation (leads with the strongest uncoined mechanism + the measurable object, 2 coinages in title):**
   > **Plexitrophy and the Plexiome: Testing Irreducible Multi-Partner Nutritional Dependence in Aquatic Microbial Assemblages**
2. **Registered-question alternative (leads with the Doolittle song-vs-singer / Black Queen puzzle, 1 coinage in title):**
   > **When One Helper Is Not Enough: Testing Plexitrophy and Compositional Recurrence Across Freshwater and Marine Metagenomes**

*(Terminology budget in main text: introduce **three** core terms — **plexitrophy** [plus subtype *anabolic plexitrophy*], **plexiome**, and descriptive **plexitype**. Keep **plexibiont** and **plexon** confined to a concise Discussion box / Glossary unless an empirical test crosses both out-of-event recurrence and non-substitutable interdependence or uses individually isolated particles; keep **holobiont** in its historical host+microbiota meaning and omit **dictyobiont/catenabiont** from main text.)*

---

## 2. Abstract Skeleton (150 words; bracketed slots filled strictly after protocol freeze)

> Microbial communities routinely exchange metabolites and electrons, yet existing concepts conflate facultative cross-feeding, taxon-blind public goods, and pairwise syntrophy. Here we define *plexitrophy*—obligate nutritional dependence, under stated external supplies, on metabolites or electrons transferred between coexisting living partners such that no single partner suffices ($d_{\text{total}} \ge 2$ simultaneous suppliers or an irreducible cycle $c \ge 3$, excluding ingestion)—and formalize the *plexiome* as an evidence-graded AND/OR dependency hypergraph separating required functional roles from substitutable donor guilds and recurring taxon subsets (*plexitypes*). Motivated by two-input chemical rescue in pseudocobalamin-remodeling microalgae and obligate bacterial tricultures combining $\text{H}_2$ and corrinoid suppliers, we evaluate *anabolic plexitrophy* ($d_{\text{anab}} \ge 2$ or $c_{\text{anab}} \ge 3$) across dereplicated aquatic metagenomes. Using completeness-calibrated support enumeration with unresolved-donor bounds, we test whether multi-partner auxotrophy exceeds fragmentation and helper-stacking nulls [P2], enriches on coarse particle fractions [P4], and separates temporal compositional transfer [P3] from geographic role recurrence across supplier replacement [P1].

*(Exact word count of blockquote: **150 whitespace-delimited words**.)*

---

## 3. Section-by-Section Outline (Outcome-Contingent)

### 1. Introduction: Three Confounded Axes in Microbial Interdependence
- **The gap between dyadic syntrophy and taxon-blind Black Queen public goods:**
  - *Syntrophy* captures obligate thermodynamic/metabolic coupling, historically dyadic ($m = 2$).
  - *Black Queen Hypothesis* (Morris et al. 2012, `PMC3315703`) captures adaptive gene loss of leaky functions, but its "Shooting the Moon" dynamic predicts function-stacking into single rare helpers ($d = 1$ hubs/nested hierarchies) and an additive, donor-blind pool.
  - *Nonspecific facilitation* (Goldford et al. 2018, `PMC6405290`) produces dense cross-feeding meshes among generalists that all grow solo ($d = 0, m = 1$).
- **Motivating empirical cases (keeping proven vs. open links explicit):**
  1. *Mixed catabolic + anabolic bacterial triculture:* Men et al. (2014, `PMC3993147`) — *Dehalococcoides mccartyi* 195 fails in dyadic coculture with either *Desulfovibrio vulgaris* ($\text{H}_2$, $d_{\text{catab}} = 1$) or *Pelosinus fermentans* (phenolic corrinoid, $d_{\text{anab}} = 1$) alone, growing sustainably only in the triculture ($d_{\text{total}} = 2$, with DMB supplied externally in medium).
  2. *Candidate two-input anabolic remodeling in plankton:* Cyanobacteria synthesize pseudocobalamin while lacking `bluB`/`bza` (Helliwell et al. 2016 `PMC4850488`; Heal et al. 2017 `10.1073/pnas.1608462114`), whereas remodeling microalgae (*Pavlova lutheri*) rescue growth only when purified pseudocobalamin and DMB are supplied together (in situ living donor count $d_X$ still open), and marine bacterial $\alpha$-ribazole/cobinamide partial salvagers co-assemble B12 (`c = 2` pair in Wienhausen et al. 2024, `PMID:38720071`).
- **The organizing empirical question (Doolittle's song vs. singer made testable):** When an aquatic assemblage reappears across stations or years, what recurs — the taxon subset (*plexitype*), the constrained AND/OR dependency structure (*plexiome*), both, or neither — and when does blocked helper stacking force $d_{\text{anab}} \ge 2$ or $c_{\text{anab}} \ge 3$?

### 2. Framework & Definitions (Table 1 + Figure 1)
- **The four orthogonal axes:** *plexitype* (descriptive conditionally predictive taxon subset), *plexiome* (evidence-graded AND/OR dependency hypergraph preserving cross-role compatibility constraints, cf. Frioux et al. 2018 `PMC6129287`), *plexitrophy* ($\forall (S, E_S) \in \mathcal{M}(X, E): d_{\text{total}} \ge 2 \lor c \ge 3$, decomposed into $d_{\text{anab}}, d_{\text{catab}}, d_{\text{total}}$, strictly excluding serial chains $d=1, c=1, m \ge 3$ and dietary ingestion), and *architecture vs. host association* (orthogonal axes; $n=3$ collider vs. $n \ge 4$ separator/treewidth identifiability rules).
- **Nine methodological traps in metagenomic interaction inference:** Draft-MAG fragmentation false gaps, multi-binner duplication in `per_mag` tables (Trap 9), coarse/any-subunit KEGG module tables (Trap 8) & archaeal *Thaumarchaeota* `bluB`-without-`cobT` exceptions, unbinned donor pools, locus-tag drift in paired RNA/DNA (`ISSUE_ko_misalignment.md`), ambient dissolved supply vs. donor specificity, facilitation vs. obligacy, genome-size streamlining confounds across filter fractions, and bulk-filter pooling vs. single-particle `plexon` boundaries (Stephens et al. 2024 `PMC10811735`).

### 3. Results (Required Empirical Modules & Proposed Decision Branches)
- **Result 1 — Calibrated prevalence and mechanism of prokaryotic anabolic plexitrophy (`P2`, `protocols/P2-plexitrophy-recipient.md`):**
  - Denominator flow from consensus (`dastool-*`) $\to$ 95%-ANI `skani` dereplicated HQ prokaryotic MAGs across freshwater and marine cohorts, scored with `+`/`,` subunit-aware KO rules verified on unit fixtures and Nanopore contig-interior synteny.
  - **P2-Q1 outcome:** Joint bootstrap-identification bounds $[a - d, b - c]$ for `anabolic plexitrophy` excess over Null 1 (contig-block downsampling) and Null 2 (genome-size/prevalence-preserving pathway permutation) against proposed margin $\delta_p = 0.05$ (pending calibration and freeze), keeping `d=0` ambient rescues and `unresolved` calls inside fixed $N_{\text{cand\_aux}}$.
  - **P2-Q2 outcome (stratified by external supply $E$):** Whether cross-regime incompatible recipient pairs (`R-I3` aerobic B12 + anoxic reduced-S; `R-I2` under DIN-depleted $E_{\text{minimal}}$) show $\ge 2\times$ higher multi-partner odds (proposed $\delta_{\text{OR}} = \ln 2$, pending freeze) than matched redox-neutral $\mathcal{C}_0 \times \mathcal{C}_0$ pairs, alongside the `R-I6` / `C_THI_SPLIT` single-regime wiring calibration stratum and single-cell compartmentalization/respiratory-protection counts ($S_{\text{spatial}}, S_{\text{temp/resp}}, S_{\text{unexplained}}$).
  - *If falsified ($U < \delta$):* Demonstrates that single-helper "Shooting the Moon" stacking, complete-pathway donors, or draft-MAG fragmentation explain apparent multi-partner auxotrophy in aquatic metagenomes — an equally high-impact corrective finding.
- **Result 2 — Conditional predictability and membership specificity of recurring taxon subsets (`P3`, `protocols/P3-plexitype-recurrence.md`):**
  - Forward-year split transfer (2024 $\to$ 2025 across distinct Lake Winnipeg station-events — testing cross-year spatial transfer per Trap 6 — plus multi-year marine series where gated $\ge 10$ independent events pass): regularized JSDM conditional Bernoulli log-predictive density gain of designated members ($M_m$) over abiotic + `V \ (C1 ∪ N_k)` background gradient ($M_e$) **and** Mahalanobis-matched outside-clique nonmembers ($M_n$, including strong outside predictors) against proposed $\delta = 0.05$ nats (pending calibration and freeze).
- **Result 3 — Microhabitat enrichment across size fractions (`P4`, `protocols/P4-size-fraction-pilot.md`, extended to $\ge 10$ verified multi-fraction events if HS22-003 pilot is identifiable):**
  - Paired within-event comparison of the **exact same dereplicated recipient MAG $\times$ endpoint** across verified coarse (`3–20 µm`) vs. fine (`0.22–3 µm`) fractions (with equal uniquely prokaryotic-mapped read-pair subsampling and `U` partitioned into ambiguous vs. undetected donor support; evaluated against proposed margin pending freeze).
- **Result 4 — Geographic role recurrence across supplier replacement (`P1`):**
  - Worked descriptive illustration from CMO mesocosms (`Cycloclasticus metabat_9342` in 2024–25 vs. `lorbin_184` [`sp001566825`] in 2025–26 aromatic dioxygenase transcription, framed explicitly as genus-level degradation role turnover with $n=1$ pool/arm), paired with the primary multi-station spatial DNA+RNA evaluation on `Tara Oceans Table W6` (`protocols/P1-song-singer-tara.md`, proposed $\pm 0.1\text{ nat}$ equivalence margin pending calibration and freeze) and benchmarked against experimental `btuB` riboswitch treatment response in Wienhausen et al. (2022, `PMC9296465`).

### 4. Discussion: When Does a Consortium Need a Name?
- Synthesizing P2–P4 (and P1): separating descriptive compositional recurrence (*plexitype*) from constrained exchange roles (*plexiome*) and multi-partner nutritional dependence (*plexitrophy*), without invoking holobiont/hologenome group-selection assumptions (Douglas & Werren 2016 `PMC4817262`).
- Concrete roadmap for moving candidates from MAG-modeled `anabolic plexitrophy` to confirmed `plexibionts` (targeted single-particle genomics for `plexon` boundaries, vitamer/precursor amendments, and co-culture; testing `P5` exchange-type-to-architecture predictions once edge sets pass resolution gates).

---

## 4. Figure and Table Plan

- **Table 1 — Operational lexicon and orthogonal axes:** One row each for `plexitrophy` (with `d_anab`, `d_catab`, `d_total`, `c`), `anabolic plexitrophy`, `plexiome`, `plexitype`, plus Glossary scope for `plexibiont` and `plexon`, stating prior art, required evidence grade, and explicit falsifier.
- **Figure 1 — Conceptual & mathematical distinction of plexitrophy and the AND/OR plexiome:**
  - *Panel A (Wiring taxonomy with $(d_{\text{anab}}, d_{\text{catab}}, d_{\text{total}}, m, c)$):* Contrasts facultative $K_4$ facilitation $(0,0,0,1,1)$ (Goldford `PMC6405290`), 1-helper BQH hub $(1,0,1,2,1)$ (Morris `PMC3315703`), serial chain `D1 -> D2 -> X` $(1,0,1,3,1)$ (**not** plexitrophy — includes donor-side remodeling Wiring 3a!), fan-in plexitrophy `D1 -> X <- D2` ($d \ge 2$), and cyclic plexitrophy `A -> B -> C -> A` ($d=1, m=3, c=3$).
  - *Panel B (Motivating empirical cases with exact provenance labels):* Left: Men et al. (`PMC3993147`) groundwater triculture with `DvH` ($\text{H}_2$, $d_{\text{catab}}=1$), `PfR7` (phenolic corrinoid, $d_{\text{anab}}=1$), and medium DMB labelled separately ($d_{\text{total}}=2$). Right: Helliwell et al. (`PMC4850488`) two-input chemical rescue (pseudocobalamin + DMB $\to$ algal `COBT/S/C`, with field supplier count marked open) vs. Wienhausen et al. (2024) *Colwellia*–*Roseovarius* $c=2$ pair.
  - *Panel C (Constrained AND/OR hypergraph):* Illustrates how `AND` across required inputs coexists with constrained `OR` donor combinations (song vs. singer) and why $n=3$ non-colliders vs. $n \ge 4$ trees/meshes differ in conditional-independence identifiability.
- **Figure 2 — Metagenomic pitfalls & calibration diagnostics (`P2` Methods / Traps 1, 8, 9):**
  - *Panel A:* Effect of Trap 9 (raw multi-binner `per_mag` vs. `dastool-*` + 95%-ANI `skani` dereplication) and Trap 8 (`+`/`,` complex grammar vs. any-subunit `evaluate_module`) on false donor/auxotrophy calls.
  - *Panel B:* Empirical false-gap curves from contig-block downsampling of $\ge 95\%$ complete Nanopore MAGs + example Nanopore contig synteny plot contrasting a true internal locus deletion ($\ge 5\text{ kb}$ conserved flanking genes at uniform read depth) against a contig-edge truncation (`unresolved`).
- **Figure 3 — Calibrated recipient-level anabolic plexitrophy and blocked-stacking test (`P2` Results):**
  - *Panel A:* Sankey/flow denominator diagram from consensus MAGs $\to$ dereplicated HQ MAGs $\to$ $N_{\text{auto}}$ vs. $N_{\text{cand\_aux}}$ partitioned into `d=0 (ambient E)`, `d=1`, `plexitrophic`, `U_ambig`, and `U_no_donor`.
  - *Panel B:* P2-Q1 observed-minus-null identification intervals $[a-d, b-c]$ (with 95% event-cluster bootstrap bars and proposed $\delta_p = 0.05$ line, frozen prior to execution) against Null 1 (fragmentation) and Null 2 (Shooting-the-Moon permutation).
  - *Panel C:* P2-Q2 recipient log-odds contrast $\Delta \theta_{\text{recip}}(E)$ across `R-I3` and `R-I2` (stratified by $E$) vs. matched $\mathcal{C}_0 \times \mathcal{C}_0$ pairs against proposed $\delta_{\text{OR}} = \ln 2$ (frozen prior to execution), alongside the `R-I6`/`C_THI_SPLIT` calibration stratum and single-genome bypass bar (`S_spatial`, `S_temp/resp`, `S_unexplained`).
- **Figure 4 — Cross-year compositional transfer (`P3`) and size-fraction contrast (`P4`):**
  - *Panel A–B (`P3`):* Held-out 2025 Bernoulli log-predictive density differences ($M_m - M_e$ and $M_m - M_n$) per Lake Winnipeg station-event against proposed $\delta = 0.05$ nats (frozen prior to execution), showing whether candidate *plexitype* boundaries beat continuous background gradients and strong outside-clique predictors.
  - *Panel C–D (`P4`):* Paired recipient MAG $\times$ endpoint missing-data bounds `[K/N, (K+U)/N]` across `3–20 µm` vs. `0.22–3 µm` fractions (at 50% breadth, 70% breadth, and equal uniquely prokaryotic-mapped read pairs, displaying prokaryotic vs. eukaryotic read-share bars).
- **Figure 5 — Decoupling functional role from taxonomic singer (`P1` descriptive, spatial cohort, and experimental benchmark):**
  - *Panel A (Descriptive time-series illustration):* CMO two-winter trajectory showing genus-level aromatic oxygenase transcription persisting across dominant MAG replacement (`metabat_9342` in 2024–25 vs. `lorbin_184` in 2025–26, coordinate-remapped KOs, plotted with explicit pre-dosing baseline and $n=1$ pool/arm caveat).
  - *Panel B (Primary multi-station geographic test — `Tara Oceans Table W6`, `protocols/P1-song-singer-tara.md`):* Station-level $\Delta\text{ELPD}$ ($M_0 - M_e$ and $M_1 - M_0$ against proposed $\pm 0.1\text{ nat}$ equivalence margin, frozen prior to execution) and adjusted `btuB` supply slope across dominant cobalamin-supplier turnover (including archaeal *Thaumarchaeota* producers) within depth and basin strata.
  - *Panel C (Independent treatment-response benchmark — Wienhausen et al. 2022, `PMC9296465`):* Controlled co-culture / precursor-amendment positive control verifying directional repression of riboswitch-regulated `btuB` transcription under known cobamide supply vs. B12 limitation.
