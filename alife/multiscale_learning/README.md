# multiscale_learning: a layered ecosystem simulator for plexitrophy timescale tests

Single writer: this directory's author (Claude, Opus 5.5, an alife sub-agent for claude-plexi).
Brief: `../BRIEF.md`, including "Additions after Astra's check" and later coordinator notes
(power-law and lognormal tails, an abundant core, phototrophic eukaryotes).
**Status: draft for cross-family review.** Under AGENTS.md it counts as one vote. Nothing was
posted to the relay. No file under `data/` was opened.

What I read: AGENTS.md; BRIEF.md (all of it); the top of `discussion/plexi-concept-review.md`
(summary, definition, P5 statement); `discussion/deep/manifold-timescales.md` (all of it, as on
disk 2026-11-18 ~20:00; its §7 still held `SIM_RESULTS_PLACEHOLDER`); `protocols/P1-*.md` and
`protocols/P3-*.md` (full drafts); `notes/plexitype-recurrence.md` (first 40 lines).

**What this is.** A world in which a community is literally a sum of layers that adapt at
different rates. The fast layers are the "activations" (expression). The slow layers are the
"weights" (genome content, SNVs). Dependencies are not coded in: they **emerge** through
Black-Queen gene loss in the genome layer and through regulated repression of a lineage's own
synthesis in the expression layer. Ground truth is logged. Synthetic observations are emitted in
the team's table shapes, the team's inference ideas are run on them, and the results are scored
against the truth.

**What this is not.** A simulation shows that a signature is *possible* under its assumptions.
It does not show that the signature occurs in lakes or oceans. Every rate below was chosen by me
and was not fitted to data.

---

## 1. Mechanisms (what is simulated)

Code: `msl/model.py`. Time is in days. Integration is explicit Euler with dt = 0.1 d. The slow
layers (genome, SNV) are stepped every 0.5 d, alternating between midnight and noon.

| layer | state | dynamics | characteristic time (set by parameters) |
|---|---|---|---|
| **environment** | temperature T, resource inflow S_in (3 resources), storms, diel light | annual sinusoid + AR(1) anomaly (τ 5 d, sd 1.5 °C); seasonal inflows with their own phases; Poisson storm pulses on resource 1 (~1 per 30 d, ×3, decaying with τ 3 d; the storm *inflow* is not recorded, but its effect on resource concentrations is measured in `env.tsv`); diel light | forcing |
| **phototrophic eukaryotes** | algae A, grazers Z | light × temperature × **B12 (metabolite 0) limitation** × logistic growth; grazing by Z produces bloom–crash cycles (spring and autumn blooms). Algae release DOC into resource 0 (exudate follows the diel light; sloppy feeding), release metabolites 1–3 on grazing, and **consume** usable B12 | weeks to seasonal |
| **expression** (activations) | per lineage: E_syn (own synthesis of 5 metabolites), E_up (uptake transporters), E_cat (allocation over 3 resources) | relaxes with τ_e = 0.2 d toward regulatory targets set by *local* pools. Regulated lineages (50%) repress their own synthesis when the external pool is high. Uptake is riboswitch-like, repressed by supply. Catabolic allocation follows substrate availability | 0.2 d (own); slaved to pools |
| **population** | biomass B_i of 1200 lineages per station (from a regional pool of 1500) | Liebig growth: min(energy-limited growth − synthesis and gene costs, each essential metabolite's supply (own leaky synthesis + uptake)); chemostat dilution D = 0.1 d⁻¹; immigration from a regional pool whose abundances follow a **power-law (Zipf 1.6)** or **lognormal (sd 2.3)** tail, so a long rare tail persists by mass effect and lineages can **bloom from rarity**. A gleaner–opportunist axis (V_max·v, K_S·v^2.5) makes opportunists bloom after storms | ~days–weeks |
| **genome** (weights) | F_i,m: fraction of lineage i's cells carrying biosynthesis gene m | selection s = growth(carrier) − growth(non-carrier) (Black Queen: loss pays when leaked metabolites are available), loss mutation 2e-5 d⁻¹, HGT gain 2e-5 d⁻¹ × community carrier frequency, drift (Ne 1e5) | months–years |
| **nucleotide** | 8 SNV loci per lineage | drift with Ne ∝ abundance (bottlenecks speed drift), weak mean reversion (1/1500 d), 2 loci under temperature-fluctuating selection, locus 0 hitchhikes on changes in gene-0 carrier frequency | years |

Metabolite pools are well mixed. Every lineage leaks a lineage-specific fraction λ ∈ [0.05, 0.45]
of its own synthesis. Lineages start prototrophic (F ≈ 1). **Nothing about dependence is coded
in.** Auxotrophy, d ≥ 2 and cycles arise only from gene loss under leaky supply (the genome layer)
and from repression of a lineage's own synthesis (the expression layer).

Variants used: `leak=False` (no exchange: the negative-control world), `euk='flat'` (algae held at
their mean, with mean light and temperature, so mean fluxes are kept but the bloom cycle is
removed), `vitamer=True` (40% of metabolite-0 producers make a "pseudo" form that only they can
use, like pseudocobalamin), and `perturb` (CMO-like dosing: a toxicant reduces growth of sensitive
lineages for 10 d and adds a substrate usable by 25% of lineages).

**Approximations to keep in mind.** The genome layer is mean-field: a lineage with F = 0.4 still
synthesises at 40% capacity, i.e. its non-carriers are fed by its carriers within the lineage.
Pools are quasi-steady and spatially uniform, so there are no particles or size fractions, and by
the gauge argument in manifold-timescales.md §4.3 identity can matter only through the vitamer
switch. Grazers eat only algae. There is no explicit viral lysis. Every number was chosen, not fitted.

## 2. Ground truth (what is logged)

Code: `msl/truth.py`. Outputs go in `results/truth/<world>/`.
- **Per-layer state time series** (`layer_series_daily_top300.npz`; W_powerlaw): T, S, pools P,
  A, Z, inflows, storms, and B, F, E_syn, E_up, E_cat and SNV frequencies, daily at noon for the
  300 most abundant lineages. The simulator records every lineage at every design sampling time
  in memory; the analyses use those records.
- **Realised transfers versus required partners**, logged separately (Astra's point). For each
  recipient, `edges.tsv` gives `realised` edges (flux attributed to a donor's share of the
  well-mixed pool) and `required` edges (donor ∈ some minimal sufficient set; the value is 1 if
  the donor is in **every** minimal set, i.e. necessary). "EUK" is a donor node like any lineage.
- **Necessity by counterfactual removal**, under the same medium, endpoint and horizon. A
  recipient is scored only if it is viable now (endpoint: net growth μ − D − mortality ≥ 0 in
  the recorded medium). Pools are rescaled to the inflow from a donor subset S plus the recipient
  itself (quasi-steady pool, linear in inflow). Subsets are searched up to size 4 among the 7 top
  donors + EUK. Two horizons:
  - `h0` (expression frozen; minutes);
  - `h1` (expression re-set to its regulatory targets under the counterfactual pools; hours).
    The genome is fixed under both.
  - d_h = the minimal |S| that keeps the recipient viable.
  - **Genome-layer dependence**: d_h1 ≥ 1.
  - **Expression-layer (facultative) dependence**: d_h0 ≥ 1 but d_h1 = 0.
  - Plexitrophic: d ≥ 2, or a directed cycle of length ≥ 3 among support edges.
- **d / cycle per lineage per time** (`d_cycle_per_lineage.tsv`): every 14 days, for the 300 most
  abundant lineages at that time. It gives d_h0, d_h1, layer, cycle3, plexitrophic_h0/h1,
  n_required, n_necessary, n_realised, n_realised_not_required, euk_required and genes_lost.
- **Dynamic knock-out check** (`knockout_validation.tsv`): from full-state snapshots every
  ~300 d, the top necessary donors (plus donors that matter only at h0) have their leak removed.
  The model is then re-run for 20 days with the same forcing, and recipient biomass is compared
  with an unperturbed control.

## 3. Synthetic observations (what an observer sees)

Code: `msl/observe.py`. Tables in `results/observations/W_powerlaw/<design>/` (TSV, same shapes as
the team's tables):
- `events.tsv`: a shared **event_id** for each station × arm × sampling time. The DNA, RNA, SNV,
  amplicon and environment tables use the same key (sample columns `DNA_<event>`,
  `AMP_<event>`), so held-out events cannot leak into training.
- `amplicon_counts.tsv`: an ASV table over all 1200 lineages plus a **chloroplast** row from
  algae. Depth 4e4, Dirichlet-multinomial (α 300).
- `dna_mag_counts.tsv`: MAG read counts, plus `unbinned` (lineages without a MAG) and
  `eukaryote` reads (algal DNA takes read share in blooms: closure). Depth ~2e6 × lognormal
  (sd 0.3), DM α 400, reads ∝ biomass × genome size.
- `rna_counts_long.tsv`: transcripts per MAG × KO (housekeeping, bio, upt, cat), plus unbinned
  and **eukaryote** transcripts (which dominate in blooms). Depth ~2e6, per-gene lognormal noise
  (sd 0.3), DM α 150.
- `gene_cov_ratio_long.tsv`: gene coverage / MAG coverage (Poisson), i.e. the observable of the
  genome layer (carrier fraction F).
- `snv_freq_long.tsv`: alt-allele frequency with Poisson depth from MAG coverage, error 0.5%,
  and calls only at depth ≥ 5.
- `env.tsv`: measured temperature, the three resource concentrations (S0–S2) and a
  chlorophyll proxy, with sensor noise. The storm inflow itself is not a column, but storms
  raise S1, which **is** measured, so storm pulses are largely visible to Me. (Correction after
  Carbon 249001105845305344: the original text said storms were not measured.) Metabolite pools
  are not measured.
- Catalogue: `per_mag_cds.tsv` in the column names `scratch/p2_scorer.parse_per_mag_rows`
  reads (mag_id, contig_id, **protein_id** = `<contig_id>_<n>`, feature_index, start, end,
  strand, **KO**, EC). It includes **CDSs with blank KO**: 40 hypotheticals per MAG and, 30% of
  the time, a remnant of a lost biosynthesis gene. `genomes.tsv` gives mag_id, **domain**,
  completeness_est from housekeeping recovery, contamination, n_contigs, genome_size_bp and
  **contig_lengths**. Contigs are dropped until a completeness target of U(0.6, 1.0) is reached.
  One **Eukaryota** bin (EUK0001, the algae) appears in genomes, CDS, DNA and RNA tables. The
  legacy two-column `per_mag.tsv` (mag_id, KO) is also written.
- KO IDs are synthetic, in the unused range **K90001–K90023** (`synthetic_ko_map.tsv` maps each
  to its role), so p2_scorer's K+5-digit parser accepts them. They are **not** KEGG
  orthologs. I checked that `parse_per_mag_rows` parses MAG0001 (50 CDS, edge status assigned).
  The `SKO_` prefix of the first sweep was dropped in the revision. The MAG → lineage key stays in
  memory (truth side) and is not written to the observation tables.

**Sampling designs** (`observe.make_designs`; day 0 = start of the 6-year observation window):

| design | cadence | n | purpose |
|---|---|---|---|
| daily_ideal | daily at noon, 6 y | 2190 | upper bound (not realistic) |
| weekly_2y | weekly, years 1–2, clock time 07:12–14:24 at random | 105 | primary realistic series |
| weekly_2y_noon | weekly at exactly noon | 105 | aliasing control |
| weekly_2y_drift | weekly, clock time 4.8 h later at each visit | 105 | aliasing demonstration (diel → ~35 d) |
| monthly_4y_gaps | monthly May–Oct (ice gap), years 2–5, 15% missed | ~20 | Lake-Winnipeg-like |
| two_campaigns | two winters (years 3, 4), every 3–4 d for 75 d, control and treated arms | 2 × 43 | CMO-like mesocosms, dosing on campaign day 14 |
| subdaily_burst | every 4.8 h for 14 d in summer | 70 | expression-layer reference |

## 4. Worlds, seeds and runtime

| world | tail | exchange | eukaryotes | purpose |
|---|---|---|---|---|
| W_powerlaw | power law | yes | dynamic blooms | main world |
| W_lognormal | lognormal | yes | dynamic | same seed, depth and forcing (same storm and bloom schedule); only the tail differs |
| W_flat_euk | power law | yes | flat (mean A, Z) | does the bloom cycle dominate prokaryotic spectra? |
| W_no_exchange | power law | **no leak** | dynamic | negative control for P3 and the circuit probes; same lineages as W_powerlaw |
| W_powerlaw_treated | power law | yes | dynamic | perturbed arm for the two campaigns (identical forcing and noise streams) |
| P1 stations | power law | yes, vitamer off vs on | dynamic | 10 stations × 2 worlds; 400 lineages each, drawn from the shared pool; station temperature offset U(−3, 3) °C |

Seeds: master seed 101 (`run_all.py: SEED`). The environment stream is `[seed, 1]`, dynamics
`[seed, 2]`, the lineage draw `[seed, 3]`, the initial state `[seed, 4]`, the species pool
`pool_seed=7` (tail `[7, 99]`), the designs `11`, the P1 stations `5000+k`, and observation
noise `[101, len(world name)]`. `results/run_meta.json` records the runtime;
`results/provenance.tsv` records the code sha256 (12 characters).

How to run (numpy and scipy only, no installs):
```
cd alife/multiscale_learning
python3 -I run_all.py            # full sweep -> results/
python3 -I run_all.py --quick    # smoke test (300 lineages, 2-y burn-in) -> results_quick/
python3 -I run_all.py --p1-only  # P1 stations only
python3 -I run_all.py --revision # rerun truth/KO/P3/bloomers/mechinterp/P1 + observer tables (§6)
python3 -I summarize.py          # results/summary_tables.md
```

Runtime of the full sweep: 1426 s wall (≈24 min) on a shared, heavily loaded 16-core machine,
single process. The longest single simulation (10 years, 1200 lineages) took about 30 s.

---

> **Revised 2026-11-18 after cross-family review** (Astra 248998526021189632, 249001223084490752;
> Carbon 249001105845305344). The truth, P3, P1, bloomer and mechinterp results in §5.0, §5.2
> and §5.4 below are from the **first sweep and are superseded** by §6. The headline "43–71%
> plexitrophic" is **withdrawn**: almost all d ≥ 2 is donor-count dilution of one limiting
> metabolite, not functional fan-in (§6.1). §5.1 (timescales, AR(1) null, perturbation, blooms)
> and §5.3 (aliasing, closure) were not rerun and still stand as written.

## 5. Results

Full tables: `results/summary_tables.md`, built by `python3 -I summarize.py` from the TSVs in
`results/`. Every number below comes from those files (run of 2026-11-18; code hashes in
`results/provenance.tsv`: run_all `2966259d8966`, model `ee6aa412d30d`, truth `e433452be041`,
observe `e7b0d0682e43`, analysis `db18e6254a6e`). Each number is from **one seed per world**, so
treat differences between worlds as single draws, not as estimates with error bars.

### 5.0 What emerged (ground truth, `truth_summary.tsv`)
Denominator: viable lineage-times among the 300 most abundant lineages, every 14 d over 6 years
(16 534 in W_powerlaw).

| world | genome-layer dependence (d_h1 ≥ 1) | plexitrophic, regulated (d_h1 ≥ 2) | expression-only dependence | cycle ≥ 3 | needs EUK | realised donors / required donors (mean) |
|---|---|---|---|---|---|---|
| W_powerlaw | 0.668 | 0.432 | 0.062 | 0.0005 | 0.018 | 6.4 / 1.25 |
| W_lognormal | 0.854 | 0.706 | 0.016 | 0.001 | 0.009 | 6.6 / 1.66 |
| W_flat_euk | 0.669 | 0.405 | 0.060 | 0.001 | 0.012 | 6.7 / 1.17 |
| W_no_exchange | 0 | 0 | 0 | 0 | 0 | 0 / 0 |

- Starting from all-prototroph lineages, Black-Queen loss produced widespread auxotrophy in four
  burn-in years. d ≥ 2 (plexitrophy by the d criterion) arose in 43–71% of viable lineage-times.
  Cycles of length ≥ 3 almost never arose (≤ 0.1%): in a well-mixed pool, fan-in from substitutable
  donors is easy, but irreducible cycles need specific partners.
- **Realised is not required.** Recipients drew on ~6.4 realised donors but needed ~1.25. In 97%
  of viable lineage-times at least one realised donor was not required. Flux-based edges would
  overstate the plexiome about five-fold.
- **Expression-layer (facultative) dependence was rare** (2–6%). Under h1, derepressing a
  lineage's own synthesis rescues it. Most dependence sat in the genome layer.
- **Static necessity versus dynamic knock-out** (`knockout_validation_summary.tsv`, 132
  necessary pairs vs 1258 non-required pairs): after 20 days with a necessary donor's leak
  removed, the median recipient biomass ratio was 0.876 (necessary) vs 0.991 (not required).
  The ratio was < 0.8 in 34% vs 7% of pairs. The static call points the right way, but most
  "necessary" recipients do not crash within 20 days: other donors and immigration compensate. A
  necessity label is horizon-specific, which supports Astra's point about stating the horizon.

### 5.1 Predictions from manifold-timescales.md

**(a) Do layer signals separate in frequency?** (`timescales_by_design_layer.tsv`; τ̂ from an
OU + white-noise likelihood after removing annual harmonics, 25 series per layer; daily series
use the r(2)/r(1) estimator.)

| design | expression τ̂ | population τ̂ | genome τ̂ | SNV τ̂ | adjacent-layer ordering AUC (pop>expr / genome>pop / SNV>genome) |
|---|---|---|---|---|---|
| truth, daily, noise-free | 15.2 d | 150 d | 297 d | 319 d | 1.00 / 0.76 / 0.38 |
| daily_ideal (noisy) | at bounds (80%) | 52 d | at bounds (68%) | 224 d | 0.77 / 0.23 / 0.75 |
| weekly_2y | at bounds (68%) | 9.7 d | at bounds (76%) | 148 d (52% at bounds) | 0.62 / 0.33 / 0.67 |
| monthly_4y_gaps | 60% at bounds | 25 d (52% at bounds) | 14 d (10 series, 60% at bounds) | no callable SNVs | — |
| two_campaigns | at bounds (92%) | 22 d | at bounds (80%) | 68 d | 0.40 / 0.54 / 0.51 |

Medians are from the TSV, and "at bounds" gives the fraction of τ̂ below Δ/2 or above the span.
- In the noise-free truth, expression's **own** τ (0.2 d) is invisible at daily resolution. The
  daily expression series has τ ≈ 15 d because it is **slaved** to pools and abundances, exactly
  as manifold-timescales.md §2.1 and §6.0 predict.
- Population and genome are cleanly ordered in the truth (AUC 1.00 and 0.76). The genome and SNV
  layers overlap (AUC 0.38): drift with abundance-scaled Ne makes SNVs no slower than gene
  content here.
- With sequencing noise and realistic cadence, ordering collapses to AUC 0.3–0.7. Most
  expression and genome τ̂ hit the bounds: gene-coverage ratios are near-binary plus Poisson
  noise, so the genome layer is mostly measurement noise at any cadence.
- Population τ is long in this world (months, not days), because the tail is held up by
  immigration and seasonal succession dominates. "Population = days" is therefore not guaranteed.

**(b) Can an AR(1) / red-noise null mimic the spectra?** Mostly yes. Fitting OU + white noise per
series and comparing the low/high spectral ratio with 49 surrogates at the same times, the null
was rejected in only 4–16% of population series, 8–32% of expression, 12% of genome and
6–16% of SNV series (W_powerlaw: weekly_2y, weekly_2y_noon, two_campaigns; 0–4% for
monthly_4y_gaps). At
realistic lengths, a one-timescale null fits most series even though the truth has several
layers. That is the expected false-negative side of the note's warning ("not red noise is not
evidence for layers").

**(c) Do perturbations propagate fast → slow?** (`perturbation_propagation.tsv`; CMO-like dosing,
treated vs control arm with identical forcing.)
- Truth, daily: the first response (t10) came at 1.5–2.5 d for expression and 1.5 d for
  population, then 2.5–5.5 d for genome and 5.5 d for SNVs. The half-response (t50) ran
  population 4.5 d < genome 5.5–9.5 d < expression 10.5–13.5 d < SNV 16.5–22.5 d.
- So onset follows fast → slow, but the bulk expression divergence comes *after* the community
  reorganises, because expression is slaved to the new community state. "Expression first" holds
  for onset only.
- Observed (every 3–4 d, sequencing noise, one pool per arm): population responded at t50 =
  8.5 d in both campaigns. Expression, genome and SNV were unresolved (NaN, or a t10 later than
  t50 driven by noise). With one pool per arm, the ordering was **not recoverable** from the
  observations. This matches manifold-timescales.md §6.2 (descriptive only).

**(d) Does the eukaryotic bloom cycle dominate prokaryotic frequency structure?**
(`eukaryote_bloom_dominance.tsv`, daily truth.)
- Algae supplied 25% of organic C, with spring and autumn blooms (CV 1.6).
- Over season harmonics + temperature, lagged algal biomass added a median R² of only 0.01
  (population), 0.07 (expression), 0.01 (genome) and 0.02 (SNV).
- With the bloom cycle removed (W_flat_euk), band-power fractions stayed nearly the same
  (population slow-band 0.979 vs 0.978).
- Reading: in this parameterisation, blooms recur at the same season every year, so their
  effect is **collinear with season**. Their *incremental* frequency signature is small.
- Limits: this does not show that blooms are unimportant. It shows that a seasonally locked
  bloom cannot be separated from season by spectra alone. One confound: the flat control
  delivered more mean algal C (0.42 vs 0.25), because flux is nonlinear in A.

### 5.2 P3 recurrence and P1 song-vs-singer against the truth

**P3** (`p3_groups.tsv`): up to 5 groups per world × design × universe, chosen deterministically
from training events (86 groups in total). Support requires the lower bound of (Mm − Me) > 0.05 nat/event and of
(Mm − Mn) > 0.
- **Support: 0 of 86 groups**, in every world and design, including the exchange worlds and the
  daily 6-year upper-bound design (train years 0–2, test years 3–5).
- In W_no_exchange (truth: no edges), there were 0 false supports. Decisions were "against" or
  "inconclusive". The procedure is therefore conservative.
- Where the truth had a direct required edge between the target and a predictor (18 groups across
  the exchange worlds), decisions were 0 support, 8 against and 10 inconclusive. The gains were
  real but small. Example: W_powerlaw weekly_2y, target 511: Mm − Me = 0.038 [0.015, 0.061] nat,
  but members did not beat matched nonmembers (Mm − Mn lower bound −0.008).
- Reading: in a world where most taxa share a few substitutable donors in a well-mixed pool,
  members carry little information about a masked member beyond environment and nonmembers.
  **Recurrence of exchange partners is not a strong compositional signal**, even with daily data.
  A P3 null on real data is therefore weak evidence against exchange.

**Abundant core vs all taxa.** Ranking was frozen on training events with equal weights
(`core_ranking.tsv`).
- In W_powerlaw weekly_2y, the top 200 held 78% of training read share. Their median detection
  frequency was 0.30 and median seasonal occupancy 0.5, and **42% of core taxa were detected in
  under half the months**: a top-ranked bloomer can be absent most of the year. Under the
  lognormal tail, the core was more persistent (occupancy 0.875, only 8% below 0.5).
- Restricting P3 to the core **changed nothing material**. The same groups were mostly
  re-selected, and decisions were identical or shifted between against and inconclusive. The
  20–80% eligibility filter already removes the rare tail.
- **Bloomers outside the core** (rare at some time, < 1e-5, and > 1e-3 at another;
  `bloomers_vs_core.tsv`): 54 in W_powerlaw and 51 in W_lognormal. Of these, 93% and 98% had a
  genome-layer dependence (d_h1 ≥ 1, or ≥ 1 required donor) at one or more evaluated times. The
  rates for core non-bloomers were 89% and 100%, and 0% in W_no_exchange. A core-only analysis
  therefore drops dependent taxa at about the background rate. Dependence is not concentrated in
  the core.

**P1** (`p1_song_vs_singer.tsv`; 10 stations per world, 400 lineages each, monthly × 12; LOSO;
margin 0.05 nat/obs):

| world | truth | guild measure | guild gain [95%] | identity gain [95%] | decision role / identity |
|---|---|---|---|---|---|
| vitamer off (identity irrelevant) | uptake RNA tracks the usable pool: r = −0.98 | guild_all (gene-presence × DNA) | −0.007 [−0.024, 0.005] | −0.010 [−0.029, 0.003] | against / against |
| vitamer off | ″ | oracle: true pool | 0.057 [0.004, 0.127] | −0.004 [−0.015, 0.003] | inconclusive / against |
| vitamer on (identity matters) | r = −0.96; usable fraction sd 0.25 | guild_all | −0.005 [−0.010, −0.002] | −0.008 [−0.022, 0.002] | against / **against (wrong)** |
| vitamer on | ″ | oracle: true pool | 0.010 [−0.006, 0.021] | −0.005 [−0.013, 0.002] | against / against |

- The song exists in the truth: recipient uptake expression is almost a function of the usable
  pool (r ≈ −0.97). The observable "guild" (summed relative abundance of MAGs carrying the
  synthesis KO) correlated only 0.21 (vitamer off) and −0.07 (vitamer on) with that pool. Most
  taxa carry the gene, and algal B12 demand and pool dynamics decouple supply from capacity.
- The P1 decision rule therefore returns **"against" for a song that is present**, and misses the
  identity effect in the vitamer world.
- Even an oracle pool gives only a marginal gain through RNA measurement noise (per-gene
  lognormal sd 0.3 + Dirichlet-multinomial).
- Proposed P1 change (a proposal only; P1 is not mine to edit): validate the guild proxy against
  a supply measurement, or a supply-responsive marker, before treating "no guild gain" as
  evidence against role recurrence.

### 5.3 Which designs identify which layer, and which signals are artefacts

| layer | daily (ideal) | weekly 2 y | monthly 4 y with gaps | two campaigns (3–4 d) | 4.8-h burst 14 d |
|---|---|---|---|---|---|
| expression (own τ 0.2 d) | slaved component only (τ ≈ 15 d) | slaved, noisy; mostly at bounds | no | no | own dynamics barely; span too short |
| population | yes (τ̂ biased low by noise: 52 vs 150 d) | partially (τ̂ 9 d, biased low) | no (n = 16) | onset after perturbation only | no |
| genome (gene content) | no: near-binary, noise-dominated | no | no | no | no |
| SNV | partially (224 vs 319 d) | no | no callable sites at monthly depth | no | no |
| perturbation ordering | yes (truth) | — | — | **not** with one pool per arm | — |

Artefacts that look like signals:
- **Diel aliasing** (`aliasing.tsv`). With the clock time drifting 4.8 h per weekly visit, the
  prokaryotic share of the RNA library (which algal diel transcription modulates through closure)
  shows a spurious peak at **35.2 d**. 22% of its spectral power falls in the 28–42 d band (vs
  2.4% at fixed noon), and sampling hour explains R² = 0.33. Normalising per housekeeping gene
  within a genome removes the artefact (R² ≈ 0.02; no 35-d peak). **Record sampling time and
  normalise within genomes.**
- **Closure by eukaryotic reads** (`closure.tsv`). In relative abundance with euk reads in the
  denominator, the median correlation with the euk read fraction is −0.31 (absolute biomass:
  +0.08), the mean pairwise correlation rises from 0.058 to 0.090, and 15% of taxon pairs flip
  the sign of their correlation. Excluding euk reads from the denominator reduces the shared
  signal (−0.14) but does not remove it. Blooms create a shared negative mode in prokaryote
  proportions.
- **Seasonal collinearity**: the bloom, temperature and succession effects are inseparable by
  spectra (§5.1d).
- **Red noise**: an AR(1) null fits most series at realistic lengths (§5.1b). A spectrum that
  does not beat red noise says nothing against layers, and one that does is not evidence for
  them.

### 5.4 Mechinterp probe: features vs modules, circuits vs plexiome edges
(`mechinterp_probes.tsv`; daily_ideal amplicon data on the 200-taxon core; module alignment = best
|cosine| between a component's loadings and a module indicator, against 300 random modules of the
same size.)
- **Features.** PCA components align strongly with environmental modules: the season groups
  (0.70 vs a null of 0.13) and algal-bloom followers (0.60). They also align with resource
  guilds (0.26–0.57) and with auxotroph sets (0.24–0.34; all at or above the 99.7th percentile of
  random modules). NMF aligns more weakly (season 0.21–0.27; auxotroph sets 0.15–0.17, 88–100th
  percentile). Auxotroph sets overlap the season and guild structure, so this alignment does not
  show that the components encode dependence.
- **Circuits: the probe fails its negative control.** The clr correlation between taxa ranks
  true required donor→recipient edges at AUROC 0.73 (W_powerlaw); a donor-preserving permutation
  null gives 0.49. But in **W_no_exchange** (same lineages, no leak, no exchange), the same
  statistic scored against W_powerlaw's edges gives **0.73**. The apparent circuit recovery is
  therefore niche similarity: required donors share season and guild with their recipients. A
  ridge "growth" probe gives 0.53–0.57 in both worlds.
- **Hidden state.** Leave-one-year-out linear probes decode the unmeasured metabolite pools:
  - environment alone: R² 0.38–0.42;
  - + uptake-gene RNA: R² 0.44–0.53;
  - + DNA abundances: R² 0.50–0.70 (W_powerlaw).
  The pools are linearly readable from the community state, but readable is not causal (the
  mechinterp caution in manifold-timescales.md §4.1).

### 5.5 What is shown vs speculation
Shown, under these assumptions:
- Plexitrophy by d ≥ 2 emerges readily from leaky public goods plus gene loss. Cycles ≥ 3 do not
  emerge in well-mixed pools.
- Realised flux overstates required partners about five-fold.
- Expression's own timescale is invisible at daily or slower cadence; only its slaved part
  shows.
- Red noise mimics most series.
- Diel aliasing through RNA closure produces spurious ~monthly cycles.
- P3 is conservative and low-powered. The P1 guild proxy can miss a real supply response.
- Correlation-based "circuits" are confounded by niches.

Speculation (not tested here):
- that real communities sit in the same regime;
- that particles or spatial structure would make cycles ≥ 3 common (manifold-timescales.md §4.3
  predicts this);
- that P3 would have more power for exchange through specific (non-substitutable) partners.

### 5.6 Limits and next steps
- One seed per world. A seed sweep (≥ 5 per world) is the first thing to add before quoting any
  number as a rate.
- The perturbation arms and the P3 gains need replicates.
- The genome layer is mean-field.
- There is no spatial structure or size fraction.
- Observation noise levels were chosen, not calibrated.
- The P1 guild proxy and the oracle comparison should be repeated with a lower-noise RNA model,
  to separate measurement noise from proxy failure.


---

## 6. Revisions after cross-family review (Astra 248998526021189632)

Also covers Carbon (Gemini) 249001105845305344 and Astra's follow-up 249001223084490752.
- **Rerun:** `python3 -I run_all.py --revision`. It reruns the 4 worlds' truth, knock-out
  validation, P3, bloomers, mechinterp and P1, and re-emits the observer tables. It does not
  redo timescales, aliasing, closure, blooms or perturbation.
- **Runtime:** 1826 s wall (~30 min) on the shared machine. Log in
  `results/run_log_revision.txt`, meta in `results/run_meta_revision.json`.
- **Code hashes** (`results/provenance_revision.tsv`): run_all `b9f898ed2f96`, model
  `81abc276214f`, truth `2fd224df3188`, observe `a0bda80f36e0`, analysis `429ad430b813`.
- **Same seeds** as the first sweep. The simulations are identical; only the truth and analysis
  code changed, plus the knock-out mechanics below.
- Still **one seed per world**.

### 6.0 What changed in the code
1. **Exhaustive singleton test** (Astra 1). `necessity_at` tests the empty set and *every*
   supplier (every lineage with any release, plus EUK) alone, so d = 0, d = 1 and the claim
   d ≥ 2 are exact under the snapshot model.
   - Larger sets are searched only among the top 30 (pairs) or top 8 (triples, quadruples)
     suppliers, ranked by their share of the metabolites the recipient cannot cover itself.
   - A pair found gives d = 2 exactly. A larger set found gives only an upper bound `d_hi_h1`.
   - No rescue up to size 4 is recorded as **unresolved** (`d_hi_h1 = -1`,
     `d_unresolved_h1 = 1`), never as d = 5.
   - Realised flux is now computed over all suppliers.
2. **Snapshot sensitivity, named as such** (Astra 2). Pools are rescaled from the recorded stock
   by instantaneous release shares, with competitor uptake fixed. Two checks validate this
   against re-simulated 20-d removal trajectories (same forcing; knocked-out leak is now *lost*,
   no longer returned to the donor):
   - (a) single-donor removal;
   - (b) keep-only-one-supplier: remove every lineage's leak except the recipient's and one
     supplier's, and switch off EUK release unless EUK is that supplier.
3. **Cycles per support** (Astra 3). κ is evaluated within each rescuing support S ∪ {X} as "X on
   a simple directed cycle of length ≥ 3" (plexitrophy-computation.md §1.3).
   plexitrophic_h1 = d_lo ≥ 2, or κ = 1 on every support found (P-all over the supports found,
   which is weaker than P-all over all irredundant supports). A singleton support cannot carry a
   ≥3-cycle, so κ never adds plexitrophy when d = 1. The old union-graph `cycle3` is gone.
4. **Functional fan-in vs dilution** (Carbon 1, Astra 249001223084490752). For each d ≥ 2
   recipient:
   - the *deficient* metabolites are those whose pool, cut to the recipient's own release,
     leaves it non-viable;
   - **dilution** = some deficient metabolite cannot be rescued by any single supplier with the
     other pools intact. It is `dilution_witnessed` if 2–4 producers of that one metabolite do
     rescue it, and `dilution_unresolved` if not;
   - **functional fan-in** = every deficient metabolite is rescuable by a single supplier, ≥ 2
     metabolites are deficient, and no single supplier covers them all (so different suppliers
     are needed across metabolites);
   - `other` = the remainder.
   Also logged: deficient set, `supports_unknown` (d ≥ 1 with no witnessed support, so donor
   identities are unknown for edge scoring) and genes_lost.
5. **P3 aligned with protocols/P3-plexitype-recurrence.md where feasible** (Astra 5):
   - rotated masked targets (every member in turn);
   - δ = 0.02;
   - 20 matched control sets (10 for the 2190-event daily design), with Mn as the mean of log
     scores;
   - background PC1 from V \ (C ∪ N_k), using proportions whose depth denominator **excludes
     the masked target** (fixes the leak), fit on train and projected onto test;
   - shrinkage λ = 5 only on added taxon slopes;
   - Bernoulli detection with log depth;
   - cluster bootstrap over 28-d event blocks;
   - verdicts supported / falsified / inconclusive / N.
   It is labelled the **"P3 pilot"**. It is not the frozen P3: no JSDM conditional, and
   discovery is by residual-correlation cliques.
6. **P3 ground truth = the estimand** (Astra 5). The same pipeline runs on noise-free daily
   truth in the same train/test years:
   - target = expected detection probability 1 − exp(−expected reads), scored as an expected
     log score;
   - predictors = true relative biomass with the target masked;
   - true environment.
   This gives the **oracle out-of-year transfer gain** and an oracle verdict. A required
   exchange edge within the group is reported only as mechanism.
7. **P1 relabelled "P1 proxy stress test"** (Astra 4, Carbon 4):
   - 20 stations per world;
   - margin 0.1 nat;
   - the protocol's replication gate is reported (≥ 20 stations and ≥ 2 dominant genera at ≥ 5
     stations; gate failure → inconclusive under the protocol);
   - dominant-genus transfer folds;
   - bootstrap slope at the same ridge penalty (λ 5) as the predictive model;
   - hk ≥ 20 gated **per recipient-event**, **recipient intercepts**, and log chl added
     (eukaryotic closure).
   Still not the protocol: one stratum, no basins, the unresolved pool is not in the dominance
   denominator, and nats are per-observation means per station.
8. **Observer format** (Carbon 5):
   - `per_mag_cds.tsv` uses p2_scorer's `protein_id`/`KO`/`EC` columns;
   - `genomes.tsv` has `domain`, `genome_size_bp` and `contig_lengths`;
   - one Eukaryota bin (EUK0001);
   - synthetic KOs K90001–K90023 with `synthetic_ko_map.tsv`;
   - verified by parsing MAG0001 with `scratch/p2_scorer.parse_per_mag_rows`.
   The README is corrected: storm pulses raise the measured S1, so they are largely visible to
   Me.
   `two_campaigns_treated` observations were not re-emitted, because the perturbation was not
   rerun, and were removed.
9. **Bloomers** (Carbon 5): `genes_lost ≥ 1` and functional fan-in are now reported alongside
   "any d ≥ 1".

### 6.1 Ground truth, revised (`truth_summary.tsv`; viable lineage-times among the top 300, every 14 d)

| world | d ≥ 1 (genome layer) | d ≥ 2 (exact) | of all viable rows: functional fan-in | dilution, witnessed | dilution, unresolved | other | d_hi unresolved | κ = 1 on all supports (certified) | supports unknown / rows with d ≥ 1 | realised / mean distinct donors in witnessed supports |
|---|---|---|---|---|---|---|---|---|---|---|
| W_powerlaw | 0.668 | 0.432 | **0.024** | 0.181 | 0.191 | 0.036 | 0.216 | 0 | 3 577 / 11 039 (32.4%) | 12.8 / 1.43 |
| W_lognormal | 0.854 | 0.706 | **0.029** | 0.339 | 0.298 | 0.039 | 0.355 | 0 | 4 011 / 9 659 (41.5%) | 22.3 / 2.00 |
| W_flat_euk | 0.669 | 0.405 | **0.012** | 0.181 | 0.178 | 0.035 | 0.207 | 0 | 3 610 / 11 695 (30.9%) | 14.0 / 1.36 |
| W_no_exchange | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 / 0 | 0 / 0 |

- **Withdrawn:** "plexitrophy (d ≥ 2) arose in 43–71%". The d ≥ 2 fractions are exact under the
  snapshot model, now that singletons are exhaustive. But **86–90% of d ≥ 2 rows are dilution**
  (6 151, 7 205 and 6 262 rows). Dilution here means: for at least one deficient metabolite, no
  single producer can rescue the recipient even with the other pools left intact.
  - **Mixed deficiencies are common.** Metabolite 0 (B12-like) is the *sole* deficient product
    in only 2 574 / 6 151, 1 762 / 7 205 and 2 381 / 6 262 dilution rows (W_powerlaw,
    W_lognormal, W_flat_euk).
  - The share of dilution attributable *solely* to B12 dilution is not measured.
  - **Functional fan-in** (different suppliers for different metabolites) is **1.2–2.9% of
    viable lineage-times** (399 / 16 534, 333 / 11 309 and 211 / 17 476; verified by Astra
    249011970736832512). Ten-fold rarer than the first headline, it is the only part that matches
  plexitrophy's intended meaning.
- **Cycles: no certified all-support cycle.** All 22 205 d ≥ 2 rows have
  `supports_complete = 0`, because the search beyond singletons covers only the top suppliers.
  truth.py requires completeness for `kappa_all_supports = 1`, so this field cannot be 1 by
  construction and is not evidence that cycles fail.
  - Witnessed cyclic supports exist in **4 W_powerlaw rows and 11 W_lognormal rows** (0 in
    W_flat_euk): `kappa_any_support = 1`, i.e. at least one witnessed support with X on a
    simple ≥3-cycle.
- **Unknown supports.** Some dependent rows (d ≥ 1) enter edge scoring with **no witnessed
  support** (d_hi unresolved), so their donor identities are unknown: 3 577 / 11 039, 4 011 / 9 659
  and 3 610 / 11 695 (W_powerlaw, W_lognormal, W_flat_euk), i.e. **30.9–41.5%**. Required-edge
  truth used in P3 and mechinterp is incomplete for those rows.
- **Realised donors vs mean distinct donors in witnessed supports** (renamed from
  `mean_n_required`). Recipients drew on 12.8–22.3 realised suppliers.
  - The second number (1.36–2.00) averages, over all viable rows, the union of donors across
    witnessed supports. It puts zeros for unknown supports and counts interchangeable donors
    separately.
  - It **does not estimate how many partners a recipient needs**, so no ratio is claimed. The
    earlier "five-fold" figure is withdrawn.
- **Liebig masking** (Carbon 2). genes_lost is bimodal: 0 genes in 44% of rows and all 5 in 29%
  (W_powerlaw). Metabolite 0 is among the deficient products in **90–92% of d ≥ 2 rows**, but it is the sole
  deficient product in only a minority of them (see above).
  - Mechanism: B12 is costly (×3) and algae draw it down. When a lineage's growth is capped by
    metabolite 0, the other synthesis genes carry cost with zero marginal benefit, because
    Liebig's minimum masks them. So the other four are lost too, a domino loss.
  - Dependence that follows: d ≥ 2 given genes_lost = 5 is 0.68 (0.88 lognormal), vs 0.13
    (0.26) given genes_lost = 0.
  - Contribution: lineages with all 5 genes lost make up 46% (W_powerlaw) and 55% (W_lognormal)
    of the d ≥ 2 rows. I did not cross-tabulate fan-in class by genes_lost, so how much of the
    fan-in comes from the domino-loss group is not measured.
- **Snapshot sensitivity vs trajectories** (W_powerlaw, `knockout_validation_summary.tsv`):

  | test | snapshot class | n | median 20-d ratio to control | ratio < 0.5 |
  |---|---|---|---|---|
  | single-donor removal | necessary | 131 | 0.876 | 3% |
  | single-donor removal | not required | 982 | 0.985 | 2% |
  | keep only one supplier | d = 1 (one supplier suffices) | 25 | 0.695 | 28% |
  | keep only one supplier | d ≥ 2 (no supplier suffices) | 43 | 0.585 | 42% |

  The snapshot classes order the outcomes the right way, but separate them **weakly**. Within
  20 days, community-wide feedbacks dominate: other lineages' decline, pools drawn down, and
  immigration. Treat snapshot necessity as a **sensitivity ranking, not a necessity verdict**.

**Hashes of the files these counts come from (sha256, 12 characters):** `results/summary_tables.md` b6ef8c7db32a; `results/truth_summary.tsv` 9d48e044ecf6; `results/truth/W_powerlaw/d_cycle_per_lineage.tsv` 2855b3b43959; `results/truth/W_lognormal/d_cycle_per_lineage.tsv` db001b64b66c; `results/truth/W_flat_euk/d_cycle_per_lineage.tsv` 6dbe3121dfcd. The §6.1 reporting fixes (Astra 249011970736832512) were made from these TSVs, with no rerun.

### 6.2 P3 pilot vs its own estimand (`p3_groups.tsv`, 90 groups; 3 worlds × 3 designs × 2 universes × 5)

| pilot verdict / oracle verdict | n |
|---|---|
| falsified / falsified | 31 |
| inconclusive / falsified | 51 |
| inconclusive / inconclusive | 5 |
| falsified / inconclusive | 2 |
| inconclusive / supported | 1 |

- Pilot: **0 supported of 90**. Oracle (noise-free, same years): **1 supported of 90**. Oracle
  gains were mostly below δ = 0.02. The incremental, out-of-year, compositional predictability
  that P3 estimates is **small in this world even without measurement noise**.
- The pilot's non-support is therefore largely **correct for its estimand**. It is not evidence
  about exchange. The one oracle-supported group (W_lognormal daily_ideal core200) was
  inconclusive in the pilot, a miss.
- "0/86" in §5.2 described the earlier algorithm and is replaced by the above.
- **Mm − Mn does not separate the worlds** under the aligned pilot (Carbon 3 had suggested it
  might). Mean Mm − Mn per world × design ranged −0.004 to +0.015 (W_powerlaw), +0.001 to +0.023
  (W_lognormal) and +0.001 to +0.012 (W_no_exchange). The first sweep's +0.035/+0.022 vs −0.004
  contrast does not survive masking, rotation and averaging over 20 control sets.
- **Reframe** (Carbon 3). The world has 5 well-mixed global pools and no partner-specific
  interactions, so it tests **diffuse guild substitution**, not partner-specific plexitypes. Me
  already conditions on measured resources (which carry the storm signal), chl and a background
  gradient that decodes much of the pool state. Little residual is left for members to explain.
  A partner-specific world (spatial pairing or specific vitamers) is needed to test P3's intended
  claim.
- Abundant core vs all: 32 groups were discovered with identical members in both universes, and
  all 32 got identical verdicts. The core restriction still changes nothing material.

### 6.3 P1 proxy stress test (`p1_song_vs_singer.tsv`; 20 stations per world)

- **The protocol's replication gate fails in both worlds.** 0 and 1 dominant genera reach ≥ 5
  stations: with 400 lineages per station from a 1500-lineage pool, no genus dominates guild DNA
  repeatedly. Under the protocol, every comparison is **inconclusive**. The decisions printed in
  the table ("against | gate_failed") are proxy-stress-test outputs only.
- Guild gain [95%]: +0.002 [−0.003, 0.007] (vitamer off) and +0.004 [−0.003, 0.010] (on). With
  the true pool as an oracle: +0.019 [0.003, 0.034] and +0.002 [−0.008, 0.012]. All are below
  the 0.1 margin.
- Dominant-genus transfer: 10 and 9 genera tested, 0 with a lower bound above the margin.
- **Why:** after per-event gating and recipient intercepts, the observed uptake log-ratio still
  correlates only r = 0.28 / 0.23 with the true uptake expression, and its sd is **7.2–7.6×**
  the true sd. RNA measurement noise dominates. The guild proxy correlates only 0.34 / 0.16 with
  the true usable pool.
- **Narrowed claim:** in this simulator, P1-style guild models cannot detect a supply response
  that exists in the truth (r(pool, true uptake) = −0.81 / −0.66). Two causes: the per-gene RNA
  noise I chose, and a gene-presence guild proxy that tracks supply poorly. The **guild-proxy
  validation point stands**. The noise level is my assumption, so this is a stress test of the
  proxy, not a power statement for Tara.

### 6.4 Bloomers and mechinterp, revised
- **Bloomers outside the core** (W_powerlaw / W_lognormal):
  - any d ≥ 1: 95% / 98%;
  - d ≥ 1 *with* genes_lost ≥ 1: **30% / 50%**;
  - functional fan-in at any evaluated time: **21% / 17%**.
  - Core non-bloomers: 61% / 85% (genes lost) and 50% / 45% (fan-in).
  - Revised reading: by the stricter criteria, bloomers outside the core carry genome-layer
    dependence **less often** than core taxa. The earlier "about the background rate" is
    withdrawn.
- **Mechinterp:** unchanged in substance. Correlation AUROC for required edges is 0.73 / 0.76,
  and 0.74 in the no-exchange control scored against W_powerlaw's edges. Required-edge truth is
  now built from witnessed supports only, so it is incomplete for 32–42% of dependent rows.

### 6.5 Headlines after revision (narrowed)
1. Leaky public goods + gene loss produce many d ≥ 2 recipients in a well-mixed chemostat.
   86–90% are **dilution**: at least one deficient metabolite cannot be rescued by any single
   producer. Mixed deficiencies are common, and B12-only dilution is not separated.
   Functional fan-in is 1.2–2.9% of viable lineage-times. There is no certified all-support
   cycle; 4 + 11 rows have a witnessed cyclic support.
2. Liebig masking by one costly, algae-drawn-down metabolite drives domino loss of all synthesis
   genes.
3. Snapshot necessity ranks removal effects correctly but weakly against 20-d trajectories.
4. The P3 pilot gives no support. Its noise-free estimand is also mostly below δ, so the null is
   largely correct here. This world cannot test partner-specific plexitypes.
5. The P1 proxy stress test fails the protocol gate. Where it runs, the gene-presence guild proxy
   and the RNA noise mask a real supply response.
6. Unchanged from the first sweep (not rerun): the timescale, red-noise, perturbation, bloom,
   aliasing and closure results in §5.1 and §5.3.
