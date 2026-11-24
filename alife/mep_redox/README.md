# mep_redox: a thermodynamically explicit redox column with assembled metabolic types

Author: Claude (Opus 5.5) sub-agent for claude-plexi, 2026-11-18. Single writer for this directory.
Status: **v3, revised after two cross-family reviews** (Astra/Codex 248992826624753664, 248994794025959424;
Carbon/Gemini 248994691424894976). See "Revisions after cross-family review" at the end.
The v1 numbers previously in this file are superseded; every number below is from the v3 sweep.
A simulation shows what is possible under its assumptions. It does not show that anything is
true of real water columns.

## What I read
- `alife/BRIEF.md`, including "Additions after Astra's check" (realised vs required partners,
  rich KO records, shared event IDs).
- `AGENTS.md`; the top of `discussion/plexi-concept-review.md` (the summary, the standing
  perspective, the P5 statement); `discussion/deep/manifold-timescales.md` §2.2 and §4.3 (the
  gauge-symmetry proposition and Astra's constraint); and
  `discussion/deep/architecture-identifiability.md` §5 (S5, the lagged-succession spurious chain).
- **Vallino & Huber 2018**, *Front. Environ. Sci.* 6:100, doi:10.3389/fenvs.2018.00100.
  - Version read: the publisher's full-text HTML, downloaded 2026-11-18 into an empty directory
    (`/tmp/mep_dl_vallino/page.html`, sha256 `e7f861579216`, the same file the
    manifold-timescales note used) and converted with a stdlib script.
  - I read the abstract, the introduction, the methods (the catalyst/ε/Ω construction, the
    photon free energy, MEP over Δt) and the results on SIO and LIO.
  - Displayed equations did not render, and I did not read the Supplementary Material.
  - This model borrows ideas, not code or equations: catalysts with allocations; a rate law
    with a thermodynamic factor; entropy production from photon capture; finite-horizon
    optimisation. My rate law is the brief's form, not Vallino's adaptive-Monod law.
- No file under `data/` was opened. No chat posts were made and no cluster jobs were run.

## Model (model.py)
**Domain.** A 1-D column of 30 m in 1 m layers.
- Mixed layer 0–5 m with Kz = 20 m² d⁻¹; Kz = 0.3 m² d⁻¹ below (meromictic-like).
- Top boundary (Dirichlet): O2 280, NO3 5, SO4 250 µM, and B12 0.002 nM (an ambient abiotic
  trickle). Other species have zero flux.
- Lateral exchange 0.002 d⁻¹ towards SO4 = 250 µM. POC sinks at 1 m d⁻¹ and hydrolyses at
  0.1 d⁻¹. POC reaching the bottom returns as DOC (a sediment source).
- Light: I(z) = I0 exp(−0.3 z − self-shading by algae). It is optionally sinusoidally forced
  (`forced`: amplitude 0.8, period 30 d).

**Chemistry.** 12 tracked species: DOC (CH2O), POC, O2, NO3, NO2, NH4, SO4, H2S, H2, CH4, DMSP
and B12. CO2, N2, H⁺ (pH 7.5) and H2O are held fixed.
- Every reaction uses ΔG = ΔG° + RT ln Q, with local activities (µM × 10⁻⁶) and standard
  ΔfG° values at 298 K applied at 288 K (an approximation). Heterotrophic steps release NH4 at
  Redfield N:C.
- 14 catabolic half-reaction pairs: AER (CH2O+O2), DMP (DMSP+O2), NAR (CH2O+NO3→NO2),
  NIR (CH2O+NO2→N2), SRB (CH2O+SO4), FER (CH2O→H2), MET (H2→CH4), HSR (H2+SO4), AMO (NH4→NO2),
  NXR (NO2→NO3), AMX (anammox), SOX (H2S+O2), SNR (H2S+NO3→NO2) and MOX (CH4+O2).

**Kinetics (per unit enzyme).** q_r = k0 · Π_s S/(K_s+S) · F_T · I_O2, in mol e⁻ per mol C per day.
- F_T = max(0, 1 − exp((ΔG_r/n_e + G_min)/(χRT))), with G_min = 1 kJ per mol e⁻ and χ = 1.
- k0 = 10 is the same for every reaction. Thermodynamics and Monod constants are the only
  differences between reactions.
- Anaerobic steps are inhibited by O2 (K = 1 µM).
- Flux J_{j,r} = (q_r/n_e) · e_{j,r} · B_j · κ_j. Here κ_j ~ lognormal(0, 0.15) is a per-type,
  per-seed kinetic jitter that breaks exact ties.

**Metabolic types (genomes).** 22 types, including single steps and stacks:
- denitrification: NAR, NIR, DEN = NAR+NIR, AERNAR, AERDEN;
- nitrification: AOB, NOB, CMX = AMO+NXR (comammox);
- H2 and S: FER, MET, FERMET, SRB, SRBH, SRBX, SOB, SNR, SOBX;
- others: AMX, MOB, AER, ROS (AER+DMP, a Roseobacter-like DMSP degrader);
- ALG.

**Growth.** dB/dt = Y_E Σ_r J(−ΔG) − maintenance − mortality (linear + quadratic) + immigration
(10⁻⁴) + mixing.
- Y_E = 10⁻³ mol C per kJ. Growth is energy capture.
- The `electron_yield` control replaces this with yield per electron, which ignores ΔG.
- **Chemotroph biomass is catalyst only (v3).** It grows on captured energy without drawing C,
  N or electrons from the pools, and its necromass is *not* returned to POC. Only algal
  necromass (half of algal mortality) becomes POC; the rest is exported.
- Catabolic ΔG covers the C/N/S redox couple only. Organic-N mineralisation (+RN NH4 per DOC
  C) is treated as non-energetic, since organic N is valued at the NH4 potential.
- **MEP is not imposed anywhere in the population model.**

**Algae (phototrophic eukaryotes).** Added at grid's request, relayed by the coordinator.
- ALG fixes C at a rate set by light × B12 limitation (K = 0.01 nM) × DIN limitation
  (K = 0.5 µM; v3).
- v3 element balance: algal biomass and exudate carry Redfield organic N, taken up as NH4
  first, then NO3 (NO3 assimilation releases 2 O2 per N). DMSP (C H2 O0.4 S0.2 per C) is made
  from CO2 + 0.2 SO4 and releases 1.6 O2 per C.
- 50% of gross fixation becomes algal biomass. The rest splits as 20% exudate DOC, 5% DMSP, and
  25% respired by the alga itself, which uses its own O2.
- The alga also contributes POC through mortality. Algae are the only source of organic C and
  of O2 below the surface.
- **B12** is made only by ROS. It is a leaky public good: release is proportional to ROS growth,
  ROS pays a 2% growth tax, and B12 decays at 0.05 d⁻¹. Algae consume B12 as a quota.
- Algal entropy production = gross fixation × (10 photons × 132 kJ mol⁻¹ − the local chemical
  potential stored per C) / T. Light absorbed by water is not counted.

**Switchable kinetic symmetry breakers** (the team's critical constraint). Splitting a net
reaction leaves ΣΔG unchanged, so any advantage of splitting must come through realised flux.
| switch | parameter | scenario |
|---|---|---|
| enzyme crowding: effective enzyme = Ω · n^(1−α) | α = 1.3 (crowding penalises stacks) / 0.8 (stacks favoured) | `crowd`, `stackadv` |
| allocation regulation (Ω ∝ power^β over the genome's steps) | β = 8 (regulated, default) / 0 (constitutive) | `constitutive`: co-expressed steps waste enzyme across the O2 micro-gradient |
| intermediate channelling (leakage) inside stacked genomes: internal consumer flux shares the consumer enzyme's capacity and obeys its co-substrate Monod terms, F_T (incl. G_min) and O2 inhibition at bulk activities (v3) | φ = 0.8 (default 0 = all intermediates leak to the pool) | `channel` |
| free-nitrite toxicity on growth | K_I = 2 µM (+φ = 0.8) | `toxic_channel` |
| dispersal (cell mixing × 30; algae excluded) | cell_mix | `dispersal` (mass-effect, non-MEP assembly) |

**Gauge check (manifold-timescales §4.3).** `FixedEnzyme` drives the same chemistry with the
enzyme field E_r(z) implied by a population state.
- Without channelling it reproduces the emergent σ to within 0.01–0.08% (T1, "gauge rel. diff").
  The genome packaging of enzyme is invisible to dissipation, as the proposition says.
- With channelling (φ > 0) the difference is 0.2%. Channelled flux is not representable as a
  pool-level enzyme field.

## Dissipation yardsticks (mep.py)
All use the same medium, with algal biomass, necromass flux and B12 held at the emergent values.
None represents channelling.
- **Best-found reference (steady):** a bounded Powell search (250 evaluations plus a 60-evaluation
  restart from the best point) over 28 controls, a log-amount and a depth shift for each
  reaction's enzyme profile.
  - The total catalyst is fixed at the emergent amount (budget error < 10⁻¹⁵), and the
    objective is steady-state chemotrophic σ.
  - **This is not an MEP optimum.** It is a best-found value under those constraints.
    Powell never reports convergence (it always hits maxfev); the restart's gain and the
    number of parameters at a bound are reported. Restart gains are < 0.1% as scenario means;
    the largest single gain is 0.148% (electron_yield seed 0, of 24 runs).
  - "Fixed-field σ within 0.5%" reports only that the emergent enzyme field, replayed without
    genomes, gives σ within 0.5% of the emergent run. One scalar cannot show that channelled
    flux is representable (the reference has no channelling), so the best-found ratio is
    descriptive only.
- **Greedy (H → 0):** at each depth and instant, the emergent total catalyst is reallocated
  to the reactions with the highest power (softmax, β = 20). This is Vallino's
  short-horizon limit.
- **Random:** 6 Dirichlet allocations per depth at the same per-depth catalyst total.
- **Horizon test:** from a sterile start (300 d of chemistry with algae but no chemotroph
  enzyme), static enzyme fields are optimised for the mean σ over [0, H], with H = 1, 10, 100
  and 1000 d, and each is then evaluated at every H. This is compared with the greedy policy
  and with the population assembling from seed biomass.

## Ground truth (analysis.py; truth/<run_id>/)
- **Realised transfers.** Each species is attributed to its producers by solving the
  steady-state transport–reaction operator once per tagged source (diffusion, sinking and
  first-order uptake).
  - Organic C is traced through POC and hydrolysis to its ultimate producer (algae or bacterial
    necromass).
  - Channelled intra-genome flows are logged separately.
  - Abiotic sources are the nodes ATM (top boundary) and LAT (lateral inflow).
  - Files: `results/realised_transfers.tsv`, `truth/<run>/realised_transfers.tsv`,
    `truth/<run>/flux_by_depth.tsv`.
- **Chemical-potential budget (v3):** `analysis.free_energy_budget`, species-resolved, closes to
  within 0.035% of T·σ_chem for all 27/27 emergent rows. Over all 123 rows (including the
  fixed-field configurations), 19 exceed 0.2%, with a maximum of 0.90% (electron_yield seed 1,
  random1) (T8). Its terms are algal capture, boundary inputs and outputs, net
  transfer to the CO2/N2/H2O/H+ reservoirs, storage, and dissipation by chemotrophs,
  diffusion, sinking, hydrolysis, mortality conversion and sediment return.
- **Whole-column C, N, S and electron budgets (v3):** `Model.column_budget` with algae on and
  open boundaries; tested in `tests/test_conservation.py`.
- **Realised partner counts.** For each recipient: shares per input, d_real (living non-self
  donors with ≥ 10% of some input), Hill numbers, and c_real (SCC size of the thresholded
  realised graph). These count **use**, not need.
- **Necessity (v3).** Counterfactual knock-outs under the same medium and horizon (600 d).
  - The endpoint is recipient column biomass ≥ 10% of an **unperturbed control** integrated
    over the same horizon from the same state.
  - Singles: every living donor is removed alone. OR-groups: all living donors of one input of
    one recipient are removed together.
  - Failures are split into direct (the donor has a realised transfer to the recipient),
    indirect (upstream) and solver failures.
  - d_lb_direct = the minimum hitting set of the failing OR-groups and failing direct donors.
    d_lb_upstream also includes indirect failures.
  - Both are **lower bounds conditional on monotone rescue**, not demonstrated viable supports.
    Only d_lb_direct speaks to fan-in; d_lb_upstream also counts serial chains.
- **Division of labour.** For six handoffs: the stacked share of the consumer step, the
  flux-weighted depth separation of the steps, and same-genome self-supply.

## Synthetic observations (observe.py; synthetic/Dnnn/)
- **Blinded (v3):** each dataset is `synthetic/Dnnn/`, with run/event/sample ids
  `Dnnn_e1_z01`, …. The scenario names are only in `truth/unblinding_key.tsv`, and the
  permutation uses OS entropy.
- Each run has 3 events (casts at days 1480, 1490, 1500) × 11 depths. DNA, RNA and chemistry
  share `event_id` and `sample_id`.
- `chemistry.tsv`: 5% lognormal noise and detection limits.
- `abundance.tsv`: reads per MAG and sample, plus an `unbinned` row; multinomial with 10⁶ reads.
- `per_mag_cds.tsv` (headers as `scratch/p2_scorer.py` reads them): mag_id, contig_id,
  protein_id (= `<contig_id>_<index>`), feature_index, start, end, strand and KO.
  - About 45% of housekeeping CDSs have a blank KO.
  - Pathway steps are contiguous operons.
  - Realistic KO ambiguities are kept: narGH = nxrAB, amoABC = pmoABC, and sqr is shared by
    SOX and SNR.
- `mag_completeness.tsv`: per-run contig loss gives 60–100% completeness, reported with ±3%
  noise. Also domain (ALG = Eukaryota), genome_size_bp, contig_ids and contig_lengths.
- `rna_counts.tsv`: per-CDS counts; pathway CDSs ∝ step flux, housekeeping ∝ biomass, algal
  psbA/rbcL ∝ fixation; lognormal noise; 10⁶ reads.
- MAG ids are anonymised. The key is `truth/<run>/mag_key.tsv`.
- Fixed yardstick configurations (`__mep_opt`, `__greedy`, `__random1`, `__random2`) emit
  the same files. Each reaction's enzyme is carried by its specialist type.

## Observer-side inference (infer.py)
- P6-like test: per sample, Spearman between RNA marker-class shares and power predicted from
  the observed chemistry, against a Monod-only prediction and the DNA share.
- Unexploited power: the share of predicted power in classes absent from DNA.
- KO-rule exchange graph (NO2, H2, H2S, CH4, B12). Since v3, edges are drawn only between MAGs
  present (≥ 0.1%) in the **same sample**, then aggregated; scored against the realised truth
  edges.
- Mean redox path length on the SCC-condensed DAG, inferred vs true (the P5 statistic).
- Co-occurrence graph (Spearman ≥ 0.8). Its false edges between depth-adjacent MAGs are the
  spatial form of the lagged-succession artefact.

## How to run
```
cd /data/dev/plexityper/alife/mep_redox
python3 -I tests/test_conservation.py  # stoichiometry, channel and whole-column budgets
python3 -I run_all.py --workers 4      # full sweep: 2026 s wall on 4 cores (v3)
python3 -I summarise.py                # writes results/summary.md
```
- Seeds 0, 1 and 2 set the initial biomass, the per-type κ (rng 1000+seed) and the observation
  noise (seed+5 / seed+9).
- The genome catalogue uses seed 12345.
- Single population run: 1500 d in 6–141 s on one core under machine load
  (`results/population_runs.tsv`).
- Sweep:
  - 81 population runs: 6 mechanisms × 3 pools × 3 seeds, plus 3 controls × 3 seeds;
  - 27 analyses;
  - 7 knock-out sets, each with an unperturbed control plus 9–22 removal runs;
  - 1 horizon test.
- numpy 2.3.5 and scipy 1.17.0. Nothing was installed.

## Results (v3; full tables: results/summary.md, sha256 5e43d9eb9f38)
Values are mean ± sd over 3 seeds unless stated. σ is in J K⁻¹ m⁻² d⁻¹ and is averaged over
the last 300 d.

**R1. At steady state, σ per unit primary production is consistent with similar net conversion
across all arrangements; it is not an identity** (T1, T8).
- Emergent / best-found reference = 0.962–0.968. The reference is not an MEP optimum
  (convergence notes in T1).
- Random allocations (same per-depth catalyst) reach 97–101% of emergent σ, and greedy reaches
  98–99%.
- σ_chem/NPP is 1.01–1.06 J K⁻¹ per mmol C in every scenario.
- The species-resolved budget (T8) closes within 0.035% for 27/27 emergent rows. 19/123 rows
  overall exceed 0.2% (maximum 0.90%, electron_yield seed 1, random1). Algal capture is ≈ 509–511 kJ per mmol C.
  T·σ_chem takes ≈ 290–309, the net transfer to the CO2/N2/H2O/H+ reservoirs ≈ 248–262,
  boundary in − out ≈ +35 to +53, and transport/hydrolysis/mortality/sediment dissipation 4–7.
- What differs between arrangements is mainly the boundary term: exports of partially oxidised
  products (NH4, H2S, NO2) and SO4 exchange, plus transport dissipation. Fixed NPP alone does not
  fix these (Astra).
- So MEP has slack only (a) over finite transient horizons (R7) and (b) where catalysts feed
  back on primary resource capture (R3). Steady-state σ is a weak discriminator here.

**R2. Algae dominate dissipation and supply.**
- Algal photon dissipation is 79–80% of total σ.
- Algae supply 88–100% of the intake of aerobic heterotrophs, 31–97% for sulfate reducers, and
  3–41% for lithotrophs.
- ALG is a failing single removal for 71 of 74 non-algal recipient instances (T4).

**R3. Growth selection does not maximise total dissipation, because of a leaky public good**
(`pool_contrast.tsv`).
- Total σ with the full pool is lower than with at least one restricted pool in 17/18
  mechanism × seed cases, and lower than both in 15/18. Neutral means: 56 full, 88 split-only,
  72 stack-only.
- The mechanism, checked on the v3 runs: generalist heterotrophs crowd the sole B12 supplier
  ROS. Removing them raises ROS's share (24% → 34%, neutral seed 0), B12 (0.004 → 0.010 nM)
  and algae (302 → 435 mmol C m⁻²).

**R4. Division of labour responds to the kinetic switches, but the dissipation consequences are
≤ 2% and only partly survive v3** (T2, T3).
- Stacked share of nitrification (comammox): crowd 0.12, neutral 0.81, dispersal 0.87,
  constitutive 0.92, stackadv 0.94, channel and toxic_channel 0.96.
- Under constitutive expression the stack AERDEN disappears (95 → 0 mmol C m⁻²). Stacked
  denitrification falls from 0.99 to 0.25: the O2 micro-gradient mechanism.
- σ per NPP, split-only minus stack-only (%), paired by seed and expressed relative to neutral:
  - crowd +0.35 to +1.21;
  - constitutive +1.2 to +2.1;
  - stackadv −0.5 to −0.9;
  - channel −0.05 to −0.21;
  - toxic_channel −0.1 to −0.5.
- Signs match the mechanisms (crowding and constitutive expression favour splitting; stacking
  advantage favours stacks). With the corrected, capacity-bounded channelling, though, the
  channel effect on dissipation is within ~0.2% of neutral, even though channelling still
  raises the comammox share.

**R5. Redox handoffs become spatially separated chains.**
- Sulfate reduction → sulfide oxidation are separated by 11–15 m, and methanogenesis →
  methanotrophy by 7–13 m. Nitrification and denitrification steps stay within 0.5–2 m of each
  other.
- True redox path lengths are 1.25–1.52, not ≥ 3.

**R6. Realised partners overstate counterfactual necessity** (T4b).
- Realised d_real ≥ 2: 6–14 recipients per run.
- d_lb_direct ≥ 2: 2–6. d_lb_upstream ≥ 2: 2–12.
- Solver failures: 0.
- The d_lb_direct = 2 cases are catabolic fan-in at the interface: sulfide oxidisers need O2
  from ALG plus H2S from ROS (via DMSP) or sulfate reducers; nitrifiers need O2 from ALG plus
  NH4 from heterotrophs. They are lower bounds conditional on monotone rescue.
- ROS is often required only indirectly, through B12 → ALG (counted in d_lb_upstream only).
- ALG itself has no failing single removal: ambient B12 keeps it above the endpoint.

**R7. Horizons** (T6, neutral seed 0, sterile start; `horizon_convergence.tsv`).
- The H = 1 d static best does best at 1 d (1878), but is the worst static field at 1000 d
  (78.2, against 94.8 for the H = 1000 d best: −17.5%). That is Vallino's "ditch".
- The greedy policy is best at 10 d.
- The population assembly is worst at 1 d (11) and best at 100 and 1000 d (253, 97). That
  comparison is not catalyst-matched, so it does not show that assembly beats MEP.
- The searches hit maxfev; restart gains are 0–0.9 units.

**R8. Observable signatures** (T5).
- RNA~power correlation is 0.76–0.80 under energy-yield selection vs 0.69 under electron-yield
  selection. The best-found reference gives 0.69 and random allocations 0.73. That is a weak
  discriminant of selection currency, with overlapping ranges, and none of closeness to an
  optimum.
- Dispersal-swamped assembly shows high unexploited power (0.19 vs 0.02–0.07).
- **No chemistry + gene signature separates near-optimal from far-from-optimal arrangements.**
- The KO-rule graph, even with per-sample co-occurrence (v3), has precision 0.07–0.22 and
  recall 0.74–0.94. Its inferred path length stays ≈ 1.45, because most MAGs exceed 0.1% in
  most samples (mixing plus immigration floors). Carbon's pooling diagnosis was therefore not
  the main cause here.
- Co-occurrence precision is 0.04–0.26, with 8–16 false depth-adjacent edges per run: the
  spatial twin of the lagged-succession artefact.

## Relation to P5 and the lagged-succession artefact
- The true redox graphs are short (path 1.25–1.52). The long handoffs (SRB → SOB, MET → MOB)
  are separated by 7–15 m, so real partners rarely co-occur, while depth-adjacent non-partners
  do.
- A depth-resolved P5 test therefore needs genome-based edges, a transport-aware null, and a
  presence threshold well above immigration and mixing floors.

## Speculation (not shown by the simulation)
- Where eukaryote–bacterium cofactor links control primary production, "MEP-organised" may be
  indistinguishable from "public-good-limited". MEP tests in nature should condition on
  primary production and boundary exports.

## Limitations
- Parameters are set by hand. Only 3 seeds.
- ΔfG° values are at 298 K. DMSP's ΔfG is a placeholder (its stoichiometry is balanced).
- Chemotroph biomass is catalyst only.
- The reference search is local and parametric. Algae and B12 are not controls, and channelling
  is not representable.
- Knock-outs use one horizon and endpoint.
- No size fractions are emitted.
- The v3 sweep took 34 min wall, a little over the 30-min target.

## Revisions after cross-family review (Astra 248992826624753664)
Also covers Astra 248994794025959424 and Carbon 248994691424894976.

1. **Conservation (Astra 1).** v1/v2 added channelled intermediates back to the pool on top of
   NU·(J + J_int), which created NO2 (+1 per AMO→internal NXR turnover).
   - Removed. Tests: Astra's counterexample, closed-box "full tendency = balanced reactions",
     and N/C/electron balance over 50 random channelled states.
2. **Channel bounds (Astra 2).** Internal consumer flux now shares the consumer enzyme's
   capacity and obeys co-substrate Monod terms, F_T (with G_min) and O2 inhibition.
   - Tested: external + internal ≤ capacity; no internal flux where F_T = 0; internal ≤ own
     production.
3. **Reference wording (Astra 3).** "MEP-opt" is renamed "best-found reference".
   - Convergence (Powell success, restart gain, parameters at bound), budget error and the
     common-feasible-set check are reported.
   - The reference lacks channelling and holds B12 and algae fixed.
4. **Knock-outs (Astra 4).** Compared with an unperturbed control at the same horizon.
   - Direct, indirect and solver failures are reported separately.
   - The hitting sets are called conditional lower bounds and split into d_lb_direct and
     d_lb_upstream (Carbon 2).
5. **Blinding.** Synthetic ids are anonymous; the key is in `truth/unblinding_key.tsv`.
6. **Whole-column conservation (Carbon 1).**
   - Algae take up Redfield N (NH4, then NO3 with O2 release).
   - DMSP consumes SO4, with a balanced formula and O2 release. DMP was rebalanced
     (1.2 O2, 4.8 e⁻ per C).
   - DMSP carries no organic N, and chemotroph biomass is catalyst only.
   - Whole-column C, N, S and electron budgets with algae on are tested to ≤ 10⁻¹⁵ relative.
7. **Identity wording (Carbon 3, Astra's qualification).** The T1 agreement is described as
   "consistent with similar net conversion", and the species-resolved chemical-potential budget
   is added (T8). It closes within 0.035% for 27/27 emergent rows; 19/123 rows overall exceed
   0.2% (maximum 0.90%).
8. **Per-sample rule graph (Carbon 4)**, and **p2_scorer headers, domain and contig lengths
   (Carbon 5).**

9. **Reporting fixes (Astra 249009418569957376).** The feasible-set column is renamed
   "fixed-field σ within 0.5%". Budget closure is given with denominators (27/27 emergent rows
   within 0.035%; 19/123 rows overall above 0.2%, maximum 0.90%). The restart-gain bound is
   stated as a scenario mean, with the 0.148% maximum. No new sweep was run.

Claims withdrawn or changed from v1:
- "Channelling gives a −2.8% split-vs-stack dissipation effect": **withdrawn**. In v3 it is
  within 0.2% of neutral.
- "Emergent sits 3% below an MEP optimum": now "0.962–0.968 of a best-found reference under
  stated constraints".
- v1 absolute σ values (≈ 100) are superseded: v3 neutral σ_chem is 56. Algal N limitation and
  catalyst-only chemotrophs lower production.
- "ROS required by many recipients" is now classed as an indirect / upstream requirement.

Claims that survived:
- σ/NPP near-invariance at steady state.
- The B12 public-good loss of dissipation: 17/18, previously 18/18.
- Algal dominance.
- The sign of the crowding, constitutive and stacking-advantage packaging effects (≤ 2%).
- Spatially separated redox chains.
- Realised ≫ required partners.
- The horizon "ditch".
- The co-occurrence depth artefact.
- The absence of an MEP-specific observable signature.

## Files and hashes (v3; sha256, first 12 hex)
| file | sha256 |
|---|---|
| model.py | 2715d0fb202a |
| analysis.py | 2c6ff545c94b |
| mep.py | e7863947a588 |
| observe.py | 336970ac8353 |
| infer.py | 4162e0ccf2b4 |
| run_all.py | 0d87b62a8407 |
| summarise.py | 4dfdb5c134c1 |
| tests/test_conservation.py | 47ea3251e7a1 |
| results/sigma.tsv | 410e13777dbb |
| results/pool_contrast.tsv | 89bace44d8d2 |
| results/division_of_labour.tsv | c4828a43649e |
| results/necessity_knockouts.tsv | df5a97203ace |
| results/partners_realised.tsv | 038dff886462 |
| results/realised_transfers.tsv | fe0795ddeff1 |
| results/mep_yardsticks.tsv | cd789d294025 |
| results/free_energy_budget.tsv | 4c7926ed4358 |
| results/horizon.tsv | afe0d7d0ba9e |
| results/horizon_convergence.tsv | af077bb9e3a9 |
| results/signatures_inference.tsv | 7d8f2a3d31f1 |
| results/biomass.tsv | 87ad83b03a4d |
| results/summary.md | 5e43d9eb9f38 |

`results/runs/*.npz` holds the final states. `synthetic/` is 98 MB and `truth/` is 1.7 MB.
