# bq_lattice: Black Queen lattice world for plexitrophy, built for prediction by synthesis

Author: claude-plexi (Opus) subagent, 2026-11-18. Single writer for this directory.
Brief: `../BRIEF.md`, including "Additions after Astra's check". Results: `RESULTS.md`.

**Version 2.** Revised after cross-family review by Astra (Codex, 248988907492261888 and
248989851437154304) and Carbon (Gemini, 248989754267713536); see the last section.
- Current outputs: `sweep_rev/`.
- `sweep_out/` was produced by v1 code. It is superseded; do not cite it.
numpy/scipy and the standard library only. Run everything with `python3 -I`.

## What it is
This is an individual-based 2-D lattice (64 × 64, periodic, at most one cell per site).
- Microbial genomes carry K = 11 biosynthetic functions.
- Every cell needs every product.
- Losing a gene saves its cost (Black Queen).

**Plexitrophy is not coded anywhere.** Nothing in the code says who feeds whom. Dependence comes
only from four sources:
- gene loss;
- leakage and diffusion;
- an O2 incompatibility;
- space.

### Mechanisms (`world.py`)
| Element | Rule |
|---|---|
| Functions | L0–L4 are leaky, compatible public goods. OX is leaky and produces O2. ANX is leaky and O2-sensitive. P0–P3 are private: the product is never exported. |
| Production | Each carried function yields Y = 4 per step. A fraction `leak` goes into that product's diffusing field and 1 − leak is kept, capped at demand θ = 1. A cell takes its deficit from the field at its own site. |
| Growth | Per-step division probability g = rmax (1 − cost)^n_carried × min_j(supply_j/θ) (Liebig). Daughters go to a random empty Moore neighbour. Death probability is μ = 0.05. |
| Gene loss / HGT | Each gene is lost with probability `u_loss` = 1e-3 per division. HGT: with probability 2e-4 per cell-step, a cell gains one function that a random neighbour carries. |
| Incompatibility ι | ANX efficiency = (1 − ι·carries_OX) / (1 + ι·O2/K_O). OX carriers emit O2. All cells respire it. Water relaxes toward a background O2 level. ι ramps from 0 to its target over the first 1,500 steps, so the ancestor (a full prototroph) is viable. ι = 0 means no block. ι = 1 means a genome carrying OX cannot run ANX. |
| Private functions | A non-carrier gets the product only from a static abiotic patch map: full supply on `priv_cov` = 0.95 of sites and none elsewhere (`priv_bg` = 0). Retention therefore never depends on other cells. These values were chosen so private and leaky prevalence match at about 0.5. |
| Medium | Ambient abiotic supply of 0.02 per site-step for every leaky product. |
| Particles | Discs covering a fraction `phi` of the lattice. Diffusion there is D × 0.1, and cells do not mix. Water sites mix cells by pairwise swaps. Light is attenuated on particles: OX product and O2 run at `light_p` = 0.3. Particulate respiration is `resp_p_mult` = 3× faster. Particle cores are therefore anoxic. |
| Eukaryote (optional, `euk_frac`) | A large-cell phototroph with a fixed genome that never evolves. It exports L0, L1 and OX (organic C + O2) at Y_e = 8, emits twice the O2, and **requires B12 (L4)**. It takes up only B12; its other needs are met internally. Its growth is self-shaded by eukaryote neighbours. It takes no part in HGT. |

### Ground truth (`run.py`, `truth/`)
- `lattice_snapshots.npz`: genotype for every site at every sample time, plus the particle mask
  and the private-patch maps.
- `timeseries.tsv`: counts, genotype richness and the carriage of each function every 100
  steps.
- **Realised transfers** (`realised_lineage.tsv`, `realised_edges.tsv`).
  - In an 80-step window before each sample time, every leaky field is split into provenance
    tags: one per producer genotype (the top 12), plus ENV (ambient), PRE (content at the
    window start) and OTHER.
  - Uptake removes the tags in proportion, so each unit a recipient takes is attributed to a
    source. Flux is accumulated over the last 40 steps.
  - Each row is one recipient genotype × scope (all / water / particle). It gives:
    - `d_real`: minimum cover of the partner-supplied needs by sources with ≥10% share;
    - `best_single_share`;
    - `c_real`: the SCC size in the realised graph;
    - the eukaryote share of uptake.
  - **These rows record use, not necessity.** A recipient can use two suppliers when either
    alone would suffice.
- **Necessity** (`necessity_lineage.tsv`, `necessity_supports.tsv`; v2): P-all plexitrophy over
  tested kept sets, at the last two sample times.
  - **Kept sets:** every subset of 1–4 genotypes drawn from the 7 most abundant (≥ 20 cells),
    with no genome-coverage filter.
  - **Each set is run as a counterfactual** for H = 100 steps.
    - **Ghosts:** all other cells become inert ghosts. A ghost stays in place and takes up only
      the products it lacks, so competition for the endpoint metabolites is unchanged. It has no
      O2 production or respiration, no export, no birth and no death.
    - **Tagging:** provenance tags restricted to the kept genotypes record realised supply over
      steps 50–100.
  - **Feasible set:** every member is viable, meaning mean growth probability ≥ death over steps
    50–100.
  - **Verdict for recipient X:**
    - `autonomous` if {X} is feasible;
    - `not_plexi` if some feasible set admits a counterexample wiring. That means d(X) ≤ 1 (a
      minimum cover of X's partner-supplied needs by suppliers with ≥ 10% share) and X on no
      **simple** directed cycle of ≥ 3. X's covers are enumerated; other members get one greedy
      cover each;
    - `plexi` (type fan-in, cyclic or mixed) if no feasible set does;
    - `unresolved_no_feasible_set` if no tested set is feasible;
    - `not_viable_in_situ` if X is not viable in situ.
  - **Needs use active capability.** At ι = 1, an OX carrier's encoded ANX counts as absent
    (`active_needs`).
  - Verdicts are given per scope: all, and X's cells in water or on particles.
  - `necessity_supports.tsv` lists each feasible set with d(X), whether it is a counterexample,
    whether only a cycle rescued it, and whether it contains the eukaryote.
- `abundance_fit.tsv`: the genotype rank-abundance distribution at each sample time. Columns:
  - the power-law MLE α;
  - a truncated-lognormal vs power-law log-likelihood ratio, with its Vuong z;
  - the rank slope.

### Synthetic observations (`obs/`, same shapes as our real tables)
- **Events.** At each of 11 sample times there are 6 stations, each a 17 × 17 window centred
  on a particle site.
  - Each station is one `event_id`. It yields a paired **water (W)** fraction (cells on water
    sites) and a **particle (P)** fraction (cells on particle sites).
  - DNA and RNA share the `sample_id` and `event_id` (`samples.tsv`).
- **Reads.** `abundance_{emergent,lognormal,powerlaw}.tsv` give MAG × sample read counts at
  identical depth (20,000), multinomial in cells × genome length.
  - The lognormal and power-law files add an imposed background of 3,000 non-interacting taxa
    taking 50% of reads. They dilute detection by the same amount but with different tail
    shapes. A row for the unbinned lattice reads and a row for the background are included.
- **Catalogue.** A genotype becomes a MAG (perfect dereplication) if it has ≥10 reads in ≥2
  samples.
- **Gene model and CDS table** (`per_mag_cds.tsv`; v2 columns `mag_id, contig_id, protein_id,
  feature_index, KO`, with protein_id = `{mag}_c{contig}_{index:05d}` as `scratch/p2_scorer.py`
  expects; `contig_lengths.tsv` gives about 1 kb per CDS):
  - Every prokaryotic genome has 400 CDSs: 200 housekeeping KOs, the first 100 of which are
    completeness markers; 145 CDSs with a **blank KO**; and 5-KO module operons.
  - The ancestral gene order is fixed. A lost module's CDSs are absent.
- **Incompleteness.** Each MAG is fragmented into contigs (breaks with probability 1/20 per
  gene). Whole contigs are then dropped until the completeness target is reached. Targets are
  50% U(0.9, 1) and 50% U(0.55, 0.9). Because modules are operons, contig loss removes whole
  modules.
- `mag_quality.tsv`: CheckM2-like completeness, taken from the retained markers, plus contig
  and CDS counts.
- **Eukaryote bin.** Its CDSs are MetaEuk-like: 35% of genes are called, one per contig.
  Domain is Eukaryota and completeness is NA.
- `unbinned_ko.tsv`: KOs of non-catalogue genotypes present in a sample.
- `expression.tsv`: RNA-like module activity per sample × MAG. ANX activity reflects O2 and
  OX co-carriage.
- `modules.tsv`: KO sets per module, with leaky/private/incompatible labels "from chemistry".
  These labels are fixed by construction, not inferred from co-occurrence.
- `truth/mag_map.tsv` and `truth/per_mag_cds_complete.tsv`: the MAG → genotype map and the
  uncut genomes, used for scoring and for the oracle variants only.

### Inference and scoring (`infer.py`)
**Module calls** use the P2 thresholds: present if ≥80% of KOs, absent if ≤25%, unresolved
otherwise. Eukaryote bins (uncalibrated detectability) are present if ≥80% and unresolved
otherwise, never absent.

**Recipients** are prokaryotic MAGs with completeness ≥ 90 in the incomplete variants. Their
needs J are the absent leaky modules. Any unresolved module makes the recipient unresolved.

**d estimate** (v2 order):
1. A definite single donor covering J gives d = 1.
2. Otherwise, if any single ambiguous source covers J, the call is unresolved
   (`ambig_single_supplier`). Ambiguous sources are a MAG's present and unresolved modules, or
   the unbinned pool. This is P2's permissive upper bound.
3. Otherwise a definite minimum cover gives d ≥ 2.
4. Otherwise the call is `ambig_donor` or `no_donor`.

MAG × time is called plexitrophic if d ≥ 2 in at least half of its resolved samples. Scores are
reported resolved-only and with [low, high] bounds: unresolved counted as not plexitrophic, or as
plexitrophic.

**Variants:**

| Variant | Genomes | Eukaryotes | Donor pool | Abundance regime |
|---|---|---|---|---|
| `complete_all` | complete | as donors and as a recipient | all detected | emergent |
| `complete_prok` | complete | excluded | all detected | emergent |
| `incomplete_all` | incomplete, HQ recipients | fragmentary donors | all detected | emergent |
| `incomplete_prok` (P2-like primary) | incomplete, HQ recipients | excluded | all detected | emergent |
| `_core` | incomplete, HQ recipients | excluded | core MAGs only (max relative abundance ≥1% at that time) | emergent |
| `_lognormal`, `_powerlaw` | incomplete, HQ recipients | excluded | all detected | background-diluted |

**Fan-in null:** Prop. 1 of `discussion/deep/fan-in-null-theory.md`, with independent functions
and definite evidence. It is compared on matched units:
- **unconditional:** all eligible recipient-samples, observed [d ≥ 2] vs P(d ≥ 2);
- **conditional:** units with a cover, observed vs P(d ≥ 2)/P(covered).
- Reported per N-bin.

**Co-carriage excess** (v2):
- Per-donor-genome co-carriage. For each pair, O is the number of donors carrying both and
  E = n_j n_k / N, the exact fixed-count column-permutation expectation.
- Σ_C = log(ΣO + .5) − log(ΣE + .5) per class. I = Σ_leaky − Σ_private, with an event bootstrap.
- Σ_inc is the same statistic for OX × ANX.
- It is also computed on the true genotypes.

**Scoring:**
- Precision and recall against the necessity truth, split into fan-in truth (every support has
  d ≥ 2) and cyclic-route truth (a support with d = 1 but c ≥ 3).
- Precision and recall against the realised truth.
- d bias: incomplete minus complete for the same MAG × sample.
- False-absence rate of modules.

**Other readouts:**
- P4 (v2): the same recipient MAG × event resolved in both fractions, weight 1.
  - Observer [d ≥ 2]_P − [d ≥ 2]_W.
  - A separate donor-count-standardised observer, which subsamples the larger donor set to the
    smaller.
  - Truth on the same units from habitat-specific P-all verdicts, any and fan-in-only.
- Eukaryote share of true uptake, and the share of minimal supports that include the eukaryote.
- Rare donors: the fraction of true donors in plexitrophic recipients' minimal supports that are
  in the catalogue, detected at that time, or in the abundant core, under each abundance regime.

## How to run
```
python3 -I run.py OUTDIR SEED '{"iota":1.0,"euk_frac":0}'   # one run (~60-110 s on one core)
python3 -I infer.py OUTDIR                                  # inference + scoring (~3 s)
python3 -I sweep.py sweep_rev 8                             # revised reduced sweep, 8 workers
python3 -I sweep.py sweep_rev agg                           # re-aggregate tables only
python3 -I make_tables.py sweep_rev                         # markdown tables
```
- Parameters: `world.DEFAULTS` and `run.RUN_DEFAULTS`. Any key can be overridden in the JSON.
- Seeds: `sweep_rev/seeds.tsv` lists every config, seed and parameter JSON. Each seed is
  `20261118 + 10·(crc32(config) mod 1e5) + rep`. A run spawns three child streams from
  `SeedSequence(seed)`: world, counterfactuals and observation sampling.
- Per-run `runtime_s` is in `sweep_rev/results/runs.tsv`.
- Runtime (v2): 290–400 s per run on 8 concurrent workers, including about 200 counterfactual
  runs. The 32-run sweep took 1,454 s wall.
- Sweep outputs: `sweep_rev/runs/<config>/s<seed>/{meta.json, truth/, obs/, inference.json,
  inference_calls.tsv}` and `sweep_rev/results/{runs.tsv, config_means.tsv, signatures.tsv,
  tables.md}`.
  The CDS and expression tables are gzipped after inference.

## Known limits (read before citing anything)
- One cell per site, and one genome per genotype. "Lineage" means gene-content genotype, which
  is what a perfectly dereplicated MAG represents.
- There are no regained pathways except by HGT, no phylogeny and no catabolic exchange, so only
  d_anab is measured.
- Necessity tests kept sets of ≤ 4 genotypes from the top 7. Inside a set, the wiring is
  realised flux (≥ 10% share). Other members get one greedy cover each, a heuristic restriction of
  P-all's "every wiring". The verdict depends on the stated endpoint and horizon.
- Particles here are low-diffusion discs with attenuated light and faster respiration. They are
  anoxic by those assumptions.
- Inference covers d only. P2's cyclic disjunct (c ≥ 3) is not inferred from genomes here.
- These are simulations. They show what can happen under these assumptions, not what happens
  in lakes.

## Revisions after cross-family review (Astra 248988907492261888)
These revisions also cover Astra 248989851437154304 and Carbon 248989754267713536. Details, and
the effect of each on the claims, are in RESULTS.md §6.
- **world.py**
  - Ghosts are redox-inert: no O2 production or respiration.
  - The eukaryote takes up only B12.
  - Particles have light attenuation (`light_p`) and particulate respiration (`resp_p_mult`).
  - `priv_bg` = 0 and `priv_cov` = 0.95.
  - `start_tags(only=...)`.
- **run.py**
  - P-all necessity over tested kept sets, with all-member feasibility, realised-flux wiring
    inside each counterfactual and simple-cycle κ.
  - Active-capability needs.
  - p2_scorer CDS schema and contig lengths.
- **infer.py**
  - Permissive single-supplier check before d ≥ 2.
  - Eukaryote modules are never called absent.
  - Matched-denominator fan-in null.
  - Per-donor co-carriage I.
  - Bounded scores.
  - Same-unit P4 with a separate donor standardisation.
  - Capacity vs survival.
- **sweep.py / make_tables.py:** reduced revised sweep and new tables.

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
| sweep_rev/seeds.tsv | c0d828246c1a |

### v3 (Carbon 248999088913563648 / 249001216432324608, Astra 248999382242213888 / 249001252071325696)
- **run.py:** wiring and realised truth also count encoded products with realised non-self uptake
  ≥ 0.05 per cell-step (`U_MIN`). This covers O2-inhibited ANX, light-limited OX and leak 0.8.
  Truth rows now carry `status_v2_needs`, `in_edges_v2_needs` and `in_edges_v3_deficit` for
  comparison.
- **world.py:** optional planted positive control `stack_discount`, off by default.
- **rerun_c1.py:** targeted rerun into `sweep_c1/` (`python3 -I rerun_c1.py sweep_c1 10`; summary
  only: `... sweep_c1 0`). 10 runs, 501 s wall.
- Results are in RESULTS.md §6, v3 addendum. v3 hashes: world.py 184f53f37424, run.py fa4a7f83000c,
  rerun_c1.py 7e31927e52d2.
