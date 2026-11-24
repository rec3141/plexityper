# Plexi concept review: group discussion and literature analysis

Started 2026-11-18 15:52 at grid's request (relay 248913233255514112). Single writer:
claude-plexi (Opus). Others write their own section notes and propose edits in #plexityper.

**Ground rule from grid:** the `plexi-` prefix stays; every other term and claim in
`plexichat.txt` is up for revision.

## Summary for grid (closing, 17:22; read this first)
About 90 minutes of three-family discussion (15:52–17:20). Everything below is agreed by Opus,
Astra and Carbon unless marked. Nothing is frozen and no compute has been run.

**Definition (one sentence):**
> *Plexitrophy* is obligate dependence of an organism, under a stated medium and growth
> endpoint, on metabolites or electrons transferred between coexisting living partners
> (including across a host–symbiont interface), such that no single partner suffices. Either at
> least two distinct partners must supply simultaneously (d_total ≥ 2), or the organism belongs
> to an irreducible cycle of at least three (κ = 1). Inputs obtained solely by ingesting
> another organism are excluded.
>
> *Anabolic plexitrophy* (d_anab ≥ 2 or c_anab ≥ 3) is the subtype for which genomes can
> identify candidates. Pathway calls do not by themselves show obligate exchange.

Serial chains (d = 1) are not plexitrophy. Counts are d_anab, d_catab and d_total; a partner
with both roles counts once.

**Vocabulary for Paper 1:**
- **plexitrophy** and **plexiome** lead. The plexiome is an evidence-graded AND/OR dependency
  hypergraph: OR = substitutable donors (the "singers"), AND = jointly required roles.
- **plexitype** is a descriptive, possibly overlapping recurring *subset* of taxa that must
  beat environment, gradients and matched nonmembers out of sample. It carries no mechanism
  claim.
- **plexon** applies only with an independent physical boundary (an isolated particle, or a
  tracked parcel).
- **plexibiont** is glossary only, until an interdependence test passes.
- Holobiont keeps its historical meaning, and host association is a separate axis.
- Dictyobiont and catenabiont are dropped; topology uses plain words.

**Organising question:** when an assemblage recurs, does its composition recur, its exchange
roles, or both? (Song vs singer.)

**Motivating cases (checked by abstract and/or full text):**
- Men et al. 2014: Dehalococcoides grows only in triculture, with H2 from Desulfovibrio and
  corrinoid from Pelosinus, plus external DMB. That gives d_anab = 1, d_catab = 1,
  d_total = 2: mixed plexitrophy.
- Helliwell et al. 2016: algae remodel cyanobacterial pseudocobalamin only when DMB is also
  supplied. This is a *candidate* two-input dependence; the number of living donors in situ is
  untested.
- Wienhausen et al. 2024: two bacteria co-assemble B12, forming a c = 2 pair.

**Protocols (drafts in protocols/, all cross-family reviewed, margins "proposed" until frozen):**
| | Test | Data | Status |
|---|---|---|---|
| P1 | Exchange role recurs across stations with different dominant suppliers (cobalamin) | Tara Oceans W6 (external) | Astra's draft (812d46c677e2); both cross-family reviews applied (station holdout with ≥ 4 training stations per genus; separate genus holdout; exact size interval; marker-normalised DNA) |
| P2 | Anabolic plexitrophy exceeds fragmentation and stacking nulls; enriched where stacking is blocked | Dereplicated prokaryotic MAGs; Stage 0 = Ebb & Flow (43 HQ, pipeline check only) | Carbon (828c0c3e0dde); reviews closed |
| P3 | Plexitype members predict a masked member in a new year better than environment and outside taxa | Lake Winnipeg 2024 → 2025 | Opus (33de07d6ebf6); reviews closed |
| P4 | The same recipient MAG × endpoint differs between particle and free-living fractions | Amundsen HS22-003 (pilot, n = 1 event) | Astra (0ed5e0034020); reviews closed |
| P5 | Exchange chemistry predicts topology (public goods nested, redox chains, amino acids reciprocal) | Same MAGs as P2 | Statement only (in the synthesis); reviewed by both |

**P6 (new, drafted 2026-11-19 at grid's request):**
protocols/P6-thermodynamic-allocation.md tests whether thermodynamic driving force
(capacity·F_K·F_T·(−ΔG)) predicts *expressed* catabolic allocation beyond capacity × kinetics.
It is the testable form of grid's free-energy (MEP / Vallino) perspective. It is gated on paired
RNA plus in-situ chemistry existing. The protocol and its reaction table (6c26362681b8) are
cross-family reviewed and closed. Remaining: grid's go, the outcome-blind data-gate check, a
frozen δ, and the count-generating calibration.

**Paper 1 outline:** discussion/paper1-outline.md (Carbon, b25922221a7d).
- Working title: "Plexitrophy and the Plexiome: Testing Irreducible Multi-Partner Nutritional
  Dependence in Aquatic Microbial Assemblages".
- A 150-word conditional abstract, outcome-contingent results for P1–P4, and a figure list.
- It is written so that null results are publishable.

**Problems found in grid's pipelines (verified in code by Opus):**
- **Trap 8:** `mag_analysis/bin/kegg_module_completeness.py` (55 modules).
  - Only 3 amino acid modules.
  - Thiamine is HMP-only.
  - M00120 mixes in aerobic cobG.
  - Any single subunit earns full step credit, so partial complexes are overcalled.
- **Trap 9:** `mag_analysis/modules/metabolism.nf` builds per_mag from every raw binner plus
  the consensus bins, which duplicates genomes 4–6×.
- **Trap 3:** the illumina_rna GFF choice causes KO locus-tag drift. ko_remap.tsv.gz is a
  workaround.

**Decisions only you can make:**
1. Fix kegg_module_completeness.py, metabolism.nf (MAP_TO_BINS) and illumina_rna's GFF choice
   in danaSeq? Or keep the workarounds in the protocols?
2. P1: is Tara W6 acceptable as the external cohort, or do you have, or have collaborators
   with, an independent paired DNA+RNA aquatic series?
3. P4: confirm whether the HS22-003 fractions (0.9 / 2.5 / 0.8 L on 0.22 / 3 / 20 µm) were one
   serial filtration.
4. Save the paywalled PDFs by hand into discussion/papers/: Venrick 1971 and Wienhausen 2024.
   Heal 2017 is open access. Chrome automation couldn't save them to disk.
5. Freeze order and compute.
   - Verified as staged locally and table-scale: **only P2 Stage 0** (Ebb & Flow per_mag and
     depths).
   - P2's primary cohorts (Freshwater 191 HQ, metaMDBG 160 HQ, Chesterfield 10 HQ) have only
     checkm2 reports locally. per_mag and breadth tables must be staged from Grex, or
     metabolism/mapping run.
   - P3 needs its feature table (coverage or ASV) located; not yet verified.
   - P4 (HS22-003) has 6 raw FASTQs on /matika and needs mapping unless matching mag_analysis
     outputs exist.
   - Fresh MAG building costs 200–500 core-h (Carbon 248933478447890432).
   - **P1 also needs competitive DNA/RNA mapping of the Tara pairs**, unless compatible
     precomputed mappings exist. W6 gives only pair metadata, and eligibility and mapping
     availability are unverified (Astra 248933325003472896).
6. Manuscript owner (manuscript/ is untouched).

## Standing perspective from grid (2026-11-18, keep open)
Microbial communities can be seen as **learned representations of an underlying manifold,
shaped by free-energy minimisation**. Several layers each learn at their own timescale and
occupy their own frequency domain:
- population
- genome
- pathway
- gene
- expression
- nucleotide

The current state of the system is roughly the sum of these. In effect, we are doing
**mechanistic interpretability on microbial ecosystems**.

**grid on energy:** the energy part is hand-wavy so far, but Joe Vallino's work is "on the
right track":
- [doi:10.3389/fenvs.2018.00100](https://www.frontiersin.org/journals/environmental-science/articles/10.3389/fenvs.2018.00100/full);
  maximum entropy production over time and space horizons; to be verified.
- So read "free energy" thermodynamically, as MEP and dissipation-maximising distributed
  metabolism, rather than only as a variational objective.

**grid's follow-up:** "so rigid rules and deterministic-like descriptions are not going to fit
well." So the crisp definitions above (d ≥ 2, κ = 1, discrete topology classes, plexitype
membership) should be read as **thresholded summaries of graded, probabilistic,
layer-dependent quantities**, not as the objects themselves. Preferred forms:
- distributions over supports or partners, and expected partner counts;
- soft feature loadings, not class membership;
- dependence measured per timescale.

Pre-registered tests stay, but their estimands become continuous.

Graded estimands per protocol (Astra 248970428366110720):
- **P1:** supplier shares, unresolved share, adjusted btuB slope, and both predictive gains,
  each with an interval. The ">50% dominance" stratum remains a design choice.
- **P2 and P4:** propagate uncertainty in pathway presence, donor detection and ambient supply
  into lower and upper bounds on dependence and on paired fraction differences.
- **P3:** held-out predictive gain against soft loadings learned in training only.
- **P5:** edge uncertainty, and class stability across blocks.

Caveats:
- A posterior "P(no single supplier suffices)" needs an explicit model over biological states.
  That is not the frequency of multi-supplier routes, and weighting wirings uniformly lets
  duplicate routes change the answer.
- Bootstrap frequencies describe stability. Posterior topology weights need a stated likelihood
  and prior.
- Detection and precision gates remain documented measurement decisions.

How this frames the work: plexitype, plexiome and plexitrophy are candidate *features* or
*circuits* in that representation. Tests should be read for what they say about which layer or
timescale carries the dependence. (Exploration: discussion/deep/manifold-timescales.md.)

## grid on habitat power laws (2026-11-18)
- Lake Winnipeg may hold ~20k ASVs, but rank abundance follows a power law. Most ASVs are so
  rare they can probably be ignored. **Only a couple of hundred taxa are abundant enough to
  matter, and most of those yield good MAGs.**
- Caveat: exponential growth means the gap from undetectable to abundant can be short
  (blooms, conditionally rare taxa).
- But time series show that most seasonal action is concentrated in the abundant lineages.

Implications to fold into the protocols (Opus's proposal):
- **P3:** restrict eligible taxa to the abundant core, about 200 lineages, defined in training
  data by a frozen abundance or prevalence rule. This cuts pseudo-cliques and the multiplicity
  burden.
  - For core taxa, detection is near-universal, so the primary likelihood should model
    abundance (log-ratio or count with offset), not Bernoulli detection.
  - The power simulation must be re-run under this design.
- **P2:** focal recipients and donor pools are drawn from the well-assembled core MAGs. Rare
  lineages are kept as an explicit "rare or unbinned donor" pool feeding the upper bound, and
  `unresolved_no_donor` is read as a probable rare or bloom donor.
- **Simulations:** use power-law rank-abundance, with rare-to-bloom transitions, so inference
  is tested under realistic dominance and detection.
- **Astra 248976328636153856 (refinements):**
  - **P3 abundant-core definition:**
    - Freeze the ranking statistic and taxonomic unit on 2024 data only, weighting
      station-events equally.
    - Report read share, detection frequency and seasonal occupancy separately: a top-ranked
      bloomer can be absent most of the year.
  - **P3 abundance model:**
    - Use a count model with library depth and compositional uncertainty.
    - Mask the target before any predictor-side normalisation, to avoid leakage.
    - Recalibrate δ and the discovery step at this scale: 0.02 nats for presence has no
      automatic abundance meaning.
  - **P2:** keep individually resolved rare MAGs as possible single-supplier witnesses; pool
    uncertain evidence only in the permissive bound.
  - **Simulations:** run both lognormal and power-law tails, with the same depth and bloom
    schedules.
- **Timescales:** the seasonal signal sits mainly in the abundant layer; bloom onsets are where
  fast dynamics cross detection limits.

## grid on eukaryotes (2026-11-18)
"We already know that euks are going to be a big source of metabolites, so can't ignore their
contributions." Phytoplankton and other protists are major donors (DOM, amino acids, DMSP,
vitamin precursors, phycosphere exudates) as well as recipients (B12-requiring algae).

Implications (Opus's proposal):
- **P2:** eukaryotes enter the donor pool even though recipients stay prokaryotic for now.
  - The freshwater run already has eukaryotic outputs at
    grex:~/scratch/project_nanopore/mag_freshwater/eukaryotic/: tiara, whokaryote, MetaEuk
    proteins and GFF, and MarFERReT contig and protein taxonomy.
  - Eukaryotic donor capacity comes from MetaEuk gene calls annotated to KO, with abundance
    from coverage and, where available, 18S/ASV or RNA.
  - Without eukaryotic donors, d ≥ 2 calls are biased upward: missing donors push toward false
    plexitrophy, as the computation note shows.
- **P4:** coarse fractions (3–20 µm) are eukaryote-rich, so particle and phycosphere enrichment
  must account for eukaryotic donors.
- **P1:** B12-requiring algae are major cobalamin recipients and competitors; include them in
  the recipient and demand terms.
- **Simulations:** phototrophic eukaryotes as primary organic-carbon and metabolite sources.
- **Astra 248976870271795200 (refinements):**
  - **P1:** keep the reviewed prokaryotic btuB endpoint and add a *separate* algal-demand
    analysis.
    - W6 pairs are prokaryote-enriched fractions, so larger-fraction algal samples need their
      own station, time, depth and filter crosswalk.
    - Unsampled algal demand is "unknown".
  - **P2 and P4:** count a resolved eukaryotic donor only when the supported pathway steps
    belong to the same organism in that event. Pooling MetaEuk hits across taxa can
    manufacture a single supplier.
    - Unresolved multi-taxon hits go in the permissive bound.
    - Distinguish biosynthetic capacity from release.
    - 18S is occurrence evidence unless copy number and coverage are calibrated.
    - Report how many recipient calls change when eukaryotic evidence is included.
- **Paper 1:** Helliwell (alga as remodelling recipient) and the eukaryote–bacteria B12
  literature become central rather than side cases.

## Workstreams (one owner each)
| Concept | Owner | Working notes |
|---|---|---|
| plexitype, plexon (recurrence, lifetime, compositional class) | codex-plexi (Astra) | [notes/plexitype-recurrence.md](../notes/plexitype-recurrence.md) |
| plexitrophy, plexiome (and architecture) | agy-plexi (Carbon) | [notes/plexitrophy-and-plexiome.md](../notes/plexitrophy-and-plexiome.md) |
| plexibiont superclass, architecture subclasses | claude-plexi (Opus) | (in vision file) |

grid 248913654321692672: duplication is fine; each family writes a full vision and Opus merges.
Visions: [Opus](../notes/vision-claude-plexi.md) · [Astra](../notes/vision-codex-plexi.md) · [Carbon](../notes/vision-agy-plexi.md)

Each section ends with: prior terms that already cover it, what survives, a proposed
revision, and testable hypotheses (question, data, pipeline, falsifier, rough compute).

## Datasets and pipelines in scope
- Local: `/matika/seqs/{metagenome,rnaseq,assembly,amplicon}`; microscape.app.
- Pipelines: github.com/rec3141 (danaSeq `nanopore_assembly`, `mag_analysis`, ...).
- HPC: Grex (via concentration), Fir. Discussion phase only: **no jobs** without grid's go.

## Findings log
Times in parentheses are approximate. The relay message ids are authoritative (check with the relay DB).
- (15:52) Astra 248913372233777152: treat *lifetime* as estimated per occurrence, with a
  separate observed date range per type; nondetections cannot date extinction (L4: one deep
  date recovered 95.4% of OTUs pooled over 72 shallow dates). Seasonal compositional
  recurrence is already established (Fuhrman et al. 2006 PNAS, doi:10.1073/pnas.0602399103;
  Gilbert et al. 2012 ISME J 6:1089). Proposed added prediction: knowing some members improves
  held-out prediction of the others beyond environment, season and depth.
  Open question: must a plexitype require interaction evidence, or can it be descriptive?

- (15:55) Opus 248913641759752192 proposed hub = asymmetric conditional occurrence.
  **Withdrawn** (248914093553401856) after Carbon 248914020543152128: asymmetry measures
  nestedness, which environmental specialist gradients and chains also produce.
- (15:58) Carbon 248914020543152128: architecture is called from conditional-independence
  structure given environment E. Star: spokes ⫫ | hub, E. Chain: ends ⫫ | middle, E. Mesh: no
  single node separates the core. Donor-specific vs ambient supply: test
  I(X; specific donors | total community pathway supply G_m, E) > 0; vitamer remodelling (bluB/cobT
  pseudocobalamin→cobalamin) and siderophore/TonB pairs as specific edges; k > 2 cycles predicted
  enriched in particle/aggregate size fractions.
- (15:58) Astra 248913854914281472: conditional tests need habitat overlap (positivity),
  since nested occurrence arises from a narrower salinity tolerance. Plexon boundary
  counterexample: rank abundances on individual chitin particles were poorly predicted by the bulk
  average (PMC9335236), so bulk series can support recurring composition but leave the instance
  unresolved.
- (16:00) Astra 248914056547057664: the RNA sign test is weak (storage, lags, feedback repression
  and abundance can reverse covariance during real exchange); diffusivity alone gives no
  dilution length without release rate and boundaries. **Fork and chain with the same middle node
  are Markov-equivalent**, so orientation needs non-covariance evidence. Goldford et al. (PMC6405290)
  Fig. 3: four isolates used one another's byproducts yet each grew alone on glucose, so dense
  facilitation can coexist with zero obligate partners.
  Opus note: in the hypergraph framing, genome annotation orients edges (donor encodes,
  recipient lacks); covariance supplies only the skeleton.

- (16:08) Opus 248914484261208064 reviewed Astra's vision (d012b71910c8): Astra's H2 (an
  exchange role predicts recipient activity despite donor replacement) is proposed as the first
  pilot. Holobiont keeps its historical meaning; host association and graph shape are separate axes.
- (16:09) Carbon 248914507304714240:
  - **Markov fix.** At n = 3, star and chain are equivalent, but a collider (sink hub) is
    identifiable. At n ≥ 4, undirected conditional independence separates star (one separator
    into ≥ 3 singletons), chain (tree diameter ≥ 3) and mesh (chordless cycle).
  - **AND/OR dependency hypergraph.** OR-edges are substitutable donors; AND-hyperedges are
    conjunctive complementation (plexitrophy). Each recipient gets a boolean support formula,
    with CheckM2-weighted gap probabilities kept as "unresolved" literals.
  - **Integer standardised:** plexitrophy means a minimal AND-clause of d ≥ 2 donors for one
    recipient (≥ 3 organisms), or an irreducible cycle of length ≥ 3.
  - **Pipeline.** mag_analysis --run_metabolism --run_mge gives KofamScan KOs, KEGG module
    completeness, MinPath, KEGG-Decoder, antiSMASH (siderophores) and MacSyFinder
    (secretion). illumina_rna --from_mag gives gene_counts and bin_tpm, so within-MAG
    normalised expression is log(gene TPM / MAG housekeeping TPM).
- (16:10) Astra 248914574040285184, H2 eligibility from metadata only:
  - Source: genice-sample-metadata cmo_rna_metadata.tsv (a2aabfeb1d76) joined to
    genohub_metagenome_metadata.tsv (ca36581041af).
  - 130 RNA rows → 112 exact DNA matches → 104 after QC exclusions → 97 dated, spanning
    39 dates.
  - These are serial observations from 2 CMO campaigns, **not 39 independent donor
    replacements**. First gate: vessel and inoculum independence, and donor turnover.

- (16:13) Astra 248914851…/ref 3373d50ac0f7 (vision now 85ed65cf38c5):
  - CMO has one pool per treatment per winter (rnaseq-exploration/audit/AUDIT.md C1), so H2 on
    CMO is only a descriptive substitution pilot until independent events exist.
  - **Pipeline caveat:** rnaseq-exploration/ISSUE_ko_misalignment.md documents a Bakta
    basic/extra locus-tag mismatch, and illumina_rna/main.nf still lists the extra GFF first.
    Check the coordinate remap against the count matrix before any KO analysis.
  - H1 candidate: Namao/Lake Winnipeg multi-year station revisits; season and depth coverage
    still need reconciling.
  - Ruling on my question: same members with different edges = same plexitype, distinct
    plexiomes.

- (16:15) Astra 248914907328069632, review of Opus's vision:
  - A degree-preserving null fixes the hub statistic, so it can't test hub concentration. Use
    it only for wiring motifs conditional on degree. For hubs, preserve detectability, resource
    prevalence, richness and edge count, and let the degrees vary. **Accepted.**
  - A serial three-species loop can have one immediate donor per recipient yet need the whole
    loop to persist. Report **local simultaneous-donor count d** separately from the
    **smallest self-maintaining consortium m**. **Accepted.**
  - Lagrangian example: Ottesen et al. PNAS doi:10.1073/pnas.1222099110, Fig. 1.

- (16:18) Astra 248915045740101632, Lake Winnipeg H1 screen:
  - 336/1906 DNA manifest rows are Namao; 327 pass the filters.
  - 22 of 130 station labels occur in all four years (2022–25), but none has a calendar month
    sampled in all four.
  - For 2024–25, 17 stations share ≥ 1 month (21 station–month cells before depth and protocol
    reconciliation).
  - So the data support **two-year conditional spatial replication, not annual recurrence**.
- (16:20) grid 248914999900553216: approved downloading PDFs/XMLs. Venrick 1971 is reachable
  through the UM proxy (full access), but Chrome automation could not save it to disk
  (blob download and localhost POST both blocked). It needs a manual save; see paywalled.tsv.
- (16:21) Astra 248915154431295488: whole-sample classes or recurring subsets? Masked-member
  prediction tests a subset, not discrete community states; claiming states needs mixture vs
  smooth models on the same taxa. Prior art: Wilkinson et al. (JSDM conditional prediction,
  doi:10.1111/2041-210X.13518). Novelty = transfer across independent aquatic events + a linked
  mechanism test. **Opus: subset**, see below.

- (16:24) Astra 248915594451533824 (vision 212d7ff23f0b): plexitype falsifier = membership
  specificity. If matched nonmembers predict the masked taxa equally well within the
  prespecified margin, the membership boundary is unsupported.

- (16:27) Carbon 248915719529873408 found a counterexample to "d ≥ 2 or m ≥ 3". The serial
  chain D1 → D2 → X has d(X) = 1 but m(X) = 3, so every two-hop chain would count. **Fix
  adopted:** report the triple (d, m, c), where c = size of X's strongly connected component. *(Later superseded: c → κ, a simple cycle of length ≥ 3 through X; see the synthesis.)*
  - Fan-in plexitrophy: d ≥ 2.
  - Cyclic plexitrophy: c ≥ 3.
  - Serial cascade (d = 1, c = 1, m ≥ 3) is **not** plexitrophy; it belongs to chain
    architecture.
- (16:27) Carbon: in CMO (diesel-dosed mesocosm), a **song-vs-singer substitution** shows up
  descriptively.
  - Opus checked the source, rnaseq-exploration/audit/AUDIT.md rows F1a and F2a–b.
    - F1a: in the diesel pool, Cycloclasticus carries 63–95% of the tmo, todC, nahAc and
      meta-cleavage transcripts.
    - F2a–b: the bloom is metabat_9342 in 2024-25 and lorbin_184 (sp001566825) in 2025-26.
  - Carbon also reports that canonical CYP153 shifted between binned and unbinned
    Alcanivorax/Sphingorhabdus.
  - This is a degradation role, not an exchange edge, with one pool per arm (C1).
  - It can still serve as a worked example of "role recurs, singer changes".
  - Coordinate-remapped KOs exist: rnaseq-exploration/data/ko_remap.tsv.gz (verified present,
    90 MB), addressing the Bakta locus-tag issue.

- (16:29) Astra 248915827713556480, refining Carbon's triple:
  - **c(X) is defined directly** as the size of X's strongly connected component in each
    feasible minimal support. Carbon's "min m(Y)" alternative is wrong in this case: an obligate
    cycle A, B, C, where A also needs self-sufficient D, gives c(A) = 3 and m(A) = 4 but
    min m(Y) = 1.
  - Alternative donor choices are kept separate before computing c.
  - Correction: AUDIT.md F1a pools binned and unbinned Cycloclasticus. lorbin_184 alone carries
    63% of the second-winter diesel-pool tmo transcripts. Evidence = a recurring **genus-level
    role** with a change of dominant MAG; any exchange edge still needs its own evidence.

- (16:39) Carbon 248925922304376832: literature pass on cobamide remodelling (details in
  notes/plexitrophy-and-plexiome.md §4, 0619824656b2). Opus checked the PubMed abstracts for
  Helliwell 2016 and Wienhausen 2024.
  - **Helliwell et al. 2016 Curr Biol** (doi:10.1016/j.cub.2016.02.041; Carbon read the full
    text): **evidence for a *possible* two-donor (fan-in) mechanism.** Per Astra's correction
    248926045801463808, d_X is left unassigned: Fig. 5 shows that two chemical inputs together
    rescue algal growth, but the number and identity of living suppliers are unresolved.
    - Synechococcus (and cyanobacteria broadly) make only pseudocobalamin. Pseudocobalamin is
      orders of magnitude less bioavailable to B12-dependent microalgae.
    - Some algae can use it if DMB is supplied, by remodelling it themselves.
    - In the field, a remodelling alga *could* need a pseudocobalamin source **and** a DMB
      source. Whether these are two living donors (d = 2) is untested.
    - Caveat (Opus): this was shown in culture with DMB supplied, not with a field DMB donor;
      and the recipient is a **eukaryote**, which our prokaryote MAG pipeline does not bin.
  - **Wienhausen et al. 2024 Nature** (doi:10.1038/s41586-024-07396-y; abstract read by Opus
    and Carbon). Two B12 auxotrophs jointly make B12: Colwellia releases α-ribazole, which
    Roseovarius uses to build B12, so they form a c = 2 pair. A recipient fed by that pair is a
    serial chain (d = 1). In the tropical Atlantic, predicted salvagers outnumber B12 producers.
  - **Shelton et al. 2019 ISME J** (full text, per Carbon): KO-scoring refinement. Use
    cbiL/cobI + cbiF/cobM + cbiC/cobH for ring synthesis (97.6% sensitivity / 99.9%
    specificity), not the chelatases; add cbiZ (K08260) for remodelling.
  - **Wienhausen et al. 2022 ISME J:** shipboard mesocosms with paired metagenomes and
    metatranscriptomes. Five ENA projects: PRJEB43936, 43939, 43941, 43944, 43946. A
    perturbation candidate for P1, pending same-bottle and same-time-point DNA/RNA matching.
    Per Astra 248926130518016000, RNA was sequenced in only 2 of 3 provinces (M1, M2), and
    triplicate bottles replicate treatments, not sites. **P1 plan: Tara W6 for geographic
    replication first; Wienhausen 2022 as a treatment-response check.** B12 repressed
    btuB without changing biosynthesis or metH mRNA, which supports the warning that the
    sign of RNA covariance is not diagnostic.
  - Implication: Paper 1 can use the B12 pseudocobalamin/DMB case as the motivating example of
    a *candidate* fan-in plexitrophy, which our field test (P2) then quantifies. P2 should add a eukaryotic-recipient arm (or flag it as out of
    scope).

- (16:47) Carbon 248926971509522432: P2 (828c0c3e0dde) restricted to prokaryotes, with
  eukaryotic remodellers reserved for an extension. Prokaryotic counterpart found: **Men et al.
  2014 AEM** (doi:10.1128/AEM.03477-13; abstract checked by Opus).
  - Dehalococcoides mccartyi 195 grew sustainably only in **triculture**, not in either
    coculture. DvH supplied H2 and Pelosinus PfR7 supplied phenolic corrinoids, which Dhc195
    remodelled to cobalamin (cbiZ induced).
  - DMB was supplied externally.
  - **Definitional point raised by this case (Opus):** H2 is an electron donor, not a
    nutrient. Counting only nutritional inputs gives d = 1 (corrinoid) plus external DMB;
    counting energy substrates too gives d = 2.
  - **The plexitrophy definition must state whether energy and electron-donor inputs count.**
    Opus proposes reporting both: d_nutr and d_all.
  - Astra 248927038870044672, independently: Men et al. supports a two-partner requirement for
    growth (H2 electron donor + corrinoid) but only one corrinoid-supplying partner. It
    supports the broad ecological claim and leaves the narrow two-donor cobamide claim open.
    The two contributions should be labelled separately in any figure. **Two families agree on
    separating energy from nutritional inputs; Carbon pending.**
  - Astra 248927119698477056 endorses d_nutr as primary (plexichat.txt line 135, "obligate
    nutritional dependence").
    - For each count, state the medium and the endpoint, and count distinct required **living
      suppliers**. Classify each input by its demonstrated role.
    - A substrate supplying both biomass and energy can count for d_nutr, and its supplier
      counts once.
    - **Men triculture worked example:** d_nutr = 1 (PfR7 corrinoids) and d_all = 2 (plus DvH
      H2), with external DMB shown as a medium input.
  - Carbon 248927611782610944, wording fix: don't call H2 "not a nutrient", because electron
    donors are catabolic nutrients (Brock; syntrophy, hydrogenotrophy). Rename to **d_anab /
    d_catab / d_total**. Opus accepts.

- (16:56) Carbon 248929042493259776: Heal et al. 2017 (full text) DOI fixed; paywalled row
  resolved.
  - Predicted remodelling capacity (no corrin-ring synthesis; bluB + cobT and salvage) occurs in
    73 of 3,408 phylogenetically representative reference genomes, including marine Methylophaga
    and Marinobacter. Its prevalence and activity in North Pacific communities are not measured
    by this count (Astra 248929136080764928).
  - **Nitrosopumilus (Thaumarchaeota) makes cobalamin by pairing the O2-independent ring
    pathway with O2-dependent bluB in one genome, and lacks bacterial cobT.** P_BLUB is now
    scored without cobT in Archaea.
  - Opus note: this is a natural **single-genome stacking of an O2-independent and an
    O2-dependent step**. It is a counterexample to treating route-level oxygen labels as
    incompatibility, and supports P2's choice to pre-register pairs by actual enzyme O2
    chemistry and to count single-genome bypasses.

### Converging (as of 16:56)
0. **A plexitype is a recurring *subset* of taxa** within samples. Subsets can overlap and sit
   on continuous gradients; it is not a whole-sample state. (Opus 16:21; Carbon agreed 248917314325889024; all three.)
1. **Plexitype is descriptive** (Astra, Carbon; Opus conceded): a compositional class that
   predicts held-out members beyond environment, season and depth, and beats a
   continuous-gradient model (the enterotype lesson). The mechanism claim lives in
   **plexibiont** and **plexitrophy**.
2. **Architecture from conditional-independence structure**, per edge set, called only when
   stable; "unresolved" is a valid outcome. Direction comes from genomes or thermodynamics, not
   covariance.
3. **Recurrence and mechanism are recorded separately**, each with graded evidence.
5. **Plexiome = AND/OR dependency hypergraph** with graded evidence, alternate and
   unresolved edges kept, and direction from genomes. Plexitrophy is reported as (d, m, κ): d = the minimal AND-clause of simultaneous donors; m = the smallest self-maintaining consortium; κ = 1 if X lies on a simple supply cycle of length ≥ 3, evaluated in every feasible wiring (the earlier |SCC| ≥ 3 rule is SUPERSEDED: see discussion/deep/plexitrophy-computation.md A3). Fan-in plexitrophy if d ≥ 2; cyclic if κ = 1. A serial cascade is not plexitrophy.
4. **Facilitation ≠ obligacy**: the Goldford control means exchange evidence alone does not
   establish plexitrophy; obligacy needs the missing-function test (minimal donor set ≥ 2).

- (16:03) Astra 248914143092326400: **plexon = an independently delimited physical or
  hydrographic occurrence**, e.g. an individually isolated particle, or a Lagrangian parcel
  tracked by hydrography and sampled against its surroundings, plus evidence the candidate occurs
  there. A pooled >3 µm filter is not a plexon (it holds thousands of particles). Carbon's
  size-fraction prediction is an enrichment test at a different scale.
- (16:03) grid 248914167570284544: for paywalled papers, list them (discussion/paywalled.tsv) or
  read them through Claude's Chrome session via the UM library.

### Open
- ~~Keep or drop holobiont/dictyobiont/catenabiont?~~ Resolved: see the merged synthesis
  (holobiont historical; the other two are dropped from Paper 1).

## Merged synthesis (draft 1, 16:09; Opus merging the three visions)
Sources: [Opus](../notes/vision-claude-plexi.md) (abstracts) · [Astra](../notes/vision-codex-plexi.md)
212d7ff23f0b · [Carbon](../notes/vision-agy-plexi.md) dc7f12d4a152 (deep dive 3ab35e89cee5) (full texts of Morris 2012, Frioux
2018, Douglas & Werren 2016, Goldford 2018, Stephens 2024). Every row below is agreed by all three
families unless marked.

### The organising question (Astra)
When an aquatic assemblage reappears, what recurs: the same populations, the same exchange roles,
both, or neither? Doolittle's "song vs singer" is the prior framing. Our contribution is a
measurement scheme that lets composition-recurrence and exchange-recurrence come apart empirically.

### Four orthogonal axes
| Axis | Term | Definition (operational) |
|---|---|---|
| Composition | **plexitype** | A reproducible, possibly overlapping **subset of taxa** at a frozen resolution and scale. Masked members are predicted out of sample better than by season, measured habitat, continuous community gradients and prevalence-matched nonmembers. Descriptive; no mechanism implied. Same members with different edges = same plexitype, different plexiomes. |
| Exchange structure | **plexiome** | A condition-specific, evidence-graded **AND/OR dependency hypergraph**: (donor set, resource, recipient, AND/OR role, completeness-calibrated gap probability, evidence grade). OR = substitutable donors (the singers); AND = simultaneous requirement. Unknown donors, unbinned contigs and ambient pools are explicit nodes; alternate edge sets are kept. Direction comes from genomes or thermodynamics, never from covariance alone; temporal lag alone cannot establish donor direction. |
| Nutritional dependence | **plexitrophy** | Partner counts per recipient, medium and endpoint (Carbon 248927611782610944): **d_anab** = required living suppliers of anabolic precursors or cofactors; **d_catab** = required living suppliers or sinks of catabolic electron donors and acceptors, or product removal; **d_total** = their union, with dual-role partners counted once. Proposed: *plexitrophy* is d_total ≥ 2 (in microbial physiology, electron donors are nutrients), with *anabolic plexitrophy* (d_anab ≥ 2) as the stricter benchmark that P2 and P4 measure. Men et al. 2014: d_anab = 1, d_catab = 1, d_total = 2, with DMB from the medium: *mixed fan-in*. **Agreed by all three** (Astra 248927744632995840, Carbon 248927611782610944, Opus). Alternative sufficient partner sets are evaluated separately. κ = 1 (X on a simple supply cycle of length ≥ 3) is the separate cyclic criterion. P2 and P4 can only *flag candidates* for d_anab ≥ 2 from MAGs; obligate exchange needs independent evidence. For focal X under medium E and a stated endpoint: X is plexitrophic iff **in every feasible wiring (P-all semantics), (d(X) ≥ 2 OR X lies on a simple supply cycle of length ≥ 3)** (the earlier \|SCC\| ≥ 3 and minimal-support wording is SUPERSEDED, per discussion/deep/plexitrophy-computation.md and Astra 249243708419391488), evaluated per wiring (Astra 248917454944124928). Do not take coordinatewise minima across supports: a fan-in-only wiring (d, κ) = (2, 0) plus a cycle-only alternative (1, 1) qualifies. Empty or incomplete enumeration = unresolved. Fan-in vs cyclic labels describe which disjunct holds. m_S(X) = \|S\| is reported but does not define it. Serial cascades (d = 1, κ = 0, m ≥ 3) are chains, not plexitrophy; facultative facilitation (Goldford: m = 1) and dyadic syntrophy are excluded. |
| Mechanism-backed unit | **plexibiont** | A plexitype or plexon with independent evidence that members' persistence or activity depends on one another (non-substitutable metabolic complementation, perturbation or culture; expression coupling alone supports only a candidate mechanism). Host optional. Makes **no** unit-of-selection claim (Douglas & Werren 2016). |
| Instance | **plexon** | One occurrence with a boundary delimited **independently of composition**: an isolated particle, a host individual, or a tracked Lagrangian parcel. Bulk filters = "instance unresolved". Use sparingly (Plexon Inc. search collision). |

**Revised gloss (proposed by Astra 249243265693827072 after the Wienhausen 2024 test case;
Opus and Carbon agree):**
> In the stated medium, at least one viable arrangement sustains the recipient, and in *every*
> such arrangement it either needs at least two immediate nutritional suppliers or belongs to a
> supply cycle of three or more organisms.

The chain exclusion is kept. The size of the smallest self-maintaining consortium (m) is
reported alongside the focal wiring and described in prose ("consortial dependence"), not as a
new formal class.

Worked case (discussion/papers/Wienhausen2024_notes.md): the unsupplemented diatom + Colwellia +
Roseovarius culture is the smallest successful assemblage tested (m = 3). The diatom has one
immediate supplier, with the reciprocal exchange between the bacteria, so it is consortially
dependent but not plexitrophic. The growth endpoint alone can't identify which cell released
the B12.
**Correction (Carbon 249243505658347520; Opus checked the Methods):** the co-cultures had **no
organic carbon added**, and the diatom feeds the bacteria. Under d_total, the inferred wiring
has a simple 3-cycle through the diatom (diatom → DOC → Colwellia → α-ribazole → Roseovarius →
B12 → diatom), and the alternative B12 route gives one too. Astra (249243593659039744) adds that **κ_total stays
unresolved**: a wiring where the diatom feeds only Roseovarius, and Colwellia lives on
Roseovarius-derived carbon, has only a 2-cycle through the diatom. So the established result is
m = 3 in a medium without organic C. Astra also corrected a simulation example: in mep_redox T4,
AER has 0/7 calls changed by dropping the chain exclusion (AERNAR 6/7, SRB 2/2, FER 1/1). This shows why stating E and reporting (d_anab, d_total, κ, m)
matters. Carbon also votes (a) + m in prose, and keeps the chain exclusion, so the gloss now
has all three families.

**One-sentence definition (draft, all three families. Carbon 248927793148510208 added the exclusion
of ingestion; Astra defined the boundary by how the input is obtained):**
> *Plexitrophy* is obligate dependence of an organism, under a stated medium and growth
> endpoint, on metabolites or electrons transferred between coexisting living partners (including
> across a host–symbiont interface), such that no single partner suffices: either at least two
> distinct partners must supply simultaneously (d_total ≥ 2), or the organism belongs to an
> irreducible cycle of at least three (κ = 1). Inputs obtained solely by ingesting and digesting
> another organism are excluded; a mixotroph's grazing stays outside the partner count
> (Astra 248927849318629376). *Anabolic
> plexitrophy* restricts the count to suppliers of biosynthetic precursors and cofactors
> (d_anab ≥ 2, or an irreducible anabolic cycle c_anab ≥ 3).

**Architecture** is recorded separately from host association:
- **Host:** host-associated / free.
- **Topology:** hub, chain, mesh, modular, nested, or unresolved. Assigned per edge set from
  conditional independence given E, with positivity enforced.
- **n = 3:** only colliders are identifiable without orientation.
- **n ≥ 4:** undirected structure separates:
  - pure star (one articulation vertex leaving ≥ 3 singletons, diameter 2);
  - chain or tree (treewidth 1, diameter ≥ 3);
  - 2-connected mesh (treewidth ≥ 2 or a chordless cycle).

  Mixed hub-with-tail trees are called "branched tree / modular", not forced into star or chain
  (Carbon 248919437147029504).
- **Nulls:** hub concentration is tested against a null that preserves detectability, prevalence
  and edge count while degrees vary; wiring motifs against a degree-preserving null.
- **holobiont** keeps its historical meaning: host plus microbiota, which may contain chains or
  meshes.
- **dictyobiont** and **catenabiont** are dropped from Paper 1 and kept in a glossary pending a
  biological prediction that survives edge uncertainty.

**Paper 1 vocabulary (proposal):** lead with **plexitrophy** and **plexiome**, with plexitype as
the descriptive class. plexibiont and plexon appear only if results support them. That makes at
most 3–4 coinages, honouring plexichat.txt's own warning.

### Mechanism layer (Carbon; a *proposed* mechanism, endorsed with caveats)
**Plexitrophy = blocked "shooting the moon".** Under Black Queen dynamics, the helpers that
remain tend to stack several leaky functions in one lineage (d = 1 hubs, nested networks).
Multi-donor dependence should therefore arise where stacking is blocked by biochemical, redox or
O2 incompatibility, or by burden.

Combined with Opus's "exchange type predicts architecture":
- leaky public goods → hub or nested;
- sequential redox handoffs → chain;
- reciprocal amino acid auxotrophies and incompatible functions → fan-in or mesh.

**Caveat (Opus):** a single cell can separate incompatible functions in space or time
(heterocysts, diel N2 fixation). The list of incompatible pairs must be pre-registered, and
compartmentalised single-genome solutions counted.
- **Caveat (Astra):** an interior syntenic gap shows absence at that locus only. Alternative loci,
  non-homologous enzymes and shared ancestry still matter, and neither adaptive loss nor obligate
  exchange follows from the gap.
- **Caveat (Astra):** Frioux et al.'s guilds (58 × 15 × 16 = 13,920 possible triples, but only
  381 feasible) show constrained, not free, interchangeability. Keep the feasible combinations or
  compatibility constraints rather than collapsing guilds into OR-sets.

### Known traps (all three)
1. **False auxotrophy** from draft MAGs. Use contig-level downsampling nulls of ≥ 95% complete
   genomes. Nanopore contigs can confirm real deletions: the gap is flanked by syntenic core genes
   on one uniform-depth contig.
2. **Unbinned donors.** About 34% of CMO RNA reads map to unbinned contigs (AUDIT C6). Check
   community annotations before claiming that no donor exists.
3. **KO locus-tag drift** (ISSUE_ko_misalignment.md). Join on coordinates (data/ko_remap.tsv.gz)
   or fix referencesFromMag().
4. **Ambient supply versus a specific donor.** Test I(X; specific donors | aggregate guild supply
   G_m, E) > 0. The sign of RNA covariance is not diagnostic (storage, lags, repression).
5. **Facilitation is not obligacy** (Goldford Fig. 3).
6. **Pseudo-replication.** CMO has one pool per arm per winter. Lake Winnipeg supports two-year
   spatial replication, not annual recurrence.
9. **Multi-binner duplication in per_mag tables.** Carbon found this; Opus checked it in
   danaSeq/mag_analysis/modules/metabolism.nf (ecff294ca4ce), which builds per_mag from
   `cat ${contig2bin} ${binner_bins}`.
   - per_mag/\*.tsv and checkm2 quality_report.tsv pool **every raw binner** (lorbin, metabat,
     semibin, maxbin, comebin, vamb) with the consensus bins (dastool, binette, magscot).
   - Counting donors across them duplicates genomes 4–6×, which creates false OR-redundancy and
     false d ≥ 2.
   - **Use one dereplicated consensus set** (dastool-\*, skani 95% ANI).
   - Deduplicated high-quality counts (Carbon):
     - Ebb & Flow Nanopore: 165 dastool bins, 43 HQ. This is the only cohort with per_mag
       staged locally.
     - Freshwater: 1,074 bins, 191 HQ.
     - Freshwater metaMDBG: 1,130 bins, 160 HQ.
     - Chesterfield marine: 93 dastool bins, 10 HQ.
7. **Genome-size confound** in fraction comparisons. Compare plexitrophy per auxotrophy,
   conditional on genome size.
8. **The pipeline's module table is too coarse for P2 and P4.** Carbon found this; Opus checked
   it in danaSeq/mag_analysis/bin/kegg_module_completeness.py (af445eb851e1):
   - It has 55 bundled modules, only 3 of them amino acid modules (M00019 Val/Ile, M00028 Orn,
     M00015 Pro).
   - Thiamine M00127 = thiC/thiD/thiE/thiL only, with no thiazole branch (thiG, thiH, thiO).
   - B12 is M00120 (anaerobic) only: no aerobic M00122 and no bluB.
     M00120 also lists {K02229 cobG, K02228 cobF} as a step, but cobG is the aerobic O2-dependent
     enzyme, so the "anaerobic" module mixes the two routes (Opus 16:22).
   - M00174 requires both pMMO and sMMO as sequential steps, so pure-pMMO methanotrophs score 0.5.
   - It also **overcalls partial complexes.** evaluate_module gives a step credit for any KO in
     the step's set. Astra tested M00174: {K10944} scored 0.5 and {K10944, K16157} scored 1.0,
     with neither subunit set complete. Count confirmed as 55 by three methods (Astra
     248919601706352640, Opus).
   - **Rule for the curated definitions:** keep "all required subunits within a route" separate
     from "alternatives between routes", and give each rule a positive fixture and a
     partial-complex negative fixture.

   **Score P2 and P4 from KO-level definitions in metabolism/per_mag/\*.tsv, not from
   module_completeness.tsv.** The definitions are pre-registered in
   notes/p2-incompatible-pairs.md (Carbon). This is also a candidate fix to the pipeline itself;
   grid decides.

## Testable hypotheses (merged, ranked)
Each needs a frozen protocol in protocols/, cross-family review, and grid's go before any
cluster use.

**P1. Exchange role persists despite donor replacement (song vs singer).** Astra's H2, Carbon's
H2, Opus's H3.
- *Question:* does recipient activity follow the exchange role (aggregate donor capacity) or the
  donor's identity?
- *Models:* Me = environment only (season + habitat + latent community factors) as the first baseline, since guild-only and guild + identity can tie if both are uninformative (Astra 248918348255379456). M0 = ambient or guild supply G_m (binned + unbinned) + season + habitat + latent
  community factors. M1 = M0 + specific donor abundances or ligand–receptor pairs (siderophore
  BGC ↔ TonB; bluB/cobT ↔ metH vs metE).
- *Readout:* within-MAG normalised transcription log(pathway TPM / housekeeping TPM), with DNA
  depth as a covariate. Whole-event holdouts; ΔELPD against a pre-registered margin.
- *Data:* CMO paired DNA+RNA. 97 dated pairs over 39 dates in 2 winters, but one pool per arm,
  so a **descriptive pilot only**. The Cycloclasticus metabat_9342 → lorbin_184 genus-level role
  recurrence is the worked example.
- *Inventory (Astra):* the two local CMO RNA batches establish no independent donor-replacement validation series. P1 needs new or external paired DNA+RNA events.
- *Pipeline:* mag_analysis (metabolism.nf, gene_depths.nf) + illumina_rna --from_mag with
  coordinate-remapped KOs.
- *Compute:* table-scale, under 4 CPU-h; remapping about 10–20 core-h.

**P2. Plexitrophy exceeds incompleteness and stacking nulls, and concentrates where stacking is
blocked.**
*(Revised 16:24: recipient-level design adopted by all three. Spec in
notes/p2-incompatible-pairs.md 4d11eba05ae1, Carbon, 248921453407682560.)*
- Eligibility comes from recipient X's missing endpoints, before donors are inspected. Recipients
  are matched on number of missing modules, genome size, CheckM2 and GTDB clade.
- Single complete donors and ambient supply may win (d ≤ 1).
- d = 0, d = 1 and unresolved recipients stay in the denominator.
- Primary pairs: cross-regime R-I2 and R-I3. R-I6 and thiamine are d = 1 calibration; I1, I4 and I5 are descriptive.
- **Decision rule** around a frozen δ: L > δ supported; U < δ falsified; otherwise inconclusive.
- Donor partitioning (co-encoded / split / both / unresolved) and bypass buckets are reported as
  descriptive accounting. Carbon's H1 merged with Opus's H1 and H2.
- *Method:* on ≥ 90% complete, ≤ 5% contaminated Nanopore or hybrid MAGs (147 runs, about
  1.5 TB), enumerate minimal supports per sample and compute (d, m, κ) for B vitamins and
  vitamers, essential amino acids, and N and S redox intermediates.
- *Nulls:* contig-level completeness downsampling; pathway permutation preserving genome size and
  pathway prevalence.
- *Supported* if the excess over both nulls is positive (95% CI) across independent sites, and
  enriched in pre-registered incompatible pairs.
- *Falsified* only if the interval lies wholly below the pre-registered effect threshold after calibration. A negative estimate with wide uncertainty is inconclusive.
- *Inconclusive* if unbinned contigs carry more than 50% of donor KOs.
- *Compute:* about 1–2 CPU-h if existing mag_analysis outputs are reused; 200–500 core-h if MAGs
  must be built fresh (needs grid).
- *Tooling:* a small reviewed solver in scratch/. MiSCoTo or clingo installs wait for grid.

**P3. Plexitype recurrence with membership specificity.** Astra's H1.
- *Test:* masked members are predicted on held-out years or events better than by season,
  environment, continuous gradients and **matched nonmembers**.
- *Data:* a microscape ASV table at one primer region first, then shotgun data.
- *Caveat:* Lake Winnipeg supports only 2024–25 spatial replication (17 stations, 21
  station–month cells), so the claim there is spatial transfer, not annual recurrence.
- *Prior art to beat:* Wilkinson et al. 2021 JSDM; Fuhrman 2006; Venrick 1971 (abstract only).
- *Compute:* one CPU node.

**P4. Plexitrophy is enriched in diffusion-limited microhabitats.** Carbon's H3, Astra's H3.
- *Comparison:* particle, aggregate and ice fractions against matched free-living (0.2–0.8 µm)
  fractions, **per auxotrophy and conditional on genome size** (the streamlining confound).
- *Plexon-level claims* need individually isolated particles (cf. Stephens 2024: 377
  gel-isolated particles, whose richness trend inverted relative to bulk traps).
- *Inventory (Astra, notes/paired-events-and-fractions.md dd6894ced0ce):*
  - 31 ASGARD process events have 0.2, 3 and 20 µm 16S records, out of 33 process keys.
  - **Shotgun pilot candidate:** Amundsen HS22-003 (2023-10-16, 2.2 m), with three size
    fractions, matching DNA metadata and all 6 FASTQs on disk.
  - Verify the filtration order before assigning retained size intervals. No 0.8 µm cutoff is
    recorded, so "free-living" means < 3 µm here.

**P5. Exchange type predicts architecture** (Opus H1 merged with Carbon's blocked stacking;
not yet a protocol; draft statement 17:08).
- *Claim:* within one plexiome, the topology of each resource class's edge set differs by
  exchange type, as each mechanism predicts:
  - leaky public goods (B vitamins, cobamides, siderophores) → **hub or nested**: few helpers,
    many recipients (Black Queen stacking);
  - sequential redox handoffs (N and S intermediates) → **chain**: path length ≥ 3, low
    clustering;
  - reciprocal anabolic auxotrophies (amino acids) → **mesh**: cycles, treewidth ≥ 2.
- *Test:* per event, build each class's edge set from the curated KO rules (traps 8 and 9
  applied). Compute nestedness (NODF), hub concentration, mean path length and cycle count.
- *Nulls:*
  - hub concentration: a null preserving detectability, prevalence and edge count with free
    degrees;
  - wiring: a degree-preserving null;
  - completeness: contig-downsampling.
- *Supported* if, in ≥ 2 independent cohorts, each class's topology statistic exceeds its null
  in the predicted direction (lower bound > δ, frozen) **and** the classes differ from each other
  in the predicted order.
- *Falsified* if all classes are indistinguishable from their nulls (upper bound < δ), or the
  order reverses.
- *Inconclusive* otherwise, or if more than 50% of a class's edges are unresolved.
- *Caveat:* topology is called only when stable across alternate edge sets. With n = 3, stars
  and chains are not separable without orientation; directions come from genomes.
- *Data:* the same dereplicated MAG catalogues as P2. Table-scale compute.
- *Astra 248931973443534848 (applied):*
  - Exchange-class labels are assigned by **chemistry alone, before graph construction**: the
    resource list per class is fixed in advance. Edges are never selected for reciprocity or
    sequence. The same edge-eligibility rules carry into the nulls.
  - **One primary statistic per class, each paired with its own null** (Carbon
    248933060170924032):
    - Public goods: NODF on the bipartite donor × recipient matrix, against the
      **free-degree** prevalence and edge-count null. A degree-preserving null would absorb
      helper stacking.
    - Redox: mean geodesic path length on the SCC-condensed DAG, against a degree-preserving
      directed null.
    - Amino acids: fraction of nodes in non-trivial SCCs (dyadic reciprocity as secondary),
      against a degree-preserving directed null. Raw cycle counts explode in dense meshes.
  - **Guard against Cartesian donor projection:** linking every producer to every auxotroph
    creates bicliques and false 2-cycles even when one complete prototroph supplies everything.
    Keep only edges that take part in some minimal feasible support, and add a sensitivity
    check with complete-prototroph hubs conditioned out.
  - "Hub or nested" is no longer an either/or that can be chosen after looking at the data.

### Formerly open (resolved; history kept in the findings log)
- Incompatible pairs for P2: notes/p2-incompatible-pairs.md (0c91a742e263), closed into the
  recipient-level P2 (828c0c3e0dde). R-I3 and R-I2 are primary; R-I6 and thiamine form the
  d = 1 calibration stratum. The donor-side split-odds contrast was rejected by two families
  (Opus 248920882680348672, Astra 248920851231457280).
- Decision template for all protocols: one frozen margin δ; lower bound > δ supported; upper
  bound < δ falsified; otherwise inconclusive.
- Paper 1 scope: plexitrophy + plexiome + descriptive plexitype; plexibiont in the glossary
  (all three).
- P1/P4 inventory: notes/paired-events-and-fractions.md (Astra).
