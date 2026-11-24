# Plexi Pond: a browser artificial-life toy for plexitrophy

`index.html` is one self-contained file (vanilla JS + canvas, no external dependencies, light
and dark themes, works at phone width). Open it in a browser. Nothing is fetched or published.
`headless_test.js` runs the same simulation core under node, with no browser.

## The world (rules)
- **Field.** 96 × 64 grid with 9 diffusing, slowly diluted metabolite channels:
  - six essentials: B12 (cobamide), B1, B7, Met, Trp and fixed N;
  - organic C, O2 and H2.

  Organic C and fixed N are supplied abiotically. O2 exchanges with the air at the top edge.
  Light decays with depth, so there is a redox gradient: oxic at the top, anoxic below.
- **Microbes** (up to 2,000) random-walk, take up, grow and divide. Each genome is 10 bits:
  - six biosynthesis genes, one per essential;
  - Photo (oxygenic: light → energy, makes O2, exudes C);
  - Resp (C + O2);
  - Ferm (C → H2, with yield falling as local H2 rises);
  - Hyd (H2 oxidation).

  A cell grows only when it holds all six essentials and enough energy. It makes an essential
  only if it carries the gene, and otherwise must take it up from the water.
- **Leaky or private.** Each gene carries a leak bit:
  - a leaky gene releases a fraction `leakiness` of its product into the water;
  - a private gene keeps it all but costs extra.
- **Gene cost / Black Queen.** Every gene costs energy per tick. At division, genes are lost at
  rate `gene-loss rate`, regained rarely (as if by HGT), and leak bits can flip. Losing a gene
  whose product is in the water saves cost.
- **Incompatibility.** B12 (anaerobic route), nitrogenase and hydrogenase are O2-sensitive:
  - their activity falls with local O2;
  - in a genome that also carries Photo they lose a further fraction `incompatibility`.

  So the obvious "stacked" donor (a phototroph that also makes B12 and fixed N) is penalised.
- **Algae (eukaryotes).** A fixed eukaryotic type, drawn as squares:
  - larger and slower;
  - three times the photosynthesis, and twice the biosynthetic output;
  - exudes 30–60% of its fixed C whatever its leak bit;
  - makes B1, B7, Met and Trp, and leaks them if those genes are leaky;
  - can never make B12, nitrogenase or hydrogenase.

  So algae appear in the graphs as big donors of C and vitamins, and as B12 recipients.
- **Particles** (optional) are discs with diffusion reduced to 8%, an extra organic-C source,
  and slow cell movement: diffusion-limited patches.
- **Mixing** stirs the water toward its mean and moves or teleports cells.

Nothing in the code refers to plexitrophy, partners or topology. Dependence arises only from
gene loss under leakage, cost and chemistry. The analysis is computed from the run, after the
fact.

## Which principle each mechanism embodies
| Mechanism | Principle |
|---|---|
| Leaky products, gene cost, loss at division | Black Queen hypothesis (public goods, adaptive gene loss) |
| Private genes cost extra | the cost of privatisation; the cheater / producer balance |
| O2-sensitive functions, Photo incompatibility | blocked stacking: one genome cannot supply everything, so several donors may be needed (P2 mechanism; fan-in note) |
| Light gradient, O2, H2 and syntrophic H2 removal | redox gradients and catabolic services (a sink serves a fermenter). This is a gesture toward MEP, not an implementation of it |
| Particles versus mixing | spatial structure and diffusion limitation (P4: particle versus free-living) |
| Algae | eukaryotic donors of C and metabolites; B12 auxotrophy (Helliwell-type cases) |
| Fast genome loss versus slow ecology | two timescales (genome versus population), which are coupled |

## Ground-truth panel: realised transfers versus required partners
The panel shows two different quantities. Everything in it is a **live, threshold-based
estimate**.
- **Realised transfers (who fed whom).** Every uptake, in an 8 × 8-cell block, is credited to
  the genotypes that recently released that metabolite in the block and its four neighbours, in
  proportion to their release. Fluxes are exponentially averaged, with a time scale of about
  165 ticks. This is an approximation, not tracked molecules.
  - A donor's share is a fraction of the recipient's **total cellular supply** of that channel:
    endogenous synthesis kept plus everything taken up.
  - Channels are first classified as *required* or not. Transfers ≥ 15% on a required channel
    are realised supplier edges. Transfers on other channels (a prototroph scavenging leakage,
    H2 removal) are drawn dotted as *opportunistic* and never counted.
  - `used` is the number of distinct suppliers on required channels.
  - `d_eff` gives each required resource equal weight. Each donor's weighted shares are summed,
    and `d_eff` is the inverse-Simpson effective number over donors, so 1 ≤ d_eff ≤ used (and
    d_eff = 1 when used = 1).
  - A supply from the medium or the air counts 0.
- **Requirements.** A requirement is any of:
  - an essential that the recipient lacks, or carries but has mostly switched off (activity
    < 30%), and that the medium does not mostly supply (medium share < 50%);
  - O2, for obligate respirers;
  - organic C from others, for non-phototrophs when the medium supplies < 50% of it.

  H2 removal is never a requirement. In the simulation step, fermentation yield only falls with
  H2 (1/(1 + H2/K)), and H2 also escapes and decays, so a sink helps but is never necessary.
  `analyze` and `step` now agree on this. A requirement with no supplier at ≥ 5% is
  *unresolved*.
- **Required partners, live estimate.** Every genotype supplying ≥ 5% of a requirement is an OR
  alternative, and requirements are ANDed.
  - `d_anab`, `d_catab` and `d_total` are the exact minimum set cover, computed by dynamic
    programming over requirement bitmasks with no cap on donors. Substitutable donors therefore
    count once, and so does a partner with two roles.
  - `minCover` returns "infeasible" (some requirement has no donor), which is distinct from d = 0
    (no requirement: autonomous).
  - A genotype is **plexitrophic** if d_total ≥ 2 or it lies on a cycle ≥ 3 (heuristic, below).
  - A genotype is *resolved* if it has no unresolved requirement, or if d_total ≥ 2 already
    holds. Adding requirements can only raise a minimum cover, so that call is safe.
  - Percentages use these denominators:
    - % plexitrophic and % redundant-use: resolved cells;
    - fan-in: resolved auxotroph cells;
    - % using ≥ 2 suppliers: cells of the analysed genotypes.

    The panel also reports coverage (the share of all cells in the analysed genotypes, which
    are those with ≥ 4 cells, top 40) and the resolved fraction.
  - A genotype that uses ≥ 2 suppliers when one would suffice gets a grey ring, labelled
    "redundant use". It is not plexitrophic.
- **Cycles ≥ 3 are a heuristic.** Only simple directed cycles through *sole-supplier* edges are
  searched. That is the coordinatewise shortcut which
  `discussion/deep/plexitrophy-computation.md` (§3.3, A1/F3) shows can miss cycles that exist
  only in alternative supports. A per-support evaluation is not implemented.
- **Required partners, counterfactual check (knock-out button).** This tests **all** realised
  suppliers of the chosen recipient (share ≥ 15% on a required channel), with no cap. Arms:
  - the control;
  - each supplier removed;
  - each kept as the only tested supplier (when there are ≥ 3);
  - all tested suppliers removed.

  Genotypes below the 15% threshold, and every other genotype, stay present in every arm as the
  "background". Removed genotypes are banned from re-evolving, and all arms share the same
  random stream. "Persists" means mean abundance over the last quarter of the horizon is ≥ 30%
  of control. Outcomes:
  - *inconclusive*: the control has < 1 cell, so it is extinct or nearly so;
  - *not needed*: removing all tested suppliers is harmless;
  - *one suffices*: removing all is fatal, but one tested supplier alone (with the background)
    suffices. Not plexitrophic among the tested set;
  - *at least L*: no single tested supplier suffices and removing all is fatal. L = max(2,
    number of individually fatal removals). This is a **lower bound** under monotone rescue
    (adding a supplier never hurts), relative to the tested set and its background. Exact counts
    would need subset tests;
  - *ambiguous*: removing one supplier hurts but removing all does not, or one supplier suffices
    while removing a different one is fatal (noise or non-monotone effects such as competition
    release).

  One seed per arm, so the outcome is noisy. Many suppliers mean many arms (2k + 2), and the
  page asks before running more than 12.
- **Architecture labels** are in plain words: none, chain, fan-in, hub, tree or mesh. They are
  computed over the 12 most abundant genotypes:
  - mesh = a directed cycle ≥ 3 or any undirected reticulation;
  - the required graph uses major (≥ 15%) or sole-supplier edges;
  - the realised graph uses required-channel edges only.

## Observer's view
Tap the field to take a sample (a disc), or take 8 random samples or a vertical transect.
- **What the observer gets.** Per-sample genotype counts, with a detection limit of 2 cells.
  Each detected genotype is a MAG drawn from a synthetic assembly (below). Its observed genes
  are those on surviving contigs, deterministic per genotype.
- **How it infers edges.** It uses genome complementarity (B lacks f, A has f) plus association
  across samples. Association is a Pearson correlation of log relative abundance, or a ridge
  partial correlation given all taxa. With fewer than 4 samples it falls back to shared
  presence.
- **What it cannot see:** leak bits (private producers look like donors), the abiotic medium,
  O2 inhibition and incompatibility, and fluxes.
- **Scoring against realised required-channel transfers.** The *eligible* ordered pairs are
  those where both genotypes were observed and are in the analysed truth set.
  - Conditional precision and conditional recall are computed on that same eligible set.
  - Predicted edges outside it are reported as unscorable.
  - Coverage = eligible true edges / all true edges.
  - End-to-end recall = correct edges / all true edges.
  - Plexitrophy calls are compared only for resolved truth genotypes.
- **Synthetic assembly and "Export samples".** The intact genome is laid out first: a few blank
  (unannotated) CDSs, then each pathway gene followed by 2–5 blank CDSs, with its order fixed
  per genotype. It is then cut at simulated breakpoints (p = 0.14 per CDS boundary), and whole
  contigs are lost with probability 1 − "contig survival". Only surviving features are
  exported, in coordinate order, with consecutive per-contig indices (`protein_id =
  <contig_id>_<n>`). The contig-end uncertainty is therefore real, not added afterwards. The
  TSV has sections:
  - abundance (genomes × samples; the columns are shared event IDs `EV0001…`);
  - samples (position, particle flag);
  - per_mag (mag_id, domain, contig_id, protein_id, feature_order, KO, EC, gene). Domain is
    Eukaryota for algae and Bacteria otherwise;
  - contigs (with lengths);
  - completeness;
  - **boundary truth**: for each exported pathway CDS, `end` or `internal` from the simulated
    breakpoints.

  `check_export_p2.py` parses the export with `scratch/p2_scorer.py`'s `parse_per_mag_rows`
  (`complete_cds_table=True`) and compares the parser's edge status with the boundary truth. In
  the test export the parser agreed on 122 of 122 pathway CDSs (37 end, 85 internal), across
  25 MAGs, 4 of them eukaryotic.
- **KO IDs are illustrative representatives** written from memory, and must be checked before
  real use. K00768 is cobT: the DMB phosphoribosyltransferase in lower-ligand assembly and
  salvage, *not* de novo corrin-ring synthesis. K00549 is metE.

## Challenges (prediction by synthesis)
Live challenges are scored as "completed at tick X", and earlier completion scores higher:
- *Fan-in*: ≥ 50% of auxotrophic cells have required d_anab ≥ 2, held for 400 ticks.
- *Ring of three*: a required cycle ≥ 3 (heuristic detector). **It never occurred**: 0 of 32
  sweep runs, plus about 20 test runs at 3,000–4,000 ticks. A likely reason, from Carbon's
  review: every biosynthesis gene is a parallel, input-free reaction under uniform C and N
  supply, so there is no precursor cascade or thermodynamic arrest to force a ring. Adding a
  precursor dependency or an H2 threshold could make rings reachable. Neither is implemented, so
  the ring may be unreachable under these rules, and the game says so.
- *Busy but not plexitrophic*: ≥ 40% of cells use ≥ 2 suppliers while < 10% are plexitrophic.
- *Fool the observer*: with ≥ 6 fresh samples, the observer calls a mesh while the true required
  graph is acyclic.
- *Black Queen market*: ≥ 600 cells and no prototroph among the common genotypes.
- *Predict, then synthesise*: you state the % of cells plexitrophic at tick 4,000 and the
  architecture. Three fresh seeds (101, 202, 303) run with your settings. Score =
  max(0, 100 − 2·|error in percentage points|) + 15 per seed whose architecture matches (max 145).

## What it can and cannot show
**Can show (under its assumptions only).**
- Multi-donor obligate dependence can emerge from leakage, gene cost and incompatible
  chemistry, with no hard-coded partners.
- Its frequency is graded and depends on the settings. The headless sweep is in
  `sweep_3000.tsv`, made by `sweep.js`: 32 settings, 1 seed each, 3,000 ticks, rerun after the
  review fixes. These are **live, threshold-based estimates over resolved cells**, and the
  analysed genotypes covered 74–97% of cells.
  - Plexitrophic: 0–13% of resolved cells in the 16 runs at leakiness 0.3; 0–47% in the 16 runs
    at leakiness 0.8; 21–47% in the 4 runs with leakiness 0.8, gene cost 0.9 and mixing 0.3.
  - Required cycles ≥ 3 (heuristic): 0 of 32.
  - Required-partner architecture: mesh in 26 of 32 runs, chain in 3, tree in 2, hub in 1.
    Chains occur at low leakiness with high gene cost and no mixing: that is where "Fool the
    observer" is reachable.
- Realised multi-supplier use on required channels and required d diverge.
  - % using ≥ 2 suppliers ranged from 0 to 100%. It reached 78–100% in well-mixed, high-cost
    worlds, where many genotypes share each pool.
  - Much of that use is redundant: in the default run, 14% of resolved cells were
    redundant-use only, against 13% plexitrophic.
- An observer working from co-occurrence and complementarity does poorly. In the default
  headless run (12 random samples, contig survival 0.85), on the 10 eligible true edges,
  conditional precision was 0.17 and conditional recall 0.30. But coverage was only 10 of 62
  true edges, so end-to-end recall was 0.05.

**Cannot show.**
- Anything about real organisms or real rates. All constants are arbitrary.
- That the mechanisms here are the ones operating in nature.
- Thermodynamic MEP: energy yields are fixed numbers, not free energies.
- Fluxes attributed by local release share are an approximation, not tracked molecules.

**Seed noise is large.** In the headless test, "No incompatibility" gave 0% plexitrophic in
2 of 2 seeds, against 1–13% at the default setting. That is suggestive of blocked stacking, but it is n = 2 and not
a result.

## How to run
- Browser: open `index.html`.
- Headless: `node headless_test.js [ticks=3000] [seeds=2]`. The node binary is at
  `/home/grid/.local/share/agent-relay-trial/node_modules/node/bin/node`.
  - It runs unit checks:
    - minCover OR/AND, Astra's 15-alternative case, DP against exhaustive search on a 40-donor
      instance, infeasible versus autonomous;
    - the architecture classifier and the 2-cycle exclusion;
    - knock-out interpretation, including inconclusive, lower-bound, any-two-of-three and
      ambiguous cases, and that all 7 of 7 suppliers are tested.

    It also checks same-seed determinism; field and population invariants; d_eff bounds; that
    H2 is never a requirement; that phototrophic prototrophs register no suppliers; and the
    observer's eligible-set accounting. It runs 5 contrasting settings × 2 seeds, the observer
    and one live knock-out, and writes an export to `$PLEXI_EXPORT` (default
    `/tmp/plexi_export_test.tsv`).
  - Then run `python3 -I check_export_p2.py /tmp/plexi_export_test.tsv`.
  - Result: 29 of 29 node checks pass, and the parser check gives 122 of 122.
  - The core steps at about 2.3–4.6 ms per tick at up to 2,000 cells (this run was on a shared,
    loaded machine; an earlier run measured 0.9–1.9 ms), so 30 fps holds at the default speeds.
- Seeds: the world seed is set in the UI (default 1). The RNG is a seeded mulberry32 copied
  into knock-out clones. Particles are placed from the seed.

## Files (sha256, first 12 hex characters)
These were current after the review fixes. Recompute with `sha256sum`.
- `index.html` `13265d978414`
- `headless_test.js` `7a0bd452c84b`
- `check_export_p2.py` `b3821ad99afb`
- `sweep.js` `ac6489690862`
- `sweep_3000.tsv` `f5c8acacf0f2`

**Review status:** written by Claude (Opus). Reviewed by Astra (Codex; relay
248981991110131712, refinement 248982684969648128) and Carbon (Gemini; relay
248982555596341248). All their points are addressed above, except that a per-support cycle
search was not implemented: cycles are labelled heuristic instead. The fixes themselves have
not been re-reviewed.
