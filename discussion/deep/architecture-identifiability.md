# Which plexiome architectures can be identified from observational data?

Author: Claude (Opus) deep-dive sub-agent for claude-plexi. Date 2026-11-18. **Status: draft for
cross-family review** (AGENTS.md: results and code need a reviewer from another model family).
Nothing was posted to the relay.

**What I read.** AGENTS.md. In `discussion/plexi-concept-review.md` (597 lines as read): Astra
248914056547057664 (l. 143-149), Carbon 248914507304714240 (l. 155-158), the "Converging" item 2
(l. 329-331), the merged-synthesis **Architecture** block (l. 382-397) and **P5** (l. 544-575).
No file under `data/` was opened. The literature below is cited from memory, not re-read for this
note. §9 marks how confident I am about each reference.

**Code** (standard library + numpy 2.3.5 + scipy 1.17.0, run with `python3 -I`), all in
`scratch/deep/`. First 12 hex characters of each sha256:

| file | sha256 | output |
|---|---|---|
| arch_lib.py (graph classes, generators, partial-correlation and PC tests, tetrachoric) | 2e5efe1e6f89 | n/a |
| population_demos.py (exact population results, §2-§4) | cf50a689e8e7 | population_demos.out |
| dag_orientation.py (CPDAG, Meek rules, PC vs GGM, §3) | 3c51a7cff091 | dag_orientation.out |
| confounding_demo.py (habitat filtering, closure, §5) | 2414433b79cd | confounding_demo.out |
| power_sim.py + summarise_power.py (sample size, §5.3) | 851eb0b43a45 / 6b65ee7cea69 | power_sim.tsv, power_summary.out |
| decision_rule.py (three-valued-edge rule, §6) | e7ed82eca55d | decision_rule.out |

Labels used below: **[Proof]** is a full argument given here. **[Std]** is a standard result
(reference in §9). **[Sim]** is a numerical result from the scripts, valid only for the simulated
model. **[Conj]** is a conjecture or judgement.

---

## 0. Verdict in brief

1. **Claim (a)** (n = 3: fork ≡ chain, collider identifiable) is **true**, but only for an acyclic
   causal model with no hidden common causes beyond the measured E, with faithfulness, and with
   the independence tests carried out **given E only, then given E + A**. It is **false for the
   usual microbial "network" estimator** (partial correlation given all other taxa). There the
   collider becomes a triangle, a "mesh", because B and C are dependent given A. At n = 3, star
   and chain are also the *same undirected graph*: the only thing that differs is the arrow
   pattern. The identifiable collider is exactly the **fan-in (AND-clause, d = 2) shape**, so the
   one plexitrophy-relevant shape is also the one observational data can see.
2. **Claim (b)** (n ≥ 4: undirected structure separates star, chain and mesh) is **true for the
   skeleton of a DAG** as recovered by separating-set search (PC/GES-type adjacency), and true
   trivially if the plexiome *is* an undirected Markov random field. It is **false for the
   concentration graph** (partial correlation given all others, graphical lasso, SPIEC-EASI-style
   networks): that graph is the *moral* graph. A sink hub with k ≥ 2 donors becomes a clique
   K_{k+1}. Of all orientations of an n-node tree, only n / 2^(n-1) keep the moral graph a tree
   (0.50 at n = 4, 0.06 at n = 8; enumerated exactly, §2.3). The errors run one way:
   **tree-type calls from a concentration graph are sound for the skeleton, mesh calls are not.**
3. Carbon's three classes **do not partition the graphs**. A cycle with a pendant is neither a
   2-connected mesh nor a tree. "Treewidth ≥ 2" and "contains a chordless cycle" are both just
   "has a cycle", which a triangle already satisfies. Disconnected graphs are not covered. §1
   gives a partition based on the block-cut tree. "Nested" is not a property of an undirected
   graph: it needs a donor × recipient orientation.
4. **Orientation.** The CPDAG orients a sink hub (≥ 2 donors) but never a source hub (a
   public-goods donor) or a directed cascade. Each of those is Markov-equivalent to n
   orientations. Genome-derived directions plus Meek's rules orient a cascade once one edge is
   known on the source side. However, **genome direction (metabolite flow) is not causal direction
   in abundance** (syntrophic pull, product inhibition), so background knowledge must enter as a
   testable hypothesis, not as a hard constraint (§3.3).
5. **Cycles.**
   - For linear feedback, including gLV equilibria where the noise enters through growth rates,
     d-separation still holds. The concentration graph is then the support of AᵀD⁻¹A. A receiver
     hub and a reciprocal ring of 5 both give K5; a one-way ring gives a chordless C5.
   - For nonlinear feedback only σ-separation is guaranteed. Under σ-separation, the topology
     *inside* a strongly connected component cannot be recovered from conditional independence:
     a ring and a clique look the same.
   - **A "mesh" call licenses only "not decomposable by conditioning on the measured taxa and
     E".**
   - *Revised after review (§10):* a chordless cycle in the **concentration graph** has no
     causal reading. A DAG with no feedback, selection or latents can produce an induced
     concentration C4 (Astra's counterexample).
   - The diagnostic pattern is a chordless cycle of length ≥ 4 in the **causal skeleton**
     (separating-set search), where every unshielded triple on it is a definite noncollider
     with separator evidence. Under acyclicity, no selection and MAG-faithfulness, that pattern
     is impossible (§4.3). Its explanations are then feedback, selection or unfaithfulness.
6. **Confounding** (§5, simulated, five taxa that do not interact):
   - Habitat filtering produces a spurious **star** around whichever taxon best proxies E
     (81 % of replicates at N = 400).
   - It produces a spurious **mesh** when niche responses are unimodal and E is adjusted linearly
     (100 % at N ≥ 400), or when a second gradient is unmeasured (94 %).
   - It produces a spurious **chain** under seasonal succession (taxa tracking E at different
     lags). Here conditioning on the current E does *not* help: path(4) in 81 %.
   - Conditioning on correctly measured E with the correct functional form removes the
     spurious edges (≥ 94 % empty graphs). Measurement error in E (reliability 0.7) brings the
     spurious star back (89 % at N = 400).
   - Closure turns a true path into a mesh. This happens when abundances are made relative
     within the focal taxa (100 %), and under whole-community clr when the community is only
     10 taxa (99 %). With 55 taxa the distortion vanishes (path 96 %).
7. **Sample size** (§5.3): *a planning benchmark for the naive pipeline only*
   (significance-selected full-order partial correlations, i.i.d. Gaussian graphical models). It
   is not a bound for PC/FCI or for the full rule (revised after review, §10). Independent samples for 80 % class accuracy:
   - Gaussian data: **190-330** for n = 4-8 at |partial r| = 0.3, and **450-715** at 0.2.
   - Zero-inflated data analysed naively: about 1.3-1.5× more. Accuracy then *falls* again at
     large N for stars (0.35 at N = 6400, n = 8).
   - Presence/absence data with naive phi partial correlations: about 3× more, and
     asymptotically biased (n = 6 star: 0.03 at N = 6400). Dichotomising breaks the Markov
     property.
   - Tetrachoric partial correlations remove the bias, at **600-1360** (0.3) or **1300-1950**
     (0.2).
   - The three-valued rule needs ~1000-3000 samples to call a tree in a single-spec, unblocked
     calibration. Even that is not a certification of the full rule.
   - No dataset-specific threshold is claimed. The CMO RNA series stays descriptive on
     independent-pool grounds.
8. **Decision rule** (§6). Each edge is present, absent (equivalence-tested) or unresolved. A
   class is called only if it is invariant over all completions of the unresolved edges, over
   all pre-registered specifications and over all alternative edge sets. Calls of the causal
   kind are restricted to the directions that are sound (§2.4).

---

## 1. Objects and a proper class partition

Three different graphs get called "the network". It matters which one is meant.

- **Skeleton S(G)** of a causal DAG G: u - v adjacent iff u → v or v → u in G. Under the causal
  Markov property, faithfulness and causal sufficiency (all common causes measured, here E
  included), u and v are adjacent **iff no subset of the other variables separates them**. That
  is what PC finds. [Std: Verma & Pearl 1990; Spirtes, Glymour & Scheines 2000]
- **Concentration graph** (pairwise Markov graph, Gaussian graphical model, Markov random field):
  u - v iff u and v are *not* independent given *all* other variables. For Gaussian data this is
  the support of the precision matrix, so a partial-correlation network, a graphical lasso or
  SPIEC-EASI estimates it. For a DAG model it equals the **moral graph** M(G): the skeleton plus
  an edge between every two parents of a common child. [Std: Lauritzen 1996]
- **CPDAG** (essential graph): the skeleton plus every orientation shared by all DAGs in the
  Markov equivalence class (v-structures, closed under Meek's rules R1-R3). [Std: Verma & Pearl
  1990; Meek 1995; Andersson, Madigan & Perlman 1997; Chickering 2002]

**Class partition used in all code** (`arch_lib.classify`). For a connected undirected graph on
n ≥ 4 vertices:

| class | definition | Carbon's wording |
|---|---|---|
| star (hub) | tree, diameter 2 (one centre, n-1 ≥ 3 leaves) | "one separator into ≥ 3 singletons" ✓ |
| path (chain) | tree, maximum degree ≤ 2 | "tree, diameter ≥ 3" also admits branched trees ✗ |
| branched tree | any other tree (spider, caterpillar) | "branched tree / modular" (248919437147029504) ✓ |
| mesh | 2-connected (no articulation vertex), so every vertex lies on a cycle | "2-connected mesh" ✓ |
| modular | connected, has a cycle, has an articulation vertex (≥ 2 blocks) | **missing**: e.g. a triangle with a pendant |
| disconnected | > 1 component | **missing** |

Notes. (i) A graph has treewidth ≥ 2 iff it has a cycle, and every cycle contains a chordless one
(take the shortest cycle). So "treewidth ≥ 2 or a chordless cycle" just means "has a cycle". It
does **not** mean 2-connected. K4 is chordal and has treewidth 3. Only a chordless cycle of length
**≥ 4** carries extra information (§4.3). (ii) At n = 3 star = path = P3. (iii) Nestedness
(NODF) is defined on a bipartite donor × recipient matrix. It needs orientation and cannot be read
from any undirected conditional-independence graph.

---

## 2. Claims (a) and (b)

### 2.1 Claim (a), n = 3

**Proposition 1** [Std; Verma & Pearl 1990]. Two DAGs are Markov-equivalent iff they have the
same skeleton and the same v-structures. Assume acyclicity, causal sufficiency given E, and
faithfulness. Then for the three DAGs on {A, B, C} with skeleton B - A - C:

- the fork B ← A → C and both chains B → A → C and C → A → B imply B ⫫ C | A, E and
  B not⫫ C | E;
- the collider B → A ← C implies B ⫫ C | E and B not⫫ C | A, E.

The first three are equivalent; the collider is identified. Exact population values
(`population_demos.out`, coefficients 0.6, unit noise):

| DAG | corr(B, C) | pcor(B, C \| A) |
|---|---|---|
| fork B ← A → C | 0.265 | 0.000 |
| chain B → A → C | 0.295 | 0.000 |
| collider B → A ← C | 0.000 | −0.265 |

**Corrections to the record.**

1. In Astra's message "fork and chain with the same middle node" is right. The P5 caveat
   (l. 566), "with n = 3, stars and chains are not separable", is right in substance but
   misleading: at n = 3 there is only one connected tree. What is unresolvable is the arrow
   pattern, not the topology.
2. The collider is identified **only by testing the marginal B ⫫ C | E together with the
   conditional dependence given A**. A partial-correlation network conditions on everything, so
   for a collider it reports the triangle A-B-C.
3. Ecological reading: a collider at A is two donors feeding one recipient, which is fan-in
   plexitrophy with d = 2 if both are required. The signature is *explaining away*. Across
   samples the donors' abundances are independent given E, but once the recipient's abundance is
   known they become dependent (negatively, if both raise A). That is a concrete,
   pre-registrable test of a genome-predicted AND-clause (§3.3).

### 2.2 Claim (b): which undirected graph?

**Proposition 2** [Proof]. Let G be a DAG, S its skeleton and M its moral graph.

- (i) S ⊆ M, and S and M have the same connected components.
- (ii) S = M iff G has no v-structure.
- (iii) If M is a forest, then S = M and G has no v-structure.
- (iv) If S is 2-connected, M is 2-connected.

*Proof.* (i) Moralisation only adds edges between two parents of a common child c. Those two
parents are already joined through c in S, so components do not change. (ii) The added edges are
exactly the non-adjacent co-parent pairs, that is, the v-structures. (iii) By (i), S is a spanning
subgraph of the forest M with the same components. A forest has no proper spanning subgraph with
the same components, so S = M, and by (ii) there are no v-structures. (iv) Adding edges cannot
create an articulation vertex. ∎

**Consequences for claim (b).**

- **(b) via skeleton search: true.** S is identified by separating-set search (PC: adjacency iff
  no separating set) under sufficiency and faithfulness, and the class is a function of S. So
  star, path, branched tree, mesh and modular are all identified at n ≥ 4. The conditioning sets
  range over *all subsets*, not only "all other taxa".
- **(b) via the concentration graph: one-way.**
  - A tree or star call is sound for S (iii), and it also certifies that **no node has two
    parents**: no fan-in hub anywhere, which is the right reading for P5's "public goods → hub".
  - A mesh or modular call is not sound. S could be any tree with at least one collider.
  - If S is a mesh, M is a mesh too (iv), so true meshes are not hidden; they simply cannot be
    told apart from moralised trees.
- **(b) as a statement about MRFs: true but empty.** If the plexiome is *defined* as an undirected
  graphical model, its graph is identified under positivity and faithfulness (pairwise Markov =
  global Markov for positive densities, Hammersley-Clifford) [Std: Lauritzen 1996]. That
  definition, however, discards the direction the team wants to read from genomes.

### 2.3 How often moralisation changes the class (exact enumeration)

**Proposition 3** [Proof]. Of the 2^(n-1) orientations of a tree on n vertices, exactly n have no
v-structure, so in exactly n of them M = S.

*Proof.* With no collider, every vertex has ≤ 1 parent. A tree with n-1 edges in which every vertex
has ≤ 1 parent has exactly one parentless vertex: an arborescence, fixed by its root. There are n
roots. ∎

Enumeration (`population_demos.out`; the class of M over all orientations):

| n | star skeleton: M is star / modular / mesh | path skeleton: M is path / modular | n / 2^(n-1) |
|---|---|---|---|
| 4 | 4 / 3 / 1 of 8 | 4 / 4 of 8 | 0.500 |
| 5 | 5 / 10 / 1 of 16 | 5 / 11 of 16 | 0.312 |
| 6 | 6 / 25 / 1 of 32 | 6 / 26 of 32 | 0.188 |
| 8 | 8 / 119 / 1 of 128 | 8 / 120 of 128 | 0.062 |

These counts are a combinatorial fact about orientations, not a prior over biology. They show
that the hub/chain/mesh contrast in a concentration graph depends mostly on how many colliders
exist, not on the skeleton.

Named examples (oracle, `population_demos.out`):

| true DAG | concentration graph | PC skeleton (oracle) |
|---|---|---|
| sink star B,C,D,E → A (fan-in) | **K5: mesh** | star, with all 6 v-structures at A |
| source star A → B..E (public-good donor) | star | star, no v-structure |
| mixed B,C → A → D,E | modular (star + B-C) | star |
| directed chain A → … → E | path | path |
| zig-zag A → B ← C → D ← E | modular (2 triangles) | path |
| DAG on chordless C4 | mesh with a chord | C4 |

**Finite samples** (`dag_orientation.out`, n = 6, linear SEM on uniformly random orientations,
|β| ~ U(0.4, 0.8), 100 reps). Class accuracy of the concentration graph vs the PC-stable skeleton
(α = 0.01):

| skeleton | N = 200 | N = 1000 | N = 5000 |
|---|---|---|---|
| star | 0.26 vs 0.79 | 0.30 vs 0.95 | 0.29 vs 0.91 |
| path | 0.28 vs 0.95 | 0.04 vs 1.00 | 0.04 vs 0.98 |
| branched tree | 0.25 vs 0.95 | 0.09 vs 0.96 | 0.09 vs 0.97 |
| mesh (C6) | 0.98 vs 0.95 | 1.00 vs 1.00 | 1.00 vs 1.00 |

For the concentration graph, more data makes things *worse*: the moral edges become detectable.
When the true star was a sink (≥ 2 parents at the hub), PC found a v-structure at the hub in
0.97-1.00 of runs.

### 2.4 What each source of evidence identifies

Assume acyclicity, sufficiency given E and faithfulness.

| target | (i) concentration graph | (ii) CPDAG (PC/GES) | (iii) CPDAG + genome orientations (Meek) |
|---|---|---|---|
| undirected class of S (star / path / tree / mesh / modular) | only when M is a forest (then = S); otherwise unidentified | **identified** | identified |
| sink hub (≥ 2 donors → recipient; fan-in) | seen as a clique; not identified | **identified** (compelled v-structures) | identified |
| source hub (public-goods donor) vs relay hub | identified *as "a hub without fan-in"* only | **not** identified (n equivalent orientations) | identified once one hub edge is known to point **into** the hub (then R1 orients the rest outward), or all edges are known |
| directed cascade (serial chain) | path only | not identified (n roots) | identified once the source end is known; one known edge u → v orients everything downstream of v (R1) |
| nested (bipartite donor × recipient) | no | only where compelled | as far as Meek closure reaches |
| within a cycle (mesh) direction | no | partially (v-structures on the cycle) | partially (R2-R4) |

Population averages over orientations (`dag_orientation.out`). Uniform random orientation is
used only as a reference distribution:

| skeleton | share of edges oriented: CPDAG | + 1 genome edge | + 2 genome edges |
|---|---|---|---|
| star n = 5 | 0.69 | 0.80 | 0.91 |
| path n = 6 | 0.40 | 0.58 | 0.76 |
| branched tree n = 6 | 0.50 | 0.67 | 0.78 |
| mesh C5 | 0.53 | 0.67 | 0.79 |

Meek's rules with background knowledge are sound and complete for orientation when the knowledge
is a set of directed edges consistent with some member of the class [Std: Meek 1995; see also
Perković, Kalisch & Maathuis 2017 on CPDAGs with background knowledge].

---

## 3. Genomes as background knowledge: what (c) gives and what it costs

### 3.1 What it gives
A genome-derived direction "donor encodes the pathway, recipient lacks it" is information that
covariance cannot provide. Combined with a CPDAG it orients cascades and hubs as shown in §2.4.

### 3.2 Caveat: metabolite flow is not causal direction in abundance [Conj, with standard ecological reasons]
The DAG over abundances needs arrows of the form "intervening on X's abundance changes Y's
abundance". Exchange can reverse this:

- syntrophy, where the H₂ or formate consumer drives the donor's growth thermodynamically;
- product inhibition relieved by the recipient;
- a recipient that consumes a public good and thereby draws down a donor-favouring resource;
- feedback over time (§4).

So the genome arrow donor → recipient is a hypothesis about the abundance arrow, not a fact.
Imposing it as hard background knowledge can produce an orientation that the data contradict,
and Meek closure will then propagate the error.

### 3.3 Use it as a testable prediction
Treat each alternative edge set G_k (with its genome orientations) as a hypothesis with
**testable implications** in the abundance data:

- **Fan-in (AND-clause) at recipient R with donors D1, D2, no other path between them:**
  predicts D1 ⫫ D2 | E **and** D1 not⫫ D2 | E, R (a v-structure).
  - If D1 not⫫ D2 | E, there is another D1-D2 path or a shared driver, and the prediction is
    *unresolved*.
  - If D1 ⫫ D2 | E, R holds with equivalence bounds, then abundance-mediated fan-in is
    *falsified* for this edge set. The genome inference is not falsified: the exchange may be
    real but not limiting.
- **Cascade D → M → R:** predicts D ⫫ R | E, M and no v-structure at M.
- A conflict with a **compelled** CPDAG edge is recorded as a finding. It is not resolved by
  overriding the data.

Caveat for AND-clauses [Conj]: if growth is ∝ min(D1, D2) (Liebig-like), the conditional
dependence given R is nonlinear and may be weak in partial correlation. A rank- or kernel-based
conditional-independence test is then needed, which costs more samples.

---

## 4. Cycles (claim (d)): what is identifiable when exchange is reciprocal

### 4.1 Linear feedback and gLV equilibria [Std + Proof]

For a linear SEM x = Bx + e with I − B invertible and independent errors (D = cov(e)), even with
directed cycles:

- cov(x) = (I−B)⁻¹ D (I−B)⁻ᵀ, so precision = (I−B)ᵀ D⁻¹ (I−B);
- generically, entry (i, j) is non-zero iff i and j are adjacent or share a child;
- d-separation is a valid global Markov property for such models [Std: Spirtes 1995].

**gLV equilibria** [Proof]. At an interior equilibrium, r + A x* = 0, so x* = −A⁻¹ r. If the
samples differ through growth rates r = r₀ + ε with independent ε (D diagonal), then
cov(x*) = A⁻¹ D A⁻ᵀ and **precision = Aᵀ D⁻¹ A**: (i, j) is *generically* non-zero iff some k is
affected by both i and j (column supports overlap). Signed contributions Σ_k A_ki A_kj / d_k can
cancel on a measure-zero set (Astra 248969617523261440). Exact results (`population_demos.out`, A_kk = −1):

| gLV interaction structure | concentration graph |
|---|---|
| donor hub (A benefits B..E only) | star |
| receiver hub (B..E benefit A) | **K5** |
| reciprocal ring of 5 (mutual cross-feeding) | **K5** (every pair at distance 2 shares a neighbour) |
| one-way ring of 5 | **chordless C5** |

So **under the most standard community model, a receiver hub and a reciprocal ring are
indistinguishable from a clique** in a partial-correlation network.

Caveats:

- Noise in A itself, as opposed to in r, makes x* nonlinear in the noise.
- Requiring x* > 0 (feasibility) is a selection step.

### 4.2 Nonlinear feedback: σ-separation [Std]

For nonlinear cyclic structural causal models, d-separation can fail [Std: Spirtes 1995; Neal
2000 gives counterexamples, flagged in §9]. Under unique solvability, the **σ-separation** Markov
property holds [Std: Forré & Mooij 2017; Bongers, Forré, Peters & Mooij 2021]. It treats a
conditioned non-collider as blocking only if it points to a node outside its own strongly
connected component (SCC).

Consequence:

- σ-separation in G is d-separation in the *acyclification* of G, which joins every pair of
  nodes in the same SCC and gives each node all of its SCC's parents.
- **Within an SCC, no conditional independence among members is guaranteed.** A ring of
  cross-feeding and an all-pairs clique are σ-Markov-equivalent.
- Constraint-based search with σ-faithfulness identifies at most a PAG-type object over SCCs
  [Std: Mooij & Claassen 2020, FCI under σ-separation; Richardson 1996 for the linear/d-sep
  cyclic case (CCD)].

Numerics, offered as honest negative evidence (`population_demos.out`). A tanh-saturating ring
with multiplicative noise (N = 2×10⁵), testing a vs c given {b, d}:

- cyclic ring: partial correlation of standardised residuals +0.0040;
- acyclic control: +0.0036;
- linear cyclic ring: +0.0005;
- null s.e. 0.0022.

The cyclic ring is no more dependent than the acyclic control. So the d-separation violation that
σ-separation allows **need not be large**. Observational data can then neither confirm nor rule
out within-SCC structure. "Mesh" calls for nonlinear feedback are not identified in either
direction [Conj, about magnitude].

### 4.3 The one diagnostic mesh pattern [Proof]

**Proposition 4.** Suppose the conditional-independence model of the observed variables (taxa +
E) is faithful to an undirected graph U that contains a chordless cycle of length k ≥ 4. Then no
DAG with latent variables, without selection, whose observed independence model is faithful to
its MAG, generates it.

*Proof.* Under U-faithfulness, u and v have a separating set iff they are non-adjacent in U. So
the MAG's adjacencies (pairs with no m-separating set) coincide with U, and the cycle is chordless
in the MAG.

1. Take consecutive vertices a, b, c on the cycle. Since k ≥ 4, a and c are non-adjacent, so they
   have a separating set. In U, every separating set of a and c contains b, because a - b - c is
   a path.
2. In a MAG, an unshielded triple whose middle vertex is in a separating set is a non-collider.
   Otherwise a collider at b with b conditioned would leave a - b - c m-connecting.
3. Without selection, a MAG has only → and ↔ edges [Std: Richardson & Spirtes 2002]. A
   non-collider at b therefore needs a tail at b on one of its two cycle edges.
4. Each of the k cycle edges carries at most one tail, and each of the k vertices needs one. So
   every cycle edge is directed and each vertex is the tail of exactly one of them. That is a
   directed cycle, which contradicts ancestrality. ∎

**Scope of Proposition 4 (revised after review, §10).** The premise is faithfulness of the
**full** independence model to U. The premise is not "U is the concentration graph". The
concentration graph only records the full-order (in)dependences, so a chordless cycle in it says
nothing about the lower-order statements the proof uses.

Astra's counterexample (verified in `c4_counterexample.out`):

- Setup: four independent roots R_i and children C_i = 0.6 R_i + 0.7 R_(i+1) + e_i (indices mod 4,
  unit variances).
- The root block of the 8-variable precision has diagonal 1.85, adjacent entries 0.42 and both
  chords exactly 0. The induced concentration graph on the roots is a chordless C4.
- Yet the roots are marginally independent, which U does not imply. The PC oracle finds **no**
  R-R adjacency (every R-R separating set is ∅) and a v-structure at every C_i.

The operational version therefore works on the causal skeleton:

**Proposition 4′** [Proof; same steps 2-4]. Assume acyclicity with latent variables, no
selection, and faithfulness to the MAG. Let K be the skeleton obtained by separating-set search
(u, v adjacent iff no subset of the measured variables separates them; this is the MAG skeleton).
Then K contains no chordless cycle of length ≥ 4 on which every unshielded triple a - b - c has b
in a separating set of a and c (a definite noncollider).

*Proof.* Each such b is a non-collider in the MAG (step 2). Steps 3-4 then force a directed cycle.
∎

Hence a K-cycle with separator evidence at every triple is evidence of feedback, selection or
unfaithfulness. A cycle whose triples are colliders, or whose separator evidence is missing, is
not.

Remark: linear feedback produces such a pattern (a one-way ring gives a chordless C_k; §4.1).
Selection produces it too: ancestral graphs with undirected edges represent undirected models.

**What "mesh" licenses** (ranked from cheap to expensive):

1. **Concentration-graph mesh** licenses "not a tree of pairwise-mediated dependence among the
   measured variables". It is compatible with fan-in (moralisation), latent E, closure,
   dichotomisation and feedback. It is *no* evidence of cyclic exchange.
2. **Skeleton-level mesh (PC/FCI over all conditioning sets, flexible E)** licenses "at least one
   pair cannot be separated by any measured set". Latent confounding still explains it.
3. **A chordless cycle of length ≥ 4 in the causal skeleton K, with separator evidence that
   every unshielded triple on it is a noncollider** (Prop. 4′) is evidence for feedback **or**
   selection. Bounding the chords in the concentration graph does *not* qualify: that is only a
   "concentration cycle", with no causal reading (counterexample above). Before calling feedback, rule out the selection
   mechanisms: detection thresholds, presence-conditioned sampling (positivity filters),
   relative-abundance normalisation.
4. **Reciprocal exchange** (a ring of mutual cross-feeding) needs orientation evidence in
   *both* directions (genomes) plus perturbation or culture. Observational conditional
   independence cannot license it, particularly in the nonlinear case (§4.2).

---

## 5. Confounding by E, compositionality, and sample size

### 5.1 Habitat filtering (`confounding_demo.out`, 5 taxa with **no interactions**, 200 reps)

Label = class of the subgraph on the non-isolated taxa / number of such taxa. Bonferroni α = 0.05
over the 10 pairs.

| scenario (generator) | adjusted for | N = 100 | N = 400 | N = 2000 |
|---|---|---|---|---|
| S1 linear responses, taxon 0 tracks E tightly (loading 0.95, others 0.6) | nothing | star/5 0.47, star/4 0.35 | **star/5 0.81** | modular 0.71, star 0.19 |
| S1 | E | none 0.96 | none 0.94 | none 0.96 |
| S1, E measured with reliability 0.7 | E_obs | 2-edge 0.36, star/4 0.14 | **star/5 0.89** | modular 0.61, star 0.35 |
| S2 unimodal (Gaussian) niches, optima −2..2 | nothing | modular 0.70, mesh 0.22 | **mesh 1.00** | mesh 1.00 |
| S2 | E linear (wrong form) | path/5 0.21, mixed | **mesh 1.00** | mesh 1.00 |
| S2 | E, E² (correct form) | none 0.95 | none 0.95 | none 0.96 |
| S3 two gradients, E2 unmeasured | E1 | mixed | **mesh 0.94** | mesh 1.00 |
| S5 succession: taxon i tracks an AR(1) E at lag 2i | nothing | **path/5 0.80** | path 0.73 | mesh 0.51, modular 0.41 |
| S5 | current E(t) | **path/4 0.76** | **path/4 0.81** | modular 0.43, mesh 0.40 |
| S4 presence/absence of S1 (threshold) | E linear | none 0.86 | none 0.76 | 1-edge 0.40, 2-edge 0.25 |

What conditioning on measured E does and does not fix:

- **Fixes:** common linear responses to a correctly measured E (S1). Nonlinear responses are
  fixed only if the adjustment has the right functional form (S2 with E²): a spline or
  generalised additive basis, not a linear term.
- **Does not fix:**
  - measurement error in E. The taxon that best proxies E becomes a spurious hub (S1, E_obs).
    This is the star mechanism, and the biggest practical danger, because one abundant
    generalist often tracks temperature or salinity better than the sensor log does.
  - unmeasured gradients (S3: mesh);
  - lagged tracking of a time-varying E (S5). Conditioning on E(t) leaves a spurious chain
    *ordered by phenology*, which looks exactly like a "sequential handoff" (the P5 chain
    prediction).
  - linear residualisation of binary data on E (S4 at large N), because the probit link is
    nonlinear.
- The spurious class often **changes with N** (star → modular, path → mesh). Instability across N
  or subsamples is itself diagnostic of confounding. That supports the team's stability rule.
- A covariate that is a *product* of the taxa (nutrients, O₂, pH, chlorophyll) is a collider or
  mediator. Conditioning on it opens paths and *creates* edges. E must be restricted to
  upstream drivers [Std: collider bias, Spirtes et al. 2000]. This needs a per-covariate ruling
  before any protocol is frozen.

### 5.2 Compositionality (true network = path on 5 focal taxa, |pcor| 0.3, N = 1000, 200 reps)

| background taxa K | transform | class distribution |
|---|---|---|
| 0 | log absolute | path 0.98 |
| 0 | log relative *within the focal set* | **mesh 1.00** |
| 5 | log relative (whole community) | mesh 0.36, path 0.34, modular 0.29 |
| 5 | clr over the whole community (10 taxa) | **mesh 0.99** |
| 50 | log relative / clr (whole community) | path 0.96 / 0.96 |
| 50 | log relative within the focal set | **mesh 1.00** |

Closure acts as conditioning on a common effect (the total), which is a selection mechanism in
the sense of §4.3. With few focal taxa, **never renormalise within the plexitype**. Use
whole-community clr (in our simulation it was safe at ≈ 55 taxa) or an absolute-abundance proxy
(spike-ins, cell counts, qPCR). [Std background: Aitchison 1986; Friedman & Alm 2012; Kurtz et
al. 2015]

### 5.3 Sample size for 80 % class accuracy (`power_summary.out`)

> **Scope (revised after review, §10).** These numbers score class recovery from
> significance-selected edges of full-order partial correlations in i.i.d. Gaussian graphical
> models (and their censored/binary versions). They do not exercise the equivalence tests,
> unresolved-edge completions, blocked resampling or PC/FCI of §6. Read them as a **planning benchmark for the naive pipeline only**. No inequality makes them a
> lower bound for PC/FCI or for the full rule of §6. They do not certify that rule, and they give no
> threshold for any particular dataset.

Truth: a Gaussian graphical model whose concentration graph is the target. Effect size =
|partial correlation| on every edge. Estimator: full-order partial correlation, Fisher z,
Bonferroni α = 0.05 over n(n−1)/2 pairs. Accuracy = classify(estimate) == true class on all n
taxa. 200 reps per cell, N grid 25-6400 log-interpolated. "falls" = reached 0.8, then dropped
below it at larger N (asymptotic bias). Variants:

- **zi**: 40 % zeros from left-truncating the latent Gaussian at prevalence 0.6, present values
  shifted to leave a log1p-like gap, naive Pearson;
- **pa_naive**: presence/absence with phi partial correlations;
- **pa_tet**: presence/absence with tetrachoric partial correlations and an approximate
  variance inflation.

Graphs: star K_{1,n−1}; path P_n; branched tree = spider; mesh = chordless ring C_n (the hardest
2-connected case); modular = two cycles sharing a cut vertex.

**gaussian, |partial corr| = 0.2, N for 80 % class accuracy [accuracy at N = 6400]**
| n | star | path | tree | mesh | modular |
|---|---|---|---|---|---|
| 4 | 450 [0.99] | 494 [0.97] | - | 515 [1.00] | - |
| 5 | 603 [0.96] | 595 [0.96] | 554 [0.96] | 584 [1.00] | 370 [0.96] |
| 6 | 649 [0.96] | 660 [0.96] | 660 [0.97] | 636 [1.00] | 471 [0.98] |
| 7 | 684 [0.97] | 716 [0.95] | 692 [0.96] | 680 [1.00] | 566 [0.97] |
| 8 | 695 [0.97] | 713 [0.97] | 712 [0.96] | 695 [1.00] | 634 [0.99] |

**gaussian, |partial corr| = 0.3, N for 80 % class accuracy [accuracy at N = 6400]**
| n | star | path | tree | mesh | modular |
|---|---|---|---|---|---|
| 4 | 190 [0.99] | 192 [0.98] | - | 190 [1.00] | - |
| 5 | 218 [0.96] | 232 [0.96] | 247 [0.98] | 276 [1.00] | 173 [0.98] |
| 6 | 278 [0.97] | 297 [0.95] | 284 [0.98] | 286 [1.00] | 187 [0.97] |
| 7 | 324 [0.98] | 312 [0.96] | 316 [0.94] | 316 [1.00] | 224 [0.98] |
| 8 | 330 [0.97] | 328 [0.95] | 326 [0.97] | 328 [1.00] | 241 [0.96] |

**zi, |partial corr| = 0.2, N for 80 % class accuracy [accuracy at N = 6400]**
| n | star | path | tree | mesh | modular |
|---|---|---|---|---|---|
| 4 | 630 [0.96] | 624 [0.98] | - | 648 [1.00] | - |
| 6 | 744 [0.91] | 763 [0.97] | 738 [0.96] | 748 [1.00] | 661 [0.97] |
| 8 | 971 [0.92] | 883 [0.95] | 940 [0.95] | 957 [1.00] | 699 [0.98] |

**zi, |partial corr| = 0.3, N for 80 % class accuracy [accuracy at N = 6400]**
| n | star | path | tree | mesh | modular |
|---|---|---|---|---|---|
| 4 | 287 [0.80] | 271 [0.91] | - | 297 [1.00] | - |
| 6 | 359 [0.54] falls | 352 [0.88] | 342 [0.90] | 343 [1.00] | 291 [0.71] falls |
| 8 | 398 [0.35] falls | 371 [0.92] | 381 [0.83] | 379 [1.00] | 340 [0.83] |

**pa_naive, |partial corr| = 0.2, N for 80 % class accuracy [accuracy at N = 6400]**
| n | star | path | tree | mesh | modular |
|---|---|---|---|---|---|
| 4 | 1349 [0.90] | 1231 [0.93] | - | 1329 [1.00] | - |
| 6 | 1568 [0.74] falls | 1506 [0.89] | 1444 [0.89] | 1452 [1.00] | 1293 [0.94] |
| 8 | 1780 [0.72] falls | 1789 [0.94] | 1592 [0.89] | 1539 [1.00] | 1398 [0.97] |

**pa_naive, |partial corr| = 0.3, N for 80 % class accuracy [accuracy at N = 6400]**
| n | star | path | tree | mesh | modular |
|---|---|---|---|---|---|
| 4 | 589 [0.41] falls | 528 [0.69] falls | - | 561 [1.00] | - |
| 6 | 774 [0.03] falls | 697 [0.51] falls | 728 [0.36] falls | 670 [1.00] | 636 [0.53] falls |
| 8 | >6400 [0.00] | 761 [0.56] falls | 760 [0.42] falls | 710 [1.00] | 664 [0.33] falls |

**pa_tet, |partial corr| = 0.2, N for 80 % class accuracy [accuracy at N = 6400]**
| n | star | path | tree | mesh | modular |
|---|---|---|---|---|---|
| 4 | 1296 [0.96] | 1310 [0.96] | - | 1315 [1.00] | - |
| 6 | 1632 [0.97] | 1470 [0.95] | 1526 [0.94] | 1474 [1.00] | 1328 [0.96] |
| 8 | 1945 [0.94] | 1755 [0.94] | 1859 [0.95] | 1822 [1.00] | 1425 [0.96] |

**pa_tet, |partial corr| = 0.3, N for 80 % class accuracy [accuracy at N = 6400]**
| n | star | path | tree | mesh | modular |
|---|---|---|---|---|---|
| 4 | 615 [0.96] | 633 [0.94] | - | 607 [1.00] | - |
| 6 | 795 [0.91] | 758 [0.88] | 735 [0.90] | 738 [1.00] | 668 [0.94] |
| 8 | 1357 [0.84] | 978 [0.91] | 985 [0.87] | 755 [1.00] | 714 [0.91] |

power_sim.tsv sha256 966ed4e510f9.

Reading the table:

- **Edges must be detected**, so every class except the empty graph needs every true edge at
  power ≈ (0.8)^(1/m). False edges are controlled at about 5 % familywise. Accuracy therefore
  saturates near 0.95, not 1.
- These are *independent* samples. Serial samples from one campaign count for much less (S5
  shows what autocorrelation does). No dataset-specific threshold is claimed. The CMO RNA series
  stays descriptive on independent-pool grounds anyway (Astra 248914574040285184).
- Zero inflation and dichotomisation are not just a loss of power. Naive Pearson or phi partial
  correlations on censored or binary data are **asymptotically biased**: X ⫫ Y | Z does not imply
  X ⫫ Y | 1[Z > τ]. Accuracy rises, then falls, and the spurious edges cluster around hubs. A
  latent-Gaussian (tetrachoric, truncated-copula) estimator removes the bias, at roughly 3-4×
  the Gaussian sample size.
- Under the three-valued rule (§6), which certifies *absent* edges with a margin δ, calling a tree
  needs roughly 1000 (δ = 0.15) to 3000 (δ = 0.10) samples at n = 6-8. Mesh calls come at ≈ 300.
  That asymmetry is exactly backwards from what the biology needs (§2.4), so a cheap mesh call
  must never be read causally.

---

## 6. Decision rule: supported / falsified / unresolved

The rule is consistent with "call only when stable across edge sets" (Converging item 2 and P5
caveat). Everything in **bold** is a pre-registered constant to freeze in the protocol.

**Inputs.**

- Taxon set V with n ≥ 4 for undirected class calls; at n = 3 only the collider test of §2.1 is
  available.
- Alternative edge sets G_1..G_K from the curated KO rules, each with genome orientations and
  AND/OR roles.
- Specifications s ∈ S: conditional-independence test family; E basis, **restricted to upstream
  drivers, with a spline basis of degree ≥ 2**; transform (whole-community clr, absolute proxy if
  available); α ∈ **{0.01, 0.05}**.
- Block bootstrap over **sampling dates/campaigns** (B = **200**).

**Step 0. Gates.** Any failure means *unresolved*, with the reason recorded.

1. **Power gate (design-specific calibration only):** N_eff counts independent blocks, not rows.
   - Before unblinding, the protocol freezes a simulation of the **full** rule under its own
     design: n, block structure, data type, ρ_min = **0.2**, equivalence tests, unresolved-edge
     completions, blocked resampling, PC/FCI and all edge sets.
   - The gate passes iff, in that simulation, the rule's correct-call rate for the contrast is at
     least **0.8** and its wrong-call rate is at most **0.05**.
   - The §5.3 table is a planning benchmark only. It belongs to a different estimator, no
     inequality links the two, and it never triggers a verdict by itself.
2. **Positivity gate:** each taxon is present and absent in ≥ **10** samples within each E
   tercile.
3. **Proxy gate:** no taxon has R² > **0.8** against the E spline. Otherwise conditioning on that
   taxon is near-conditioning on E, and a star centred on it is uninterpretable (S1).
4. **Closure gate:** the focal taxa make up < **20 %** of the clr reference, or absolute
   abundance is available.

**Step 1. Three-valued edges** (per edge set k, spec s, bootstrap b), `decision_rule.py`:

- *present* if the Bonferroni test rejects zero;
- *absent* if the Bonferroni two-one-sided 90 % CI of the partial correlation lies within ±δ,
  with **δ = 0.15**;
- *unresolved* otherwise.

The graph object is the **PC-stable skeleton** (adjacency iff no separating set among all subsets
of V ∪ E). The concentration graph is reported alongside but is used only for tree-type calls,
because of Proposition 2(iii).

**Step 2. Invariant class.** c(k, s, b) = the unique class over all completions of the unresolved
edges. Undefined if completions disagree or more than **12** edges are unresolved.

**Step 3. Stability.** π(c | k, s) = the share of bootstraps with c(k, s, b) = c.

**Supported(c)** iff all of the following hold:

- (a) all gates pass;
- (b) for **every** k and every s, π(c | k, s) ≥ **0.8**;
- (c) any *directed* qualifier ("sink hub", "source hub", "cascade", "nested") is compelled in the
  CPDAG, or in the Meek closure of genome orientations that is **identical in every G_k** and
  does not conflict with a compelled data edge;
- (d) a cyclic qualifier is at most "mesh (non-decomposable)". "Mesh with feedback evidence"
  additionally needs, in every (k, s), all three of:
  - a chordless cycle of length ≥ 4 in the PC/FCI skeleton K, with every chord *absent* by
    separating-set evidence, not only by a full-order equivalence test;
  - for each unshielded triple a - b - c on the cycle, a recorded separating set of a and c
    containing b, stable over bootstraps (conservative-PC style: b in every separating set
    found, not merely in one);
  - a written exclusion of selection mechanisms (§4.3).

  A chordless cycle seen only in the concentration graph is reported as a "concentration
  cycle", with no causal reading.

**Falsified(c_pred)** for a mechanism's predicted class (P5) iff the gates pass, the power gate
passes for the specific contrast, π(c_pred | k, s) ≤ **0.05** for every (k, s), and the
excluding evidence is of the **sound** kind:

- predicted *hub without fan-in* (public-goods source hub) or *cascade*: both say no node has
  ≥ 2 parents. They are falsified by a stable collider at some node of the edge set: the two
  parents independent given E, dependent given E + the node, and the arrowheads confirmed by an
  FCI/PAG run, because latent E can mimic arrowheads. A skeleton-level cycle alone makes the
  verdict *unresolved*, not falsified, because a latent driver can create an adjacency;
- predicted *fan-in hub* (AND-clause) is falsified by a stable *absent* conditional dependence
  between donors given the recipient (no v-structure, §3.3), or by a stable concentration-graph
  tree (Proposition 2(iii) excludes any node with two parents);
- predicted *mesh / reciprocal* is falsified by a stable concentration-graph forest with all
  non-edges *absent* (Proposition 2(iii); for linear feedback also §4.1). Under nonlinear
  feedback this cannot be done observationally (§4.2), and the verdict is *unresolved*.

**Unresolved** otherwise, including: a class that flips across edge sets or specifications; any
gate failure; a genome orientation that conflicts with a compelled CPDAG edge (reported as a
separate finding); or an observed class that is only explainable by moralisation, latent E or
closure.

**Calibration of Steps 1-2 on Gaussian truth** (`decision_rule.out`, one spec, no bootstrap,
|pcor| 0.3, 200 reps). Wrong calls stay at 0-5 % (the Bonferroni false-edge rate) in every cell.
Correct calls need the following N:

| n | class | δ = 0.10 | δ = 0.15 |
|---|---|---|---|
| 4 | star / path | 3000 (0.96 / 0.98); 1000 gives 0.17 | 1000 (0.97 / 0.98) |
| 6 | star / path | 3000 (0.93 / 0.94) | 1000 (0.62 / 0.69); 3000 (0.97 / 0.98) |
| 6 | mesh (C6) | 300 (0.96) | 300 (0.95) |
| 8 | path | 3000 (0.88) | 3000 (0.95); 1000 gives 0.28 |

---

## 7. Proven vs conjectured

| statement | status |
|---|---|
| n = 3: fork ≡ chain; collider identified (given sufficiency, acyclicity, faithfulness, tests given E then E + A) | Std (Verma & Pearl) |
| Concentration graph of a DAG = moral graph | Std (Lauritzen) |
| Prop. 2: M a forest ⇒ M = S and no v-structures; S 2-connected ⇒ M 2-connected | Proof (§2.2) |
| Prop. 3: n of the 2^(n-1) tree orientations have M = S | Proof + enumeration |
| CPDAG orients sink hubs but not source hubs or cascades; one genome edge orients the downstream subtree | Std (Meek) + enumeration |
| gLV equilibrium with growth-rate noise: precision = AᵀD⁻¹A | Proof (§4.1) |
| Linear cyclic SEM obeys d-separation | Std (Spirtes 1995) |
| Nonlinear cyclic: σ-separation; within-SCC structure not identifiable from conditional independence | Std (Forré & Mooij; Bongers et al.) |
| Prop. 4: a chordless cycle ≥ 4 in U, with the *full* independence model faithful to U, is impossible for acyclic + latent + no selection | Proof (§4.3), relying on standard MAG facts; confirmed by Astra. The premise is not checkable from the concentration graph (counterexample) |
| Prop. 4′: no chordless ≥ 4 cycle in the separating-set skeleton with noncollider separator evidence at every triple (same assumptions) | Proof (§4.3) |
| σ-vs-d violations are small in realistic nonlinear rings | Conj (one negative simulation) |
| Habitat-filtering, closure and sample-size numbers | Sim (model-specific) |
| Genome flow direction ≠ abundance causal direction in syntrophy | Conj (mechanistic argument) |
| AND-clause (min-type) donors give weak linear explaining-away | Conj |

---

## 8. Concrete recommendations to the team

1. Reword the Architecture block (l. 386-394):
   - "n = 3: only colliders are identifiable without orientation, *and only by marginal-vs-
     conditional tests, not by a partial-correlation network*."
   - "n ≥ 4: the **skeleton**, recovered by separating-set search, separates the classes. A
     partial-correlation network separates them only in one direction (tree calls are sound,
     mesh calls are not)."
   - Replace the three classes with the six-class partition of §1.
2. P5 caveat (l. 566): "at n = 3 the topology is fixed (P3); only the arrow pattern is
   unresolved".
3. P5's chain prediction for N/S redox handoffs is the class most easily faked by phenological
   succession (S5) and by moralisation of zig-zags. If abundance data are used, require a genome
   orientation of the source end *and* the no-v-structure test.
4. The fan-in (AND) prediction is the most testable observationally (§3.3). The sink star is the
   one compelled shape, so lead the abundance-side test with it.
5. Freeze the E covariate list with a per-covariate "upstream driver vs product" ruling. Never
   renormalise within the focal set.

---

## 9. References (from memory; confidence flagged)

High confidence that the work exists with this venue and year:

- Verma T, Pearl J (1990). Equivalence and synthesis of causal models. *UAI 1990*.
- Meek C (1995). Causal inference and causal explanation with background knowledge. *UAI 1995*.
- Andersson SA, Madigan D, Perlman MD (1997). A characterization of Markov equivalence classes for
  acyclic digraphs. *Annals of Statistics*.
- Chickering DM (2002). Optimal structure identification with greedy search. *JMLR*.
- Spirtes P, Glymour C, Scheines R (2000). *Causation, Prediction, and Search*, 2nd ed. MIT Press.
- Lauritzen SL (1996). *Graphical Models*. Oxford University Press.
- Richardson T, Spirtes P (2002). Ancestral graph Markov models. *Annals of Statistics*.
- Zhang J (2008). On the completeness of orientation rules for causal discovery in the presence
  of latent confounders and selection bias. *Artificial Intelligence*.
- Spirtes P (1995). Directed cyclic graphical representations of feedback models. *UAI 1995*.
- Bongers S, Forré P, Peters J, Mooij JM (2021). Foundations of structural causal models with
  cycles and latent variables. *Annals of Statistics*.
- Colombo D, Maathuis MH (2014). Order-independent constraint-based causal structure learning
  (PC-stable). *JMLR*.
- Uhler C, Raskutti G, Bühlmann P, Yu B (2013). Geometry of the faithfulness assumption.
  *Annals of Statistics*.
- Kalisch M, Bühlmann P (2007). Estimating high-dimensional DAGs with the PC-algorithm. *JMLR*.
- Aitchison J (1986). *The Statistical Analysis of Compositional Data*. Chapman & Hall.
- Friedman J, Alm EJ (2012). Inferring correlation networks from genomic survey data (SparCC).
  *PLoS Computational Biology*.
- Kurtz ZD et al. (2015). Sparse and compositionally robust inference of microbial ecological
  networks (SPIEC-EASI). *PLoS Computational Biology*.
- Blanchet FG, Cazelles K, Gravel D (2020). Co-occurrence is not evidence of ecological
  interactions. *Ecology Letters*.
- Liu H, Lafferty J, Wasserman L (2009). The nonparanormal. *JMLR*.

Medium confidence (exists; title or venue details may be slightly off):

- Forré P, Mooij JM (2017). Markov properties for graphical models with cycles and latent
  variables. *arXiv* 1710.08775.
- Mooij JM, Claassen T (2020). Constraint-based causal discovery using partial ancestral graphs
  in the presence of cycles. *UAI 2020*.
- Richardson T (1996). A discovery algorithm for directed cyclic graphs (CCD). *UAI 1996*.
- Neal RM (2000). On deducing conditional independence from d-separation in causal graphs with
  feedback. *JAIR*.
- Perković E, Kalisch M, Maathuis MH (2017). Interpreting and using CPDAGs with background
  knowledge. *UAI 2017*.
- Fan J, Liu H, Ning Y, Zou H (2017). High dimensional semiparametric latent graphical model for
  mixed data. *JRSS-B* (relevant to the pa_tet variant).
- Yoon G, Carroll RJ, Gaynanova I (2020). Sparse semiparametric canonical correlation analysis for
  data of mixed types. *Biometrika* (truncated latent Gaussian copula for zero-inflated data).
- Ramsey J, Zhang J, Spirtes P (2006). Adjacency-faithfulness and conservative causal inference.
  *UAI 2006*.

Not cited, because I am unsure of the exact reference: a specific paper deriving gLV-equilibrium
precision matrices. The derivation in §4.1 is self-contained.

---

## 10. Revisions after cross-family review (Astra 248969617523261440)

**Review base:** document sha256 a67a98145fa0.

**Confirmed by Astra:**

- the rooted-tree count n/2^(n−1);
- linear-Gaussian moralisation;
- the gLV precision AᵀD⁻¹A (generically non-zero; cancellation caveat added in §4.1);
- the Prop. 4 tail-count proof, under faithfulness of the **full** independence model to U and
  to the selection-free MAG.

**Fix 1: the chord-equivalence gate did not establish Prop. 4's premise.**

- Counterexample: four independent roots R_i, children C_i = 0.6 R_i + 0.7 R_(i+1) + e_i (mod 4,
  unit variances). Verified exactly in `scratch/deep/c4_counterexample.py`
  (sha256 2b1d35bc7557; output `c4_counterexample.out`, ef3f6ee73d3f):
  - root precision block: diagonal 1.85, adjacent entries 0.42, chords 0.0000 (both);
  - root marginal covariance = I;
  - induced concentration graph on the roots = chordless C4;
  - PC oracle: no R-R adjacency (all R-R separating sets = ∅) and v-structures R_i → C_i ← R_(i+1)
    at all four children.
- So a DAG with no feedback, selection or latents gives a chordless concentration C4.
- Changes:
  - §0 item 5 revised;
  - §4.3 gains a scope note and Prop. 4′ (skeleton from separating-set search, with
    noncollider separator evidence at every unshielded triple on the cycle);
  - item 3 of the "What mesh licenses" list (§4.3) revised;
  - §6 Supported (d) now requires PC/FCI-skeleton cycles with recorded, bootstrap-stable
    separating sets containing each middle vertex;
  - a cycle seen only in the concentration graph is reported as a "concentration cycle" with no
    causal reading;
  - §7 table updated.

**Fix 2: scope of the sample-size table.**

- `power_sim.py` scores class recovery from significance-selected edges in i.i.d. Gaussian
  graphical models. It does not exercise equivalence tests, unresolved edges, blocked resampling
  or PC/FCI.
- Changes:
  - §0 item 7 and §5.3 now call the table a planning benchmark for the naive pipeline only;
  - the claim that the CMO RNA series is below a universal threshold ("unresolved by
    construction") is dropped. CMO stays descriptive on independent-pool grounds;
  - §6 power gate: purely a design-specific calibration of the full rule (see follow-up below).

**Code hashes** (first 12 hex characters of sha256):

- unchanged: arch_lib.py 2e5efe1e6f89, population_demos.py cf50a689e8e7, confounding_demo.py
  2414433b79cd, dag_orientation.py 3c51a7cff091, decision_rule.py e7ed82eca55d, power_sim.py
  851eb0b43a45, summarise_power.py 6b65ee7cea69, power_sim.tsv 966ed4e510f9;
- new: c4_counterexample.py 2b1d35bc7557, c4_counterexample.out ef3f6ee73d3f.

`decision_rule.py` currently supports **concentration-graph descriptions only**: full-order
three-valued edges and invariant class. It does not implement the PC/FCI skeleton, the separator
check for Supported (d), CPDAG/Meek orientation or blocked resampling. Its calibration
(`decision_rule.out`) is therefore not a calibration of the full rule.

**Follow-up (Astra 248970131862372352, reviewed on 36364dea7eab).** Astra closed the C4 /
Prop. 4′ concern. One remaining fix: the §6 power gate's automatic "unresolved" below the §5.3
table is removed, because accuracy simulated for a different estimator is only a planning
benchmark and no inequality makes it a lower bound for PC/FCI or for the full rule. The power
gate is now purely a design-specific, pre-registered calibration of the full rule. §0 item 7 and
§5.3 were reworded to match. No code changed.

The sha256 of this revised document is reported in the hand-off message, not here, since it
cannot contain its own hash.

