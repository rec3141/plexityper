# Plexitrophy as a computational object: formalisation, complexity, exact decider, bounds

Author: Claude (Opus, delegated deep-work agent for `claude-plexi`), 2026-11-18. **Status: revised
after cross-family review by Astra (Codex), relay 248969335645061120; see §6.** Further review
still needed (Codex or Gemini). Under AGENTS.md this is one vote. Proofs are marked
**Proof**; claims I have not proved are marked **Conjecture**.

**What I read.** AGENTS.md (full). `plexichat.txt` (grep for plexitrophy, l.87–137, incl. the
dictionary entry at l.135). `discussion/plexi-concept-review.md` main text l.200–260 (Carbon
248915719529873408, Astra 248915827713556480) and l.325–380 (resolved item 5, "Four orthogonal
axes", one-sentence definition). `protocols/P2-plexitrophy-recipient.md` (full, 113 lines, draft
not frozen). `notes/vision-agy-plexi.md` l.25–36 (Astra 248917454944124928 quantifier rule) and
`notes/plexitrophy-and-plexiome.md` l.50–75. I opened nothing under `data/`. I did not re-read
Frioux et al. 2018 or the other literature cited in §2.6 this session; those citations are from
memory and need checking before they go into the manuscript.

**Code** (stdlib + numpy/scipy, run with `python3 -I`). These are the revised versions; the
original hashes are in §6:
- `scratch/deep/plexisolve.py`, sha256 `533a8bcf6245`
- `scratch/deep/test_plexisolve.py`, sha256 `3cc157681850`

Run with `cd scratch/deep && python3 -I test_plexisolve.py 250 --scale`. Result: **53/53 checks pass, 0 fail, in the revised run (log was in /tmp, not kept)**
(see §3.3 and §6).

---

## 1. Formalisation

### 1.1 Instance
An instance is I = (V, M, E, τ, {R_Y, Exp_Y, Need_Y}_{Y∈V}, X):
- V is a finite set of organisms (dereplicated genomes, plus explicit pseudo-organisms such as the
  unbinned pool, §4.4). X ∈ V is the focal recipient.
- M is a finite set of metabolites. τ : M → {anab, catab} gives each one's type. Catabolic
  *services* (an electron-acceptor supply, removal of H₂ or formate) are modelled as
  pseudo-metabolites flowing **from the service provider to the organism served**. The physical
  flow in the opposite direction (H₂ from X to the scavenger) is an input to the provider's
  service rule, which makes classic syntrophy a 2-cycle. This matches the team's exclusion of
  dyadic syntrophy.
- E ⊆ M is the external medium, available to every organism.
- R_Y is a set of Horn rules A → b (A ⊆ M, b ∈ M): given inputs A available inside Y, Y makes b.
  Several rules with the same product are **OR-routes**; the inputs of one rule are an
  **AND-clause**. A route that takes in a finished product (X imports B12) needs no rule: imports
  are wiring, not rules.
- Exp_Y ⊆ M is what Y can export. Need_Y ⊆ M is what Y requires (its growth endpoint).

### 1.2 Wiring
For a set S ⊆ V, a **wiring** W on S is a *well-founded derivation*. It gives every needed pair
(Y, m) (Y ∈ S, m "available in Y") exactly one justification:
1. `ext`: m ∈ E;
2. `rule`: a rule A → m in R_Y, with all (Y, a), a ∈ A, justified; or
3. `import from Z`: Z ∈ S∖{Y}, m ∈ Exp_Z, and (Z, m) justified **by a rule of Z**.

The pairs and justifications must form an acyclic proof DAG, so nothing justifies itself through a
loop. The organism-level graph *may* be cyclic. (S, W) is **feasible** if every Y ∈ S has Need_Y
justified (every member is alive). The **supply graph** G(W) has an arc Z → Y, labelled with the
type set {τ(m)}, whenever Y imports some m from Z.

The definitional choices, and why I made them:
- **A partner may supply several inputs.** It counts once in d. That is the "dual-role partner
  counted once" rule.
- **X may supply itself** through its own rules (remodelling, salvage), and so may any member.
  Self-supply and E-supply count 0 in d.
- **No relaying.** Only a producer can export: the exporter's *single* node (Z, m) must be
  rule-justified, which rules out relaying through an internal futile cycle p → q → p (fixture
  F5, §6). If an organism could pass on a metabolite it had imported, a hub that gathers D1's and D2's products and forwards both would turn a fan-in (d=2)
  into d=1 for free. Remodelling (cobinamide + DMB → B12 in D2, wiring 3a) is a rule, so it is
  allowed. If true uptake-and-release relays matter biologically, they must be added as explicit
  rules with evidence.
- **One justification per (Y, m).** The model is qualitative (no fluxes). Split supply of one
  metabolite is never needed, and Lemma 3 shows it would not change any decision.

### 1.3 Quantities
- d^T_W(X) = |{Z : Z → X in G(W) with a label in T}|. Here T = {anab} gives d_anab, T = {catab}
  gives d_catab, and T = both gives d_total. Each is the in-degree of X in the T-restricted supply
  graph. It counts all imports X uses in W, which is the "{D : (D→X) ∈ E_S}" definition in
  `notes/plexitrophy-and-plexiome.md` l.52.
- m_W(X) = |S|.
- c^T_W(X) = |SCC of X in G^T(W)| (the team's definition).
- κ^T_W(X) = 1 if X lies on a **simple** directed cycle of length ≥ 3 in G^T(W), else 0. This is
  the proposed replacement, see §1.5 I2.

**Counterexample property** (the negation of the per-support disjunction):
Bad^T(W) :⇔ d^T_W(X) ≤ 1 ∧ c^T_W(X) ≤ 2 (team) or ∧ κ^T_W(X) = 0 (proposed). Bad is
**closed under taking subgraphs**: removing arcs cannot raise an in-degree, enlarge an SCC or
create a simple cycle.

### 1.4 Basic lemmas (all proved)
**L1 (Horn semantics; revised, §6).** `closure(S)` is the least fixpoint of a *two-node*
Horn program with separate avail(Y, m) and prod(Y, m) nodes. That is the network-expansion
"scope" generalised to compartments (Dowling–Gallier 1984; Ebenhöh/Handorf scope). It is linear
in the ground rule set.

**L1′.** With no import restrictions, the set of members that can be kept alive is the same under
the two-node closure and under the single-node semantics of §1.2. So feasibility, Lemma L2 and
S_max are unaffected. *Proof:* build an earliest-stage derivation. Justify each node by its
earliest justification. An importer of m takes it from an exporter whose prod(·, m) stage is
minimal. That exporter's own import stage is strictly larger than its production stage, so its
node is rule-justified. Every premise has a strictly smaller stage, so the derivation is acyclic. ∎

L1′ does **not** carry over to *supply graphs* or to restricted imports: see L5 and Theorem 3′
in §6.

**L2 (union closure; unique maximal community).** If S₁ and S₂ are feasible, so is S₁ ∪ S₂.
Therefore every V′ ⊆ V has a unique maximal feasible subset S_max(V′). It is computed by
alternating `closure` with deletion of members that are not viable, and every feasible S ⊆ V′
satisfies S ⊆ S_max(V′). *Proof:* `closure` is monotone in the member set, so S₁'s and S₂'s
members stay viable in the union. Deletion never removes a member of a feasible S ⊆ V′: by
induction S ⊆ current set, so by monotonicity S's members are viable there. ∎

**L3 (wiring minimality is immaterial).** "Minimal wiring" can be read as minimal in arcs or in
justifications. Either way, a counterexample exists among all wirings on S iff one exists among
the minimal wirings on S, because every wiring contains a minimal one with a subgraph supply graph
and Bad is subgraph-closed. ∎ So the team's "minimal support-*and-wiring*" needs no separate
wiring-minimality clause.

**L4 (irredundant core).** Let (S, W) be feasible and let S* = {X} ∪ {Y : Y ⇝ X in G(W)}. Then
(S*, W|S*) is feasible and G(W|S*) ⊆ G(W). *Proof:* every supplier of a member of S* reaches X,
so it is in S*; justifications stay inside S*. ∎ Call a feasible (S, W) with S = S* **irredundant**:
every member contributes to X.

### 1.5 Ambiguities and inconsistencies in the team definitions, with fixes

**I1. "Every feasible *minimal* support" makes the call non-monotone, and it misses chain
communities X can actually live in.** In the team rule (Astra 248917454944124928, review
l.366), minimality is inclusion-minimality of S. A community can be irredundant (every member is
used) without being inclusion-minimal, and the counterexample X needs can live only there.
Adversarial case **A2** (`A2_minimal_vs_irredundant`) shows this:
- X needs p and q.
- B exports both, but B needs r from A, A needs s from C, and C is autonomous. That gives the chain
  C → A → B → X with d=1, c=1.
- A can also export p and C can export q, so {X, A, C} is a fan-in support (d=2), and it is a
  proper subset of {X, A, B, C}.
- The only inclusion-minimal support is therefore the fan-in one, and the team rule calls X
  **plexitrophic**. Yet X lives in a viable community in which it has exactly one direct partner
  and sits on no cycle. The chain is not "redundant": B does real work.
- Remove the single route "A makes p" and the team rule flips to **not plexitrophic**. *Gaining*
  a capability makes X *more* plexitrophic. That non-monotonicity breaks the P2 interval logic
  (§4) and makes the decision coNP-complete (Theorem 2).
- **Fix (P-all semantics):** X is plexitrophic under E iff X is viable and **every feasible wiring
  on every support** satisfies d ≥ 2 ∨ (cycle criterion). By L4 this is the same as quantifying
  over all *irredundant* supports.
- It is strictly stronger than the team rule: Plexi_all ⇒ Plexi_min, because minimal supports are
  a subfamily. It matches the dictionary definition (plexichat l.135: "no single pairwise
  interaction is sufficient") and the one-sentence draft ("no single partner suffices"). It is
  monotone (Theorem 8), and for d_total it is polynomial (Theorem 3).
- On the existing fixtures the two semantics agree. On 750 random instances they never disagreed
  (the A2 structure is rare at random), and Plexi_all ⇒ Plexi_min held throughout.

**I2. "|SCC(X)| ≥ 3" does not mean "X lies on an irreducible cycle of at least three".**
Adversarial case **A3** (`A3_dyad_chain`): X ↔ D ↔ Y, two dyadic mutualisms that share D.
- SCC(X) = {X, D, Y}, so c = 3, and the team rule calls X *cyclic-plexitrophic*.
- But X's only partner is D, and X lies on no cycle longer than 2.
- The same structure is the excluded serial chain Y → D → X with the excluded dyad X ↔ D overlaid.
  Two excluded architectures should not combine into an included one. The one-sentence definition
  says "belongs to an irreducible cycle of at least three", which is κ, not c.
- **Fix:** use κ(X) = [X lies on a simple directed cycle of length ≥ 3]. This is cheap: with
  d ≤ 1 it reduces to "X exchanges with nobody but its one supplier", Theorem 3(b). Keep
  c = |SCC| as a *reported* descriptor if wanted.
- Every team fixture gives the same answer under c and κ (F2: A → B → C → A is a simple
  3-cycle).

**I3. "d counted over required inputs".** Read literally, this could exclude imports X uses
only to make something it exports (X as a processing step for a partner). In an irredundant wiring
every such import is still needed for X's own viability (its product feeds a member that X
depends on). I therefore define d as the full in-degree (§1.3), which agrees with
`notes/plexitrophy-and-plexiome.md` l.52. If the team wants the narrower reading, the counterexample
search changes (some in-arcs become free, like catabolic arcs in Theorem 4) and loses
polynomiality.

**I4. m has two definitions.** Resolved item 5 (review l.334) says "m = the smallest
self-maintaining consortium" (a minimum over supports). The axes table says m_S(X) = |S| per
support. Neither enters the criterion. Report m as the size of the smallest witness: MILP
objective, or minimum support size. Note that minimising m is NP-hard (Theorem 5), unlike the
decision itself.

**I5. "Empty or incompletely enumerated support set ⇒ unresolved."** Under the deciders below
nothing is enumerated. "Unresolved" should come only from uncertain literals (§4) or a solver
time-out. An empty support family is the separate state `no_donor`.

**I6. Typed variants.** c_anab and κ_anab are computed on the anabolic-only subgraph. Because
d_total ≥ d_anab and G^anab ⊆ G, **anabolic plexitrophy ⇒ total plexitrophy**, which fits the
team's "stricter benchmark". Note that d_total ≥ 2 does not require d_anab + d_catab ≥ 2 with
different partners. A dual-role partner (A4, split=False) has d_anab = d_catab = 1 and
d_total = 1, so it is not plexitrophic, while a naive sum would say 2.

**I7. Qualitative model.** Any producer is assumed to supply enough. Quantitative insufficiency,
where two producers of one metabolite are jointly needed, is invisible. That is a modelling limit,
not an inconsistency, and should be stated in P2.

### 1.6 Team fixtures under this formalisation (all PASS; `run_fixtures`)
| Fixture | Encoding | Result |
|---|---|---|
| Serial chain D1 → D2 → X | D1 makes cbi; D2 has rule {cbi, dmb} → b12 and exports b12; X needs b12 | not plexitrophic; witness (d, c) = (1, 1), m = 3 ✓ |
| D → A plus A → B → C → A | A needs a1 (D) and a2 (C); B needs b (A); C needs c (B) | (d, m, c) = A (2,4,3), B (1,4,3), C (1,4,3); all plexitrophic under both c and κ; D not plexitrophic (d=0) ✓ |
| Alternatives (2,1) and (1,3) | {X,P,Q} fan-in, or X → B → A → X | plexitrophic under P-all and minimal semantics ✓ |
| P2 3.2 fixture 3 (3a / 3b / both) | as in P2 | 3a only: not plexitrophic; 3b only: plexitrophic; both: not plexitrophic ✓ |

---

## 2. Computational complexity

Throughout, n = |V| and the instance size is the total size of rules, exports and needs.

### 2.1 Viability is polynomial
**Theorem 1.** X has a feasible support iff X ∈ S_max(V). This can be computed with at most n
Horn fixpoints. *Proof:* L2. ∎

### 2.2 The team's rule as written (inclusion-minimal supports) is coNP-complete
**Theorem 2.** Deciding "X is plexitrophic" under inclusion-minimal semantics is coNP-complete.
*Scope (revised, §6):* "minimal" here means inclusion-minimal in **members only**, as in the
review text (l.366). I make no claim about P2's current scorer, which minimises (members, edges)
jointly; I have not examined it, and neither Theorem 2 nor A2 is attributed to it.
This holds for d_total, with every metabolite of one type, and for either cycle criterion (c or κ).

*Membership.* A certificate of non-plexitrophy is (S, W). Check in polynomial time:
- W is a well-founded feasible wiring;
- Bad(W);
- S is inclusion-minimal. By L2, S is minimal iff X ∉ S_max(S∖{Y}) for every Y ∈ S∖{X}: any
  feasible T ⊊ S with X ∈ T lies inside some S∖{Y}, hence inside its S_max.

Viability is in P (Theorem 1). Union closure is exactly what keeps this at coNP and stops it from
reaching Π₂ᵖ.

*Hardness (reduction from SAT; implemented as an exhaustive check in the A2 spirit, see
`plexi_brute_min`).* Take a CNF φ over x₁..x_n with clauses C₁..C_k. Build the instance:
- **X** needs x. X has rules {u_i, w_i} → x for every i.
- **B** exports x and needs c₁..c_k.
- **Literal organisms** L_i⁺ and L_i⁻ are autonomous. L_i^± exports c_j for every clause C_j that
  contains that literal. L_i⁺ also exports u_i, and L_i⁻ exports w_i.

The minimal supports are as follows:
- {X, L_i⁺, L_i⁻} for each i. These are fan-in supports, d = 2.
- Supports through B. Any such support must contain a set of literals that covers every clause.
  If it contains a complementary pair, then {X, L_i⁺, L_i⁻} is a proper feasible subset and the
  support is not minimal. If the literals are consistent and form a minimal cover, the support is
  minimal, and X has d = 1 (only B supplies X) and c = 1.

So X is not plexitrophic ⇔ φ is satisfiable. ∎

### 2.3 Under P-all semantics, total plexitrophy is polynomial
**Theorem 3 (revised, §6).** *Precondition:* the global metabolite dependency graph is
acyclic. Its arcs are a → p for every rule input a of product p, in any organism. Under that
precondition and P-all semantics, total plexitrophy (d_total) can be decided with O(n) Horn
fixpoint computations, for either cycle criterion. **Without the precondition the theorem fails:
the problem is coNP-complete (Theorem 3′, §6).**
*Running time and bounds (stated per Astra 248973173642936320):* the decider performs at most
2n + 2 calls to S_max (or n + 1 for κ). Each S_max is at most n two-node closures. **It never
enumerates import bundles or supports, so no fixed bound on endpoints or routes is assumed.**
Under the acyclicity precondition, L5 guarantees that each candidate D can jointly export
everything it produces in the two-node closure. X is **not** plexitrophic iff X is viable on E
alone, or for some D ≠ X one of the following holds:
- **(a) team criterion c ≤ 2.**
  - **Case A** (X is a pure sink): D ∈ C := S_max(V∖{X}), and X is viable on
    E ∪ (Prod_C(D) ∩ Exp_D). Here Prod_C(D) is what D produces in closure(C).
  - **Case B** (closed dyad): with U := S_max(V∖{X, D}) and Pool_U := ∪_{Z∈U} Prod_U(Z) ∩ Exp_Z,
    both X and D are viable in the two-organism system {X, D}. In that system D additionally sees
    Pool_U and X imports only from D.
- **(b) proposed criterion κ = 0** (private attachment): X ∈ S_max(V) when imports are restricted
  so that X imports only from D and only D may import from X.

*Proof.* (⇐) In each case the construction is a well-founded wiring. In Case A, X has
out-degree 0, so c = 1. In Case B, nothing flows from {X, D} into U, so SCC(X) ⊆ {X, D}. In (b),
X's only neighbour in either direction is D, so the only cycle through X is X → D → X. In every
case d ≤ 1.

(⇒) Take a counterexample. By L4 it is irredundant: every member reaches X.
- If d = 0, X's derivation uses only E and its own rules.
- Otherwise X has a unique supplier D. Every path into X ends with D → X.
- *(a) c ≤ 2.*
  - If X does not reach D, X can have no out-neighbour u: u reaches X, necessarily through D, so
    X would reach D. Then S∖{X} is feasible on its own, so S∖{X} ⊆ C (L2). By monotonicity of
    closure, D produces in C at least what it produced in W, which gives Case A.
  - If X reaches D, then D ∈ SCC(X) = {X, D}. Any u ∈ U = S∖{X, D} that received from X or D
    would be reachable from X and would reach X, so it would be in the SCC. Hence U is feasible
    alone, U ⊆ S_max(V∖{X, D}), and Pool only grows. This gives Case B.
- *(b) κ = 0.* If X had an out-neighbour u ≠ D, take u's path to X and stop at its first arrival
  at D. That path avoids X, so X → u ⇝ D → X is a simple cycle of length ≥ 3. Hence X exchanges
  only with D, which is exactly the restricted system. ∎

The threshold 2 is what makes this easy. Deciding "∃ wiring with d ≤ k" is solvable in n^O(k)
time for fixed k but W[2]-hard in k (set cover, Theorem 5). Plexitrophy only ever asks k = 1.

### 2.4 The anabolic variant is coNP-complete even under P-all semantics
**Theorem 4.** Deciding anabolic plexitrophy (d_anab ≥ 2 ∨ cycle criterion on G^anab, for every
feasible wiring) is coNP-complete, for both cycle criteria.

*Membership:* a counterexample wiring is a polynomial certificate.

*Hardness (implemented as `A5_sat_gadget`, Theorem 5 in the code comment):* take a CNF φ.
- **X** needs anabolic x and catabolic g_i for every variable i. X exports anabolic h_i.
- **D** is the only maker of x and needs anabolic c_j for every clause.
- **Literal organisms** L_i^±: L_i^± makes g_i from h_i (rule {h_i} → g_i), and exports c_j for
  the clauses it satisfies.

Every feasible wiring imports g_i from some literal, and that literal must receive h_i from X: it
is "tainted". D can take c_j only from untainted literals. Otherwise
X → L (h_i) → D (c_j) → X (x) is an anabolic simple 3-cycle, which also puts L in SCC(X). For each
i at least one of L_i⁺, L_i⁻ is tainted, so the untainted literals are consistent. Hence a
counterexample exists ⇔ φ is satisfiable.

The catabolic arcs L → X are "free" (they count in neither d_anab nor G^anab). That freedom is
the whole source of hardness, and why d_total escapes it (the same gadget always has d_total ≥ 2).

Checked: the brute-force oracle and the MILP agree on 2 hand-made gadgets (SAT and UNSAT) under
both criteria, and the MILP agrees with a brute-force SAT solver on 40 random gadgets (31 SAT,
9 UNSAT). ∎

**Proposition 4′ (tractable special case).** If nobody can use X's exports (X is a pure sink,
the usual P2 recipient), there is no cycle through X. Anabolic non-plexitrophy is then: X is
viable in S_max(V∖{X}) with anabolic imports restricted to a single D, for some D. That is
n fixpoints. *Proof:* as Case A of Theorem 3, with catabolic arcs unrestricted. ∎ Like Theorem 3, Proposition 4′
also needs an acyclic dependency graph (L5). Using this
shortcut means removing X's exports, which is a ⊑-decrease (§4). It is exact only if X's exports
really are unusable. Otherwise use the MILP, which takes seconds at P2 scale (§3.3).

### 2.5 Optimisation and enumeration
**Theorem 5 (minimum d or m is NP-hard and inapproximable).** Let X need u₁..u_n, and let every
other organism be autonomous and supply a subset of the u's. Then min over wirings of d_total(X)
equals the minimum set cover, and min m equals 1 + minimum set cover. Both are NP-hard, inherit
set cover's (1−ε) ln n inapproximability (Feige 1998; Dinur–Steurer 2014), and are W[2]-hard
parameterised by the bound. So "d = the minimal AND-clause of simultaneous donors" (resolved
item 5) is set-cover-like as an *optimisation*. The plexitrophy *decision* never needs it. ∎

**Theorem 6 (enumerating minimal supports).**
- (i) Even with all donors autonomous and X a pure sink, the minimal supports minus X are exactly
  the **minimal transversals** of the hypergraph {producers(u) : u ∈ Need_X∖E}. Enumerating them
  is therefore at least as hard as monotone dualization, for which the best known algorithm is
  quasi-polynomial (Fredman–Khachiyan 1996) and output-polynomial status is open.
- (ii) In general, minimal supports cannot be enumerated in output-polynomial time unless P = NP.
  In the reduction of Theorem 2, the n supports {X, L_i⁺, L_i⁻} are listable in polynomial time,
  and an extra minimal support exists iff φ is satisfiable. An output-polynomial enumerator, run
  for its polynomial budget on n outputs, would decide SAT. ∎

This is why Frioux et al.'s community counts (team notes: median 977 equivalent minimal
communities) are not a nuisance one can enumerate past: the support family can be exponential,
and an exact rule should never require it.

### 2.6 Relation to known results
- **Horn / AND-OR reachability, network expansion and scope.** Single-organism or pooled-community
  producibility is Horn least-fixpoint reachability: linear time (Dowling & Gallier 1984), the
  "scope" of Ebenhöh/Handorf/Heinrich. Our feasibility is a compartmentalised version wrapped in a
  greatest fixpoint (L2): still P.
- **MiSCoTo / ASP (Frioux et al. 2018, Bioinformatics 34:i934).** As I recall it (not re-read
  here), MiSCoTo selects **minimum-cardinality** communities that make host targets producible,
  and enumerates all optimal communities with ASP. Two differences from our object:
  - Minimum cardinality (Theorem 5, NP-hard) is not inclusion-minimality (Theorem 2).
  - Our rule quantifies over *all* supports, which makes it a coNP-type question rather than an
    NP optimisation.

  MiSCoTo's ASP machinery could encode our counterexample search (it is an NP query), but no
  installs without grid, so I used scipy's MILP.
- **Minimal hitting sets.** See Theorem 6(i). The coordinatewise-minimum shortcut is **sound but
  incomplete**: if min_S d ≥ 2 or min_S c ≥ 3 then every support qualifies. It only produces false
  negatives, as in A1 and F3.
- **Σ₂/Π₂.** The ∀(support) ∀(wiring) form looks Π₂ ("∀S: S not minimal ∨ ...", where "not
  minimal" is ∃). Union closure (L2) makes minimality polynomial and collapses the problem to coNP.
  **Conjecture:** if the model gains *non-monotone* constraints, union closure fails and the
  minimal-semantics rule becomes Π₂ᵖ-complete. Examples of such constraints are mutually exclusive
  O₂ regimes inside one support, toxicity, or competition for one donor's export. Viability itself
  becomes NP-complete under such constraints (easy SAT reduction), and P-all semantics becomes
  Σ₂-free but NP-hard to evaluate. Not proved here.

### 2.7 Summary
| Question | Complexity |
|---|---|
| X viable (feasible support exists) | P (≤ n Horn fixpoints) |
| Total plexitrophy, **P-all**, c or κ, acyclic dependency graph | **P**: O(n) fixpoints (Theorem 3) |
| Total plexitrophy, **P-all**, c or κ, cyclic dependencies allowed | **coNP-complete** (Theorem 3′, §6) |
| Total plexitrophy, team inclusion-minimal semantics | **coNP-complete** (Theorem 2) |
| Anabolic plexitrophy, P-all, c or κ | **coNP-complete** (Theorem 4); P if X is a pure sink (Proposition 4′) |
| min d, min m over supports | NP-hard, ln n-inapproximable, W[2]-hard (Theorem 5) |
| Enumerate minimal supports | ⊇ hypergraph dualization; not output-polynomial unless P = NP (Theorem 6) |
| With non-monotone (incompatibility) constraints | Viability NP-complete; minimal semantics Π₂ᵖ (conjecture) |

---

## 3. Exact deciders usable at P2 scale

The principle is to **decide plexitrophy by searching for a witness of non-plexitrophy**: a
feasible wiring with d ≤ 1 and no long cycle. Never enumerate supports.

### 3.1 Polynomial decider for total plexitrophy (`nonplexi_poly`)
Theorem 3 *is* the algorithm. It is valid only when `dep_acyclic(inst)` holds, and the code
refuses otherwise (use the MILP). Per recipient and per E:
1. `S_max` once for viability.
2. Then, for each candidate supplier D (only organisms that export something X uses), one or two
   `S_max` calls.

Correctness is Theorem 3. The returned witness ('d0', 'A'/'B' with D, or 'P' with D) says which
structure defeats plexitrophy.

### 3.2 Exact MILP for any typing or criterion (`nonplexi_milp`)
Variables (after pruning):
- s_Y: Y alive.
- a_{Y,m}: m available in Y.
- j_{Y,k}: rule k fires.
- x_{Z,Y,m}: Y imports m from Z.
- t_{Y,m} ∈ [0, N]: proof depth.
- h_{Z,Y}: arc in the relevant typed graph.
- r_Y ∈ [0, 1]: taint (κ) or reach (c).

Constraints:
- s_X = 1.
- s_Y ≤ a_{Y,m} for m ∈ Need_Y∖E.
- a_{Y,m} ≤ Σ(j for rules making m) + Σ_Z x_{Z,Y,m}.
- j_{Y,k} ≤ a_{Y,i} for each input i, and j ≤ s_Y.
- x_{Z,Y,m} ≤ Σ_k j_{Z,k} with product m (no relaying); x ≤ s_Z, s_Y.
- Well-foundedness, N = #nodes + 1:
  t_{Y,m} ≥ t_{Y,i} + 1 − N(1 − j) and t_{Y,m} ≥ t_{Z,m} + 1 − N(1 − x).
- h_{Z,Y} ≥ x_{Z,Y,m} for τ(m) ∈ T.
- d ≤ 1: Σ_Z h_{Z,X} ≤ 1.
- κ = 0: r_Y ≥ h_{X,Y} − h_{Y,X}; r_Y ≥ r_Z + h_{Z,Y} − 1 (Z ≠ X); r_D + h_{D,X} ≤ 1.
- c ≤ 2: r_Y ≥ h_{X,Y}; same propagation; h_{V,D} + r_V + h_{D,X} ≤ 2 for V ∉ {X, D}.
- Objective: min Σ s_Y. That gives the smallest witness, so m is reported for free.

Pruning, both safe by L2 and subgraph closure:
- (i) restrict V to S_max(V);
- (ii) drop imports of m into Y when Y can make m from E alone (switching to self-supply deletes an
  arc).

**Theorem 7.** The MILP is feasible ⇔ a counterexample wiring exists (P-all semantics).

*Proof.*
- (⇐) Set all variables from the wiring, with t = proof depth and r = exact taint or reach sets.
  The constraints hold by the arguments of Theorem 3: a tainted D gives a simple cycle of
  length ≥ 3; a reached V with V → D puts V in SCC(X).
- (⇒) t strictly increases along every active justification, so the active justifications form a
  well-founded derivation in which all needs are met. h dominates the arcs used. r dominates the
  true taint/reach indicators, by induction along paths. So the constraints forbid every violating
  path, and Bad holds. ∎

Every returned witness is **self-certified**: it is replayed through the Horn fixpoint and the
graph predicate (`assert` in the code).

For the team's inclusion-minimal semantics I give no MILP. Minimality is a coNP side condition
(Theorem 2), so only an exponential check (`plexi_brute_min`) is provided. That is one more
reason to adopt P-all.

### 3.3 Implementation and evidence (`test_plexisolve.py`; seeds in code)
- **Three deciders (revised, §6).** `nonplexi_brute_all` is an exhaustive search over typed arc
  sets G *and* rule/import choices. It uses the single-node fixpoint `dag_smax`, not `closure`
  and not Theorem 3, plus a DFS for simple cycles. `nonplexi_milp` and
  `nonplexi_poly` complete the trio.
- **Agreement.** All agree on 750 random small instances (3 generators: 250 each; per-generator
  tallies for (total, c) include 222 plexitrophic and 373 viable overall), in all 4 variants
  ({total, anab} × {c, κ}). That is 3000/3000 comparisons. *(Revised run: 77 instances have
  cyclic dependencies; there the poly decider is skipped and the brute force and the MILP agree.)*
  The P-all and minimal semantics never disagreed on these instances, and Plexi_all ⇒ Plexi_min
  held in all of them.
- **Fixtures.** All team fixtures pass (§1.6).
- **Adversarial cases.** All PASS:
  - **A1, coordinatewise minimum is wrong.** Three minimal supports, (d, c, m) = (2,1,3),
    (1,3,3) and (1,4,4). The triangle partner also receives an export from X, so it *looks*
    dyadic. The truth is plexitrophic under both semantics. Coordinatewise min = (1,1) gives
    "not plexitrophic".
  - **A2, minimal vs irredundant.** The team semantics says plexitrophic and P-all says not.
    Deleting one route flips the team semantics (non-monotone).
  - **A3, dyad chain.** c-criterion: plexitrophic. κ-criterion: not.
  - **A4, dual-role partner.** d_total = 1, so not plexitrophic; the split version (two partners)
    is total-plexitrophic but not anabolic-plexitrophic.
  - **A5, SAT gadget (Theorem 4).** Agrees on SAT and UNSAT.
- **Scale (synthetic only, no project data).** 200 genomes × 30 endpoint metabolites; the
  recipient is missing 10 endpoints, with 1–2 alternative precursor routes on 5 of them; about 8%
  remodelling-type routes. *Revised:* the precursors are now drawn so that the dependency graph
  is acyclic and Theorem 3 applies. The earlier instance was cyclic and the revised poly decider
  refuses it. Revised run (shared, busy machine):
  - MILP (HiGHS via scipy 1.17): 16–23 s per total variant and 3.2–3.3 s per anabolic variant,
    all without a counterexample.
  - The earlier timing, 6.2–7.4 s on the old instance with a single complete helper added
    (witness = that helper), was not re-run.
  - The naive pure-Python poly decider takes 34–42 s (11–22 s in the first run), because it recomputes fixpoints and is not
  counter-based. A Dowling–Gallier implementation would make it sub-second. Either is fine for
  hundreds of recipients × a few E scenarios on one core.

---

## 4. Monotonicity and conservative interval logic

### 4.1 Order on instances
Write I ⊑ I′ if I′ has:
- a superset of the organisms,
- for each organism, a superset of its rules and exports,
- E ⊆ E′,
- for each organism, a **subset** of its needs.

Types are unchanged. "More capability, fewer requirements" is upward.

**Theorem 8 (monotone quantities).**
- Viable(X) and, under P-all semantics, NonPlexi^T(X) are monotone non-decreasing along ⊑. This
  holds for every T and both cycle criteria. *Proof:* a wiring of I is a wiring of I′ with the same
  supply graph. ∎
- The minimum over wirings of d^T, of c^T and of m is non-increasing along ⊑.
- For one fixed wiring, adding arcs can only raise d, c and κ (subgraph closure).
- **Not monotone:**
  - (i) The inclusion-minimal semantics. It is monotone when *organisms* are added (minimality of
    S is intrinsic to S), but not when routes, exports or E are added, or needs removed: see A2.
  - (ii) Splitting or merging bins. These are not ⊑-comparable.
    - Splitting one genome into two partial bins can turn d = 1 into d = 2 (P2's Trap 9).
    - Merging two genomes (a chimera, or 95%-ANI over-dereplication of complementary strains)
      can turn d = 2 into d = 1. It can also *create* a cycle: X → A → D1 and D2 → X become
      X → A → D → X when D1 and D2 are merged into D.

### 4.2 Two-point identification (proved, from Theorem 8)
Classify every uncertain literal: organism presence, a route or KO, an export, a need of X or of a
donor. I⁻ takes every uncertain capability absent and every uncertain need present. I⁺ takes the
reverse. If the truth satisfies I⁻ ⊑ I_true ⊑ I⁺, then with **two solver calls** per (X, E):

| State | Condition | Meaning |
|---|---|---|
| plexi (robust) | Viable(I⁻) ∧ ¬NonPlexi(I⁺) | plexitrophic in every completion |
| not_plexi, d=0 ambient(E) | X viable alone on E in I⁻ | autonomy |
| not_plexi, d≤1 (robust) | NonPlexi(I⁻), witness has a supplier | a single-partner wiring certainly exists |
| no_donor (robust) | ¬Viable(I⁺) | no feasible support in any completion |
| unresolved (U) | otherwise | robust state undetermined |

The states are mutually exclusive by Theorem 8. *(Revised, §6.)* If I⁻ and I⁺ are themselves
admissible completions, each robust state is exact: it holds iff **all** completions agree on
it. **U means only that completions disagree.** It does *not* mean that both "plexitrophic" and
"not plexitrophic" are attainable. Astra's example: the only donor's presence is uncertain, so
the completions are no_donor or d = 1, and plexitrophy is unattainable, yet the recipient is in
U. The bounds below remain valid; the upper one can be loose. With literals shared across recipients (one gappy donor
serves many), the Manski interval [N_plexi/N, (N_plexi + U)/N] in P2 §4.1 is valid but may not be
jointly sharp. Under the inclusion-minimal semantics there is **no** two-point bound. "Truth value
flips under unresolved literals" (Astra's rule) would then need all 2^k completions.

### 4.3 Direction of error for each data defect
| Defect | Relation of observed to truth | Error in the observed call |
|---|---|---|
| Unbinned or undetected donor (below breadth cutoff) | I_obs ⊑ I_true | **false plexitrophy** (and false no_donor). A not_plexi call is safe. |
| False gap in a donor (incomplete MAG) | I_obs ⊑ I_true | **false plexitrophy** (two partial producers look like fan-in) |
| False gap in recipient X | I_obs ⊑ I_true (X's true rules ⊇ observed, true needs ⊆) | **false plexitrophy**; at the extreme X is autonomous (d=0) and should not be a candidate auxotroph |
| Ignoring donors' own auxotrophies | I_model ⊒ I_true | conservative for plexi calls; anti-conservative for not_plexi (a "helper" that cannot itself live) |
| Contamination, KO false positive, "present" but nonfunctional pathway, organism detected but not co-active | I_obs ⊒ I_true | **false non-plexitrophy** |
| Chimeric bin or over-dereplication (merge) | not comparable | d deflated (false not_plexi); c or κ can be inflated (false cyclic) |
| One genome split across bins (Trap 9) | not comparable | d inflated (**false fan-in**) |

The dominant MAG defects (incompleteness and unbinned donors) both push toward **false-positive
plexitrophy**. So a raw observed plexitrophy rate is an upper-biased estimate unless it is
computed on I⁺. This is the formal reason P2's Null 1 (fragmentation) is necessary, and it says
which bound to report as primary: N_plexi computed with I⁺.

### 4.4 The unbinned pool
**Proposition 9.** Replace the unbinned pool by one relaxed pseudo-organism P*. P* has no needs,
makes every product any unbinned KO-route can make **without inputs**, and exports all of them.
Then NonPlexi(I with P*) ≥ NonPlexi(I with any partition of the pool into organisms, with any
needs).

*Proof:* map every pool organism to P* in a counterexample wiring.
- X's in-degree does not increase, because pool suppliers merge into one.
- P* has no in-arcs, so every cycle in the new graph avoids P* and already existed in the old one.
  SCC(X) shrinks and no new simple cycle appears.
- Non-pool derivations stay valid, because P* supplies everything without prerequisites. ∎

So P* belongs in I⁺ only. This formally justifies P2's rule that "unbinned can never be a single
d = 1 donor": a recipient that is not_plexi only through P* is `unresolved_ambig_donor`. It is
deliberately permissive. Any refinement of P* (for example contig-linked KO co-occurrence) is
valid only if it keeps P* ⊒ every plausible pool organism.

### 4.5 Remark on P2's one-level model
P2 scores donors as complete producers and does not model the donors' own needs. In that model no
organism depends on X, so **κ = 0 and c = 1 always, and the cyclic criterion can never fire**.
Remodelling chains (3a) are present but acyclic. P2 as written therefore tests **fan-in
plexitrophy only**. That should be said, or donors' genome-defined auxotrophies should be added as
literals (needs in I⁻, dropped in I⁺).

---

## 5. Recommendations for P2
1. **Change the quantifier** from "every feasible *minimal* support" to "every feasible wiring on
   any support" (P-all, equivalently irredundant supports).
   - It implies the team rule and matches "no single partner suffices".
   - It is monotone (Theorem 8) and polynomial for d_total (Theorem 3).
   - The team rule as written is coNP-complete and non-monotone (Theorem 2, A2).
2. **Replace |SCC| ≥ 3 with "X lies on a simple cycle of length ≥ 3" (κ).** Report c only as a
   descriptor. Add fixture A3 (X ↔ D ↔ Y) with expected `not plexitrophic`.
3. **Decide by witness search, never by enumerating M(X, E).** Use `nonplexi_poly` for d_total
   and `nonplexi_milp` for d_anab (coNP-complete in general, seconds at P2 scale). Drop the
   "incompletely enumerated" category: unresolved comes only from uncertain literals or solver
   time-out.
4. **Freeze the wiring semantics in the protocol**:
   - one source per (organism, metabolite);
   - no relaying (only producers export);
   - self- and E-supply count 0;
   - d = in-degree of X in the (typed) supply graph, with dual-role partners counted once;
   - catabolic services modelled as pseudo-metabolites from the provider.
5. **Use two-point identification** (I⁻, I⁺) per (X, E) to produce P2's five states (§4.2):
   - unbinned pool as relaxed P* in I⁺ only;
   - donors' genome-defined auxotrophies in I⁻ only where called with confidence;
   - contamination-flagged or "present-but-unresolved" capabilities removed in I⁻.

   This makes the [a, b] identification interval in §4.1 a valid outer bound, which may be
   loose (§6.2).
6. **State the direction of bias explicitly.** Incompleteness and unbinned donors give false
   *positive* plexitrophy. Contamination, chimeras and over-dereplication give false *negatives*.
   Bin splits (Trap 9) inflate d. Report the I⁺-based count as the conservative plexitrophy count.
7. **Say that P2's one-level model tests fan-in only** (κ ≡ 0, §4.5), or add donor-need literals
   if the cyclic arm is meant to be tested.
8. **Add A1–A5 to the mandatory fixture suite in §3.2**, and the "3a and 3b both feasible ⇒ not
   plexitrophic" case, which already passes. Keep the brute-force oracle as the reference
   implementation for fixtures.
9. **Report m as the minimum witness size** (MILP objective). Do not use any minimum of d or m in
   the decision; those are NP-hard optimisations (Theorem 5) and not needed.
10. **Review.** This analysis and code were written by one family (Claude). Before P2 is frozen,
    Codex or Gemini should review Theorems 2–4, Theorem 8, and `plexisolve.py`.

---

## 6. Revisions after cross-family review (Astra 248969335645061120)

Astra (Codex) reviewed plexisolve.py `c29070242a59` and test_plexisolve.py `c526692b8188` and
found one blocking bug and two overclaims. Astra confirmed two parts as sound: the monotonicity
of §4.1 under capability extension, and Theorem 4's gadget under its typing. Revised hashes:
plexisolve.py `533a8bcf6245`, test_plexisolve.py `3cc157681850`, and this file as reported to the
coordinator.

### 6.1 (Blocking) Relay through a futile cycle: the solvers disagreed
**Instance (now fixture F5).**
- A makes and exports p.
- D makes s, has rules p → q and q → p, and exports p and s.
- X needs p and s.

The old `nonplexi_poly` and the old brute-force oracle both used the two-node `closure`. They
found the witness "X imports p and s from D" (d = 1). The MILP found no counterexample. The
closure let D import A's p, run p → q → p, and re-export the result as "produced". Under §1.2 the
node (D, p) has one justification: if it is imported it cannot be exported, and if it is
rule-justified its proof needs (D, q), which needs (D, p), which is a cycle. The old brute force
shared the closure, so it was not an independent check on this point. Astra was right.

**Decision on semantics: strict single-node well-founded derivations (S1), as §1.2 already
said.** Re-exporting an imported metabolite is not allowed, including via a futile cycle. An
export of m by Z requires the single node (Z, m) to be rule-justified inside an acyclic proof
DAG. The MILP already implemented S1. If several justifications are active for one node, they
are harmless: an active rule justification has premises with strictly smaller proof depth, so it
never depends on (Z, m), and the import justification can simply be dropped.

I considered and did not adopt a more permissive alternative, the *lineage* semantics: Z may
export m if the proof of prod(Z, m) does not pass through Z's own pool of m. Under it, a
bidirectional interconverter that imports p and q could export both, by converting each into the
other. Under S1 it can export only one of them. Reasons for not adopting it:
- S1 is what §1.2 states and what the MILP implements.
- A both-ways interconverter that exports both is a relay of the pool in all but name.
- I have not determined the complexity of the lineage semantics (open).

The team should confirm the choice of S1.

**Lemma L5 (when the two-node closure is safe).** Suppose the global metabolite dependency graph
is acyclic. Then every two-node wiring has a single-node wiring on the same members whose supply
graph is a subgraph. This holds under any import restriction and any per-organism extra external
pool. *Proof:*
- Justify each available (Y, m) by a rule if prod(Y, m) holds in the two-node derivation.
  Otherwise import it from the same Z the two-node wiring used; that Z has prod(Z, m), so its
  node is rule-justified.
- Rank nodes by (topological position of m, 0 if rule-justified or external, 1 if imported).
  Rule premises have a smaller metabolite. An imported node depends on a rule-justified node of
  the same metabolite. So ranks strictly decrease and the DAG is acyclic.
- Only imports are dropped, so the supply graph shrinks. ∎

Theorem 3's proof used only closure monotonicity, S_max and import restrictions. **Theorem 3
therefore holds under an acyclic dependency graph, and the same goes for Proposition 4′.** The
code (`dep_acyclic`) checks the precondition and refuses otherwise. Pathway-module encodings such
as those in P2 (complete pathways, plus remodelling rules like cbi + DMB → B12) are typically
acyclic. Interconversion pairs (p ↔ q) are what break it.

**Theorem 3′ (Theorem 3 fails in general).** With cyclic dependencies allowed, deciding *total*
plexitrophy under P-all semantics and S1 is coNP-complete, for either cycle criterion.

*Membership:* a single-node wiring is a polynomial certificate.

*Hardness (fixture F6, `F6_interconverter_gadget`):* take a CNF φ over x₁..x_n.
- **A_i⁺** exports v_i⁺ and **A_i⁻** exports v_i⁻. Both are autonomous.
- **D** makes w and has rules v_i⁻ → v_i⁺ and v_i⁺ → v_i⁻. It exports w and all the v's.
- **X** needs w and c_j for every clause j, and has rules v_ℓ → c_j for each literal ℓ ∈ C_j.

Only D exports w, so a d ≤ 1 wiring must take everything from D. X exports nothing, so c = 1.
Under S1, D can export v_i⁺ only by rule-justifying (D, v_i⁺). That needs (D, v_i⁻) to be
imported, so v_i⁻ is not exportable at the same time. The set of v's D can export together is
therefore exactly the consistent literal sets. So a counterexample exists ⇔ φ is satisfiable.
X is always viable, by importing from A's and D, with d ≥ 2. ∎

Checked:
- Two hand gadgets (SAT and UNSAT) agree under the brute force and the MILP.
- 40 random interconverter gadgets (3 variables; 34 SAT, 6 UNSAT) agree under the MILP and
  a brute-force SAT check.

The earlier gadgets (Theorems 2, 4, 5) have acyclic dependency graphs, so their statements are
unaffected by the change of semantics.

**Code changes.**
- `nonplexi_brute_all` is now independent of `closure`. It enumerates allowed typed arc sets
  *and*, for every node that could be either produced or imported, whether it is rule-justified,
  then evaluates the single-node fixpoint `dag_smax`.
- `feasible_set`, `minimal_supports` and `plexi_brute_min` use the same single-node oracle.
- `nonplexi_poly` raises an error on cyclic dependencies.

Test results for the revised suite: 53/53 checks pass.
- **Fixtures:** all old ones, plus F5 (Astra's p ↔ q relay), F6 (the interconverter gadget,
  SAT and UNSAT), F7 (Carbon's {p1} → p2, {p2} → p1) and F8 (Astra's cross-organism cycle). F8
  without the rule r → p is acyclic, and all three deciders agree on it (`plexi`). In F5, F7 and
  F8 the poly decider correctly refuses; the brute force and the MILP give `plexi`.
- **Randomized:** 750 random instances × 4 variants, all decider pairs agree (3000/3000).
  - 77 of the instances have cyclic dependencies; there only the brute force and the MILP run.
  - On the 673 acyclic ones, the poly decider also agrees in the total variants.
- **SAT gadgets:** 40/40 anabolic gadgets (31 SAT) and 40/40 interconverter gadgets (34 SAT)
  agree.

The previous "3000/3000 agreement" was weaker evidence than I claimed. Two of the three deciders
shared `closure`, and the random generators rarely produced the relay pattern.

### 6.2 §4.2: the sharp-bound claim is withdrawn
Astra's counterexample: X needs p, and its only possible donor either exists and supplies p
(d = 1) or is absent (no support). No completion is plexitrophic, yet the old table put the
recipient in U and claimed "both answers attainable". That claim was false and is removed (§4.2
edited). What remains proved:
- Each robust state (plexi, not_plexi, no_donor) holds **iff all completions agree on it**,
  provided I⁻ and I⁺ are admissible completions (Theorem 8).
- U means only that completions disagree on the three-way state.

The interval [N_plexi/N, (N_plexi + U)/N] stays valid as an outer bound. Its upper end can be
loose, because U includes recipients for which plexitrophy is unattainable. Tightening it needs
the recipient's intermediate completions: 2^k for k relevant uncertain literals, which is cheap
for k ≲ 10. That is optional and not needed for validity.

### 6.3 Theorem 2 and A2: scope
Both concern **member-only** inclusion-minimality, which is the reading of the team text at
review l.366. Astra reports that P2's current scorer minimises (members, edges) jointly. I have
not read that implementation and attribute neither result to it. Whether joint minimisation has
the same non-monotonicity is an open question for its author.

### 6.5 Carbon's export-bundle condition (relay 248973006869020672) and Astra's counterexample (248973173642936320)
**Carbon's proposal.** For the bundle K that D exports to X, require
K ⊆ rule-products and K ∪ Need_D ⊆ HornClosure(R_D, (E ∪ Pool_to_D) ∖ K). It is implemented as
`can_jointly_produce_and_live` in `scratch/p2_scorer.py` (`6cfc29c2c5da`). I read that function
and `decide_nonplexi_p_all`, and ran the latter through an adapter on F5, F7 and F8. I did not
edit the file.

**Local claim: correct.** *Proposition 6.5a.* Fix D and fix a pool P of metabolites that D can
import, where P's availability does not depend on any node of D. Then D has a single-node
derivation in which every k ∈ K is rule-justified and Need_D is met iff Carbon's condition holds.
*Proof:*
- (⇐) Import every m ∈ P∖K. Closure stages order the rule firings, so the derivation is acyclic,
  and the K nodes are rule-justified.
- (⇒) In any such derivation D's imports lie in P∖K, because K nodes are not import-justified.
  Every rule-justified node is in the closure from those seeds. ∎

This correctly blocks the intra-organism cases: Astra's p → q → p (F5) and Carbon's
{p1} → p2, {p2} → p1 (F7). Carbon's scorer returns `plexi` on both, which is correct.

**Global use: incorrect, because the premise "P does not depend on D" fails.** Astra's
counterexample is now fixture **F8**:
- A makes p.
- D makes s and has rules p → q and r → p.
- B has q → r.
- X needs p and s.

D's pool contains r, but r exists only because B converts D's q, and D made q from the p it
imported. Removing K = {p} from D's final pool does not remove the products upstream of it that
were derived from D's own imported p. The only "own" p D could export would come from
D:p → D:q → B:r → D:p, which is not well-founded under S1.

I confirmed with the adapter that `decide_nonplexi_p_all` returns `not_plexi` with witness
(sink, D) under both c and κ. My MILP and the independent single-node brute force both return
`plexi` (d = 2, from A and D). The global dependency graph p → q → r → p is cyclic, so
`nonplexi_poly` refuses.

**Two further limits of the bundle approach.**
1. *Enumerating bundles is exponential.* `_minimal_import_bundles_for_org` lists the minimal
   bundles, which are minimal transversals in general (Theorem 6(i)). It is polynomial only under
   a stated bound b on the number of minimal bundles, at O(b · n) checks. Even with K fixed, the
   rule-or-import choices of *other* organisms (interconverters elsewhere, as in F8) remain a
   global choice. Theorem 3′ shows that this residual problem is coNP-hard. So no local
   per-donor condition can be exact in polynomial time unless P = NP.
2. *Under acyclic dependencies the condition is redundant.* By L5 every organism can export its
   whole two-node production simultaneously.

**Recommendation.** Do not adopt the local condition as an exact test. For the P2 scorer:
- run `dep_acyclic` first;
- if the dependency graph is acyclic, use the Theorem 3 decider: two-node closures, no bundles;
- otherwise use the MILP. It encodes global well-foundedness through the proof-depth variables
  t, so it handles F5, F7 and F8.

Carbon's function can be kept as a fast *necessary* pre-filter: if it fails for every bundle, D
cannot be X's sole supplier. This holds provided the pool passed in is a superset of D's possible
imports, because the closure is monotone in the seed set and in K ⊆ K′. It must not be used as a sufficient one.

### 6.6 Documentation fixes (Astra 248974912022560768)
Astra closed the review of the relay counterexample and the complexity results on plexisolve.py
`533a8bcf6245`, and endorsed the convention of one justification per node with producer-only
exports for P2. Two stale phrases in this file (`299f7978bb5e`) were corrected:
1. The U row of the §4.2 table now reads "robust state undetermined". It no longer says "both
   answers are attainable".
2. §5 recommendation 5 now calls the [a, b] interval "a valid outer bound, which may be loose".
   It no longer calls it "exact per recipient".

The code is unchanged (plexisolve.py `533a8bcf6245`, test_plexisolve.py `3cc157681850`).

### 6.4 Effect on the recommendations (§5)
- Recommendation 3: use `nonplexi_poly` only when `dep_acyclic` holds, and otherwise the MILP
  (which is exact under S1 in all cases).
- Recommendation 4: add "single justification per node; no relaying, including via futile
  cycles".
- Recommendation 5: two-point identification gives valid outer bounds and exact robust states,
  not a sharp U.
- Add F5–F8 to the fixture suite. F5 is Astra's p↔q relay, F6 the interconverter SAT gadget, F7
  Carbon's p1↔p2 pair, and F8 Astra's cross-organism cycle. Also add F8 without r → p, which is
  acyclic and must still give `plexi` from all three deciders.
