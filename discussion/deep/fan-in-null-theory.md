# Fan-in (d >= 2) null theory for P2

Author: claude-plexi (Opus 5.5), written as a subagent deep dive on 2026-11-18. This is a discussion
note, not a protocol, and has not been posted to any chat. **Revised the same day after
cross-family reviews by Astra (Codex, relay 248967796767506432) and Carbon (Gemini, relay
248968598093479936)**; see §7. All numbers below are from the revised scripts.

Inputs read (full text): `protocols/P2-plexitrophy-recipient.md` (draft, not frozen; §1, §3.3,
§3.4, §4), `notes/p2-incompatible-pairs.md` §§1–4, `discussion/plexi-concept-review.md`
("Mechanism layer", Known traps, P2, P5). No file under `data/` was opened. Everything numeric
below is either a proof or comes from synthetic simulations in `scratch/deep/`.

| Script | Purpose | Seed |
|---|---|---|
| `scratch/deep/fanin_closed_form.py` | exact formula, MC check, Tables A–E | `default_rng(20261118)` |
| `scratch/deep/fanin_curveball.py` | Null 2 (curveball) vs other nulls on a synthetic world with no recipient mechanism; nulls run on donor-only matrices | metacommunity `default_rng(20261119)`; workers `SeedSequence(20261121).spawn(8)` |
| `scratch/deep/bq_stacking_sim.py` (output `bq_stacking_sim.out`) | Black Queen loss model with stacking, incompatibility and a private-function control | `SeedSequence(20261123).spawn(40)` |
| `scratch/deep/check_carbon_inc.py` | exact enumeration of P(d = 1 \| inc, eligible X) and the endpoint co-supply example | none (exact) |

All run with `python3 -I`, numpy only. sha256 values are in §6.

**Scope.** P2 defines plexitrophy as d >= 2 OR c >= 3 over all minimal supports, under a
supply scenario E. This note covers only the fan-in part: d is the minimum number of co-occurring
donors whose union covers X's missing set J, under E_minimal, with no cycles, no unbinned donors
and no ambient supply. So d >= 2 means J is covered but no single donor carries all of J.
Unresolved calls are ignored here. Their effect is additive and is handled by P2 §4.1.

---

## 1. Closed-form null (Task 1). [PROOFS]

### 1.1 Exact expression

Setup:
- X misses the set J, with |J| = k.
- There are N donors, independent of each other but not necessarily identical.
- Donor i carries the whole of J with probability a_i.
- Donor i carries no function in S ⊆ J with probability z_i(S).
- Functions within one genome may be dependent in any way. Stacking is a within-genome dependence.

Define:
- A = {every j in J is carried by some donor} (J is covered);
- B = {some donor carries all of J}.

Then d = 1 ⇔ B and d >= 2 ⇔ A \ B.

**Proposition 1.**

    P(d >= 2) = ∏_i (1 − a_i)  −  Σ_{∅≠S⊆J} (−1)^{|S|+1} ∏_i z_i(S)
    P(d = 1)  = 1 − ∏_i (1 − a_i)
    P(no donor) = Σ_{∅≠S⊆J} (−1)^{|S|+1} ∏_i z_i(S)

*Proof.*
1. B ⊆ A, so A^c ⊆ B^c. Hence P(A \ B) = P(B^c) − P(A^c).
2. Donors are independent, so P(B^c) = ∏_i (1 − a_i).
3. A^c = ∪_j U_j, where U_j = {no donor carries j}. Donors are independent, so
   P(∩_{j∈S} U_j) = ∏_i z_i(S).
4. Inclusion–exclusion over the U_j gives the result. ∎

No assumption about within-genome dependence is used. For identical (iid) donors with
independent functions of prevalence q_j:
- a = Q = ∏_j q_j;
- z(S) = ∏_{j∈S} (1 − q_j).

For the symmetric case q_j = q:

    P(d >= 2) = (1 − q^k)^N − Σ_{s=1}^{k} (−1)^{s+1} C(k,s) (1 − q)^{sN}.

For a pair (k = 2) with joint cell probabilities p11, p10, p01, p00:

    P(d = 2) = (1 − p11)^N − (1 − q1)^N − (1 − q2)^N + p00^N.                     (★)

A Monte Carlo check with 200,000 replicates (Table "MC verification" in the script output) agrees
with the formula to within 0.002–0.003. It was run for independent cases with k = 2, 3, 4 and
for a latent-helper (stacking) case.

### 1.2 Asymptotics in N and dependence on k

Behaviour at the ends of the N range:
- **Small N.** P(d >= 2) = 0 at N = 1, because one donor cannot be a fan-in. At N = 2 with
  k = 2 it is 2 q1 q2 (1 − q1)(1 − q2).
- **Large N.**
  P(d >= 2) = (1 − Q)^N − Σ_j (1 − q_j)^N + O(max_{j<l} z({j,l})^N) ~ e^{−NQ}, provided Q < min_j q_j.
  - The decay rate is −ln(1 − Q) ≈ ∏_j q_j, which shrinks geometrically in k.
  - The characteristic N at which fan-in disappears is about 1/∏ q_j. For q = 0.3 this is
    11 (k = 2), 37 (k = 3) and 123 (k = 4).

So **P(d >= 2) is unimodal in N**. It rises while coverage is being filled in and falls once a
single all-carrier becomes likely. The peak, from e^{−NQ} vs k e^{−rN} with r = −ln(1 − q), is at

    N* ≈ ln(k r / Q) / (r − Q),

which agrees with the exact argmax to within 1 in every row of Table A.

**Table A** (independent functions, symmetric q; exact values; from `fanin_closed_form.py`):

| q | k | N=5 | N=10 | N=20 | N=40 | N=80 | N* | max |
|---|---|---|---|---|---|---|---|---|
| 0.1 | 2 | .119 | .329 | .590 | .640 | .447 | 32 | .657 |
| 0.1 | 3 | .064 | .266 | .658 | .917 | .922 | 55 | .937 |
| 0.2 | 2 | .267 | .462 | .419 | .195 | .038 | 13 | .481 |
| 0.2 | 3 | .265 | .634 | .817 | .725 | .526 | 20 | .817 |
| 0.3 | 2 | .316 | .334 | .150 | .023 | .001 | 7 | .359 |
| 0.3 | 3 | .448 | .678 | .576 | .335 | .112 | 11 | .682 |
| 0.3 | 4 | .439 | .814 | .847 | .722 | .522 | 15 | .866 |
| 0.5 | 2 | .176 | .054 | .003 | .000 | .000 | 4 | .195 |
| 0.8 | 2 | .005 | .000 | .000 | .000 | .000 | 2 | .051 |

**How large the confound is.** The null alone produces differences far bigger than P2's
δ_OR = ln 2 = 0.69:
- **Richness (q = 0.3, k = 2).** P(d >= 2) = 0.334 at N = 10 but 0.023 at N = 40. That is a
  logit difference of 3.06, about 4.4 × δ_OR, from richness alone.
- **Missing-set size (q = 0.3, N = 20).** k = 2, 3, 4 give 0.15, 0.58, 0.85, a logit difference
  of 3.4 between k = 2 and k = 4.
- **Prevalence (Table E, phi = 0).** Take a pair with B12/dsr-like prevalences
  (q = 0.2, 0.1) and a C0-like amino-acid pair (q = 0.85, 0.85):
  - the rare pair gives P(d = 2) = 0.40–0.54 for N = 10–40;
  - the C0-like pair gives P(d = 2) < 3·10⁻⁶;
  - the log-OR is > 12 (the script clips probabilities at 1e-9).

  **Table E assumes independent functions in both pairs** (phi = 0, p11 = q1 q2, so the
  "incompatible-like" pair here *does* have co-carriers). Its log-OR > 12 therefore comes from
  **unequal marginal prevalence alone**. Structural exclusion (p11 = 0 for a truly incompatible
  pair) is a separate effect, treated in §1.3 (Prop. 2). It pushes the incompatible side further
  toward d = 2, and the two effects add on the logit scale.

  So **matching recipients on missing-module count, genome size and completeness (P2 §3.3) does
  not match the donor-side quantities that set P(d >= 2): event richness N, k, and the donor
  prevalences q_j of the specific missing functions.** C0 modules (his, trp, leu, ...) are
  near-universal among donors, so C0 × C0 recipients are d = 1 almost surely under any null.

The direction of the richness effect also changes. P(d >= 2) rises with N below N* and falls
above it. Rich events (deep Illumina co-assemblies) and poor events (Nanopore, few HQ MAGs) can
therefore give opposite-signed differences even with identical biology. The pooled log-OR in
P2 §4.3 will depend on the cohort's mix of N.

### 1.3 Effect of within-genome correlation (stacking)

**Proposition 2 (k = 2).**
- Fix the marginals q1 and q2, and let c = p11 − q1 q2 be the covariance (phi = c / sd).
- Then ∂P(d = 2)/∂c = N [p00^{N−1} − (1 − p11)^{N−1}] < 0 for N >= 2 whenever p10 + p01 > 0.
- So P(d = 2) is strictly decreasing in phi.
- At the mutual-exclusion bound (p11 = 0): P(d = 2) = P(A), so every covered recipient is a
  fan-in.

*Proof.*
1. Differentiate (★). Both p11 and p00 increase by c.
2. 1 − p11 = p00 + p10 + p01 > p00, so the bracket is negative. ∎

**Corollary (pathway-level exclusion only).** Suppose the two missing functions are *pathways*
that cannot share a cytosol, so p11 ≈ 0 apart from compartmentalisers and chimeric bins. Then
d = 2 holds for every covered eligible recipient, by construction. The only route to d = 1 is a
co-encoding genome, which P2 already tallies in its S_spatial / S_temp/resp / S_unexplained
buckets.

**Endpoint-level co-supply (Carbon, 248968598093479936).** P2's recipient strata are defined by
missing *endpoints*, and endpoints can have OR-routes. In R-I2 and R-I3, the cobalamin endpoint is
supplied by the aerobic route (cobG/cobNST + bluB) **or** the anaerobic route (cbi + bza).
An anaerobe encoding cbi + bza together with dsr (R-I3) or nif (R-I2) therefore supplies both
endpoints at d = 1 with no compartmentalisation. So, at the endpoint level:
- p11 = P(one genome supplies both endpoints) > 0;
- the O2 argument excludes only the aerobic-B12 × anaerobic-partner *route* combination.

Prop. 2 still applies with the endpoint-level p11. The fan-in rate is very sensitive to it at
moderate and high N (`check_carbon_inc.py`, (★), q1 = 0.3 for B12 any route, q2 = 0.1):

| p11 | N = 10 | N = 20 | N = 40 | N = 80 |
|---|---|---|---|---|
| 0 | 0.629 | 0.878 | 0.985 | 1.000 |
| 0.01 | 0.535 | 0.696 | 0.654 | 0.447 |
| 0.02 | 0.449 | 0.545 | 0.431 | 0.198 |
| 0.05 | 0.235 | 0.236 | 0.114 | 0.016 |

Consequences:
- For the real R-I2 and R-I3 pairs the inc-stratum d = 2 rate is **not** tautological. It
  measures how rare anaerobic-route all-in-one suppliers are.
- It is still a donor-composition quantity, (q1, q2, p11) at N, and is fully predicted by
  the donors without reference to the recipient. The main conclusion (Q2 contrast = donor
  composition) is unchanged.
- The inc-pair null must be built at the **endpoint** level, scoring each donor's endpoint supply
  through every route. The synthetic O± modules in §2.2 are pathway-level, so §2.2 is the
  p11 = 0 limiting case.

Table B gives magnitudes (q1 = q2 = 0.3):

| phi | N = 10 | N = 20 | N = 40 |
|---|---|---|---|
| −0.43 (exclusion) | 0.94 | 1.00 | 1.00 |
| 0 | 0.33 | 0.15 | 0.02 |
| +0.3 | 0.14 | 0.035 | 0.001 |

For k = 3 and 4 the stacking direction is numerical, not proved. Table C uses a latent-helper
model (a fraction h of donors are helpers carrying each function with probability 0.9, with the
marginal fixed at q = 0.3, N = 20). P(d >= 2) falls monotonically with h:

| | h = 0 | h = 0.05 | h = 0.1 | h = 0.2 |
|---|---|---|---|---|
| k = 3 | 0.58 | 0.32 | 0.17 | 0.04 |
| k = 4 | 0.85 | 0.46 | 0.24 | 0.06 |

Positive within-genome association raises coverage slightly (by FKG, since "some donor has j" is
an increasing event). It raises single-donor coverage much more, so the net effect is fewer
fan-ins.

### 1.4 Donor-side incompleteness (a gap in Null 1)

If each donor route is detected with sensitivity s, then q → sq, and P(d = 1) falls by roughly
a factor s^k while coverage falls only by s. Table D shows the inflation of P(d >= 2) with
q = 0.3 and k = 2:

| s | N = 20 | N = 40 | N = 80 |
|---|---|---|---|
| 0.9 | ×1.44 | ×2.1 | ×4.4 |
| 0.8 | ×2.0 | ×4.1 | ×16 |

P2's Null 1 fragments only the *recipient*, to catch false auxotrophy. Donor fragmentation and
the 80% step threshold for donors manufacture fan-ins, and the effect grows with N. **Null 1 must
also be applied to donors.** Draft-genome donors are mainly the non-HQ MAGs that P2 admits to
the donor pool.

---

## 2. What the curveball Null 2 holds fixed (Task 2)

### 2.1 Fixed and free

| Quantity | Held by curveball within an event |
|---|---|
| event richness N, the module set | yes |
| each module's prevalence in the event (column sums) | yes |
| each genome's module count (row sums) | yes |
| pairwise co-occurrence / phi between modules | **no**: driven to its degree-implied value |
| biochemical feasibility (no aerobe with nifHDK + bluB, ...) | **no**: impossible genomes are generated |
| lifestyle or phylogenetic block structure | **no** |
| the recipient's own row | eligibility is frozen, so it should be held fixed (P2 says "permute donor module presences") |

Three consequences follow.

**(i) For incompatible pairs, observed minus Null 2 is positive automatically** (for any
endpoint-level p11 below the curveball's margin-implied value; at pathway level p11 = 0).
- Observed pathway-level p11 = 0 by biochemistry.
- Curveball creates co-carriers at the margin-implied rate, so null d = 1 is common.
- By Prop. 2, the observed d = 2 rate exceeds the null for *any* recipient-level biology.
- P2-Q1's "excess over Null 2" in the R-I2 and R-I3 strata is therefore largely a restatement
  of the incompatibility, plus the rarity of anaerobic-route all-in-one suppliers (endpoint p11),
  not evidence of plexitrophy.

**(ii) For compatible pairs, the sign of the excess is minus the sign of phi_obs − phi_CB.**
Genuine Black Queen stacking (phi > 0) gives a *negative* excess over Null 2.

**(iii) Null 2 absorbs stacking.** Its row sums are counts *of the modules being permuted*. A
helper that "shoots the moon" is exactly a row with a large sum, and curveball preserves it. In
the limit, one helper carries all leaky functions and the rest carry none; the matrix is then
the unique one with those margins and the null equals the observation.
- §3 shows this numerically. Under strong stacking (S = 10, m = 1, β = 1), P(d = 1 | comp) is
  0.884 observed and 0.884 under curveball, but 0.111 under the exact-count independent null.
- Null 2 therefore cannot detect the mechanism it is named after.
- "Genome size" should be preserved through an *exogenous* covariate (assembly bp, or the count
  of modules outside the tested set), not through the row sum of the tested modules.

Curveball and the binomial independent model (IND, (★)) differ in two ways at the same
prevalences, with opposite effects:
- **Exact column counts.** Curveball fixes the count of each column, so a function present in
  the donors can never become uncovered. This raises d >= 2 relative to IND at small N, where
  IND puts mass on "no donor".
- **Row-sum heterogeneity.** Under a max-entropy approximation p_ij ∝ r_i c_j, the expected
  number of donors carrying both functions is Σ_i p_i1 p_i2. By Chebyshev's sum inequality this
  grows with Var(r_i), which lowers d >= 2.

In §2.2 the first effect dominates at N = 10, where R×R gives N2 0.370 vs IND 0.254. The second
dominates at N >= 20 (0.336 vs 0.412, 0.137 vs 0.259, 0.016 vs 0.065). This is an approximation
argument, not a proof.

### 2.2 Simulation: a world with no recipient mechanism

`fanin_curveball.py` builds the world and nulls as follows:
- **Metacommunity.** 3,000 genomes; 55% aerobes; log genome size g ~ N(0, 1).
- **Modules.**
  - 10 C0 modules, prevalence about 0.85;
  - 2 O+ modules (aerobes only) and 2 O− modules (anaerobes only);
  - 4 rare compatible modules R, prevalence about 0.18, structured only by g. These
    prevalence-match the O± modules.
- **Events and recipients.** Events are random N-subsets. Recipients' missing sets are
  independent of the donors.
- **Nulls.** 30 events per N. **Each recipient X is removed before any trading** (Astra fix 1),
  so each null keeps the exact row and column sums of X's donors. Each recipient gets one chain:
  5N trades of burn-in, then 10 draws thinned by N trades.
  - N2 = curveball;
  - N2c = curveball that forbids O+/O− co-carriage;
  - IND = **i.i.d. plug-in** null: closed form (★) with each function's donor prevalence plugged in
    as an independent Bernoulli probability (counts are binomial, so they are not fixed). This differs
    from §3's fixed-count IND.
  - GI = genome-intact null: 20 fresh random events per observed event, recipients selected
    as in the design.
- **Variant "filter".** Recipients are kept only where both of their missing functions are
  covered. This is a stand-in for supply-dependent persistence, the one recipient-level
  mechanism.

Results are in `scratch/deep/fanin_curveball.out` (revised run). Each cell is P(d = 2) among
eligible recipients, summed over 30 events. The world has no recipient mechanism. Metacommunity
mean phi is −0.21 for inc pairs, +0.10 for rare compatible pairs (R × R) and +0.08 for C0 × C0.

| N | stratum | observed | N2 curveball | N2c (exclusion kept) | IND | GI |
|---|---|---|---|---|---|---|
| 10 | inc | 0.774 | 0.462 | 0.774 | 0.282 | 0.642 |
| 10 | R×R | 0.364 | 0.370 | 0.334 | 0.254 | 0.351 |
| 20 | inc | 0.919 | 0.404 | 0.919 | 0.376 | 0.925 |
| 20 | R×R | 0.217 | 0.336 | 0.284 | 0.412 | 0.323 |
| 40 | inc | 1.000 | 0.209 | 1.000 | 0.342 | 0.997 |
| 40 | R×R | 0.103 | 0.137 | 0.098 | 0.259 | 0.104 |
| 80 | inc | 1.000 | 0.041 | 1.000 | 0.120 | 1.000 |
| 80 | R×R | 0.006 | 0.016 | 0.008 | 0.065 | 0.013 |
| all | C0×C0 | 0.000 | 0.000 | 0.000 | ≤0.002 | 0.000 |

Readings, in a world with **no** recipient mechanism and **no** Black Queen dynamics (and with
pathway-level exclusion, p11 = 0):
1. **P2-Q1 against Null 2.** The excess in the inc stratum is +0.31, +0.52, +0.79 and +0.96 at
   N = 10, 20, 40, 80. Q1 would be "Supported" by biochemistry alone.
2. **P2-Q2's observed inc vs R×R log-OR** is 1.79, 3.71, 11.4 and 14.3.
   - It is above ln 2 at every N, even against the prevalence-matched compatible control.
   - Against C0 × C0, the drafted control, it is unbounded: the C0 rate is 0 throughout.
   - The contrast grows steeply with N, so the pooled value depends on the cohort's richness mix.
3. **Null 2's own inc vs R×R log-OR** is 0.30 to 0.95. So "observed minus Null 2" is also large:
   Null 2 manufactures a mechanism-free effect.
4. **Keeping the exclusion (N2c) reproduces the observation.** The N2c log-OR is 1.92, 3.36,
   11.4 and 14.0; the inc rates are identical to the observed ones.
5. **GI also tracks the observation**, with log-OR 1.20, 3.26, 7.95 and 13.5. The residuals
   are sampling noise: there are only 30 observed events per N and they are event-clustered,
   and logits near 0 or 1 are clipped at 1e-4. The largest gap is inc at N = 10 (0.774 vs 0.642).
   This shows the GI residual needs an event-level bootstrap before it is read.

With supply-dependent recipient persistence ("filter" variant), the observed inc rate becomes
1.000 at every N:

| N | GI inc | GI R×R | observed R×R |
|---|---|---|---|
| 10 | 0.671 | 0.359 | 0.531 |
| 20 | 0.943 | 0.317 | 0.303 |
| 40 | 0.998 | 0.111 | 0.105 |

So the GI residual detects the recipient mechanism at N = 10, with an inc residual of +0.33 and
an R×R residual of +0.17. By N = 20 the R×R residual has vanished. At N >= 40 coverage is near
certain, and no null can detect the mechanism from d alone. There the signal is only in the d = 0
residual, and only where coverage is incomplete.

**Answer to "is the Q2 contrast already explained by marginal co-occurrence structure?":** yes.
Given the donor genomes, d is determined, and the inc vs comp contrast follows from:
- the prevalences of the missing functions;
- the donor-level p11 of the pair: ≈ 0 for pathway-level exclusion, small but > 0 at endpoint
  level through OR-routes (Prop. 2, §1.3);
- event richness N.

None of these involves the recipient. The only recipient-level information is where recipients
occur relative to supply. The GI null isolates that, and it is weak once N is large.

### 2.3 The right null for P2-Q2, and a statistic that isolates the mechanism

For a fixed event and a fixed recipient, d is a *deterministic* function of the donor rows.
Given the donor pool, the Q2 contrast between incompatible (inc) and compatible (comp) pairs
decomposes completely into three parts:
- (a) the prevalences of the missing functions (Table E);
- (b) between-function co-occurrence within donor genomes, phi (Prop. 2, including biochemical
  exclusion);
- (c) which events the recipients occur in, relative to their donors.

A Black Queen recipient mechanism can act only through (c), and helper stacking only through
(b) for *compatible leaky* functions. P2-Q2 as drafted (inc vs C0 × C0, δ_OR = ln 2) is
dominated by (a) and by the forced part of (b). It can be "Supported" in a world with no
mechanism, as §§2.2 and 3 show.

Recommended nulls and statistics:

1. **Genome-intact (GI) null for recipient-level claims.**
   - Keep every genome's full content, so all phi, exclusions and prevalences are exact.
   - Randomise the genome × event incidence matrix by curveball, preserving each event's
     richness N_e and each genome's occupancy. Do this within habitat/redox strata, so that
     anoxic events keep anaerobes.
   - Recompute eligibility and d.
   - The statistic is the difference-in-differences

         Δ_mech = [logit P_obs(d>=2 | inc) − logit P_GI(d>=2 | inc)] − [logit P_obs(d>=2 | comp) − logit P_GI(d>=2 | comp)],

     computed within (N-bin, k) strata.
   - This removes (a) and (b) exactly and leaves the recipient–donor co-assembly signal.
   - Report the GI residual of P(no donor) alongside it. Supply-dependent persistence shows up
     first as a *deficit of d = 0* relative to GI, in both strata.

2. **Stacking test for the "shooting the moon" claim (donor side, no recipients needed).**
   - For compatible *leaky* pairs, define

         Σ_stack = logit P_obs(d = 1 | comp) − logit P_IND,g(d = 1 | comp),

     where IND,g is an independent-carriage null conditional on exogenous genome size:
     logit p_ij = α_j + β log bp_i, fitted per cohort, with event-level prevalences.
   - The mechanism predicts Σ_stack > 0 for leaky public goods (B vitamins, cobamide precursors,
     siderophores). It predicts the interaction I = Σ_stack(leaky) − Σ_stack(private) > 0
     against non-leaky controls of matched prevalence.
   - The IND null redistributes each function's donor count exactly (rows free), on the
     donor-only matrix.
   - A co-occurrence structure that comes only from lifestyle or phylogeny should affect both
     classes alike and give I ≈ 0. That is an assumption for real data: it fails if private
     functions are more clade-concentrated than leaky ones. This is why conditioning on clade is
     needed (§4).
   - **The test requires leakiness labels that are independent of the co-occurrence data**, for
     example from chemistry (membrane permeability, secretion) and transport (presence of
     dedicated uptake systems in recipients). Labels inferred from the same co-occurrence
     patterns make I circular.
   - §3.3 calibrates I in the model: β = 0 gives I ≈ 0, and stacking gives I > 0.

3. **Blocking.**
   - At the pathway level the blocking is the definition, so it is not testable as an
     enrichment. What is testable is its *exceptions*: the per-event rate of co-encoders, which
     P2 §4 already tallies as descriptive buckets.
   - At the endpoint level, the informative quantity is the donor-level rate of all-in-one
     endpoint suppliers through alternative routes (anaerobic cbi + bza with dsr or nif).
   - Report P(d = 1 | inc) and that rate, not a log-OR.

---

## 3. A minimal Black Queen model with helper stacking (Task 3)

### 3.1 Model

The model has S species and Kc = 6 leaky compatible functions, Kp = 6 **private** compatible
functions (added in revision, Astra fix 3), one leaky O+ and one leaky O−. Every species needs
every leaky product.
- **Types.** A species is an aerobe (probability ρ = 0.5, allowed O+), an anaerobe (allowed O−),
  or, with probability κ, a compartmentaliser allowed both.
- **Ancestor.** Every species is a full prototroph within what its type allows.
- **Loss.** A carried (i, j) is lost at rate w_i = e^{−β n_i}, where n_i is the number of
  functions i still carries. The mutation supply scales with population size, and population
  size falls with biosynthetic burden. β >= 0 sets the strength of this "burdened helpers lose
  slowly" feedback, i.e. rich get richer.
- **Black Queen rule.** A loss fixes only if at least m other holders of j remain. Losing the last
  holder makes j community-limiting and is selected against. m > 1 stands for low leakiness or
  high demand.
- **Private functions.** Species i has a fixed niche need r_ij ~ Bernoulli(m/S), chosen so
  prevalence matches the leaky class. A private loss fixes if and only if i does not need j,
  whatever other species carry. Private functions count toward the burden n_i.
- **Absorption.** The run ends when no admissible loss remains. Each leaky function then has
  exactly m holders. Each private function is held by exactly the species that need it.
- The O± functions are **pathway-level**: there is no alternative route, so endpoint co-supply
  (§1.3) is absent except through compartmentalisers.

### 3.2 Analytic results (β = 0) [PROOFS]

With β = 0 the loss process of each function is a uniform random elimination among its allowed
carriers, and functions are independent. So the final holder set of each function is a uniform
m-subset of its allowed carriers, independently across functions.

**Compatible pair (m = 1).** The recipient X must not hold either function, so condition on
X ∉ {h1, h2}. Then h1 and h2 are independent and uniform on the other S − 1 species, and

    P(d = 2 | comp) = (S − 2)/(S − 1).

**Compatible pair (general m).** Each holder set is a uniform m-subset of the S − 1 non-X
species, so

    P(d = 2 | comp) = C(S − 1 − m, m) / C(S − 1, m).

**Incompatible pair, κ = 0.** h+ is an aerobe and h− is an anaerobe, so P(d = 2 | inc) = 1.

**Incompatible pair, n_a aerobes, n_b anaerobes, n_k compartmentalisers** (corrected per
Carbon, 248968598093479936). d = 1 needs h+ = h− = one compartmentaliser. Unconditionally this
has probability n_k / AB, where A = n_a + n_k and B = n_b + n_k. But eligibility conditions on X
holding neither function:
- if h+ = h−, there are S − 1 eligible recipients, all with d = 1;
- otherwise there are S − 2, all with d = 2.

So the pooled rate over eligible recipients, which is what `tally()` computes, is

    P(d = 1 | inc, X eligible) = n_k (S − 1) / [ n_k (S − 1) + (AB − n_k)(S − 2) ].

Exact enumeration in `check_carbon_inc.py` agrees to 5 decimals:

| (n_a, n_b, n_k) | conditional | unconditional |
|---|---|---|
| (4, 5, 1) | 0.03734 | 0.03333 |
| (2, 2, 2) | 0.15152 | 0.12500 |
| (9, 10, 1) | 0.00959 | 0.00909 |

The same pooled-ratio argument reproduces the compatible formula, P(d = 1 | comp) = 1/(S − 1).

**Consequence.**
- Take S = 20 with (n_a, n_b, n_k) = (9, 10, 1). Then P(d = 2 | comp) = 0.947 and
  P(d = 2 | inc) = 0.990, giving log-OR ≈ 4.64 − 2.89 = 1.75 > ln 2.
- That is "Supported" with **no stacking at all**. The log-odds scale amplifies differences
  between probabilities that are both near 1.
- In absolute terms the contrast is 4 percentage points.

### 3.3 Simulation (β > 0)

Revised run: 200 replicates per cell, all 40 cells. The nulls are computed per recipient on the
**donor-only** matrix, for every replicate:
- IND here is a **fixed-count permutation** null, unlike §2's i.i.d. plug-in IND. It
  redistributes each function's **exact** donor count over uniformly random donor rows
  (Astra fix 2);
- CB is curveball on the leaky-compatible block (5S burn-in, 5 draws thinned by S).

Readouts:
- **Full community:** P(d >= 2) for compatible leaky pairs (15 pairs) and for the inc pair.
- **Detected subset:** the event observes N_det = S/3 species, sampled in proportion to
  abundance w_i. Helpers are rare when β > 0.
- **Stacking statistics:**
  - Σ_C = logit P_obs(d = 1 | class C) − logit P_IND(d = 1 | class C), for C = leaky or private;
  - I = Σ_leaky − Σ_private, with a 95% percentile bootstrap over replicates (300 resamples).

Results are in `scratch/deep/bq_stacking_sim.out` (40 cells). In the table, d2c and d2i are
P(d = 2 | comp) and P(d = 2 | inc), and d1c is P(d = 1 | comp).

| S | m | κ | β | d2c | d2i | log-OR | d1c obs | IND | CB | Σ_leaky | Σ_priv | I [95%] | analytic d2c (β = 0) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 10 | 1 | 0 | 0 | 0.895 | 1.000 | 7.1 | 0.105 | 0.110 | 0.105 | −0.06 | 0.01 | −0.07 [−0.27, 0.07] | 0.889 |
| 10 | 1 | 0 | 1 | 0.116 | 1.000 | 11.2 | 0.884 | 0.111 | 0.884 | 4.11 | −0.05 | 4.16 [3.94, 4.45] | |
| 10 | 2 | 0.05 | 0 | 0.592 | 0.954 | 2.7 | 0.408 | 0.415 | 0.407 | −0.03 | −0.06 | 0.03 [−0.06, 0.13] | 0.583 |
| 30 | 1 | 0 | 0 | 0.968 | 1.000 | 5.8 | 0.032 | 0.035 | 0.032 | −0.09 | −0.04 | −0.05 [−0.39, 0.24] | 0.966 |
| 30 | 1 | 0 | 0.25 | 0.728 | 1.000 | 8.2 | 0.272 | 0.035 | 0.272 | 2.33 | 0.06 | 2.27 [2.06, 2.50] | |
| 30 | 1 | 0 | 0.5 | 0.367 | 1.000 | 9.8 | 0.633 | 0.034 | 0.633 | 3.88 | −0.05 | 3.92 [3.68, 4.18] | |
| 30 | 1 | 0 | 1 | 0.139 | 1.000 | 11.0 | 0.861 | 0.034 | 0.861 | 5.16 | −0.07 | 5.23 [4.95, 5.55] | |
| 30 | 1 | 0.05 | 0 | 0.966 | 0.995 | 1.9 | 0.034 | 0.034 | 0.034 | −0.02 | 0.13 | −0.15 [−0.41, 0.10] | |
| 30 | 1 | 0.05 | 1 | 0.109 | 0.652 | 2.7 | 0.891 | 0.035 | 0.891 | 5.42 | 0.04 | 5.38 [5.09, 5.66] | |
| 30 | 1 | 0.05 | 2 | 0.044 | 0.392 | 2.6 | 0.956 | 0.035 | 0.956 | 6.39 | −0.06 | 6.45 [6.07, 6.94] | |
| 30 | 2 | 0 | 0 | 0.858 | 1.000 | 7.4 | 0.142 | 0.136 | 0.141 | 0.05 | −0.00 | 0.05 [−0.09, 0.21] | 0.865 |
| 30 | 2 | 0.05 | 0.5 | 0.059 | 0.707 | 3.7 | 0.941 | 0.135 | 0.945 | 4.64 | −0.09 | 4.72 [4.46, 5.05] | |

The β = 0 rows match the §3.2 analytic values to within Monte Carlo error (0.895 vs 0.889;
0.968 vs 0.966; 0.858 vs 0.865; 0.592 vs 0.583).

What the model says:
1. **The model predicts the P2-Q2 enrichment, but so does every variant of it, including
   no stacking (β = 0).**
   - Over the 8 β = 0 cells, the inc vs comp log-OR is 5.8–8.9 without compartmentalisers and
     1.6–2.7 with κ = 0.05.
   - Stacking (β > 0) adds to the contrast by lowering d2c, but the contrast is above ln 2
     without it. **P2-Q2 cannot distinguish blocked stacking from no stacking.**
2. **Under stacking, compartmentalisers become the moon-shooters.** With κ = 0.05 and β >= 1,
   the compartmentaliser keeps both O functions *and* the compatible ones. P(d = 2 | inc) falls to
   0.30–0.65 at S = 30, so stacking *reduces* fan-in for incompatible pairs when any bypass exists.
   - Prediction: co-encoders of an incompatible pair carry more leaky compatible functions than
     other donors.
3. **Stacking statistic, leaky vs private (revised; replaces the earlier "1.4–4.8" claim).**
   - In the earlier version, "1.4–4.8" was the observed compatible-pair d = 1 rate against a
     Bernoulli-count independent null. It was **not** a leaky-minus-private interaction (Astra
     fix 3).
   - With the private class and exact-count nulls:
     - **β = 0 (null calibration):** all 8 cells have I in −0.15 to 0.12, and every 95% CI
       contains 0.
     - **β = 0.25:** I = 1.36–2.80.
     - **β = 0.5:** I = 2.71–4.72.
     - **β = 1:** I = 4.14–7.56.
     - **β = 2:** I = 5.38–11.1.
     - The lower bounds exceed ln 2 in every β > 0 cell.
     - These are bootstrap intervals for the **simulated mean contrast** across 200 replicate
       communities per cell. They do not establish cohort-level power.
   - Σ_private stays within ±0.25 throughout. This is **by construction** in the model, because
     private retention ignores other species. The model shows that the statistic responds to
     stacking and not to burden heterogeneity alone. It does not show that real private
     functions lack lineage structure.
4. **Curveball (Null 2) tracks the observed P(d = 1 | comp) at every β** (0.884 vs 0.884;
   0.941 vs 0.945). This confirms §2.1(iii): Null 2 conditions on the stacking.
5. **Detection.**
   - Helpers are rare when β > 0, because w = e^{−βn}.
   - With N_det = S/3 detected by abundance, P(d = 0 | inc) is 0.72–1.0. Most eligible
     recipients become `unresolved_no_donor`, and P2 gate 3 (U/N > 0.5) trips.
   - Prediction: the d = 0 fraction falls with sequencing depth, and the d = 1 donors that are
     found have lower abundance than recipients.
6. **Realistic magnitudes.** β in about [0.1, 1] is an assumption, not a measurement. The
   biosynthetic cost per pathway is O(10⁻³–10⁻²) of growth, and the abundance response is
   unknown. Over β = 0.25–1, I is about 1.4–7.6.

**An observable that separates the mechanism from the §2 co-occurrence null.** The primary
observable is I > 0: leaky compatible functions co-occur within single donors more than private
functions of matched prevalence. The comparison is against exact-count, donor-only, genome-size-
and clade-conditioned nulls, with leakiness labelled independently of co-occurrence. Supporting
evidence:
- d = 1 donors are low-abundance;
- incompatible-pair co-encoders, where present, are enriched for leaky functions.

A pure lifestyle or phylogeny co-occurrence null predicts I = 0, *after clade conditioning*, and
no abundance gradient.

---

## 4. Recommendations for P2 (Task 4)

1. **Separate the claims.** P2-Q1/Q2 currently fold three claims into one log-OR:
   - fan-in exists beyond artefacts;
   - stacking occurs where it can;
   - stacking is blocked by O2.

   Under the theory above, the inc vs comp contrast is determined by donor composition
   (Prop. 2, §2.2) and is above ln 2 even with no stacking (§3). Recommendations:
   - **Demote the Q2 contrast as drafted to descriptive.**
   - **Make the primary mechanism test the interaction I = Σ_leaky − Σ_private** (§2.3 item 2,
     §3.3 item 3). Build the nulls as exact-count, donor-only nulls conditional on exogenous
     genome size and clade, and pre-register the leakiness labels.
     - Proposed margin: ln 2. Supported if the lower bound exceeds ln 2; falsified if the upper
       bound is below ln 2.
     - In the model, every β >= 0.25 cell clears this, and every β = 0 cell has a 95% CI
       containing 0. These intervals describe the simulated mean contrast across 200 replicate
       communities, not the power of a real cohort, which must be assessed separately at the
       cohort's event and recipient counts.
   - **Report the incompatible-pair blocking** as the per-event co-encoder rate
     P(d = 1 | inc), with the bypass buckets.
2. **Replace or demote Null 2.**
   - Plain curveball on the tested modules (a) generates biochemically impossible genomes,
     inflating every inc excess, and (b) preserves stacking through row sums.
   - Use **N2c**, a curveball that forbids co-carriage of pre-registered incompatible pairs,
     with row sums taken over modules *outside* the tested set or replaced by assembly length.
   - Use the **GI null** (genome-intact, event incidence curveball within habitat/redox strata)
     for any recipient-level claim. The P2-Q1 excess should be stated against N2c and GI, not
     against N2.
3. **Stratify by N and k, and never pool across them.**
   - Bin event richness N (detected dereplicated MAGs per event), for example 5–9, 10–19,
     20–39, 40–79 and 80+. Report each bin with its denominator.
   - Analyse k = 2 pairs separately from k >= 3.
   - Within a bin, also match comp controls on the donor prevalences q_j of the missing
     functions, not only on recipient genome size. A q-matched rare compatible set (like R
     here) is needed, because C0 × C0 has P(d >= 2) ≈ 0 under every null.
   - Combine bins by Mantel–Haenszel or a stratified model. Show the per-bin estimates, since the
     sign can differ across N* (§1.2).
4. **Report observed-over-expected on the probability scale.**
   - Report O − E with E from the closed form (Prop. 1) at each event's N and q, and from N2c and
     GI.
   - Logits near 0 or 1 make log-OR margins unstable: §3.2 gives log-OR 1.75 for a
     4-percentage-point difference.
   - If a log-OR margin is kept, add an absolute-risk margin, for example ≥ 10 pp, that must
     also be cleared.
5. **Apply Null 1 to donors as well** (§1.4), i.e. joint recipient + donor downsampling. Donor sensitivity s = 0.9 already doubles
   P(d >= 2) at N = 40. Either restrict the donor pool to calibrated-complete routes, or propagate
   donor-side fragmentation into E.
6. **Expect gate 3 to trip** where helpers are rare (§3.3 item 4). Pre-register a
   depth-sensitivity readout: d = 0 and unresolved fractions vs sequencing depth or detection
   breadth threshold. A falling d = 0 fraction with depth is a mechanism-consistent pattern, not
   only a quality problem.
7. **Score incompatible strata at the endpoint level** (§1.3). Count all-in-one endpoint
   suppliers through alternative routes, for example anaerobic cbi + bza + dsr for R-I3, as
   d = 1. Report their donor-level rate per event. Build N2c's forbidden set from route
   combinations, not endpoints.
8. **The pre-registration should name the expected null values, not just the nulls.** For each
   stratum, before any data are touched, freeze the closed-form (independent) and N2c predictions
   computed from donor prevalences. The test then reads as a deviation from a number written down
   in advance.

---

### 4.1 Does the revised analysis support Carbon's proposed design (248968598093479936)?

**Q1 nulls: joint recipient + donor contig downsampling, plus a habitat/redox-stratified
genome-intact incidence null, stratified by N-bin and k.** Yes, supported.
- Donor-side downsampling is required by §1.4: s = 0.9 doubles P(d >= 2) at N = 40.
- The GI null reproduces the no-mechanism world in §2.2. N2 does not; it gives +0.31 to +0.96
  spurious inc excess.
- N and k stratification follows from Table A.

Three cautions:
- The GI residual is noisy with about 30 events per stratum: up to 0.13 at N = 10 with no
  mechanism. It needs an event-cluster bootstrap.
- GI has almost no power at N >= 40, where coverage saturates. Pre-register the d = 0 residual as
  the readout there.
- Q1 against GI measures recipient–supply co-assembly, not the existence of fan-in. Fan-in
  existence should be stated as O − E against the endpoint-level closed form or N2c, per N-bin.

**Q2 mechanism: donor-side leaky-vs-private co-carriage excess, conditional on log genome bp and
GTDB clade, with leaky and private lists pre-registered from chemistry and transport alone.**
Yes, supported. This is the interaction I in §3.3, and it is the design I would make primary.

Five caveats:
1. Use exact-count, donor-only nulls; Bernoulli counts bias Σ.
2. Match the leaky and private classes on prevalence. Unmatched classes saturate the logit: in a
   first revised run with private prevalence 0.5, Σ_private hit the clip at S = 30.
3. **The clade rank must be pre-registered and not too fine.** *This is reasoning, not a
   simulation result: the simulator has no phylogeny.* Stacking is a lineage-level
   outcome: the helper lineage is the stacker. Conditioning on genus or species can absorb the
   very signal tested. Suggested: family or order as a random effect, with a sensitivity
   analysis at class.
4. Leakiness labels must not be inferred from co-occurrence (§2.3 item 2).
5. The model only shows that I responds to stacking. Whether real private functions are
   phylogenetically more clustered than leaky ones is an empirical risk. Report Σ_private on its
   own; a large |Σ_private| flags residual clade structure.

## 5. Limits of this note

- Donors are treated as independent draws for Props. 1–2. Real events have donor–donor
  structure, such as redox-gradient co-assembly. The GI null handles this empirically, and the
  formulas do not.
- The Black Queen model has no population dynamics beyond the abundance weight e^{−βn}. It has no
  competitive exclusion, no horizontal gene transfer, no regain of functions and no phylogeny.
  It assumes every species needs every leaky function. m is a crude stand-in for leakiness.
  Private retention is community-independent by construction.
- The model's O± functions are pathway-level, with no OR-routes. The endpoint co-supply case is
  handled only analytically (§1.3).
- c >= 3 (cycles), ambient supply E and unbinned donors are outside the scope of the formulas.
- Cross-family reviews so far: Astra (Codex), who checked the closed form and the pair
  derivative by independent enumeration, and Carbon (Gemini). Neither has yet re-reviewed the
  revised simulations.

## 6. Hashes

sha256, first 12 characters (recompute with `sha256sum`), all in `scratch/deep/`, revised run:

| File | sha256 | Seed |
|---|---|---|
| `fanin_closed_form.py` | 4338c5c47ff1 | `default_rng(20261118)` (unchanged) |
| `fanin_closed_form.out` | 2d9cecca2e76 | |
| `fanin_curveball.py` | f8213b45f159 | metacommunity `default_rng(20261119)`; workers `SeedSequence(20261121).spawn(8)` |
| `fanin_curveball.out` | bccdb96a7282 | |
| `bq_stacking_sim.py` | 47212b49099d | `SeedSequence(20261123).spawn(40)`; bootstrap uses each cell's seed |
| `bq_stacking_sim.out` | 40f59da07497 | |
| `check_carbon_inc.py` | b665f27ba201 | none (exact enumeration) |
| `check_carbon_inc.out` | 3be341a92a8b | |

Superseded first-version hashes, kept for the record:

| File | First-version sha256 |
|---|---|
| `fanin_curveball.py` | 74e8a5722ef2 |
| `fanin_curveball.out` | b0cbbee13660 |
| `bq_stacking_sim.py` | 08ee467ef884 |
| `bq_stacking_sim.out` | 518cc9f76e29 |
| this note | 26f8788d925d |

Wall times: under 3 s for the closed form; about 6 min for the curveball run (8 processes);
about 3 min for the Black Queen simulation (14 processes).

## 7. Revisions after cross-family review (Astra 248967796767506432)

This section also covers the second review, by Carbon (Gemini, relay 248968598093479936).

**Astra (Codex), 248967796767506432.** Astra independently confirmed the closed form and the pair
derivative by enumeration at N = 2 and 3, and the values q = 0.3: 0.333719 at N = 10 and 0.022995
at N = 40. Fixes applied:
1. **Recipient row in the curveball null.** `fanin_curveball.py` previously shuffled the full event
   matrix, including X, and then deleted X. That changed the donors' margins. Now X is removed
   first, and every null runs on the donor-only matrix with exact donor row and column sums, using
   one chain per recipient. The §2.2 numbers are updated:
   - N2 inc excess: +0.31 to +0.96, was +0.25 to +0.96;
   - observed inc vs R×R log-OR: 1.79–14.3, was 0.94–13.7.

   Conclusions are unchanged.
2. **Bernoulli independent null.** `bq_stacking_sim.py` drew the IND null with Bernoulli cells, so
   column sums matched only in expectation. IND now places each function's **exact donor count**
   on random donor rows. It is computed on the donor-only matrix per recipient, for all
   replicates, not just the first 60.
3. **The 1.4–4.8 figure.** It was the observed compatible-pair d = 1 rate against IND, not a
   leaky-minus-private interaction. The model now has a private-function class (Kp = 6,
   community-independent retention, prevalence matched to the leaky class), and the interaction
   I = Σ_leaky − Σ_private is calibrated against its own null:
   - β = 0: I between −0.15 and 0.12, every 95% CI containing 0;
   - β >= 0.25: I = 1.36–11.1, every lower bound above ln 2.

   §2.3 now states that the test needs leakiness labels independent of the co-occurrence data.

   A first revised run used private-need probability 0.5. That saturated P(d = 1 | priv) at
   S = 30 (Σ_private hit the logit clip), so the private class was prevalence-matched (need
   probability m/S) and rerun with a new seed. Only the matched run is reported.
4. **Table E** now states that it assumes independent functions throughout. Its log-OR > 12 comes
   from unequal marginal prevalence alone, and structural exclusion is a separate argument
   (Prop. 2).

**Carbon (Gemini), 248968598093479936.**
1. **§3.2 P(d = 1 | inc) was the unconditional probability.** It is now the eligible-X pooled
   rate, n_k(S − 1) / [n_k(S − 1) + (AB − n_k)(S − 2)]. This was verified by exact enumeration,
   which matches Carbon's 0.03734 and 0.15152. The S = 20 example becomes log-OR 1.75; it was 1.8.
2. **Pathway-level vs endpoint-level exclusion.** The §1.3 corollary and the §2.2 world assumed
   p11 = 0. For P2's R-I2 and R-I3 *endpoints*, OR-routes (anaerobic cbi + bza with dsr or nif)
   give p11 > 0 at d = 1 with no compartmentalisation. A new subsection in §1.3 quantifies this:
   p11 = 0.01 already cuts P(d = 2) from 1.00 to 0.45 at N = 80. Text in §§2.1, 2.2, 2.3, 3.1 and
   4 now distinguishes the two levels, and a new recommendation 7 covers endpoint-level scoring.

   The main conclusion holds, and Carbon agrees with it: the inc vs comp contrast is a
   donor-composition quantity (q1, q2, p11, N), not a recipient-level mechanism.
3. **Carbon's proposed P2 design** is assessed in §4.1. Both parts are supported, with caveats:
   - the GI null needs an event bootstrap and has low power at N >= 40;
   - the leaky and private classes must be prevalence-matched;
   - the clade-conditioning rank must not be so fine that it absorbs the stacking lineage;
   - the leakiness labels must be independent of co-occurrence.

**Astra (Codex) closing review, 248971150990163968, of version 248fac81a0d1.** Two wording
edits, with no numbers changed:
1. **Interval scope.** §3.3 and §4 now state that the bootstrap intervals are for simulated mean
   contrasts across replicate communities and do not establish cohort-level power. §4.1 caveat 3,
   on clade rank no finer than family or order, is labelled as reasoning, not a simulation
   result, because the simulator has no phylogeny.
2. **Two different IND nulls.** IND is now labelled separately in each section:
   - §2 (`fanin_curveball.py`) uses an i.i.d. plug-in probability;
   - §3 (`bq_stacking_sim.py`) uses fixed-count permutations.
