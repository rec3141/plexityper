# P3 design and power: is δ = 0.05 nats resolvable with ~10–20 held-out events?

Author: Claude (Opus) subagent, for claude-plexi. 2026-11-18. **Design-only and outcome-blind.**
Nothing under `data/` was opened and no abundance table was read. The only inputs were the
current `protocols/P3-plexitype-recurrence.md` (read in full) and Astra's metadata screen in
`discussion/plexi-concept-review.md` (msg 248915045740101632: 327 eligible Namao rows; for
2024–25, 17 shared stations and 21 station–month cells before reconciliation). The correction
msg 248921820019212288 that the protocol cites does not appear in the review file, so it was
not read.

> **Revised after cross-family review (Astra 248969908310163456; Carbon 248973006869020672).**
> See the final section. It supersedes the verdict rule in §3(b), the coverage reading in §2.5,
> and every oracle-discovery power figure in §2–3. Sections 1–3 are kept as reviewed
> (c7c81057a4cf) for the record.

Code (numpy/scipy only, run with `python3 -I`), in `scratch/deep/`:

| file | sha256 (12) | role |
|---|---|---|
| `p3_analytic.py` | afc7915a18ac | §1 derivations and tables (deterministic); output in `out/analytic.txt` (09c7a22df837) |
| `out/p3_sim_v1_used_for_main_grid.py` | 08ddd66b67f1 | §2 main grid and sensitivity runs (the v1 code that produced them) |
| `p3_sim.py` | 8393aeed6e2f | v2: adds the "rival" outsider factor and a best-outsider control (§2.5). The extra RNG draw means v2 does **not** reproduce v1 numbers bit-for-bit |
| `p3_tables.py` | 949be92c95ea | builds `out/tables.md` (d7277e89afd3) from `out/*.json` |

Seeds: each scenario has a fixed `seed` (listed in `SCENARIOS`, e.g. `rho0.5_r0.2` → 1520). Each
of its 1000 replicates uses `SeedSequence(seed).spawn(1000)`. The bootstrap (B = 999) draws
from the replicate's own generator.

---

## 1. Derivations (analytic, no simulation)

### 1.1 What the score measures
Write the outcome as Y, the member predictors as M, and the covariates as X. Let q be the true
conditional, p_e the reduced model (env only) and p_m the full model (env + members). The
expected held-out log-score difference per collection is

  E[log p_m(Y|M,X) − log p_e(Y|X)] = I(Y; M | X) − [R(p_m) − R(p_e)],

where R(p) = E_X,M KL(q ‖ p) is the excess risk of a fitted model. With correct specification
and MLE, the excess out-of-sample log loss is about dim/(2n) (the AIC/Takeuchi result).
Adding k member predictors therefore gives

  **G(Mm − Me) ≈ I(Y; M | X) − k / (2 n_eff)**,  n_eff = effective number of independent
  training collections (≈ events × design-effect correction).

- Ridge or JSDM shrinkage lowers the penalty, at the cost of bias toward 0.
- For Mm − Mn the predictor counts match, so the penalties cancel to first order:
  **G(Mm − Mn) ≈ I(Y; M | X) − I(Y; N | X)**, where N is the matched nonmember set.
- Consequence: **Mm − Me is the harder of the two contrasts.** It is the only one that pays the
  overfitting penalty net. All the simulations below confirm this.

Scale: H(Y) ≤ ln 2 = 0.693 nats. So δ = 0.05 nats means removing 7.2% of the outcome's
entropy at p0 = 0.5, and 10% at p0 = 0.2 or 0.8.

### 1.2 Model A: one binary summary predictor
Take S ~ Bern(1/2) and logit P(Y=1|S) = a ± β/2, with a set so that the marginal equals p0.
The log-odds shift β needed for I(Y;S|X) to reach a target:

| p0 | 0.2 | 0.3 | 0.4 | 0.5 | 0.6 | 0.7 | 0.8 |
|---|---|---|---|---|---|---|---|
| β for I = 0.02 (OR) | 1.04 (2.8) | 0.89 (2.4) | 0.83 (2.3) | 0.81 (2.2) | 0.83 | 0.89 | 1.04 |
| **β for I = 0.05 (OR)** | 1.74 (5.7) | 1.45 (4.3) | 1.33 (3.8) | **1.30 (3.7)** | 1.33 | 1.45 | 1.74 |
| β for I = 0.08 (OR) | 2.38 (10.8) | 1.90 (6.7) | 1.72 (5.6) | 1.67 (5.3) | 1.72 | 1.90 | 2.38 |
| sd of per-collection difference d at I = 0.05 | 0.30 | 0.31 | 0.31 | 0.31 | 0.31 | 0.31 | 0.30 |

### 1.3 Model B: a factor-probit module, linked to the P3 clique rule
Here the outcome and k predictors share an equicorrelated residual latent correlation ρ, given
X, with every marginal prevalence equal to p0. ρ is the quantity that P3 thresholds at > 0.3.
I(Y; M₁..M_k | X) is computed exactly by Gauss–Hermite quadrature over the shared factor:

| k | ρ | p0 = 0.2 | 0.3 | 0.4 | 0.5 |
|---|---|---|---|---|---|
| 3 | 0.2 | 0.013 | 0.017 | 0.019 | 0.020 |
| 3 | **0.3** | 0.028 | 0.036 | 0.040 | **0.042** |
| 3 | 0.4 | 0.048 | 0.061 | 0.068 | 0.071 |
| 3 | 0.5 | 0.074 | 0.092 | 0.103 | 0.106 |
| 3 | 0.6 | 0.106 | 0.131 | 0.146 | 0.151 |
| 5 | 0.3 | 0.038 | 0.048 | 0.053 | 0.055 |
| 5 | 0.5 | 0.091 | 0.113 | 0.125 | 0.129 |

(The table is symmetric in p0 ↔ 1 − p0. k = 2 is in `out/analytic.txt`.)

ρ needed for a **net** gain of 0.05, i.e. I = 0.05 + k/(2n) with MLE-level penalty:

| k, p0 | oracle (n = ∞) | n = 40 | n = 60 | n = 100 |
|---|---|---|---|---|
| 3, 0.5 | 0.33 (pairwise OR 2.4) | 0.45 | 0.41 | 0.38 |
| 3, 0.3 | 0.36 | 0.49 | 0.45 | 0.41 |
| 5, 0.5 | 0.28 | 0.46 | 0.41 | 0.36 |
| 5, 0.2 | 0.35 | 0.56 | 0.50 | 0.45 |

Precision: under Model B the per-collection score difference has sd 0.28–0.51 nats. The number
of **independent** test collections needed so that the lower 95% bound exceeds δ = 0.05 with
80% power is ((1.96 + 0.84)·sd/(I − δ))²:

| k | ρ | p0 | I | sd(d) | n for 80% power (oracle, no penalty) |
|---|---|---|---|---|---|
| 3 | 0.4 | 0.5 | 0.071 | 0.36 | ≈ 2,400 |
| 3 | 0.5 | 0.5 | 0.106 | 0.43 | ≈ 450 |
| 3 | 0.6 | 0.5 | 0.151 | 0.49 | ≈ 185 |
| 5 | 0.5 | 0.5 | 0.129 | 0.46 | ≈ 260 |
| 5 | 0.6 | 0.5 | 0.177 | 0.51 | ≈ 127 |

### 1.4 Reading of §1
δ = 0.05 nats per collection is **neither trivial nor unreachable in principle, but it is
large for this design.**
- It equals the full oracle information of a k = 3 clique sitting just above P3's own selection
  threshold (ρ ≈ 0.33).
- With 40–60 training collections, the estimation penalty k/(2n) ≈ 0.025–0.04 nats eats most of
  a ρ ≈ 0.3–0.4 clique's information. Netting 0.05 then needs ρ ≈ 0.4–0.5.
- For comparison, a single-taxon odds ratio of 3.7 at p0 = 0.5 is a strong residual
  co-occurrence for microbes after environment, though not an unheard-of one.
- The binding constraint is precision, not the size of the margin. Even an oracle model needs
  hundreds of independent test collections to put a lower bound above 0.05. The 2025 test year
  has about 20 events.

---

## 2. Simulation (design-only, synthetic data)

### 2.1 Generative model
Probit latent per taxon (full spec in the docstring of `p3_sim.py`):
- **Environment.** A measured gradient g1, plus a latent background gradient g2 seen only as
  g2hat with corr c = 0.85. These stand in for the background-taxon gradient projected from
  2024 loadings.
- **Season and depth.** Season is a cosine niche over May–Oct. Library depth enters as logL.
- **Clustering.** Each event (cast) has D = 2 depths. Env is shared within the event, with
  sd 0.3 between depths. The module factor and the taxon residuals share within-event
  correlation r.
- **Cruise and year effects.** Three cruises per year, with a per-taxon random effect (sd 0.25).
  Per-taxon year intercept shift sd 0.2. A test-year env shift of +0.3.
- **Module.** K = 4 taxa (1 outcome at prevalence 0.5, plus 3 predictors at 0.3–0.7) with
  residual latent corr ρ. Members partly share the outcome's env niche.
- **Distractors.** 32 ordinary taxa plus 8 strongly env-correlated outsiders: they load heavily
  on g2, the axis the analyst sees only noisily, so they carry real leaked env information about
  the outcome. Prevalence 0.2–0.8.
- **v2 "rival" scenario.** The 8 outsiders additionally share a residual factor with the
  outcome only (outcome–rival residual corr 0.35). This models Astra's "strong predictors of the
  outcome that fail a clique edge".
- **Training size.** 30 events (60 collections); sensitivities at 20 and 45 events. Test year:
  30 events, nested subsets E = 5, 10, 15, 20, 30.

### 2.2 Analysis model (mirrors P3 §4–5 as closely as numpy allows)
- **Models.** Ridge logistic, i.e. a N(0, 1/λ) prior on standardised slopes, with λ = 1 at
  baseline and λ = 5 as a sensitivity.
  - Me = g1, g2hat, sin/cos(month), logL.
  - Mm = Me + member detections.
  - Mn = Me + matched nonmember detections.
- **Matching.** Nearest Mahalanobis neighbour on (logit 2024 prevalence, circular peak-month
  distance), without replacement. 20 random predictor orders; Mn is the mean of the 20 sets'
  log scores. Abundance was not simulated, so it was not matched. The control pool is all
  eligible nonmembers, outsiders included. "Screened" is the same pool with the outsiders
  removed, an oracle version of P3's ρ ≤ 0.3 screen.
- **Scoring.** Per-collection log score, averaged per target, then per event, then across
  events. Event-level percentile bootstrap, B = 999.
- **Verdicts.**
  - S (supported): LB(Mm−Me) > δ and LB(Mm−Mn) > δ.
  - F (falsified): UB of either contrast < δ.
  - I (inconclusive): otherwise.
- **Truth.** The conditional-on-training gain, scored on 3,000 fresh test-year collections
  drawn from new cruises.
- **Aggregations.** "prim" = one outcome member, as drafted. "rot" = every member masked in
  turn, with the others as predictors, averaged over targets per collection (§3d).
- **Not modelled.** Discovery error. Members are known to the analyst, so Mm carries no clique
  winner's curse. Every power figure below is therefore **optimistic**.

### 2.3 What the fitted models actually gain (truth, mean ± sd over 1000 training draws; nats/collection)

| scenario | Mm−Me | Mm−Mn | Mm−Mn (screened pool) |
|---|---|---|---|
| null ρ = 0, r = 0.2 | −0.034 ± 0.040 | −0.006 ± 0.051 | −0.005 |
| null ρ = 0, r = 0.8 | −0.047 ± 0.057 | −0.008 ± 0.076 | −0.007 |
| small ρ = 0.2, r = 0.2 | −0.018 ± 0.045 | +0.011 ± 0.056 | +0.012 |
| moderate ρ = 0.35, r = 0.2 | +0.006 ± 0.048 | +0.035 ± 0.058 | +0.035 |
| large ρ = 0.5, r = 0.2 | +0.036 ± 0.060 | +0.064 ± 0.067 | +0.066 |
| very large ρ = 0.7, r = 0.2 | +0.095 ± 0.067 | +0.126 ± 0.076 | +0.127 |
| ρ = 0.5, ridge λ = 5 | +0.061 ± 0.034 | +0.070 ± 0.041 | +0.072 |
| null, ridge λ = 5 | −0.012 ± 0.022 | −0.004 ± 0.029 | −0.001 |
| ρ = 0.5, no cruise/year effects | +0.058 ± 0.036 | +0.087 ± 0.047 | +0.088 |
| ρ = 0.5, 45 training events | +0.042 ± 0.049 | +0.063 ± 0.057 | +0.064 |
| v2 rival, null | −0.040 ± 0.052 | **−0.018** ± 0.067 | −0.006 |
| v2 rival, ρ = 0.5 | +0.036 ± 0.062 | **+0.056** ± 0.074 | +0.069 |

- **Penalty.** The null Mm−Me of −0.034 at 60 training collections matches k/(2n) = 3/120 =
  0.025 plus within-event correlation. At r = 0.8 the effective n falls, and the penalty grows
  to −0.047.
- **Size.** At realistic training size, a "moderate" module at the clique threshold (ρ = 0.35)
  nets about 0 on Mm−Me. A large one (ρ = 0.5) nets about 0.036, which is below δ = 0.05. Only
  ρ ≈ 0.7 clears δ on average.
- **Year-to-year spread.** The truth itself has sd 0.04–0.08 across training/test-year draws.
  One forward split is one year-level sample, so the event bootstrap cannot represent this
  spread.

### 2.4 Verdict probabilities S/F/I as drafted (one outcome member, δ = 0.05, primary pool, ridge λ = 1)

| scenario | E=5 | E=10 | E=15 | E=20 | E=30 |
|---|---|---|---|---|---|
| null, r = 0.2 | .00/.44/.55 | .00/.50/.49 | .00/.57/.43 | .00/.63/.37 | .00/.75/.25 |
| null, r = 0.8 | .01/.36/.62 | .00/.40/.60 | .00/.46/.54 | .00/.51/.49 | .00/.61/.39 |
| small ρ = 0.2, r = 0.2 | .01/.29/.70 | .01/.28/.71 | .01/.32/.68 | .00/.36/.64 | .00/.44/.56 |
| small ρ = 0.2, r = 0.8 | .02/.28/.70 | .01/.27/.73 | .01/.30/.69 | .01/.32/.68 | .00/.40/.59 |
| moderate ρ = 0.35, r = 0.2 | .06/.22/.72 | .03/.16/.80 | .02/.15/.83 | .01/.17/.82 | .01/.20/.79 |
| moderate ρ = 0.35, r = 0.8 | .04/.21/.75 | .02/.17/.81 | .02/.17/.81 | .01/.20/.79 | .01/.22/.77 |
| large ρ = 0.5, r = 0.2 | .06/.10/.84 | .06/.07/.88 | .06/.06/.89 | .06/.06/.88 | .06/.07/.87 |
| large ρ = 0.5, r = 0.8 | .07/.13/.80 | .06/.09/.85 | .05/.07/.88 | .05/.08/.87 | .04/.09/.87 |
| v. large ρ = 0.7, r = 0.2 | .15/.07/.78 | .14/.03/.83 | .14/.03/.83 | .14/.03/.83 | .15/.02/.82 |
| v. large ρ = 0.7, r = 0.8 | .14/.07/.79 | .11/.03/.86 | .12/.03/.85 | .13/.03/.85 | .15/.03/.82 |

- **Support is flat in E.** P(S) does not rise with E up to 30. The mean 95% half-width of
  Mm−Me at E = 30 is 0.13–0.18 nats (primary aggregation), against a δ of 0.05.
- **False support is ~0.** At δ = 0.05, false support under the null is ≤ 1% in every scenario.
- **Falsification works.** The null is falsified 40–75% of the time at E = 10–30.
- **Real modules end inconclusive.** At E = 10 the design is inconclusive 83–88% of the time
  for real modules (ρ ≥ 0.35). It also "falsifies" a ρ = 0.5 module 7–9% of the time, which is
  correct on its own terms, since that module's true net gain is below 0.05.

Sensitivities: K = 6 members, 20 or 45 training events, and c = 0.6 all change P(S) at E = 10
by ≤ 0.02. Removing cruise/year effects changes it by ≤ 0.01. Full tables are in
`scratch/deep/out/tables.md`.

### 2.5 Interval calibration
Percentile event-bootstrap coverage of the truth is **0.84–0.89** (nominal 0.95) when cruise
and year effects are present. It is 0.89 at E = 10 and 0.93 at E = 30 without them.

Under the null, P(LB(Mm−Mn) > 0) is 0.06–0.10 against a one-sided nominal 0.025. The intervals
are therefore anti-conservative. Two causes:
- shared cruise effects that event resampling cannot see, with only about 3 cruises per year;
- percentile intervals at small E.

The joint requirement (both contrasts) keeps false support at δ = 0 to 1–4%.

### 2.6 Improved design: rotate all members as targets, stronger shrinkage

P(S) / P(F) / P(I), "rot" aggregation, ridge λ = 5, primary pool:

| δ | scenario | E=10 | E=20 | E=30 |
|---|---|---|---|---|
| 0.05 | null | .00/.69/.30 | .00/.86/.14 | .00/.93/.07 |
| 0.05 | ρ = 0.5 | .05/.07/.88 | .05/.05/.90 | .06/.04/.90 |
| 0.05 | ρ = 0.7 | .18/.02/.80 | .26/.01/.73 | .37/.01/.62 |
| 0.02 | ρ = 0.5 | .14/.02/.84 | .18/.01/.80 | .24/.01/.75 |
| 0.02 | ρ = 0.7 | .31/.01/.68 | .48/.00/.51 | .61/.00/.39 |
| 0 | null | .02/.14/.84 | .01/.15/.84 | .01/.17/.83 |
| 0 | ρ = 0.5 | .22/.01/.78 | .33/.01/.67 | .44/.00/.55 |
| 0 | ρ = 0.7 | .40/.00/.60 | .62/.00/.38 | .75/.00/.25 |

- **Rotation alone** (λ = 1) cuts the CI half-width by about 35–40%. For example, null r = 0.2
  at E = 30 goes from 0.07/0.10 to 0.05/0.06; ρ = 0.5 goes from 0.13/0.15 to 0.09/0.10. The
  mean estimate is unchanged under the null (−0.036 vs −0.036), so rotation adds no bias there.
- **λ = 5 shrinkage** cuts the null penalty from −0.034 to −0.012. It raises the ρ = 0.5 truth
  from 0.036 to 0.061.

### 2.7 Specificity contrast (Mm−Mn) with strong outsiders in the control pool
- **Env-correlated outsiders** (base scenarios) shift Mm−Mn by ≤ 0.006 nats relative to the
  screened pool (§2.3).
- **Rival outsiders** (v2) that share residual correlation 0.35 with the outcome lower Mm−Mn by
  0.012–0.013. For example, at ρ = 0.5 it falls from 0.069 to 0.056.
- **Verdicts barely move** in either case: at E = 10 and δ = 0.05, the rival scenario gives S = .05 and
  F = .09, against S = .07 and F = .08 for the matching v2 base scenario. The reason is that Mm−Me, not Mm−Mn, is the binding contrast.
- **Matching dilutes the outsiders.** Matching on prevalence and season picks outsiders only in
  proportion to their share of the pool (8/40). Averaging over 20 sets dilutes them further. The
  primary-pool Mn therefore tests "better than prevalence/season-matched typical taxa", not
  "better than the strongest outsiders".
- **A naive "best outsiders" control fails.** Taking the top-k nonmembers by 2024 in-sample
  gain over Me overfits badly (winner's curse). Mm−Mbest = **+0.08 under the null** in the base
  scenario and +0.019 in the rival scenario. That control is biased in favour of membership and
  must not be used unless outsiders are selected by the same cross-validated procedure as the
  clique.

---

## 3. Recommendations

**(a) The ≥ 10-event gate is not adequate as a power criterion.**
- **As drafted.** With δ = 0.05 and E = 10:
  - a very strong module (ρ = 0.7) is supported 11–14% of the time;
  - a module at the clique threshold (ρ = 0.35) is supported 2–3% of the time;
  - a null module is falsified only 40–50% of the time.
  - Even in the best variant (rotation + λ = 5), E = 30 gives P(S) = 0.37 at ρ = 0.7.
- **Keep 10 as a feasibility floor**, as the protocol already says, but add a **precision gate**:
  run this design simulation with the reconciled event counts and the frozen model, and
  pre-register the verdict probabilities.
- **Rule for the likely outcome.** If P(S | ρ = 0.5) < 0.5, state in advance that the realistic
  outcomes are "falsified (no gain ≥ δ)" or "inconclusive for lack of design". Do not reach for
  "supported".
- **Two-sided reading.** With about 20 test events (Lake Winnipeg 2025), the test can rule out
  large gains, but it cannot establish moderate ones.
- **More events.** More test events (other held-out years or lakes, combined at the event level
  with a frozen clique) are the only real lever. Bootstrap tricks are not.

**(b) A calibrated δ.**
- **Split the margins.**
  - **Falsification margin δ_F = 0.05.** "Upper bound < 0.05" rules out a gain equal to the
    full oracle information of a k = 3, ρ ≈ 0.33 clique. That is substantive, and it is
    reachable: null falsified 60–93% at E = 10–30 with rotation, with or without λ = 5.
  - **Support threshold δ_S = 0.** Both lower bounds > 0, plus a pre-registered minimum point
    estimate of Mm−Me ≥ 0.02 so that "supported" is not trivial. δ_S = 0.02 on the bound is
    the alternative if the team prefers bounds only.
- **Same δ in both rules.** Keeping 0.05 for both makes "supported" nearly unreachable
  (≤ 0.15 at E ≤ 30 as drafted).
- **Inconclusive band.** The band 0 < effect < 0.05 stays explicitly inconclusive.
- **Fix the interval method.** δ should be frozen together with it: studentized or BCa
  intervals, plus cruise as a fixed effect or a leave-one-cruise-out check. The percentile
  bootstrap here under-covers (0.84–0.89).

**(c) Mm−Mn power with strong outsiders retained.**
- Mm−Mn has usable power relative to Mm−Me, because both models pay the same overfitting
  penalty. It is **never the binding contrast** in these runs. Retaining env-correlated or rival
  outsiders costs about 0.006 and 0.013 nats respectively.
- The real problem is validity, not power. Prevalence/season matching plus averaging over 20
  sets rarely confronts the clique with the outsiders that matter.
- If membership specificity against "strong predictors that fail a clique edge" is the claim:
  - add a pre-registered **secondary contrast** against the k outsiders with the highest 2024
    residual correlation with the outcome;
  - select them by the **same leave-one-event-out procedure** used to rank cliques;
  - do **not** select them by in-sample gain, which gave a +0.08 bias toward Mm under the null.

**(d) Aggregation changes that help without bias.**
1. **Rotate the masked member.** Freeze all clique members as outcomes in turn, with the others
   as predictors. Average per collection over targets, then per event.
   - It cuts the CI half-width by about 35–40%, and under the null the mean is unchanged.
   - The estimand becomes the module-average gain, which arguably fits the "plexitype" claim
     better than one hand-picked outcome. Keep the single 50%-prevalence member as a reported
     secondary.
2. **Stronger, pre-specified shrinkage on member coefficients.** In the JSDM this means a
   tighter prior on the residual covariance or factor loadings. It removes about two-thirds of
   the overfitting penalty that makes Mm−Me the binding contrast. The prior must be frozen, not
   tuned on 2024 ranking.
3. **Keep event-level averaging** (equal weight per event) and resample events. With about 3
   cruises, a cruise-level bootstrap is degenerate. Model cruise in Me and report a
   leave-one-cruise-out sensitivity instead.
4. **Freeze how Mn averages over the 20 sets.** Mean of log scores and log of mean probability
   differ by Jensen's inequality, and the protocol does not say which.
5. **Not recommended.** Pooling collections without event averaging (it understates variance
   when r is high), or choosing among aggregations after seeing 2025.

## Assumptions and limits
- **Synthetic, oracle discovery.** The model is probit, with additive-logistic analysis models.
  Discovery is oracle (members known), so there is no clique winner's curse, and the power
  figures are upper bounds.
- **Gradient.** The latent gradient is modelled as a noisy observation, not refit from
  background taxa.
- **Sizes and depths.** About 30 training events with 2 depths each. Real event counts await
  P3's reconciliation gate.
- **Within-event correlation.** r was set at 0.2 or 0.8 and is not estimated from data.
- **Oracle-to-simulation gap.** The oracle I in §1 (e.g. 0.106 for k = 3, ρ = 0.5) exceeds the
  large-training simulated gain (≈ 0.06 with 3,000 training events in a pilot). Causes:
  env-varying conditional prevalence, the additive-logistic approximation to the factor
  conditional, and year/cruise shifts. The §1 numbers are optimistic bounds.
- **Abundance and depth.** Abundance-based detectability matching and the true library-depth
  offset were reduced to a logL covariate.

---

## Revisions after cross-family review (Astra 248969908310163456)

- Reviewed version: c7c81057a4cf.
- Also folds in Carbon's input (relay 248973006869020672).
- Still design-only: nothing under `data/` was opened and nothing was posted.

### New code and outputs (`scratch/deep/`)

| file | sha256 (12) | role |
|---|---|---|
| `p3_sim_v3.py` | fb694e451f9f | discovery inside every replicate; per-year/per-cruise effects; two estimands; two verdict rules |
| `p3_tables_v3.py` | 41775b1afa3f | builds `out_v3/tables_v3.md` from `out_v3/*.json` |
| `out_v3/tables_v3.md` | fe1a6248ac99 | all v3 tables, with columns S/F/I/N (N = no clique found) |

**Seeds and replicates.**
- Scenario seeds come from the formulas in `SCENARIOS`: 11000 + 10·(100ρ) + λ, or 12000 + 100ρ
  for the rival scenarios. Each scenario runs 600 replicates via `SeedSequence(seed).spawn(600)`.
- A first 1000-replicate pass was stopped to stay within the compute budget. Its one finished
  scenario is archived as `out/v3_rho0.0_lam1.0_1000reps_aborted_run.json`. It agrees with the
  600-replicate run: P(no clique) 0.163 vs 0.153.

**Scenarios run.** Each took about 3.5 min on 16 cores.
- Null ρ = 0 at λ = 1 and λ = 5.
- ρ = 0.35, 0.5 and 0.7 at λ = 5.
- ρ = 0.5 at λ = 1.
- Rival scenario at ρ = 0.5, λ = 5.
- **Not run (budget):** λ = 1 at ρ = 0.35 and 0.7, and the rival null.

**Fixed design.** r = 0.5, 30 training events × 2 depths, K = 4 true module members, 32
ordinary taxa and 8 env-correlated outsiders.

### R1. Verdict rule (Astra point 1)
Astra is right. My §3(b) rule (support if LB > 0; falsify if UB < 0.05) can return both
verdicts for the same result, e.g. intervals [.021, .039]. **That rule is withdrawn.** Two
rules were calibrated instead:

- **Rule 1: one common frozen margin δ for both contrasts.**
  - S if LB(Mm−Me) > δ and LB(Mm−Mn) > δ.
  - F if UB(Mm−Me) < δ or UB(Mm−Mn) < δ.
  - I otherwise.
  - N if no clique is found.
  - The verdicts exclude each other by construction: LB > δ implies UB > δ.
  - If specificity is unresolved (< 2× matched controls), only Mm−Me can falsify. This never
    happened here (P(unresolved) = 0 in every scenario; on average 38 taxa were eligible).
- **Rule 2: two separate claims.**
  - C+ "positive gain": both LBs > 0.
  - C− "a gain ≥ 0.05 is ruled out": either UB < 0.05.
  - They can co-occur in principle. Across all scenarios and all E they **never did**, because
    the CI half-widths (≥ 0.06) are wider than the 0.05 window. They are reported as separate
    columns in `tables_v3.md`.

**Recommended primary:** Rule 1 at **δ = 0.02**. Rule 2 is reported as a secondary description.

Calibration of Rule 1 (S/F/I/N, rotated masks, λ = 5, primary pool). These are
**unconditional** rates: they include runs where no clique was found and runs where discovery
picked the wrong taxa.

| δ | scenario | E=10 | E=20 | E=30 |
|---|---|---|---|---|
| 0 | null ρ = 0 | .03/.10/.80/.07 | .03/.10/.80/.07 | .04/.12/.78/.07 |
| 0 | ρ = 0.5 | .07/.09/.82/.02 | .10/.08/.81/.02 | .10/.07/.80/.02 |
| 0 | ρ = 0.7 | .19/.04/.76/.01 | .32/.04/.64/.01 | .37/.04/.58/.01 |
| **0.02** | null ρ = 0 | .01/.16/.76/.07 | .01/.19/.73/.07 | .01/.22/.70/.07 |
| **0.02** | ρ = 0.35 | .03/.15/.75/.07 | .01/.16/.77/.07 | .01/.18/.74/.07 |
| **0.02** | ρ = 0.5 | .04/.12/.82/.02 | .05/.14/.79/.02 | .05/.16/.77/.02 |
| **0.02** | ρ = 0.7 | .13/.06/.80/.01 | .23/.07/.69/.01 | .27/.08/.65/.01 |
| 0.05 | null ρ = 0 | .01/.29/.64/.07 | .00/.37/.56/.07 | .00/.44/.49/.07 |
| 0.05 | ρ = 0.7 | .08/.12/.79/.01 | .12/.14/.73/.01 | .15/.17/.67/.01 |

Why δ = 0.02:
- **At δ = 0,** false support under the null is 3–4%. That exceeds the nominal ≤ 2.5% because
  of the small-E under-coverage in R2 below.
- **At δ = 0.05,** a real ρ = 0.7 module is "falsified" 12–17% of the time, about as often as it
  is supported. The discovered clique's true gain averages only about 0.06 (R3).
- **At δ = 0.02,** false support is ≤ 1% for E ≥ 10 and null falsification is 16–22%. Support
  is low everywhere: ≤ 0.27 even at ρ = 0.7 and E = 30.

Correction to §2.6. The **61%** figure there came from oracle discovery, the common-margin rule
at δ = 0.02, rotation, λ = 5 and E = 30. With discovery inside each replicate, the same cell is
**0.27 unconditional** and **0.43 conditional on recovering the true module** (n = 354/600).

### R2. Estimand and coverage (Astra point 2)
**The defect.** Astra is right. v1/v2 drew one year vector and three cruise vectors per
`gen_year` call. So the 3,000-collection "truth" was one other year with three other cruises,
and it kept between-year variation however large it got. The v1 coverage numbers (0.84–0.89)
compared against that hybrid reference and are **withdrawn**.

**The fix in v3.** `draw_year` creates the year shift, per-taxon year effects and three cruise
effects once per year. Two estimands are computed for the same fitted models and the same
discovered clique:
- **Realised-year.** 1,500 fresh events from the same 2025 year, with the same year effects,
  cruise effects and env shift.
- **Future-year.** 10 independent new years × 150 events, averaged.

**Which estimand P3 targets.** P3 should target the **realised-year estimand**: the transfer
gain in 2025 at the sampled stations, given the frozen 2024 fit. It is the only one that a
single held-out year identifies.
- A future-year (recurrence) claim needs at least two held-out years. The between-year sd of the
  gain is **0.02–0.05 nats** (table below), about the size of δ.
- That sd follows from assumed effect sizes (year sd 0.2, cruise sd 0.25 on the probit scale),
  not estimated ones.

Coverage of the event-bootstrap 95% percentile CI. Each cell lists [Mm−Me, Mm−Mn, Mm−Mn
screened]. Coverage is computed over runs where discovery found a clique.

| scenario (λ = 5, rotation) | E | realised-year | future-year | half-width | between-year sd |
|---|---|---|---|---|---|
| null | 10 | .89/.89/.90 | .89/.90/.90 | .10 | .02 |
| null | 30 | .93/.93/.94 | .91/.90/.90 | .06–.07 | .02 |
| ρ = 0.5 | 10 | .86/.87/.87 | .86/.86/.86 | .11 | .02 |
| ρ = 0.5 | 30 | .91/.92/.91 | .90/.90/.89 | .07 | .02 |
| ρ = 0.7 | 10 | .89/.89/.88 | .89/.87/.88 | .13 | .03 |
| ρ = 0.7 | 30 | .93/.93/.92 | .91/.92/.91 | .08 | .03 |

The primary (single-outcome) aggregation has the same pattern with wider intervals:
- realised-year 0.84–0.91 at E = 10 and 0.90–0.94 at E = 30;
- future-year 0.81–0.92;
- between-year sd 0.03–0.05.

Reading:
- **Realised year.** The remaining shortfall (0.86–0.90 at E = 10, 0.91–0.94 at E = 30) is
  small-sample behaviour of the percentile interval. With the year and cruise effects fixed,
  events are the right resampling unit. A studentised event-level interval is a legitimate fix
  for this estimand.
- **Future year.** Coverage is only modestly lower here, but only because the assumed
  between-year sd is small next to these wide intervals. **I withdraw §2.5's suggestion that
  BCa or studentisation helps for cross-year inference.** No event-level method can supply a
  missing year or cruise cluster.

### R3. Discovery inside each replicate (Astra point 3)
In v3, `members = np.arange(K)` is gone. Every replicate re-runs P3 step 2b on training data
only:
- **Eligibility:** prevalence 20–80%, with at least 5 detections and 5 non-detections.
- **Residual correlation:** pairwise bivariate-probit likelihood on a ρ grid with a flat prior.
  Margins come from the ridge Me fits (logit/1.7). An edge needs posterior mean > 0.3 and a
  2.5% quantile > 0.
- **Cliques:** maximal cliques of size 3–8 (Bron–Kerbosch).
- **Ranking:** CV(Mm−Me) over 10 blocked folds, each holding out 3 training events. Tie-breaks
  as in P3: smaller group, then higher mean ρ, then lexicographic ids.
- **Outcome member:** prevalence closest to 50%.
- **Matched controls:** drawn from eligible taxa outside the discovered clique.
- **Screened pool:** the estimated residual correlation with the outcome is ≤ 0.3, i.e. P3's
  real screen, no longer the oracle one.
- **Rotation:** every discovered member is masked in turn. Scores are averaged over targets
  within each collection, then within each event, **before** event resampling.

Discovery outcomes over 600 replicates (from `tables_v3.md`). "Hit" means the outcome is a true
module member and the clique contains ≥ 3 true members.

| scenario | P(no clique) | mean clique size | P(outcome is true member) | P(clique contains an env outsider) | P(hit) | realised truth Mm−Me of discovered clique, rotated |
|---|---|---|---|---|---|---|
| null, λ = 1 | .153 | 3.06 | .10 | .74 | 0 | −0.063 |
| null, λ = 5 | .070 | 3.23 | .11 | .86 | 0 | −0.003 |
| ρ = 0.35, λ = 5 | .067 | 3.21 | .20 | .81 | .055 (33 runs) | −0.001 |
| ρ = 0.5, λ = 1 | .075 | 3.11 | .41 | .54 | .243 | −0.038 |
| ρ = 0.5, λ = 5 | .018 | 3.24 | .35 | .66 | .213 (128) | +0.012 |
| ρ = 0.7, λ = 5 | .008 | 3.42 | .70 | .37 | .590 (354) | +0.062 |
| rival ρ = 0.5, λ = 5 | .002 | 4.42 | .17 | .90 | .082 | +0.079 |

Findings:
1. **Pseudo-cliques under the null.** Under the null, P3's rule returns a clique 85–93% of the
   time, almost always made of env-correlated taxa. Their residual correlation is the part of
   the latent gradient that Me's noisy g2hat misses.
   - With λ = 5 these pseudo-cliques gain about 0 out of sample. Rule 1 then mostly returns I,
     or F at a larger δ.
   - With λ = 1 they lose 0.06 nats and are falsified.
   - **"No clique found" will be rare**, so its absence is not evidence of a design problem.
2. **The true module is often missed.** Discovery recovers it only 6% of the time at ρ = 0.35,
   21–24% at ρ = 0.5 and 59% at ρ = 0.7. Pseudo-cliques with higher CV gain often outrank it.
   The discovered clique's realised gain at ρ = 0.5 is 0.012, versus an oracle 0.061.
3. **Rival scenario.** The outcome plus its residual partners form genuine, if unintended,
   cliques. Their positive gain is real conditional dependence, not an artefact.
4. **Operating characteristics at δ = 0.02** (rotation, λ = 5, E = 10/20/30):

| scenario | unconditional S/F | conditional on a clique S/F | conditional on hit S/F |
|---|---|---|---|
| null | .01/.16, .01/.19, .01/.22 | .02/.17, .01/.21, .01/.24 | — |
| ρ = 0.5 | .04/.12, .05/.14, .05/.16 | .04/.12, .05/.14, .05/.16 | .09/.02, .11/.04, .18/.04 |
| ρ = 0.7 | .13/.06, .23/.07, .27/.08 | .13/.06, .23/.07, .27/.08 | .19/.01, .37/.00, .43/.00 |

The ρ = 0.35 hit-conditional rates rest on only 33 replicates and are not tabulated. Both
contrasts are in `tables_v3.md`. The screened control pool changes S by ≤ 0.02 except in the
rival scenario, where S rises from .16 to .23 at E = 10 once rivals are screened out of the
controls.

### R4. Mn averaging (Carbon)
**Adopted.** Mn is the mean over the 20 matched sets of the per-collection log scores,
(1/20) Σ_s log p_{N_s}(y), and not the log of the mean probability. Averaging probabilities first
would give Mn an ensemble bonus that the single-model Mm does not get. Every version of the code
(v1–v3) already computed it this way, so no number changes. P3 §3 should state it explicitly.

### Revised recommendations (supersede §3 where they differ)
- **(a) Gate.** Keep the 10-event floor as **feasibility only**. Power is set by discovery and
  precision, not by the floor.
  - Even with every improvement, an E = 20 test year supports a ρ = 0.7 module 23% of the time
    unconditionally, and a ρ = 0.5 module 5%.
  - Pre-register these discovery-inclusive operating characteristics, recomputed at the
    reconciled event counts. Also state that the likely outcomes are I, or F for pseudo-cliques.
- **(b) Verdict rule.** Rule 1, one common margin, **δ = 0.02** nats per collection, mutually
  exclusive S/F/I, plus N (no clique).
  - Rule 2 claims (C+ with LB > 0; C− with UB < 0.05) are reported as secondary descriptions.
  - The estimand is the realised 2025 transfer gain. The interval is an event-level studentised
    (or percentile) bootstrap.
  - Cross-year language ("recurs") is reserved for designs with ≥ 2 held-out years.
- **(c) Specificity.** Unchanged: Mm−Mn is not the binding contrast. Discovery now produces
  env pseudo-cliques that Mm−Mn alone does not reject (their Mm−Mn truth is about 0).
  - The rule only stops them because S needs a positive lower bound on both contrasts.
  - P3 should report how many members of the selected clique load strongly on the background
    gradient.
- **(d) Aggregation.**
  - Rotate masks within each event before resampling.
  - Freeze shrinkage at λ = 5-equivalent. Carbon endorses both.
  - Use the Mn mean-of-log-scores rule from R4.
  - Shrinkage matters doubly with discovery: it moves pseudo-cliques from −0.06 (λ = 1) to
    about 0, and the ρ = 0.5 discovered gain from −0.038 to +0.012.

### Limits of v3
- **Residual correlation.** The pairwise bivariate-probit likelihood is not a JSDM posterior:
  margins are fixed at ridge-logistic fits, there is no joint shrinkage, and within-event
  clustering in training is ignored.
- **Ranking.** It refits ridge per candidate on 10 blocked folds. P3 instead specifies
  leave-one-event-out ranking with the JSDM conditional and no per-candidate refit.
- **Model assumptions.** The pseudo-clique rate depends on how well Me captures the gradient
  (c = 0.85) and on the 8 strongly loaded outsiders. Both are assumptions, but the mechanism is
  the one P3 should check in its 2024 data.
- **Fixed design.** One within-event correlation (r = 0.5), one training size (30 events) and
  K = 4.
