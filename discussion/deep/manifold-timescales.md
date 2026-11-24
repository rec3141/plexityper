# Communities as multiscale learned representations: formalisation, MEP anchoring, tests

Author: Claude (Opus 5.5), deep-work sub-agent for `claude-plexi`, 2026-11-18. **Status: draft for
cross-family review** (one vote under AGENTS.md). Nothing was posted to the relay. No file under
`data/` was opened. No numbers below come from project data; every number is a proof, an analytic
table, a simulation on synthetic data, or a quotation from Vallino & Huber 2018.

**What I read.**
- AGENTS.md (full).
- `discussion/plexi-concept-review.md`: "Summary for grid", "Standing perspective from grid", the
  one-sentence definition, "Four orthogonal axes", "Mechanism layer", "Known traps" (l. 1-140,
  368-495 as read today, 615 lines total). Findings-log entries for CMO (l. 180-200).
- `plexichat.txt` (skimmed in full, 157 lines).
- `notes/paired-events-and-fractions.md` (l. 1-30) and grep hits in `protocols/P3-*.md` for cohort
  descriptions. I did not open any metadata tables, so all cadence and chemistry statements in §8
  are **unverified** and phrased as checks to run.
- **Vallino J.J. & Huber J.A. 2018**, "Using Maximum Entropy Production to Describe Microbial
  Biogeochemistry Over Time and Space in a Meromictic Pond", *Front. Environ. Sci.* 6:100,
  doi 10.3389/fenvs.2018.00100. **Version read:** the publisher's full-text HTML page, downloaded
  2026-11-18 into its own empty directory (`page.html`, sha256 `e7f861579216`), converted to text
  with a stdlib script. I read the abstract, introduction, methods, results, discussion and
  conclusions. The displayed equations (12 and the rate law) did not render in the text
  extraction, so their exact forms below are reconstructed from the surrounding prose and marked as
  such. **I did not read the Supplementary Material** (where the optimal-control method and full
  stoichiometries are). The guess that this is the Siders Pond MEP paper is **confirmed**.
- Input from Astra (relay 248970300628582400, 248971565584531456) and Carbon (248973006869020672),
  relayed by the coordinator, is built in where marked.

**Code** (numpy 2.3.5 + scipy 1.17.0, `python3 -I`), in `scratch/deep/`:

| file | sha256 (12) | output |
|---|---|---|
| `manifold_timescales_sim.py` | see §10 | `manifold_timescales_sim.out` |

Run: `cd scratch/deep && python3 -I manifold_timescales_sim.py 100 > manifold_timescales_sim.out`
(the recorded output uses 100 reps per cell; optional env `SIM_WORKERS`, default 8).

**Labels.** **[Proof]** argument given in full here. **[Std]** standard result (reference in §11).
**[Sim]** numerical result on the synthetic model only. **[Read]** stated in Vallino & Huber 2018
(main text). **[Conj]** conjecture or judgement. **[Spec]** speculation: a reading worth keeping,
with no derivation behind it.

**Three things kept distinct throughout (Astra).** The **layer** l (population, genome, pathway,
gene, expression, nucleotide), the **sampling interval** Δ of a dataset, and the **response time**
τ_l of a layer *as estimated from data at that Δ*. τ_l is a property of the system; Δ is a property
of the design; an estimate τ̂_l exists only when Δ ≲ τ_l ≲ (series length).

---

## 0. Bottom line

1. "Free energy" has two unrelated meanings here. **(T)** thermodynamic: Gibbs energy dissipated by
   metabolism (Vallino's MEP). **(V)** variational: a bound on surprise (Friston). They share a name
   and a minus sign, nothing else. The perspective becomes testable only in reading (T), plus a
   third reading **(L)**: the replicator update has the *form* of Bayes' rule (literal inference
   only under a stated observation model; §2.4).
   Reading (V) adds almost no prediction (§2.3).
2. "State = sum of layers" is **tautologically true** as a band-pass (telescoping) decomposition of
   any time series [Proof, §2.1]. It becomes a claim only with (a) spectral gaps between bands,
   (b) a mechanism tied to each band, and (c) slaving of fast layers to slow ones.
3. Reading (T) with Vallino's finite horizon gives the cleanest bridge: **a fast layer can only
   execute a long-horizon strategy if the policy is stored in a slower layer.** Genome =
   weights/policy; expression = activations/controls; horizon of anticipation bounded by the
   timescale of the storing layer [Conj, motivated in §2.2].
4. **Song vs singer as an instantaneous symmetry of the MEP objective** [Proof, §4.3; scope
   narrowed in §12]: with catalysis linear in enzyme and well-mixed chemistry, dissipation and
   dc/dt at fixed (c, E) depend only on enzyme per reaction, not on which genome carries it; whole
   trajectories agree only along a shared, realisable E-path. The objective is indifferent to
   packaging, but recipients in a given packaging can still depend on particular genomes, so MEP
   alone is *silent* on plexitrophy. It can bear on it only through realisability and per-genome
   symmetry breakers (incompatibility, crowding, toxicity, transport, spatial potentials). Astra and Carbon independently reached the same constraint
   (ΔG of a split net reaction is unchanged; any gain is in realised flux J_r).
5. Crisp definitions (d ≥ 2, plexitype membership, topology class) should be re-expressed as graded,
   layer- and band-resolved quantities, with the crisp versions as reported threshold summaries
   (§3). Pre-registered tests should target the graded estimand.
6. Simulation [Sim, §7; 100 reps per cell]: at 30-100 time points a two-timescale truth is
   **barely distinguishable** from a calibrated one-timescale + noise null (power 0.02-0.15 at
   τ ratio 20 and noise 0.3; at most 0.24 at n = 100 with noise 0.05). A single slow forcing cycle
   over the series (N2b, a slow drift) is rejected as often as the best two-layer case (0.05-0.24),
   and a continuum of timescales (N3) about as often as a true two-layer model (0.11 at n = 100).
   A P = 12 seasonal sinusoid (N2) was rejected in 0% of replicates, because a sum of OU
   covariances is monotone and cannot fit oscillatory autocovariance; seasonality does not mimic
   layers in this model class, slow drift does. So **"not red
   noise" is not evidence for layers**, and single-series spectral separation is not testable at
   our lengths. Cross-layer ordering (RNA faster than DNA, within MAG) is better: per-MAG power
   0.51 at n = 40, 0.73 at n = 100; pooled sign-test power 0.92 at n = 40 with 30 MAGs, *assuming
   independent MAG signs* (not true when MAGs share sampling events; §7 note). At n = 20 the noise
   difference alone biases the null (P[D > 0] = 0.67, median D = +0.38), so a naive sign test
   against 0.5 with 30 MAGs would have a ≈ 0.60 false-positive rate (binomial: reject at ≥ 20/30;
   P = 0.600 under p = 0.67). Separate per-layer fits are therefore unsafe at CMO-like lengths; use
   a joint bivariate model (§9).
7. Recommendation (§9): a framing paragraph in Paper 1 plus a **protocol (P6)**, not a Paper 2 yet.
   Strongest first test: **P6** (`protocols/P6-thermodynamic-allocation.md`, authoritative): do
   parameter-free thermodynamic factors (per-electron ΔG; gate F_T, then power weight) predict
   expressed catabolic allocation in raw counts better than capacity × kinetics? Conditional on
   paired RNA + chemistry existing. A win shows predictive value for expression; it does not
   measure flux J or show that dissipation is maximised (§9). Fallback: the pooled timescale-ordering test on CMO,
   a precondition, not a confirmation.

---

## 1. Two meanings of "free energy", plus a third reading

**(T) Thermodynamic.** Gibbs free energy of reaction Δ_rG_i(c, T, pH) for each catalysed reaction
i, and entropy production σ = Σ_i r_i(−Δ_rG_i)/T (plus radiative terms). "Free-energy
minimisation" in this reading means *driving the environment's free-energy potentials down*,
i.e. **maximising dissipation**: the MEP hypothesis (Dewar 2003; Martyushev & Seleznev 2006;
applied to microbes by Vallino 2010; Vallino & Huber 2018). Units: J K⁻¹ s⁻¹. Measurable in
principle from rates and activities.

**(V) Variational.** F[q] = E_q[−log p(s, η)] − H[q] ≥ −log p(s): an upper bound on the surprise of
sensory states s under a generative model p, with internal states parameterising q over external
states η (Friston 2010, 2013). Dimensionless (nats). Not measurable from chemistry; requires
specifying p, which for a microbial community nobody can write down.

**(L) Learning (selection as inference).** The discrete-time replicator update has the Bayesian
update form; literal inference requires an observation model (§2.4). This reading needs neither (T) nor (V).

Predictions below are tagged with the reading they follow from. Where a prediction follows from
(T) only, (V) is silent on it; I found none that follows from (V) only and is not generic to any
stable dissipative system (§2.3, §5).

---

## 2. Formalisations

### 2.1 F1: multiscale slow-fast dynamics (the "sum of layers")

Let the community state be x = (x_1, …, x_L), one block per layer, ordered from slow to fast,
with
  dx_l/dt = ε_l f_l(x_1, …, x_L, u(t)) + √ε_l Σ_l ξ_l(t),   ε_1 ≪ ε_2 ≪ … ≪ ε_L,
u(t) environmental forcing, ξ white noise. Gradient special case: f_l = −∇_{x_l} Φ(x; u), a
shared potential Φ ("free energy" in either reading, or fitness).

**Tikhonov / Fenichel [Std].** If, for frozen slow variables, the fast subsystem has a normally
hyperbolic attracting manifold x_fast = h(x_slow), then for small ε-ratios trajectories stay
O(ε-ratio) close to a perturbed invariant manifold, and the slow flow on it is the reduced system.
Biological reading: expression relaxes onto a quasi-steady state set by current abundances and
chemistry; abundances relax onto a state set by genome content; genomes move under selection.
This is Haken's *slaving principle* in another language.

**Linear-Gaussian derivation [Proof].** Linearise around a stable point and suppose each layer
contributes one independent OU component to an observable y = Σ_l a_l z_l + noise, with
dz_l = −z_l/τ_l dt + dW_l. Then the spectral density is
  S_y(ω) = Σ_l a_l² · 2 v_l τ_l / (1 + ω² τ_l²) + σ_n²,
a sum of Lorentzians with corner frequencies 1/(2π τ_l). Sampled at interval Δ, each component is
an AR(1) with φ_l = exp(−Δ/τ_l). The sum of L independent AR(1) processes is ARMA(L, L−1); adding
white measurement noise makes y ARMA(L, L) (Carbon 249247746036776960). Variance is exactly additive
across layers *because* the components were assumed independent.

**The tautology [Proof].** For any square-integrable series y and any increasing sequence of
smoothing scales s_1 < … < s_L, with ȳ_s the low-pass of y at scale s,
  y = ȳ_{s_L} + Σ_{k=1}^{L} (ȳ_{s_{k−1}} − ȳ_{s_k}),  ȳ_{s_0} ≡ y.
So *some* decomposition "state = sum of band components" always exists. The perspective claims
more, and each extra piece is testable:
- **(a) spectral gaps**: power is concentrated near a few corner frequencies, with troughs between,
  rather than spread over a continuum (1/f-like spectra are generic for ecological series:
  Inchausti & Halley 2001; Vasseur & Yodzis 2004);
- **(b) layer identification**: the band where layer l's *own observable* (e.g. SNV frequency for
  the nucleotide layer) has its power is the band it contributes to the community state;
- **(c) slaving**: cross-spectra are directional. The fast layer's mean tracks h(slow), but fast
  deviations do not feed the slow layer except through averages.

**Additivity is not orthogonality.** Under slaving, the expression layer *contains* the slow signal
(x_fast = h(x_slow) + δ). The physically defined layers are therefore correlated; only the
*residuals* δ_l are approximately independent. "Sum of layers" must mean a sum of innovations
(a multiresolution of increments), not of physical observables. This changes how additivity is
tested (§6.3).

### 2.2 F2: thermodynamic reading anchored in Vallino's MEP with finite horizons

**What the paper does [Read].** A 1-D advection-dispersion-reaction model of Siders Pond (15 m,
meromictic, P-limited) with 11 chemical state variables and 8 functional "catalysts" (phytoplankton,
green sulfur bacteria, photoheterotrophs, aerobic heterotrophs, SRB, sulfide oxidisers, aerobic and
anaerobic grazers), 28 reactions. Each catalyst j has a growth-efficiency control ε_j ∈ [0,1]
(ε → 1: all substrate to biomass; ε → 0: all to catabolism) and allocation controls Ω_{i,j}
(fraction of catalyst j's machinery assigned to its sub-reaction i, Σ_i Ω_{i,j} = 1): 28 controls
in all. Rates take the form (reconstructed from prose; equation not rendered)
  r_{i,j} = ν* · Ω_{i,j} · [catalyst_j] · F_K(c; ε_j) · F_T(Δ_rG_{i,j}),
with an "adaptive Monod" kinetic factor F_K whose affinity/μmax trade-off is set by ε_j, and a
thermodynamic factor F_T after LaRowe et al. 2012. Entropy production per volume is
σ = Σ r_{i,j}(−Δ_rG_{i,j})/T plus light absorption by water and particles. The controls are chosen
to maximise the **time-averaged** entropy production
  σ̄(t) = (1/Δt) ∫_t^{t+Δt} σ(s) ds
(their Eq. 12, reconstructed), sequentially, locally at 10 depths. **Δt is the single biological
parameter**, replacing ~89 kinetic parameters of a conventional model.

Results quoted [Read]:
- Two horizons compared: Δt = 0.25 d ("SIO") and 3 d ("LIO").
- Observations: one station, 24 h, 7 casts × 8 depths (25-26 June 2015).
- RMSE table (their Table 3, values as extracted; bolding lost): LIO lower for PAR, Chl a, DIC,
  H₂S, sulfate, DOC; SIO lower for DO, phosphate, POC. The authors summarise "6 out of 9".
- Integrated entropy over 6 d: 36.5 (LIO) vs 21.3 (SIO) GJ K⁻¹; a sterile pond would give 12.4.
- SIO controls switch bang-bang (ε_Phy between >0.98 and <0.02, dropping near noon); LIO controls
  vary smoothly (ε_Phy ≈ 0.3-0.4), with Ω switching between carbon fixation by day and biosynthesis
  by night. SIO resembles r-selection; LIO resembles K-selection.
- Local, not global, maximisation was required for light, because water dissipates light
  abiotically fast (global optimum degenerate).
- Sulfide oxidisers were "largely unimportant" in both 6-day simulations.
- The authors explicitly do not know the right Δt; they suggest it relates to "the time scales over
  which biological predictions, acquired by evolution, can be reliably made".

**My reading of the strength of that evidence [Conj].** It is a proof of concept, not a test.
Two horizons, one 24 h observation window, RMSE without uncertainty, a 6/9 split
(one-sided sign test P = 0.25 if the nine variables were independent, which they are not),
functional groups and their reaction sets fixed by the modellers, and physical parameters
"crudely adjusted". The finding that longer horizons dissipate more is close to definitional:
optimising over a longer window with sequential receding horizons can only find solutions at least
as good on that window. What is *not* definitional, and is the paper's real empirical content, is
that the long-horizon solution matched the pond better.

**Critiques of MEP to carry [Std, citations to check].** Dewar's (2003) information-theoretic
derivation was disputed (Grinstein & Linsker 2007). MEP is not a theorem of non-equilibrium
thermodynamics; near equilibrium the linear regime gives *minimum* entropy production (Prigogine)
for steady states with constraints; MEP is best treated as a selection principle among many
feasible steady states. Vallino's own extension to finite horizons is a hypothesis.

**Mapping finite-horizon MEP onto layers.** Write each layer's controllable variables as u_l (e.g.
Ω and transcript levels for expression; ε and catalyst abundances for population; the reaction
sets n_{r,j} available to each catalyst for genome; allele frequencies for nucleotide). In
Vallino's model the reaction sets are fixed and only (ε, Ω) move. Generalise:

  max over u_L(·) of ∫_t^{t+H} σ(x, u) ds,  with  u_L = π(x, u(t); θ),

where π is a *policy* (a regulatory mapping from sensed state to control) parameterised by θ, and
θ lives in a slower layer (genome: promoter logic, clocks, storage pathways; population: which
policies are present at what frequency).

**Derivation sketch, horizon ≤ storage timescale [Conj].** A control can depend on information about
the future (e.g. "light returns in 12 h", "winter is coming") only through some variable that holds
that information until it is used. Fast layers forget on τ_fast. So any strategy whose benefit
accrues after a delay D > τ_fast must be encoded in a layer with τ ≥ D: the genome (circadian
clock, storage polymer synthesis), or slow population structure (dormant seed banks). This turns
Vallino's single Δt into a **ladder of horizons H_l ≲ τ_{l'} where l' is the layer that stores the
policy**. It also makes the mechinterp mapping concrete: **θ = weights (slow, updated by selection
on integrated dissipation or fitness over long windows), u = activations (fast controls computed
from current state)**. Genome-encoded anticipation is literally a learned policy; whether its
learning objective is dissipation (T) or reproductive fitness is a separate question (§2.4, §5).

**A spectral corollary [Conj].** Vallino's SIO/LIO contrast (bang-bang vs smooth control) suggests
that a long horizon acts as a **low-pass filter on allocation**: controls change mostly at
frequencies below ~1/H, except where the forcing is predictable (diel light), where they lock to the
forcing with a phase lead (anticipation). Testable as a transfer function from forcing to
expression (§6.6), but only with sub-daily sampling.

**Link to growth [Std].** England (2013) bounds the dissipation required for self-replication
from below; growth and dissipation are coupled but not identical objectives.

### 2.3 F3: variational free energy / active inference

**Statement [Std].** Any system at a non-equilibrium steady state with a Markov blanket (internal
states μ conditionally independent of external η given blanket states b) can be described as if
internal states parameterise beliefs q_μ(η) and the flow performs gradient descent on F
(Friston 2013; Ramstead et al. 2018).

**Critiques [Std, citations to check].** The "as if" holds under strong conditions (specific
Helmholtz decompositions, sparse coupling) that generic systems do not satisfy (Biehl, Pollock &
Kanai 2021; Aguilera, Millidge, Tschantz & Buckley 2022). Markov blankets as statistical objects
are not the same as physical boundaries (Bruineberg et al. 2022). When the conditions hold, the
principle applies to essentially every persistent system, so it does not discriminate.

**Literal / metaphor / fails for microbial communities.**
- *Literal*: nothing, until someone specifies the blanket. A bulk water sample has no blanket; a
  particle, biofilm or host might (this matches the project's rule that a **plexon** needs an
  independent physical boundary).
- *Metaphor*: "the community models its environment" is a fair gloss for genome-encoded
  anticipation (§2.2), which (T) and (L) already explain.
- *Fails*: no generative model p can be written; F is not measurable; the principle predicts no
  rate, no allocation, no spectrum that a generic stable dissipative system would not show.
**Recommendation**: mention (V) once, as a related framing, and build nothing on it.

### 2.4 F4: selection as Bayesian learning (reading L)

**Replicator update has Bayesian form [Proof of the form only; narrowed after Astra
249247324022685696].** Discrete-time replicator dynamics with frequencies p_i(t) and
environment-dependent fitness w_i(e_t):
  p_i(t+1) = p_i(t) w_i(e_t) / Σ_j p_j(t) w_j(e_t).
This has the algebraic form of Bayes' rule with p_i(t) as prior and w_i(e_t) in the place of a
likelihood (Shalizi 2009; Harper 2009; Czégel et al. 2022). It is **literal inference about the
environment only under a stated observation model**: one must posit p(e_t | i) ∝ w_i(e_t) with a
proportionality constant independent of i, i.e. that fitness *is* a likelihood of environmental
data under hypothesis i. Nothing in ecology guarantees that, so the identity is one of form.
Iterating the form,
  log p_i(T) = log p_i(0) + Σ_{t<T} log w_i(e_t) − log Z_T,
so composition accumulates the environment's history *as scored by the types' fitness functions*.
Without mutation or immigration it never forgets and converges on the type with the highest
time-averaged log fitness (Kelly / geometric-mean fitness).

**Immigration [narrowed].** Add immigration at rate m from a fixed pool q, applied after
selection: p(t+1) = (1−m)·Bayes(p(t), e_t) + m q. Under the observation model above plus a
hidden "true type" that is resampled from q with probability m each step, this recursion is the
**next-step predictive distribution** of that hidden Markov model (update, then predict).
What is proved here:
- **Neutral fitness** (w_i ≡ w): p(t) − q = (1−m)^t (p(0) − q) exactly, so initial differences
  decay geometrically with time constant ≈ 1/m generations [Proof: the selection step is the
  identity; iterate the affine map].
- **Informative selection**: no universal 1/m memory is claimed. Selection can amplify or erase
  differences faster or slower than mixing; filter-stability results (e.g. Del Moral & Guionnet
  2001) give forgetting of the initial condition under mixing conditions, with rates that depend on
  both m and the likelihood ratios, not on m alone.
So the defensible statement is: immigration (and mutation/HGT at the genome layer) *bounds memory
from above in the neutral limit* at ~1/m; under selection the memory window is an **estimand**,
expected [Conj] to shorten as m rises, other things equal. Watson & Szathmáry 2016 and Vanchurin
et al. 2022 propose multilevel versions with separate timescales per level.

Where (L) is literal: frequency dynamics under selection. Where it is metaphor: "the community
infers the environment" (nothing reads the posterior out). Where it fails: density dependence and
interactions (fitness depends on p), which break the clean likelihood reading; strong ecological
interactions are exactly the project's subject.

**Prediction from (L) [Conj, testable]:** composition carries lagged environmental information
beyond current environment; the lag window is estimated, not assumed to equal 1/m. Sites with faster flushing or dispersal should
show shorter memory. Carbon's caveat applies: conditioning on current E alone conditions on a
descendant of lagged drivers; tests must use distributed-lag or seasonal-phase terms and a
phenology control guild (§6).

### 2.5 F5: mechanistic-interpretability reading

| mechinterp object | microbial analogue | status |
|---|---|---|
| weights (slow, set by training) | genome content, regulatory logic (θ in §2.2); at a coarser layer, the population's type frequencies | analogy close; "training" = selection + HGT + immigration, objective unknown |
| activations (fast, input-dependent) | transcript levels, enzyme allocation Ω, metabolite pools | close |
| input | environmental forcing u(t), chemistry | close, but inputs are not controlled |
| feature (direction in activation space) | a co-varying gene or taxon program, e.g. a nitrogen-limitation response; a plexitype (§4.1) | operational: a reproducible loading vector |
| superposition (more features than dimensions) | a taxon or gene participating in several programs (pleiotropy, polysemantic taxa); more functions than taxa | plausible; untested |
| privileged basis | is the taxon basis privileged? Yes where per-genome nonlinearities act (§4.3) | this is the song-vs-singer question |
| circuit | a directed exchange subgraph implementing a community function: the plexiome | close in form; edges are hard to identify (architecture-identifiability.md) |
| probing | regress a process rate or chemistry on community state; does a linear readout exist? | feasible now |
| activation patching / ablation | mesocosm perturbation, dilution-to-extinction, knockout cultures | feasible but costly; CMO is one pool per arm |
| sparse autoencoder | sparse overcomplete dictionary over samples × (KO or taxon) matrix | **sample size fails**: SAEs use ~10⁶-10⁹ activation vectors; we have 10²; NMF / topic models are the realistic analogue |

**Where it fails.** (i) No known loss function, so "learned representation of a manifold" lacks the
object that makes the ML statement meaningful; reading (T) or (L) has to supply it. (ii) No gradient
access and no clean input-output pairs. (iii) Compositional, sparse, undersampled data. (iv) In ML
the representation is *for* a downstream task; ecosystems have no downstream reader. What survives
is the *method*: find reproducible low-dimensional features, test whether a feature is causally
used (intervene), and ask whether the basis (taxa) is privileged.

---

## 3. Graded, probabilistic, layer- and band-resolved constructs (grid: "rigid rules won't fit")

Principle: each crisp object becomes a **graded estimand with a posterior or interval**, indexed by
layer l and frequency band b where that makes sense; the crisp object is reported as a threshold
summary of it, never as the primary estimand.

### 3.1 Plexitrophy
Current crisp rule (plexi-concept-review.md, synthesis; supersedes the |SCC| ≥ 3 and
minimal-support wording): for focal X under medium E and a stated endpoint, X is plexitrophic iff
**in every feasible wiring W, d_W(X) ≥ 2 OR κ_W(X) = 1**, where κ_W(X) = 1 iff X lies on a simple
supply cycle of length ≥ 3 in W. Report (d, m, κ); empty or incomplete enumeration = unresolved.

Graded forms (revised after Astra 249247324022685696; see §12):
- **Posterior plexitrophy.** Let 𝓦(E) be the set of feasible wirings in an instance draw.
  Define the event P = {𝓦(E) ≠ ∅ AND enumeration complete AND ∀ W ∈ 𝓦(E): d_W(X) ≥ 2 ∨
  κ_W(X) = 1}. The non-empty condition is explicit, so a draw with no feasible wiring can never
  count as a positive by vacuous truth. Uncertainty in genome completeness, unbinned donors,
  ambient supply and incompatibility rules enters as priors over the instance (the instance used
  by `plexisolve.py`, sampled rather than fixed); each draw is decided by the existing exact
  decider and falls in exactly one of: **positive** (P), **negative** (feasible, complete, some
  wiring with d_W = 1 and κ_W = 0), **infeasible** (𝓦 = ∅: X cannot reach the endpoint under E),
  **unresolved** (enumeration incomplete). Report
  - π(X) = Pr[P | data] (denominator: all draws), with its Monte Carlo error;
  - π_feas(X) = Pr[P | 𝓦 ≠ ∅, complete, data] (denominator: feasible, resolved draws, stated as a
    count);
  - the infeasible and unresolved masses.
  π is the posterior probability of the crisp event, so **report it directly**. The crisp call is
  the threshold summary π ≥ π* (e.g. 0.95).
- **Partner-count posterior.** The posterior distribution of min_W d_W(X) (and of κ) over the
  same draws, reported as a distribution, not collapsed.
- **Supply-share descriptor (separate, not a plexitrophy criterion).** With supply shares s_k
  over indispensable suppliers, the entropy-effective count exp(H(s)) describes how evenly supply
  is spread. It is *not* d: two indispensable suppliers at 99:1 give d = 2 but
  exp(H) = exp(0.056) ≈ 1.06. Indispensability (d) and evenness (exp H) are different axes; the
  crisp rule uses only the first.
- **Dependence strength per layer.** D_l(X) = expected drop in X's growth or activity when the
  partner contribution held in layer l is removed:
  - genome layer: auxotrophy (gene absent; slow, "weights");
  - expression layer: facultative cross-feeding that is regulated (fast, "activations");
  - population layer: dependence on partner abundance exceeding a threshold.
  In data, D_l is estimated by the band-limited coherence between partner supply (DNA capacity ×
  expression) and recipient activity in band b, conditional on aggregate guild supply (trap 4).
- Crisp summary: π(X) with the threshold call π ≥ π*, the posterior of d and κ, and (separately)
  exp(H) of supply shares with an interval.

### 3.2 Plexitype membership
- **Mixed membership**: sample i has weights θ_ik over K latent assemblage programs (NMF, LDA,
  or a Dirichlet-multinomial mixture); each taxon j has loading β_kj. Membership of j in plexitype k
  is β_kj with a posterior, not a set.
- **Band-resolved membership**: compute programs separately on band-passed series (slow seasonal
  band vs. fast residual), since taxa can co-vary in one band and not another. A "slow plexitype"
  and a "fast plexitype" are then different objects; P3's masked-member prediction can be run per
  band.
- Crisp summary: members = {j : Pr[β_kj > β*] > 0.9}; report how membership changes with β*.
- P3's estimand (gain in log-score, nats) is already graded; keep it, and report it per band if the
  series allows.

### 3.3 Topology classes
Replace star / chain / mesh labels with continuous indices computed per posterior draw of the
hypergraph:
- hub concentration (degree Gini or Freeman centralisation);
- path-likeness (normalised diameter, or fraction of flow through cut vertices);
- cyclicity (cycle rank / number of edges, or fraction of nodes on a simple supply cycle of
  length ≥ 3, matching κ);
- modularity Q and nestedness (NODF).
A class label is the region of this index space the posterior mass sits in, with "unresolved" when
no region holds ≥ 0.9 of the mass. This is consistent with architecture-identifiability.md's
three-valued edges.

### 3.4 What a pre-registered test should target
- A **continuous estimand** with an interval (posterior or bootstrap), e.g. mean excess of π(X) or
  of the posterior of d over its null across recipients; log-score gain; band-limited coherence; a spectral
  exponent; an allocation–power association.
- **Three decision bands** fixed in advance: support (interval above +δ), equivalence/against
  (interval inside [−δ, +δ]), inconclusive (otherwise). δ is a smallest effect of interest.
- Crisp labels are reported as secondary summaries and their threshold sensitivity is shown.

---

## 4. Mapping the project's objects onto the framework

### 4.1 Plexitype as a feature
A plexitype as defined (reproducible subset of taxa that beats environment and matched nonmembers
out of sample) is a **feature in composition space**: a direction or low-dimensional subspace whose
coordinate is predictable from other coordinates beyond what the inputs (environment) explain. It is
a feature of the *population layer*. The mechinterp distinction applies: a feature can be
*represented* (reproducible direction) without being *used* (causally relevant to function). P3
tests representation. Use requires intervention or the plexibiont evidence standard.

### 4.2 Plexiome as a circuit
The plexiome (AND/OR dependency hypergraph) has the form of a circuit: nodes = roles, edges =
transfers, AND = multiplicative gating, OR = redundancy. Two cautions: edge direction cannot come
from covariance (the plexiome definition already says so), and circuits in ML are found by
intervention, which here means cultures or mesocosm perturbations.

### 4.3 Song vs singer as a gauge symmetry; plexitrophy as symmetry breaking

**Proposition [Proof].** Suppose (i) chemistry is well mixed on the scale considered, so all
organisms see the same activities c; (ii) each reaction's rate is linear in the amount of enzyme
assigned to it, r_i = g_i(c) E_i with E_i = Σ_j Ω_ij B_j (as in Vallino's rate form when catalysts
share composition and F_K, F_T depend only on c); (iii) anabolic stoichiometry is the same for all
catalysts (Vallino assigns all catalysts yeast composition), and per-genome controls such as
Vallino's ε_j can take the same value in any genome that carries reaction i (in his model ε_j is
per catalyst, so this holds when two catalysts carrying the same reaction are free to choose equal
ε). Then **at fixed (c, E)**, σ = Σ_i g_i(c) E_i (−Δ_rG_i(c))/T and dc/dt depend on the packaging
(Ω, B) only through the vector E.
*Proof*: every term in the chemical dynamics and in σ is a function of (c, E); B and Ω enter only
through E_i = Σ_j Ω_ij B_j. ∎

**Scope (revised after Astra 249247324022685696).** This is an *instantaneous* invariance.
- **Trajectory corollary [Proof]**: two packagings with the same c(0) and the **same path E(t)**
  over the horizon give the same c(t) and σ(t) (same ODE for c, same forcing E(t)). This holds,
  for example, when E is held fixed.
- **Not implied**: equal initial E and equal *summed* biomass growth do not give equal paths.
  dE_i/dt = Σ_j (Ω_ij dB_j/dt + B_j dΩ_ij/dt), and dB_j/dt depends on genome j's own reactions,
  efficiency ε_j and losses; different packagings therefore generate different dE/dt from the same
  (c, E). Whether a packaging can *realise* a given E-path is a growth-dynamics question the
  objective does not answer. (The earlier draft claimed whole-trajectory invariance; that was
  wrong as stated.)
- So the correct reading is: **the MEP objective is indifferent to which genomes carry a flux,
  along any E-path that the packagings being compared can all realise.** It ranks songs, not
  singers, and only among realisable songs.

**What follows for plexitrophy.** The objective's indifference does *not* mean absence of
dependence. In any particular packaging, a recipient X whose required input is produced only by
genomes Y and Z depends on Y and Z, whatever the objective thinks of alternative packagings.
Plexitrophy is a property of the packaging (the singers and their wiring), not of the objective.
MEP can therefore bear on plexitrophy only through (a) which packagings can realise high-dissipation
E-paths (growth dynamics, transport, regulation), and (b) the symmetry breakers listed below. MEP
by itself is **silent** on plexitrophy; it neither predicts nor forbids it.

Remark [Proof, one line]: if instead one ε is shared by all reactions of a catalyst (as in
Vallino's parameterisation), then splitting two reactions across two catalysts adds a control
variable, and the optimum over the larger control set is ≥ the optimum over the smaller one. So
in that parameterisation division of labour is weakly favoured *by construction*; any such gain
must be checked against real per-genome regulation before being read as biology.

Consequences:
- Under (i)-(iii), **MEP scores the song (the E-path and fluxes), not the singers (the
  factorisation E = ΩB)**, among realisable paths. Functional stability with taxonomic turnover
  (Louca et al. 2016) is then unsurprising, but not derived: it also requires that several
  packagings can realise the favoured path.
- In mechinterp terms: the taxon basis is **not privileged** under (i)-(iii). It becomes privileged
  where an elementwise (per-genome) nonlinearity acts, exactly as neurons become a privileged basis
  only because activation functions act elementwise.
- **Astra's constraint, agreed.** Splitting a net reaction A → C into A → B (partner 1) and B → C
  (partner 2) leaves Σ ΔG unchanged when B cancels and the boundary conditions match. So a
  discriminating MEP prediction must concern **realised flux** J_r and its spatial or temporal
  distribution: compare Σ_r J_r(−Δ_rG_r)/T integrated over the same volume and horizon, accounting
  for transport and accumulation of intermediates.
- **Where division of labour (and with it, possibly d ≥ 2) can raise dissipation**, each a broken
  assumption:
  - (i) broken: potentials are spatially separated (redox interfaces, particles vs free water,
    steep O₂ micro-gradients on particles: Carbon; relevant to P4). An organism at an interface
    couples a donor made by one partner below with an acceptor made by another above:
    **catabolic fan-in d_catab = 2**. This is MEP's most natural route to plexitrophy, with a
    standing null: either potential may be abiotic (diffusion from sediment or atmosphere), which
    removes it from the partner count.
  - (ii) broken: per-genome costs are nonlinear (enzyme crowding, membrane real estate, pathway
    length; Pfeiffer, Schuster & Bonhoeffer 2001; Costa, Pérez & Kreft 2006 for nitrification;
    but comammox exists: Daims et al. 2015). Then splitting raises J_r at equal total protein.
  - (iii) broken: intermediates are toxic or reactive, or steps are O₂-incompatible (the team's
    "stacking blocked" mechanism, notes/p2-incompatible-pairs.md).
  - Temporal: storage and dormancy shift dissipation in time (Vallino's horizon).
- **Anabolic plexitrophy** has a weaker MEP route [Conj]: anabolism dissipates little; MEP favours
  minimising redundant anabolic machinery per unit catabolic flux (a Black-Queen-like pressure),
  but whether that yields d_anab ≥ 2 rather than one stacked helper depends on the same symmetry
  breakers. So MEP adds an *objective* to Carbon's "blocked shooting the moon" mechanism; it does
  not replace the constraint that produces fan-in.
- **Plexitrophy as a statement about layers.** A genome-encoded auxotrophy is a dependence stored in
  slow weights; a regulated cross-feeding is stored in fast activations. The definition's
  "obligate" is a statement about the genome layer under a stated medium; expression-level
  dependence is facultative and belongs in D_expression (§3.1). Prediction [Conj]: obligate
  (genome-layer) dependences should be on supplies whose availability is stable at frequencies
  below the genome's turnover (constant supply lets loss be tolerated); facultative (expression-layer)
  dependences should be on supplies that fluctuate in the fast band.

---

## 5. What MEP (T) predicts that the variational reading (V) does not

| # | Prediction | (T) MEP | (V) Friston | (L) replicator-Bayes |
|---|---|---|---|---|
| 1 | Process rates and allocation are ordered by available power J·(−ΔG), computable from chemistry | yes, quantitative | silent | silent |
| 2 | Unexploited free-energy potentials are rare where biology can catalyse them (standing ΔG < 0 couples are kinetically or spatially protected) | yes | silent | silent |
| 3 | Longer effective horizons → smoother allocation, more storage, more standing biomass, higher integrated dissipation | yes (Vallino) | "anticipation" in general | no |
| 4 | Objective indifferent to which genomes carry a flux at fixed (c, E), and along shared realisable E-paths; packaging selected only via realisability and per-genome nonlinearities (§4.3) | yes (instantaneous proof under assumptions) | silent | silent |
| 5 | Catabolic fan-in concentrated at interfaces where both potentials are biogenic | yes [Conj] | silent | silent |
| 6 | Local vs global optimisation depends on abiotic decay rate of the potential (light: local; stable chemical potentials: possibly global) | yes (Vallino) | silent | silent |
| 7 | Composition carries lagged environmental memory (window ≤ ~1/m in the neutral limit; an estimand under selection) | silent | yes, loosely | yes, given an observation model |
| 8 | Anticipatory expression with phase lead to predictable forcing | yes if it raises integrated σ | yes, for any predictable variable | via selection of clocks |

The decisive asymmetry: (T) makes **quantitative predictions in physical units from measured
chemistry**; (V) does not. (L) adds the memory prediction, which (T) does not make.

---

## 6. Testable consequences on the team's data types

For each: prediction, the null that makes it trivial, and what a decisive test needs.
Notation: Δ sampling interval, n series length, τ response time.

### 6.0 What each cadence can see (analytic, S4)
Lag-one autocorrelation φ = exp(−Δ/τ). φ ≲ 0.1: the layer is aliased into white noise and
indistinguishable from measurement error. φ ≳ 0.9 over a short series: indistinguishable from a
trend. φ ∈ [0.1, 0.9] corresponds to **0.43 Δ ≲ τ ≲ 9.5 Δ**; in addition the series must span
several decorrelation times (nΔ ≳ 5τ) for τ̂ to be usable, so the resolvable window is roughly
0.4 Δ ≲ τ ≲ min(10 Δ, nΔ/5).

The τ values are illustrative orders of magnitude, **assumed, not measured**:

| layer (assumed τ, days) | τ | 6 h | daily | weekly | monthly | annual |
|---|---|---|---|---|---|---|
| expression, mRNA turnover (min-h) | 0.05 | 0.007 | 0.000 | 0.000 | 0.000 | 0.000 |
| expression program (diel) | 1 | 0.779 | 0.368 | 0.001 | 0.000 | 0.000 |
| population (doubling ~1-10 d) | 5 | 0.951 | 0.819 | 0.247 | 0.002 | 0.000 |
| community turnover (weeks) | 30 | 0.992 | 0.967 | 0.792 | 0.368 | 0.000 |
| strain sweep / SNV (months) | 120 | 0.998 | 0.992 | 0.943 | 0.779 | 0.048 |
| genome content (years+) | 1000 | 1.000 | 0.999 | 0.993 | 0.970 | 0.694 |

Consequences: **mRNA-level expression (τ minutes-hours) is invisible as dynamics in any
daily-or-slower series**; it appears only as white noise plus whatever is slaved to slower layers.
Weekly sampling resolves population dynamics; monthly sampling resolves community turnover;
annual revisits (Lake Winnipeg year-to-year) resolve only strain sweeps and genome content, and
then only with many years. Nyquist: forcing at period P needs Δ < P/2; diel needs Δ < 12 h, and
since any sub-Nyquist series aliases the diel cycle to a spurious low frequency, irregular
sampling times of day must be recorded and modelled.

### 6.1 Spectral separation (layer-specific bands)
- *Prediction (F1), qualified*: taxon abundance, gene content, transcript levels and SNV frequencies
  have power concentrated in distinct bands, ordered fast → slow. Two qualifications from the
  multiscale_learning simulator (same README, §5.1a; synthetic, one vote): (i) expression observed
  daily shows its *slaved* timescale (≈ 15 d), not its own (≈ 0.2 d); (ii) gene-content (≈ 297 d) and
  SNV (≈ 319 d) timescales overlapped even in noise-free truth (ordering AUC 0.38), so the
  gene-content vs SNV order is not a firm prediction. Under sequencing noise and realistic cadence
  all adjacent-layer AUCs fell to 0.3-0.7.
- *Nulls*: (a) generic red noise: one AR(1) per series, or a 1/f continuum; (b) shared forcing:
  everything inherits the seasonal spectrum; (c) aliasing; (d) compositional closure mixing spectra
  across taxa (S3); (e) measurement noise differing by assay (RNA noisier → looks "faster" if τ is
  estimated naively).
- *Decisive test*: noise-robust τ estimation (OU + white noise, exact likelihood; or the ratio
  r(2)/r(1), which is noise-invariant for AR(1)+noise), absolute or reference-normalised abundance
  rather than closed proportions, seasonal-phase and distributed-lag terms removed first, and
  comparison *within MAG* across layers (paired design cancels shared forcing).

### 6.2 Timescale ordering of perturbation responses
- *Prediction (F1 + slaving), qualified*: after a perturbation (CMO diesel dosing), the **onset**
  appears first in expression and abundance, then gene content and SNVs. The **half-response** need
  not follow that order: slaved expression peaks only after the community has reorganised. In the
  multiscale_learning simulator (synthetic, another Claude agent's code; `alife/multiscale_learning/
  README.md` sha256 `6ae8c001aef9`, §5.1c), t10 was 1.5-2.5 d for expression and 1.5 d for
  population, but t50 ran population 4.5 d < genome 5.5-9.5 d < expression 10.5-13.5 d < SNV
  16.5-22.5 d (Carbon 249247746036776960). Pre-register onset (t10) and half-response (t50)
  orderings as separate predictions.
- *Nulls*: different detection limits per assay; treatment confounded with season and pool
  (one pool per arm).
- *Decisive test*: needs replicated pools; with one pool per arm this is descriptive only (trap 6).

### 6.3 Additivity across layers
- *Prediction (F1)*: innovations (one-step prediction residuals) of each layer are approximately
  uncorrelated after conditioning on slower layers; variance of a community-level observable
  decomposes into layer-attributed parts with small interaction terms.
- *Null*: additivity of *physical* layers is false by construction under slaving (§2.1), so
  testing raw layer variances is a category error. Closure (proportions sum to 1) induces
  negative covariances that look like non-additivity.
- *Test*: hierarchical state-space model (slow latent → fast latent), compare with a model with a
  fast → slow feedback term (Granger-type), within MAG.

### 6.4 Low dimensionality per layer
- *Prediction (F1, Tikhonov)*: fast layers slaved to slow ones should have *lower* intrinsic
  dimension conditional on the slow state than unconditionally; the unconditional dimension of
  the expression layer can be high.
- *Nulls*: n samples ≪ features makes every estimate low-dimensional (rank ≤ n); shared forcing
  produces one dominant axis (season) in any layer; closure removes one dimension.
- *Test*: participation ratio or cross-validated PCA dimension with n held fixed across layers, after
  removing seasonal-phase terms; compare conditional (residual on slow layer) vs unconditional.

### 6.5 Where plexitrophic dependencies should sit in the frequency domain
- *Prediction [Conj, §4.3]*: genome-layer (obligate) dependences involve supplies whose variance is
  concentrated at low frequency; expression-layer (facultative) dependences track high-frequency
  supply variance. Recipient–donor coherence for obligate pairs should be highest in the slow band.
- *Null*: abundant donors are both stable and common, so slow-band coherence is expected from
  prevalence alone; use prevalence- and detectability-matched non-partner donors (as in P2/P3).

### 6.6 MEP-specific (reading T)
- *Allocation follows thermodynamic kinetics*: per sample, expressed catabolic allocation across
  redox pathways (O₂ respiration, nitrification, denitrification, sulfate reduction, sulfide
  oxidation, methanogenesis/methanotrophy, hydrocarbon oxidation) is predicted better by
  capacity × kinetics × F_T than by capacity × kinetics. Needs RNA + DNA + chemistry (activities,
  T, pH). This is P6 in §9; it tests predictive value for expression, not maximisation of
  dissipation, which needs measured fluxes.
- *Interface fan-in*: d_catab ≥ 2 candidates (organisms needing biogenic donor and biogenic acceptor
  from different partners) concentrate at redox transitions. Needs depth-resolved MAGs + chemistry
  + a way to call potentials biogenic vs abiotic.
- *Horizon / transfer function*: allocation follows forcing with low-pass behaviour above 1/H and a
  phase lead at predictable frequencies. Needs sub-daily sampling over ≥ 5-10 days. Not available
  in any dataset I know of in this project (unverified).

---

## 7. Simulations [Sim]

Synthetic only. Sum-of-OU components plus white measurement noise, exact Gaussian likelihood on
the sample times (works for irregular spacing). Time in units of the sampling interval Δ.
Signal variance 1, measurement-noise variance 0.3 unless stated. Reps per cell: see output header.

*(Filled in by the coordinator, claude-plexi-2, on 2026-11-19 from the agent's own output. The agent session ended before writing this section, so the output is pasted verbatim with no new interpretation. Script scratch/deep/manifold_timescales_sim.py 96a7e139ad6c, output 378ae83bb604.)*

```
reps per cell = 100; numpy 2.3.5
==============================================================================
S4  Analytic: lag-one autocorrelation exp(-Delta/tau). ~0 = aliased into white noise;
    ~1 = indistinguishable from a trend over a short series. Informative band ~0.1-0.9.
  tau (days)                                    6 h    daily   weekly  monthly   annual
  expression (mRNA half-life ~min-h)   0.05    0.007    0.000    0.000    0.000    0.000
  expression program (diel)            1.00    0.779    0.368    0.001    0.000    0.000
  population (doubling ~1-10 d)        5.00    0.951    0.819    0.247    0.002    0.000
  community turnover (weeks)          30.00    0.992    0.967    0.792    0.368    0.000
  strain sweep / SNV (months)        120.00    0.998    0.992    0.943    0.779    0.048
  genome content (years+)           1000.00    1.000    0.999    0.993    0.970    0.694
==============================================================================
S1  LR(M2: two OU + noise  vs  M1: one OU + noise), regular sampling Delta = 1
    total signal variance 1 in every scenario; tau in sampling intervals
  n= 30 N1a one OU tau=2, noise .3       95th pct LR = 0.06
  n= 30 N1b one OU tau=5, noise .3       95th pct LR = 0.12
  n= 30 N1c one OU tau=15, noise .3      95th pct LR = 0.27
  n= 30 calibrated threshold (max of null 95th pcts) = 0.27
  n= 50 N1a one OU tau=2, noise .3       95th pct LR = 0.47
  n= 50 N1b one OU tau=5, noise .3       95th pct LR = 0.41
  n= 50 N1c one OU tau=15, noise .3      95th pct LR = 0.70
  n= 50 calibrated threshold (max of null 95th pcts) = 0.70
  n=100 N1a one OU tau=2, noise .3       95th pct LR = 0.54
  n=100 N1b one OU tau=5, noise .3       95th pct LR = 1.16
  n=100 N1c one OU tau=15, noise .3      95th pct LR = 1.50
  n=100 calibrated threshold (max of null 95th pcts) = 1.50

  Fraction of replicates with LR > calibrated threshold (power for A*, false 'multiscale' for N*)
                                                n=30     n=50     n=100    
  A1 two layers tau 1 & 20 (ratio 20)           0.02     0.05     0.15     
  A2 two layers tau 2 & 10 (ratio 5)            0.07     0.04     0.08     
  A3 two layers tau 0.3 & 20 (fast < Delta)     0.00     0.03     0.04     
  A4 two layers tau 1 & 20, fast var .2         0.05     0.07     0.14     
  A5 two layers tau 1 & 20, noise .05           0.05     0.13     0.24     
  N2 one OU tau 3 + seasonal (P=12)             0.00     0.00     0.00     
  N2b one OU tau 3 + one cycle over series (P=n)0.05     0.13     0.24     
  N3 continuum: 8 OU, tau 0.5..50 log-spaced    0.03     0.03     0.11     
==============================================================================
S2  Cross-layer ordering: tau-hat of a 'DNA' series vs an 'RNA' series of the same MAG
    statistic D = log(tau_hat_DNA / tau_hat_RNA), each from a one-OU + noise fit
    NULL: both observe the SAME latent OU (tau=5) with different noise (DNA .2, RNA .6)
    ALT : DNA = slow latent (tau 10) + noise .2; RNA = fast (tau 1) + 0.5*slow + noise .6
  n= 20: null median D=+0.38 P[D>0]=0.67 | alt median D=+0.98 P[D>0]=0.72 | 1-MAG power=0.17 | sign-test power m=10: 0.04, m=30: 0.12
  n= 40: null median D=+0.10 P[D>0]=0.54 | alt median D=+1.23 P[D>0]=0.82 | 1-MAG power=0.51 | sign-test power m=10: 0.44, m=30: 0.92
  n=100: null median D=+0.03 P[D>0]=0.54 | alt median D=+1.59 P[D>0]=0.92 | 1-MAG power=0.73 | sign-test power m=10: 0.81, m=30: 1.00
==============================================================================
S3  Closure: slow taxon (tau 20, log-sd .3) next to a dominant fast taxon (tau 1, log-sd 1, ~5x)
    plus 8 background taxa (tau 5, log-sd .3); n = 60. tau-hat of the slow taxon, one-OU + noise fit
  log absolute : median tau-hat  11.26  IQR [6.76, 17.09]  (true 20)
  log relative : median tau-hat   1.85  IQR [1.29, 3.27]  (true 20)
  CLR          : median tau-hat   6.14  IQR [3.26, 9.33]  (true 20)
```

**Reading of §7 (author, added 2026-11-19 with the §12 revisions; no re-run).**
- *S1*: the calibrated null thresholds are small (95th percentile LR 0.06-1.50), and power against
  them stays low at every length tried. A one-cycle forcing (N2b) and a timescale continuum (N3)
  are rejected about as often as genuine two-layer truths, so a rejection would not identify
  layers. A P = 12 seasonal sinusoid (N2) is never rejected: a sum of positive OU covariances
  cannot fit the negative autocorrelation at half a period, so the two-OU model gains nothing; this
  is a property of the model class, not reassurance about seasonality.
- *S2*: per-MAG power is the fraction of alternative replicates above the null 95th percentile.
  **The pooled sign-test power (script lines 202-206) is conditional on independent MAG signs**
  and treats the null rate P[D > 0] as known from 100 reps. MAGs sampled on the same dates share
  forcing and measurement batches, so their signs are positively correlated and the effective m is
  smaller. Before use as field power, calibrate the pooled test by resampling whole sampling
  events (block bootstrap or date permutation), not MAGs.
- *S3*: closure shortens apparent timescales of a slow taxon next to a dominant fast one (median
  τ̂ 1.85 in relative abundance vs 11.3 absolute; truth 20; one-OU fits also underestimate τ in
  absolute data at n = 60). CLR recovers part of the gap (6.1). Absolute or reference-normalised
  abundance is needed for layer-timescale estimates.

---

## 8. What our datasets can and cannot support (unverified; checks to run)

| dataset | what is known from team notes | cadence / n | chemistry | supports |
|---|---|---|---|---|
| CMO mesocosms (2 winters) | 97 date-matched DNA+RNA pairs over 39 dates; one treated + one control pool per winter; ~34% of RNA reads on unbinned contigs | ≈ 20 dates per winter, spacing **unknown to me** (check) | **unknown** (diesel dosing; check for O₂, nutrients, hydrocarbons in metadata) | per-MAG layer ordering (DNA vs RNA vs SNV) pooled across MAGs (§7 S2); descriptive perturbation response; allocation test only if chemistry exists |
| Lake Winnipeg / Namao | multi-year station revisits; 2024/25: 17 stations, 21 station-month cells | a few months per year per station; annual revisits | temperature, conductivity, depth, nutrients "where available" (P3 protocol) | spatial and year-to-year composition only; **no** within-season dynamics per station; slow-band (inter-annual) questions with very few points; lagged-memory test (§2.4) only with pooling across stations |
| Amundsen HS22-003 size fractions | 3 fractions, 1 event | snapshot | CTD likely (check) | P4 particle-vs-free contrast; interface fan-in candidates only descriptively |
| ASGARD 2017 | 31 process keys × 3 fractions, 16S only | survey | check | fraction-resolved membership, no genome layer |
| MAG cohorts (Ebb & Flow, Freshwater, metaMDBG, Chesterfield) | dereplicated HQ MAGs | n/a | unknown | genome-layer (slow) dependence calls; SNV layer if read mapping per sample exists |

The only cohort with a genuine time series and paired RNA is CMO. Nothing in the notes indicates
sub-daily sampling anywhere, so the expression layer's own dynamics are not observable in our data;
only its slaved component is. Depth-stratified chemistry + RNA (needed for the strongest MEP test)
has not been shown to exist locally.

---

## 9. Recommendation

**How to carry the perspective.**
- **Paper 1: one framing paragraph** in the Discussion, not the Introduction. Content: (1) MEP with
  finite horizons (Vallino) scores fluxes, not which genomes carry them (§4.3, instantaneous
  invariance), so function conserved under carrier turnover is unsurprising but not derived; that is
  P1's question; (2) MEP alone is silent on plexitrophy, which can only enter through realisability
  and per-genome symmetry breakers (incompatibility, transport, crowding); (3) obligate vs facultative dependence corresponds to
  slow vs fast layers. Cite (V) at most once. Do not use "mechanistic interpretability" as a claim;
  use it as a description of method (find features, test their causal use).
- **Graded estimands (§3) should go into Paper 1's protocols now**, independent of the perspective:
  π(X) and the posterior of (d, κ) for P2, mixed membership for P3, posterior topology indices for P5. This is the
  part of grid's perspective with the most immediate payoff and least risk.
- **Not a Paper 2 yet.** The perspective has one quantitative anchor (MEP) and one formal
  correspondence (the replicator update has the Bayesian update form; literal inference requires an
  observation model); neither has been tested on our data. Write a **protocol P6** first; Paper 2
  follows if P6 gives a result in either direction.

**Strongest first pre-registerable test: P6, now drafted as `protocols/P6-thermodynamic-allocation.md`**
(DRAFT, not frozen; read 2026-11-19, sha256 `609046bfda52`; owner claude-plexi-2). This note no
longer carries its own version; the protocol is authoritative. Summary of how it answers the
reviews of this section:
- *What it can and cannot show* (Astra 249247324022685696): it scores **raw transcript counts**
  against parameter-free predictors per assumed unit of enzyme. A win shows that thermodynamic
  driving force carries predictive information about expressed allocation; it does **not**
  measure flux J and does **not** show that dissipation is maximised. Testing maximisation needs
  measured rates and a comparison of Σ_r J_r(−Δ_rG_r)/T over the same volume and horizon across
  realisable allocations (Astra 248971565584531456).
- *Models* (π_r ∝, normalised over pathways; cap = marker-normalised DNA capacity):
  M_cap = cap; **M_capK = cap·F_K** (the baseline to beat); **M_capKT = cap·F_K·F_T** (gate only);
  M_power = cap·F_K·F_T·(−Δ_rG/n_e). **Every ΔG enters per electron (Δ_rG_r/n_e)** in both F_T and
  the power weight, so how a reaction is written cannot rescale it (Carbon 249247746036776960).
- *Outcome and likelihood (frozen)*: raw counts, Dirichlet-multinomial, one common dispersion
  estimated on training folds under M_capK and held fixed for all models.
- *Estimands*: primary Δ = S_power − S_capK; pre-registered split Δ_gate = S_capKT − S_capK
  (thermodynamic *feasibility* information) and Δ_weight = S_power − S_capKT (energy-*proportional*
  allocation). A claim of incremental thermodynamic information rests on Δ_gate; a claim about
  power-proportional allocation rests on Δ_weight. Carbon's mep_redox T5 shows F_T gating alone can
  produce much of an RNA-power association, which is why the split is required.
- *Bands*: supported (lower bound > +δ), equivalent (interval inside [−δ, +δ]), **against (upper
  bound < −δ)**, inconclusive otherwise; δ frozen as a smallest meaningful difference before any
  calibration outcome is seen.
- *Data gate first*: metadata-only check for paired DNA + RNA + in-situ chemistry, ≥ 20 independent
  events across ≥ 2 redox regimes. If it fails, P6 is "not identifiable on local data".

**Fallback / precondition test on CMO (reading F1): layer ordering** (also in the P6 file).
- **Not** separate per-layer fits with a sign test against 0.5: §7 S2 shows unequal DNA/RNA noise
  biases separate fits at n ≈ 20 (null P[D > 0] = 0.67, ≈ 60% false positives with 30 MAGs;
  Carbon 249247746036776960).
- Instead, per HQ MAG, fit a **joint bivariate (DNA, RNA) state-space model** with separate
  measurement-noise terms and seasonal-phase terms; test τ_DNA > τ_RNA against τ_DNA = τ_RNA by
  **parametric bootstrap from the H0 joint fit**; combine across MAGs with resampling over whole
  sampling dates (MAG signs are not independent; §7 note). Add the SNV layer only for MAGs with
  enough depth (prediction τ_SNV > τ_DNA, which the multiscale_learning simulator suggests may fail
  even in truth; §6.1).
- Honest status: a positive result is expected under almost any biological model, so it confirms
  only that layers are separable at the sampling cadence, not that they "learn". A negative result
  (no ordering) would be informative against F1 at this cadence.

---

## 10. Provenance

- Script: scratch/deep/manifold_timescales_sim.py  sha256 96a7e139ad6c
- Output: scratch/deep/manifold_timescales_sim.out  sha256 378ae83bb604
- Run: python3 -I scratch/deep/manifold_timescales_sim.py (numpy 2.3.5; seeds in the script; reps per cell 100)
- Author: Opus subagent (Claude family). The session ended before §7 and §10 were filled; the coordinator (claude-plexi-2) filled both on 2026-11-19 from the agent's own files, with no new analysis.
- Status: cross-reviewed by Astra (Codex, 249247324022685696) and Carbon (Gemini,
  249247746036776960); both checked §7 against the output and the hashes above. Revisions in §12.
  The script and output are unchanged by the revision (no re-run).

---

## 11. References (confidence; ✓ = read today, otherwise from memory and to be checked)

- ✓ Vallino J.J., Huber J.A. 2018. Front. Environ. Sci. 6:100. doi 10.3389/fenvs.2018.00100 (main
  text HTML; supplement not read).
- Vallino J.J. 2010. Ecosystem biogeochemistry considered as a distributed metabolic network ordered
  by maximum entropy production. Phil. Trans. R. Soc. B 365:1417. (high)
- Vallino J.J. et al. 2014 (periodically forced methanotroph communities; cited by Vallino & Huber;
  journal not checked).
- LaRowe D.E. et al. 2012. Thermodynamic limitations on microbially catalyzed reaction rates.
  GCA. (cited by Vallino & Huber; medium)
- Dewar R.C. 2003. J. Phys. A 36:631; Grinstein G., Linsker R. 2007. J. Phys. A 40:9717 (critique).
  (medium)
- Martyushev L.M., Seleznev V.D. 2006. Phys. Rep. 426:1. (high)
- England J.L. 2013. Statistical physics of self-replication. J. Chem. Phys. 139:121923. (high)
- Friston K. 2010. Nat. Rev. Neurosci. 11:127; Friston K. 2013. Life as we know it. J. R. Soc.
  Interface 10:20130475. (high)
- Ramstead M. et al. 2018. Phys. Life Rev. 24:1. (medium)
- Biehl M., Pollock F.A., Kanai R. 2021. A technical critique of some parts of the free energy
  principle. Entropy 23:293. (medium)
- Aguilera M., Millidge B., Tschantz A., Buckley C.L. 2022. How particular is the physics of the free
  energy principle? Phys. Life Rev. 40:24. (medium)
- Bruineberg J. et al. 2022. The Emperor's New Markov Blankets. Behav. Brain Sci. 45:e183. (medium)
- Shalizi C.R. 2009. Dynamics of Bayesian updating with dependent data and misspecified models.
  Electron. J. Stat. 3:1039 (the replicator–Bayes remark appears in his notes; exact source
  **uncertain**).
- Harper M. 2009. The replicator equation as an inference dynamic. arXiv:0911.1763. (medium-high)
- Czégel D., Giaffar H., Tenenbaum J.B., Szathmáry E. 2022. Bayes and Darwin: how replicator
  populations implement Bayesian computations. BioEssays 44:2100255. (medium)
- Watson R.A., Szathmáry E. 2016. How can evolution learn? Trends Ecol. Evol. 31:147. (high)
- Vanchurin V., Wolf Y.I., Katsnelson M.I., Koonin E.V. 2022. Toward a theory of evolution as
  multilevel learning. PNAS 119:e2120037119. (high for title/venue; article number medium)
- Del Moral P., Guionnet A. 2001. On the stability of interacting processes with applications to
  filtering and genetic algorithms. Ann. Inst. H. Poincaré Probab. Statist. 37:155. (medium)
- Tikhonov A.N. 1952. Mat. Sb.; Fenichel N. 1979. Geometric singular perturbation theory for ODEs.
  J. Differ. Equ. 31:53. (high)
- Haken H. 1983. Synergetics. Springer. (high)
- Inchausti P., Halley J. 2001. Science 293:655 (reddening of population spectra). (medium)
- Vasseur D.A., Yodzis P. 2004. The color of environmental noise. Ecology 85:1146. (medium)
- Louca S., Parfrey L.W., Doebeli M. 2016. Decoupling function and taxonomy in the global ocean
  microbiome. Science 353:1272. (high)
- Pfeiffer T., Schuster S., Bonhoeffer S. 2001. Science 292:504. (high)
- Costa E., Pérez J., Kreft J.-U. 2006. Why is metabolic labour divided in nitrification? Trends
  Microbiol. 14:213. (medium)
- Daims H. et al. 2015. Complete nitrification by Nitrospira bacteria. Nature 528:504. (high)
- Elhage N. et al. 2022. Toy models of superposition (Transformer Circuits Thread); Bricken T. et al.
  2023. Towards monosemanticity; Olah C. et al. 2020. Zoom in: an introduction to circuits. Distill.
  (high)

---

## 12. Revisions after cross-family review (Astra 249247324022685696; Carbon 249247746036776960)

Made 2026-11-19 by the author, in place, **without re-running the simulation** (script
`96a7e139ad6c` and output `378ae83bb604` unchanged). Before this revision the note was
`712ca75c9367` (as filled by the coordinator).

| # | Reviewer point | Change |
|---|---|---|
| A1 | §4.3 proof gives instantaneous invariance at fixed (c, E), or trajectory invariance only along a fixed/shared E-path; equal initial E and summed growth do not give equal dE/dt | §4.3 restated as an instantaneous proposition plus a trajectory corollary for shared realisable E-paths; the earlier whole-trajectory claim is withdrawn. Added: the objective's indifference to packaging does not remove dependence within a packaging, so MEP alone is silent on plexitrophy. §0 item 4 and §5 row 4 updated. |
| A2 | §9 P6 scores counts against power per assumed enzyme: predictive, not J or maximisation; needs a capacity × kinetics baseline vs capacity × thermodynamic kinetics; frozen raw-count DM outcome with common dispersion; a band below −δ | §9 now summarises and points to `protocols/P6-thermodynamic-allocation.md` (609046bfda52), which implements all of these; §0 item 7 and §6.6 reworded. |
| A3 | §3.1 must use the current rule (every feasible wiring: d ≥ 2 OR κ = 1); exp(H) is a separate descriptor (99:1 gives ≈ 1.06 with d = 2); report π itself | §3.1 rewritten accordingly; "Pr[π > 0.95]" and d_eff-as-criterion removed; §3.4 and §9 references updated. |
| A4 | §2.4: Bayesian *form* only; literal inference needs an observation model; mixing after selection is a next-step predictive; (1−m)^t decay proved only under neutral fitness | §2.4 narrowed: form identity; observation model stated as an assumption; neutral-limit decay proved; under selection the memory window is an estimand. §0 item 1, §2.4 prediction and §5 row 7 updated. |
| A5 | S2 pooled power assumes independent MAG signs (script l. 202-206) | Labelled in §7 note, §0 item 6 and §9; calibrate by resampling sampling events before field use. |
| C1 | §0 wrongly said seasonal forcing "also beats the null"; N2 was 0%; only N2b drift and N3 continuum did. Naive sign test at n ≈ 20 has ≈ 60% false positives | §0 item 6 corrected (N2 = 0%, monotone OU covariance; N2b 0.05-0.24; N3 0.11 at n = 100); FPR checked: reject at ≥ 20/30, P = 0.600 under p = 0.67. §9 fallback now uses a joint bivariate state-space model with parametric bootstrap. |
| C2 | P6: ΔG per electron in F_T and power; gate-only model cap·F_K·F_T | Already in the P6 protocol; §9 aligned (M_cap, M_capK, M_capKT, M_power; Δ_gate / Δ_weight split). |
| C3 | §6.1/§6.2 qualified by multiscale_learning: onset vs half-response order differ; gene content vs SNV overlap (AUC 0.38) | §6.1 and §6.2 qualified, citing `alife/multiscale_learning/README.md` (6ae8c001aef9) §5.1a/c as synthetic, one-vote evidence. |
| C4 | §2.1: sum of L AR(1) is ARMA(L, L−1); with white noise ARMA(L, L) | §2.1 corrected. |

**Second round (Astra 249248487229014016, after e992457ec82e):** "exactly Bayes" (§1) and
"literal identity" (§9) replaced by "has the Bayesian update form; literal inference requires an
observation model", consistent with §2.4. §3.1 now makes non-empty, complete wiring sets explicit
in the posterior event, adds π_feas with its stated denominator, and reports infeasible and
unresolved masses separately. The P6 reference is updated to `609046bfda52` (current protocol).

Also corrected while revising: the header run line now says 100 reps (the recorded output),
and §7 gained the author's reading of S1-S3 (no new numbers).

