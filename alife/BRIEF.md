# Artificial-life simulations for plexitrophy: shared brief

Requested by grid (2026-11-18): build different kinds of artificial-life simulations or games
based on the proposed principles of plexitrophy, and think about **prediction by synthesis**.

## Prediction by synthesis
1. Build generative worlds from explicit mechanisms: leaky public goods and gene loss (Black
   Queen), incompatible chemistry, redox gradients and dissipation (MEP, after Vallino), and
   multi-timescale adaptation.
2. Let plexitrophy, plexiomes and plexitypes **emerge**. Never hard-code them.
3. Know the ground truth: who actually feeds whom, d / m / cycle per organism and time, and
   topology.
4. Emit **synthetic observations** in the same shape as our real data, then run the team's
   inference ideas (P1–P5) on them and score recovery against the truth.
5. Ask which mechanisms produce which observable signatures. A signature that only one mechanism
   produces is a prediction to test on real data.

The synthetic observation formats are:
- per-sample abundance tables (genomes × samples);
- per-genome KO presence (a per_mag-like TSV: `mag_id, ko`);
- incompleteness, by dropping genes as if from contig loss;
- optional "RNA" expression;
- size fractions or particles.

## Read first
- /data/dev/plexityper/AGENTS.md (rules)
- the top of /data/dev/plexityper/discussion/plexi-concept-review.md:
  - the standing perspective (learned representations, multi-timescale, MEP / Vallino);
  - "rigid rules won't fit": prefer graded quantities;
  - the one-sentence plexitrophy definition; d_anab / d_catab / d_total; c as a simple cycle ≥ 3.
- /data/dev/plexityper/discussion/deep/*.md (theory notes: fan-in nulls, computation,
  identifiability, P3 power)
- protocols/P1–P4 and the P5 statement in the review doc.

## Rules for every simulation
- Directory: /data/dev/plexityper/alife/<your-sim>/ (single writer per directory).
- Include a README.md: mechanisms, parameters, the ground-truth outputs, the synthetic
  observation outputs, how to run, seeds and the runtime.
- Python: numpy/scipy, standard library, `python3 -I`. No pip/conda installs, and no network
  installs.
- Browser games: one self-contained HTML/JS file, with no external dependencies.
- Keep any single run under about 15 min on one core, and the whole sweep under about 45 min.
  The machine is shared.
- Do NOT open anything under /data/dev/plexityper/data/ (blinded real data).
- No chat posts and no cluster jobs.
- Separate what the simulation shows from speculation. A simulation shows that something is
  possible under its assumptions, not that it is true.

## Additions after Astra's check (relay 248974646464397312)
- **Realised transfers are not required partners.** A recipient may *use* two suppliers when
  either alone would suffice. Log realised transfers separately. Score *necessity* against
  viable alternatives under the same medium, endpoint and time horizon, e.g. by knock-out or
  counterfactual removal in the simulation.
- **Observer inputs must be rich enough for P2's logic.** Emit contig/CDS IDs, feature order,
  CDS records with blank KO fields, and genome completeness, not just `mag_id, ko`.
- **Shared event IDs.** Keep shared event IDs across paired assays (DNA/RNA) and fractions, so
  held-out events can't leak into training (P1, P4).
