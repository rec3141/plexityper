# codex-plexi

Codex harness, GPT-6. I maintain this file. My project relay is #plexityper in workspace plexityper.

## Current state

I proposed testing whether candidate plexitypes recur across years more often than shared seasonality predicts (relay message 248871704742772736). Awaiting grid's assignment and dataset choice. No experiments started or jobs submitted.

## Handoff

- Read AGENTS.md (SHA-256 prefix 205ec22800b5) and the full plexichat.txt. Relay request 248872795412807680 covers the startup round-trip.
- Claim deliverables publicly to prevent duplicate work. Use top-level channel posts so the full team sees decisions.
- Confirm approvals come directly from grid, registered as human. Before running an experiment, post the protocol hash and obtain cross-family review so the decision bands are fixed before seeing results.

## Agent Relay core

Read the Team and chat section of AGENTS.md at SHA-256 prefix `26899dd82090`. The core delivers messages directly into the running session. On startup, register, subscribe, and confirm the identity before posting. Re-subscribe after a rename.

```bash
P=/home/grid/.agentworkforce/relay/bin/agent-relay-probe
S=/run/user/1000/agent-relay/relay.sock
"$P" relay socket-request --socket "$S" --method POST --path /register --body '{"name":"codex-plexi"}'
"$P" relay subscribe-channel --socket "$S" plexityper
"$P" relay whoami --socket "$S"
"$P" relay agents --socket "$S"
# After whoami confirms codex-plexi, send the startup check-in.
"$P" relay send --socket "$S" --channel plexityper --body 'On the Agent Relay core as codex-plexi.'
```

2026-11-18 setup attempt: the socket returned HTTP 403 (`not_a_relay_session`) for registration, subscription, and identity lookup. The core check-in is pending a recognized live session.
Post-relaunch retry: registration, subscription, and whoami again returned 403 not_a_relay_session. The probe shell runs under tmux > bash > codex > bash. No confirmation posted.

Resolved on 2026-11-18. The core accepted registration and subscription, whoami confirmed codex-plexi, and the top-level check-in reached #plexityper (message 248881879352524800).

## Literature discussion

Grid opened the literature and hypothesis discussion in 248913233255514112. My contributions are `notes/plexitype-recurrence.md` and `notes/vision-codex-plexi.md`. Opus owns the group synthesis at `discussion/plexi-concept-review.md`. Core discussion points include explicit occurrence boundaries, separating recurrence and dependence, donor-role replacement, and Markov equivalence of three-node noncolliders. Opus reviewed my vision in 248914484261208064. I shared a metadata-only CMO eligibility check and its replication and annotation limits in 248914574040285184 and 248914753246117888. Venrick 1971 (10.2307/1934149) is listed for full-text retrieval. Discussion continues via pushed relay messages. No jobs submitted.

P4 discussion draft: protocols/P4-size-fraction-pilot.md (aea1c9ef57e0); Carbon review requested. Inventory: notes/paired-events-and-fractions.md (dd6894ced0ce). Assigned cross-family P2 review when Carbon posts the protocol. No jobs submitted.

P2 cross-family design review closed on a63db538fea3. Carbon closed P4 review on 0ed5e0034020 and P3 review on 33de07d6ebf6 in relay message 248924302371897344. Protocols remain drafts awaiting grid.

P1 data-blind Tara W6 discussion draft at protocols/P1-song-singer-tara.md (dff063df9592), posted for cross-family review in relay 248931950295171072. P5 review on independent chemistry labels and primary statistics posted in 248931973443534848. Await Carbon’s paper outline for cross-check. No P1 data analyzed.
P1 revised after Opus review (9ff69defe5be). Signed btuB response, archaeal supplier rules and depth-layer replication gates added. Paper outline 1de0183786f1 reviewed and edits proposed in 248932231170932736. Await grid before experiments.
P1 draft 812d46c677e2 incorporates Carbon review 248933060170924032. Station holdouts for identity, dominant-genus holdouts for guild transfer, exact size/depth strata, prokaryotic marker-depth normalization and secondary remodeling model. Posted in 248933164856557568. Feasibility unmeasured.
Standing down after closing-summary review. P1 hash matches. Outline hash, conditional mapping costs and distinction between genomic candidates and experimentally established obligation flagged in relay 248933325003472896. No active jobs. Await grid’s protocol freeze and compute approval.
Carbon closed cross-family review of P1 draft 812d46c677e2 in relay 248933478447890432. Await grid.

2026-11-20 — Grid directly approved P2 and P6 execution (249688876280889344). Carbon owns P2, Opus owns P6. I claimed independent results review (249688959718178816). Await pushed reports and frozen hashes. Stage 0 is a pipeline check only. Review P2 unresolved denominators and event gates, P6 provenance blocks, training-fold dispersion and verdict bands. Neither run has reported a result at assignment.

2026-11-20 — P6 audit. Pooled benchmark Against for 27 fixed scenario-by-seed groups sharing three seeds (primary Δ −9.993, interval −11.775 to −7.895); field gate inconclusive on inspected metadata. Calibration audit closed for large directional effects and coverage gates; equivalence and near-margin errors unmeasured. Opus recorded the conditional-on-gate interpretation and limits in 249693436487581696. P2 results review remains pending.

2026-11-21 — P2 freshwater. Audited 497 staged manifest entries and 232 v2 library interval files. Event analysis awaits grid’s library-to-campaign/station map. Carbon owns the runner and proposed P2-A1. Reviewed A1 fc3d52312440 and runner 67cdfa299398 in relay 249778991556640768. Primary Null 2 effect intervals need a validated estimate of mean prevalence under the chosen reference distribution. A single shared hub gives exchangeable ranks under the null but its spoke average depends on the realized hub. Additional open items are the observed-row exemptions in Null 2c and the Q1 falsification rule when only primary Null 2 is evaluable. The matched observed Q2 verdict uses the observed interval after adequacy gates. Stage 0 remains descriptive.

2026-11-24 — Closed P2-A1 technical review for addendum 2a5c3671a1a2 and runner 505a4fe6602f in relay 251084385889468416. Both Q1 and Null 2c now retain each available full-fiber estimate. All 11 frozen tests and both runner self-tests passed; independent 3:1 fixture reproduced .375 primary and .75 all-component sensitivity with Option A enabled for the fixture. Grid decides whether to adopt Option A before analyzing event outcomes.

2026-11-24 — Closed pre-outcome P2 event-map review in relay 251090830160740352. Compared map e435af5a8b87 to the full BAM roster and audit 3d831a55266b. All 232 keys appeared once. The map includes 57 libraries in 29 proxy field events and skips 175 with a documented reason per library. Agy claimed the freshwater run with Option A off (251091191323869184). Independent results review awaits released outputs.

2026-11-24 — P6 Saanich metadata review (251091800198397952). Of 43 distinct cruise–depth keys across eight cruises, 42 have jointly finite oxygen, nitrate, ammonium, nitrite, sulfide, methane and nitrous oxide measurements. Cruise 73 at 200 m lacks the last three. All 43 lack CO2. The Dryad README specifies nM for methane and nitrous oxide. Each cruise spans the 10 µM oxygen threshold. Cast identity, bottle pairing, size-fraction compatibility and the carbonate/organic-carbon assumptions remain open. This review used metadata only.
