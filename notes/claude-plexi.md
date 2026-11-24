# claude-plexi (Opus) state

Updated 2026-11-20 (after grid's go to run P2/P6, 249688876280889344).

## Who I am
- Claude Code (Anthropic family) in grid's tmux session 8.
- Relay identity: **claude-plexi-2** (re-registered 2026-11-18 after a session restart; the old claude-plexi name is held by a stale session; grid approved the rename).
- Moderator and single writer of discussion/plexi-concept-review.md.

## Done this session
- Ran grid's two-hour concept discussion (relay 248913233255514112, from 15:52).
- Closing summary at the top of discussion/plexi-concept-review.md (decisions for grid listed
  there).
- Wrote protocols/P3-plexitype-recurrence.md (33de07d6ebf6); both cross-family reviews closed.
- My vision: notes/vision-claude-plexi.md (abstracts only).

## P6 execution (mine; reviewer Astra)
- Frozen protocol 9a47534487ba (freeze 0a9b44648826 + seed-family scope addendum).
- Block map: scratch/p6_block_map.py c0707beddb88 → results/p6/synthetic_block_map.tsv 84731f8fccc8.
- Scorer: scratch/p6_scorer.py 92306119ff5a (regress D001/D002 OK).
- Synthetic benchmark: AGAINST, Δ −9.99 [−11.78, −7.90] (results/p6/synthetic_benchmark.json d17de8aceeec). Secondary: scratch/p6_secondary.py 89cafe0db90c → synthetic_secondary.tsv 214c95a53ee5; Δ < 0 in 41/41 cells.
- Field gate NOT MET: results/p6/field_gate.md 1d1f0c629070.
- Calibration done: results/p6/calibration.json 050fe76fb747 (correct verdict in all scored cohorts; Equivalent band uncalibrated; φ=0.2 cells all gate-fail). Posted for Astra review.
- Waiting on grid (249690616493424640): external data? precondition test spec or skip?
- P2 is Carbon's to run (frozen e0dd8645dd32); Astra reviews.
- P2 breadth DONE per library (232; manifest_v2.tsv 6635acc1a6ee; breadth_per_bam.tsv fd9625b15747). Task 2 via depth variant 8195122_2. Per-event union waits on grid library→campaign×station map. Equal-depth pass later (BAMs on Grex).
- P2 event map (grid delegated 2026-11-24 "do your best"): scratch/p2_build_event_map.py 4d7d3291b024 -> event_map_best_effort.tsv dc0d9421aa24 (57 libs, 29 events). Carbon runs, Astra reviews. Saanich: Dryad tab opened in Chrome for grid manual save (BotStopper blocks automation).

## Staged data
- data/grex_mag_freshwater/: read-only copy of the 490 derep-representative per_mag tables,
  plus derep_clusters and the checkm2 report, from grex:~/scratch/project_nanopore/mag_freshwater
  (commit 84c3ca3eba58). Hashes are in PROVENANCE.md. KO tables are not opened (P2 blind).
- Grex is reached via concentration (it holds the ControlMaster): one check, then single
  commands.

## Open (waiting on grid)
- The 6 decisions in the summary: pipeline fixes, P1 cohort, HS22-003 filtration order,
  paywalled PDFs, freeze and compute, manuscript owner.
- waiting-on-grid.tsv: the PDF download row was approved; manual save is still needed.

## Standing rules (with reasons)
- GitHub activity log (grid 2026-11-24): keep a project repo updated with activity logs. No repo exists yet; creation and visibility wait for grid's go. Never commit .agentworkforce/, data/, paper PDFs/HTML or MetaCyc material.
- No cluster jobs or installs without grid's go (AGENTS.md, CLAUDE.md).
- Verify peers' code claims in the code before logging them. Trap 8's count was misreported
  once (53 vs 55).
- Times in logs: use the relay DB, not my estimates. I mis-stated the clock once (said 17:00
  at 16:33).
- Chrome automation cannot save PDFs to disk on grid (blob download and localhost POST were
  blocked).
- Background watcher: scratch/quiet_watch.py (read-only on the relay DB).
- Pathway definitions: metapanclust = panclustome; the danaseq-analysis session on grid owns them. Align, don't fork.
