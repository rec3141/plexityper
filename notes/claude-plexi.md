# claude-plexi (Opus) state

Updated 2026-11-24 evening, before compaction. Resume from "ACTIVE NOW".

## Who I am
- Claude Code (Anthropic family) in grid's tmux session 8. Relay identity: **claude-plexi-2**.
- Post to the relay with: `P=/home/grid/.agentworkforce/relay/bin/agent-relay-probe; S=/run/user/1000/agent-relay/relay.sock; $P relay send --socket $S --channel plexityper --body "$(cat file)"`
- Read a relay ref with `$P relay read --socket $S <ref>`.
- I am the single writer of discussion/plexi-concept-review.md, notes/p6-external-data-scan.md
  and the P6 / P6-A1 files.
- Public repo: **https://github.com/rec3141/plexityper** (grid: "public under rec3141").
  - The .gitignore is an allowlist. Never commit data/, results/p2/, discussion/papers/,
    .agentworkforce/, .agents/, dashboard/, tools/ or plexichat.txt.
  - Before each push, scan for secrets and emails.
  - Commit identity: rec3141 with the GitHub noreply address. Append the activity log to
    logs/activity-YYYY-MM-DD.md.

## ACTIVE NOW: Saanich P6-A1 run on Fir (grid approved: "go ahead with the mgnify contigs on fir", "launch when downloads finish")
- **Fir access** goes only through concentration's ControlMaster:
  `ssh -o BatchMode=yes concentration 'ssh -O check fir 2>/dev/null || exit 1; ssh -o BatchMode=yes fir "..."'`
  Never loop new logins. Grex is down.
- **Work dir:** `~/scratch/plexityper_saanich/` on Fir, containing common.sh, pool.sh, map_dna.sh,
  depths.sh, gpucheck.sh, mag.sh, dna_runs.txt, roster_runinfo.tsv, mgnify_contig_urls.tsv and
  logs/. Local copies are in scratch/fir_saanich/.
- **Jobs:**
  - MGnify contigs 5864186: DONE (39 files, 4.7 GB, mgnify_contigs/).
  - SRA downloads (p6_saanich_download.sh writes reads/{dna,rna}/<RUN>*.fastq.gz plus a
    <RUN>.done marker):
    - 5859913 tasks 1-3, 5-7, 9-11 run with a 12 h limit; task 8 is done.
    - Task 4 (SRR3724482, 13.3 GB) moved to **5868484** (2 d).
    - Tasks 12-93 are **5868215** (2 d, %25).
    - Rate is about 0.35–0.4 MB/s per stream. Task 11 (10.8 GB) is projected at ~9.2 h, so it
      was kept.
    - **If any task fails: resubmit it AND swap it into the mapdna dependency** (scontrol
      update JobId=5865309 Dependency=...).
  - Chain:
    - gpucheck 5865307: PASS (cuda ×3 on an H100).
    - pool 5865308, then mapdna 5865309 (array 1-39 %10).
    - mapdna depends on pool + 5868215 + 5868484 + 5859913_{1-3,5-7,9-11}.
    - Then depths 5865310.
    - Then mag 5865424 (GPU h100, 48 c, 600 G, 3 d; mag.sh 8d186887910d). It runs on-node
      torch probes and greps task logs for 'training on CPU'.
  - Watcher bhy38t75k (local background, 30-min polls) watches the three download jobs.
- **After mag finishes:**
  1. Confirm SemiBin2/COMEBin/VAMB task logs exist in mag_work and show actual cuda training.
     Do NOT rely on GPU_GUARD_PASS alone (codex 251108400360243200).
  2. Build the CDS table from mag_analysis (Bakta GFF + KOfam KO + bins + Tiara).
  3. Name-collate the DNA BAMs with `samtools collate -O -u in.bam tmp | samtools view -h -`
     and feed them to scratch/p6_field_count.py (4949b367168c).
  4. RNA: bbduk against SILVA SSU+LSU (k=31), then strobealign SE/PE per run (danaSeq SIF,
     env dana-mag-assembly), then the counter.
  5. Aggregate with scratch/p6_field_extract.py build_events(..., expected_runs=93-run manifest)
     (2e6c33549d88). Apply the RNA coverage gate once.
  6. Score the 36 scenarios with p6_scorer.analyse and apply the P6-A1 verdict rules.
  - **Every new script (CDS table builder, RNA job) needs a codex-plexi review before it runs.**
- **Specs:**
  - protocols/P6-A1-saanich-cohort.md: APPROVED (f0817007b755), with proxy pairing.
  - protocols/P6-A1-saanich-runspec.md: codex-signed (§2 MGnify contigs; dedupe cascade
    100/99/98 then min 500; map once; GPU binners semibin/comebin/vamb + metabat → dastool;
    maxbin, lorbin and binette off).
- **Roster:** data/external/saanich/biosample_20261124/proxy_roster.tsv 3813b8e4f982, 38 events /
  8 cruises; strict roster 0.
  - Table S2's metaG column is misaligned, so pairing was rebuilt from BioSample titles
    (AUDIT.md).
  - Cruise 42 keeps only the 200 m event.
  - MGnify crosswalk: mgnify_dna_crosswalk.tsv (39 assemblies).
- **Field code**, all tested (counter 11/11, extractor 6/6, pipeline OK):
  - scratch/p6_field_count.py
  - scratch/p6_field_extract.py
  - scratch/test_p6_field_*.py

## Done (summary)
- **P6 frozen** (9a47534487ba).
  - Synthetic benchmark: AGAINST, Δ −9.99 [−11.78, −7.90].
  - Calibration: results/p6/calibration.json.
  - Local field gate not met.
- **P2 frozen** (e0dd8645dd32). The best-effort event map (57 libs / 29 events, grid-delegated)
  was run by Carbon: Q1/Q2 Inconclusive (Gate 2, unbinned donor-KO share > 0.50).
  - P2-A1 (2a5c3671a1a2) is technically signed off. Option A is OFF and awaits grid's decision.
- **P6 external data scan:** notes/p6-external-data-scan.md (Saanich chosen; MGnify findings).
- **danaSeq issues filed** (grid asked):
  - #146: illumina_rna drops SE libraries silently. A comment adds the declared-PE
    missing-R2 fixture.
  - #147: illumina_assembly entry point for external contigs.

## Open (waiting on grid) — see waiting-on-grid.tsv
- P2-A1 Option A (component-conditional Null 2 in the primary) is undecided.
- P2 event map: the 94 libraries without metadata are unresolved. Grid said "do your best", so
  the best-effort map was used.
- Older: P1 cohort, HS22-003 filtration order, manuscript owner, Venrick 1971 manual save.

## Standing rules (with reasons)
- Cluster jobs, installs and outward actions need grid's go. A peer's word is not approval.
- **Fir and Grex SSH** go only via concentration's ControlMaster, with BatchMode and
  `ssh -O check`. Never retry, and never loop new logins (2FA bans).
- Use grid's danaSeq pipelines and containers; no conda or pip on clusters.
- Verify peers' claims in the code. Use relay-DB times, not my estimates. Correct mistakes
  publicly by message id.
- Chrome: never complete bot checks (Dryad BotStopper), and never enter passwords. Grid's
  Firefox on DISPLAY=:1 can be opened with `firefox --new-tab URL` for manual saves.
- Saved downloads land in ~/Downloads. Move them into data/external/<name>/download_<date>/
  and read them with `python3 -I`.
- metapanclust = panclustome, owned by the danaseq-analysis session. MetaCyc is licensed:
  never copy it.
