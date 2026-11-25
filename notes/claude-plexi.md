# claude-plexi (Opus) state

Updated 2026-11-25 (chain resubmitted after 3 download failures). Resume from "ACTIVE NOW".

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
- **Jobs (resubmitted 2026-11-25):**
  - MGnify contigs 5864186 DONE; pool 5865308 DONE; gpucheck 5865307 PASS (cuda x3, H100).
  - First download wave: 90/93 runs .done. 3 failed, so the old chain (5865309/10, 5865424)
    was cancelled by SLURM.
    - SRR3724482 (DNA, task 4) failed on a stale prefetch lock left by cancelled 5859913_4.
      The lock, .prf and .tmp were removed.
    - SRR3719367 and SRR3719449 (RNA, tasks 43 and 82): prefetch served **SRA Lite**
      (.sralite, simplified qualities), and the script only read .sra.
    - Fix: p6_saanich_download.sh 6c3af7d8c14c uses .sra or .sralite (explicit if/elif; the 24add63f0cd4 ls|head version failed under pipefail, codex 251217963302436864) and logs the source.
      Note the SRA Lite runs in the results (bbduk/strobealign do no quality trimming).
  - **New:** dl task 4 6015941; dl RNA 43,82 6015942 (both 2 d; replace cancelled 6015859/60).
    mapdna 6015861 (afterok 6015941), then depths 6015862, then mag 6015863.
    The RNA downloads are not in the DNA chain.
  - Watcher (local, 30-min polls) watches 6015941/42 + 6015861 and exits on any failure.
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
