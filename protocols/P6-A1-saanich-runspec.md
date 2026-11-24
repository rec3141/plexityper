# P6-A1 Saanich run spec: assembly, annotation and mapping (fixed before any RNA is opened)

- **Status.** Drafted 2026-11-24 by claude-plexi-2 for the approved P6-A1
  (`protocols/P6-A1-saanich-cohort.md`, f0817007b755).
- **Compute.** Fir (grid: "grex is still down i believe, try fir").
- **Read before running.** Roster and runs only: `proxy_roster.tsv` 3813b8e4f982 (38 events in 8
  cruises, codex-plexi review 251097403444609024) and its run list (93 SRA runs). No expression
  values have been opened.

## 1. Inputs
- **DNA.** 39 SRA runs for the 38 events, about 302 GB as .sra and about 514 Gbp.
  - Platforms: HiSeq 2000 (17 runs), HiSeq 2500 (21) and GAII (1, cruise 42 at 200 m).
- **RNA.** 54 SRA runs, about 306 GB and about 599 Gbp. These include the A/B replicate
  libraries, which are collapsed per event at the counting stage.
- **Download.** sra-toolkit 3.0.9 (`prefetch` then `fasterq-dump`, gzipped), into
  `~/scratch/plexityper_saanich/reads/{dna,rna}/`.
  - The run list and its sha256 are copied to the same directory.
  - Read counts per run are checked against SRA `spots`.

## 2. Reference: one gene catalogue shared by all events
All steps use grid's danaSeq pipelines and containers. No new installs.
1. **Per-event metaG assembly.** `illumina_assembly` with the Apptainer profile on SLURM, one
   sample per event, default assemblers, output `<sample>.dedupe.fasta`.
2. **Pooled reference.** Concatenate all 38 per-event `dedupe.fasta` files and dereplicate at the
   contig level with the pipeline's dedupe step (BBTools dedupe, defaults), giving
   `saanich_pooled.fasta`.
   - A pooled, not per-event, reference keeps one gene catalogue across events, so capacity and
     RNA counts are measured on the same features.
3. **Annotation and binning.** `mag_analysis` on `saanich_pooled.fasta`, with all 38 metaG
   read sets for depth, default settings.
   - Bakta annotation with coordinate-joined KOs (trap 3).
   - Binning (DAS Tool set), CheckM2 and GTDB-Tk.
   - KO calls follow the pipeline's merged annotations. The frozen table's
     `make_binned_ko_classifier` disambiguation then uses the MAG membership.
4. **RNA mapping.** `illumina_rna --from_mag <mag_analysis out>` with default QC: adapter trim,
   rRNA split against SILVA, strobealign and featureCounts per gene.
   - Human removal is skipped because these are marine samples.
   - The setting is fixed now (`--human_ref` not given) rather than tuned later.

## 3. P6 inputs derived from these outputs (rules frozen in P6 / table notes)
- **cap_r,e.** Gene depth × effective length per event, from each event's metaG mapping to the
  pooled reference. It is normalised by prokaryotic single-copy marker depth and computed as
  in `scratch/p6_scorer.py` `extract_run`.
- **n_r,e.** Uniquely assigned panel RNA reads per event, summed over its A/B libraries.
  Read-assignment rules follow table notes §5.2.
- **Unbinned contigs** enter through community annotation, as in frozen P6. **Eukaryotic
  contigs** (Tiara) are excluded from the panel, and their share is reported.

## 4. Pre-registered diagnostics (descriptive only; no change to approved P6-A1)
- **S-O2sulf (Carbon 251096563254214656).** Set O₂ to ½ DL for events with H₂S > 0 and
  N₂O = 0. This targets the unzeroed SBE43 floor on cruise 74 at 165/200 m.
- **S-SO4 coefficient.** Dissimilatory sulfate reduction is 1:1 in S. The approved addendum's
  ½ × ΣH₂S is kept as written, and 1 × ΣH₂S is also reported. The effect is about 0.04%.

## 5. Order and gating
1. Download and verify spot counts.
2. Run the assemblies, the pooled reference and mag_analysis.
3. Run `illumina_rna` up to the per-gene counts (its output, not yet read).
4. Apply the RNA coverage gate once, by the frozen rule.
5. Score the 36-scenario grid with the frozen scorer, extended only by the P6-A1 chemistry
   supply.

Compute is substantial: per-event assembly of about 514 Gbp of metaG, on SLURM under
def-rec3141. Job IDs and resource use are logged in `notes/claude-plexi.md`.
