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
1. **Per-event metaG contigs from MGnify (grid's go in chat, 2026-11-24: "go ahead with the mgnify
   contigs on fir").** Instead of running `illumina_assembly`, the processed contigs come from
   MGnify MGYS00005865 (pipeline v5.0). There are 39 assemblies for the 38 roster DNA BioSamples
   (crosswalk `mgnify_dna_crosswalk.tsv`, reviewed by codex-plexi 251104688602529792).
   - Each assembly is built from exactly one roster DNA run, and no assembly is shared.
   - SAMN05224393 has two assemblies, one per run.
   - No RNA run is a source.
   - Only the contig FASTAs are used. **MGnify's CDS, eggNOG/KO and other annotations are not
     used.**
   - Download: `scratch/p6_saanich_mgnify_contigs.sh`, with sha256 manifest `SHA256SUMS`.
2. **Pooled reference.** Concatenate the 39 contig sets, prefixing contig IDs with the assembly
   accession, and dereplicate at the contig level with BBTools dedupe (defaults; danaSeq
   container). This gives `saanich_pooled.fasta`.
   - A pooled, not per-event, reference keeps one gene catalogue across events, so capacity and
     RNA counts are measured on the same features.
   - **DNA mapping for binning.** All 39 DNA runs go through strobealign to the pooled
     reference. The read-order output is teed to `p6_field_count.py` (§2.4) and also
     coordinate-sorted for `jgi_summarize_bam_contig_depths` (danaSeq container). That gives
     the `--depths` file and `--bam_dir` input for mag_analysis.
3. **Annotation and binning.** `mag_analysis --assembly saanich_pooled.fasta --depths <jgi depths>
   --bam_dir <sorted DNA BAMs>`, with default settings.
   - Bakta annotation with coordinate-joined KOs (trap 3).
   - Binning (DAS Tool set), CheckM2 and GTDB-Tk.
   - KO calls follow the pipeline's merged annotations. The frozen table's
     `make_binned_ko_classifier` disambiguation then uses the MAG membership.
4. **DNA and RNA mapping and counting: one path for both (revised after codex-plexi
   251098919257686016).**
   - **Why not `illumina_rna`.** It is not used for counting, because it skips single-end
     libraries. The RNA layouts are 37 SE and 17 PE runs; the 39 DNA runs are all PE.
   - **Alignment.** Each run is aligned to `saanich_pooled.fasta` with strobealign from the
     danaSeq container, in SE or PE mode per run. Output stays in read order, so it is
     name-grouped.
   - **RNA rRNA removal.** RNA reads first go through bbduk with SILVA SSU+LSU (k=31). This is
     the same rRNA step as `illumina_rna` QC.
   - **Counting.** `scratch/p6_field_count.py` 4949b367168c (bounded name-grouping check), with fixtures
     `scratch/test_p6_field_count.py` a9189483ff99 (11/11 OK). The same rule applies to DNA and
     RNA:
     - drop flags 0x4, 0x100 and 0x800;
     - one count per SE read or PE fragment;
     - mates' aligned blocks are unioned, so a base covered by both mates counts once;
     - the fragment goes to the CDS with strictly the largest aligned-base overlap, and exact
       ties are dropped;
     - the output gives both fragments and overlap_bp per CDS;
     - MAPQ ≥ 20 is a sensitivity.
   - **Event totals.** All 54 RNA runs are summed into their 38 event totals (A/B libraries
     together); the DNA runs are summed the same way.

## 3. P6 inputs derived from these outputs (rules frozen in P6 / table notes)
- **Field extractor.** `scratch/p6_field_extract.py`, with fixtures
  `scratch/test_p6_field_extract.py` (6/6 OK). There is also an end-to-end synthetic test, `scratch/test_p6_field_pipeline.py`. It covers SE+PE SAM → counter → extractor over all 36 scenarios, with zero/nondetect chemistry, nonuniform marker depths, eukaryotic and ambiguous-unbinned exclusions, the anammox exclusions, exact A+B summation, and Cohort/analyse. It replaces the synthetic `extract_run`, which
  hard-codes Tier C plus H2.
  - **Panel.** Core + Tier B_N2O, 11 rows; h2_gated rows are off.
  - **cap_r,e.** Σ overlap_bp over the CDSs of r, divided by d_prok,e. d_prok,e is the true median
    over the 14 marker families of Σ overlap_bp / Σ CDS length, over prokaryotic marker CDSs
    (table notes §2.2).
  - **Binned CDSs.** The frozen classifier with the MAG's KO set. The Tier C nirK/norB
    exclusions apply to DEN_NIRNOR only. DEN_NOS (K00376) follows its frozen row, which has no
    anammox exclusion.
  - **Unbinned CDSs.** unambiguous_kos only.
  - **Eukaryotic CDSs** are excluded.
  - **N2.** Weiss (1970) air saturation.
  - **Chemistry.** The P6-A1 supply and the 36-scenario grid.
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
3. Align and count DNA and RNA per §2.4 (strobealign → `p6_field_count.py`), producing per-gene counts. The counts are not read yet.
4. Apply the RNA coverage gate once, by the frozen rule.
5. Score the 36-scenario grid with the frozen scorer, extended only by the P6-A1 chemistry
   supply.

Compute: no de novo assembly (MGnify contigs). The costs are the dedupe, mapping of about 1.1 Tbp of reads, and mag_analysis, on SLURM under def-rec3141. Job IDs and resource use are logged in `notes/claude-plexi.md`.
