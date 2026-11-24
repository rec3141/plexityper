# P2 freshwater breadth job: spec, pending grid's go (relay 249719437766606848)

Nothing runs until grid approves and Grex is logged in on concentration.

## Steps
1. **Check the depths header first.** Read the header of `mapping/depths.txt` once, through
   concentration's ControlMaster. If it holds only per-contig mean depth and variance, it
   cannot give breadth, and the job below is needed.
2. **Per event, take the union of covered positions** (Astra 249719554330509312, P2 l.51).
   - An event is **cruise/campaign × station**, as in frozen P2 l.51. Same-cast depths and
     technical/sequencing replicates are grouped into one event; depths are not split. The
     BAM → event map follows this definition (Astra 249719710081794048).
   - Where several BAMs in an event map to the same reference, take the union of their
     covered positions before counting bases, e.g. `samtools merge` per event, then
     `samtools coverage` or `samtools depth -a`, counting positions with depth ≥ 1.
   - Per-BAM counts alone would lose the overlap between BAMs.
3. **MAG breadth = union covered bases / full MAG length,** with zero-coverage contigs
   included in the length.
4. **Keep the ID maps:** BAM → event ID, and contig → MAG ID.
5. **Sensitivity run:** repeat with prokaryotic-mapped reads equalised across events
   (`samtools view -s` subsampling to the minimum event depth). Include it in the job
   estimate.
6. **Where it runs:** the danaSeq container (dana-mag-assembly env, which has samtools), as
   one small sbatch array, one task per event (campaign × station). Outputs come back as TSVs with hashes added
   to `data/grex_mag_freshwater/PROVENANCE.md`.

## Status 2026-11-20
- Primary pass (per-BAM merged intervals) is scratch/p2_breadth_array.sh 6698e9da8f46; test job 8194718.
- The equalized-depth sensitivity needs a **second pass over the BAMs** once event depths are known, because the merged intervals discard read identity and multiplicity (Astra 249753353072328704). The BAMs stay in place on Grex; this pass is tracked separately.
