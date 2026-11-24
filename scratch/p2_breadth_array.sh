#!/bin/bash
# P2 freshwater breadth, per-BAM covered intervals on the 490 derep-representative MAG contigs.
# Spec: plexityper scratch/p2_breadth_job_spec.md; grid's go in chat 2026-11-20 ("stage it").
# Per-BAM intervals let the event-level union (campaign x station, P2 l.51) be computed later
# by merging intervals across the BAMs of an event; the full-length MAG denominator is applied
# then. Flags match samtools coverage defaults (exclude UNMAP, SECONDARY, QCFAIL, DUP);
# deletions are not counted as covered (bamtobed -splitD).
#   sbatch --array=1-232%40 p2_breadth_array.sh
#SBATCH --account=def-rec3141
#SBATCH --partition=genlm,genoa,genoacpu-b
#SBATCH --job-name=plexi-p2-breadth
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=03:00:00
#SBATCH --output=/home/rec3141/scratch/project_nanopore/logs/%x-%A_%a.out
set -euo pipefail
W=/home/rec3141/scratch/project_nanopore/plexityper_p2_breadth
D=/home/rec3141/scratch/project_nanopore/assembly_freshwater/mapping
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$W/bams.txt")
B=${BAM%.sorted.bam}
SIF=$(readlink -f /home/rec3141/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif)
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
BED=${BED:-targets_v2.bed}; O=$W/out_${BED%.bed}; mkdir -p "$O"
echo "[INFO] $(date -Is) $(hostname) $BAM BED=$BED sha=$(sha256sum "$W/$BED" | cut -c1-12) out=$O"
apptainer exec -B /project/6043541/rec3141/scratch -B /home/rec3141 "$SIF" bash -euo pipefail -c "
source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/dana-mag-binning
export TMPDIR=$W/tmp_\${SLURM_JOB_ID}; mkdir -p \$TMPDIR
samtools view -u -F 0x704 -L $W/$BED $D/$BAM \
  | bedtools bamtobed -splitD -i stdin \
  | cut -f1-3 | sort -S 3G -T \$TMPDIR -k1,1 -k2,2n \
  | bedtools merge -i stdin | gzip -c > $O/$B.covered.bed.gz.tmp
samtools idxstats $D/$BAM > $O/$B.idxstats.tmp
"
mv "$O/$B.covered.bed.gz.tmp" "$O/$B.covered.bed.gz"
mv "$O/$B.idxstats.tmp" "$O/$B.idxstats"
echo "[INFO] $(date -Is) done $B"
