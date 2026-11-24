#!/bin/bash
# Streaming variant of p2_breadth_array.sh for BAMs too large for the sort-based pipeline
# (task 8194877_2 / 8195114_2 OOM at 8G and 32G). Same targets BED, same read filter
# (samtools depth default excludes UNMAP, SECONDARY, QCFAIL, DUP = -F 0x704), deletions not
# counted as covered (depth default without -J). Positions with depth >= 1 are merged into
# half-open BED intervals by awk; no sort needed because depth output is coordinate-ordered.
# Validation: run on task 1 as well and compare with the bamtobed-based v2 output.
#   sbatch --export=ALL,BED=targets_v2.bed --array=1,2 p2_breadth_depth.sh
#SBATCH --account=def-rec3141
#SBATCH --partition=genlm,genoa,genoacpu-b
#SBATCH --job-name=plexi-p2-breadth-depth
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=06:00:00
#SBATCH --output=/home/rec3141/scratch/project_nanopore/logs/%x-%A_%a.out
set -euo pipefail
W=/home/rec3141/scratch/project_nanopore/plexityper_p2_breadth
D=/home/rec3141/scratch/project_nanopore/assembly_freshwater/mapping
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$W/bams.txt")
B=${BAM%.sorted.bam}
SIF=$(readlink -f /home/rec3141/GENICE/danaSeq/mag_analysis/.danaseq-mag-analysis.sif)
export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
BED=${BED:-targets_v2.bed}; O=$W/out_${BED%.bed}_depth; mkdir -p "$O"
echo "[INFO] $(date -Is) $(hostname) $BAM BED=$BED sha=$(sha256sum "$W/$BED" | cut -c1-12) out=$O method=depth"
apptainer exec -B /project/6043541/rec3141/scratch -B /home/rec3141 "$SIF" bash -euo pipefail -c "
source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/dana-mag-binning
samtools depth -b $W/$BED $D/$BAM \
  | awk -F'\t' -v OFS='\t' '
      \$3 > 0 { p = \$2 - 1
                if (\$1 == c && p == e) { e = p + 1 }
                else { if (c != \"\") print c, s, e; c = \$1; s = p; e = p + 1 } }
      END { if (c != \"\") print c, s, e }' \
  | gzip -c > $O/$B.covered.bed.gz.tmp
"
mv "$O/$B.covered.bed.gz.tmp" "$O/$B.covered.bed.gz"
echo "[INFO] $(date -Is) done $B"
