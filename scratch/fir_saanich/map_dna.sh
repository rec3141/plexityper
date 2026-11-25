#!/bin/bash
# array over DNA runs (line N of dna_runs.txt)
#SBATCH --account=def-rec3141
#SBATCH --job-name=plexi-saanich-mapdna
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/home/rec3141/scratch/plexityper_saanich/logs/%x-%A_%a.out
set -euo pipefail; source /home/rec3141/scratch/plexityper_saanich/common.sh
RUN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" $W/dna_runs.txt); D=$W/reads/dna; O=$W/bam_dna; mkdir -p $O
[ -s $D/$RUN.done ] || { echo "[ERROR] $RUN not downloaded"; exit 1; }
T=${SLURM_TMPDIR:-$W/tmp_$RUN}; mkdir -p "$T"
inenv dana-mag-assembly "strobealign -t 14 $W/ref/saanich_pooled.fasta $D/${RUN}_1.fastq.gz $D/${RUN}_2.fastq.gz | samtools sort -@ 4 -m 2G -T $T/$RUN -o $O/$RUN.sorted.bam.tmp - && mv $O/$RUN.sorted.bam.tmp $O/$RUN.sorted.bam && samtools index $O/$RUN.sorted.bam && samtools flagstat $O/$RUN.sorted.bam > $O/$RUN.flagstat"
