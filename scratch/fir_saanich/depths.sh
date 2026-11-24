#!/bin/bash
#SBATCH --account=def-rec3141
#SBATCH --job-name=plexi-saanich-depths
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --output=/home/rec3141/scratch/plexityper_saanich/logs/%x-%j.out
set -euo pipefail; source /home/rec3141/scratch/plexityper_saanich/common.sh
n=$(ls $W/bam_dna/*.sorted.bam | wc -l); [ $n -eq 39 ] || { echo "[ERROR] $n BAMs"; exit 1; }
inenv dana-mag-binning "jgi_summarize_bam_contig_depths --outputDepth $W/ref/depths.txt $W/bam_dna/*.sorted.bam"
