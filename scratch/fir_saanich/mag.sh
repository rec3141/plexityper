#!/bin/bash
#SBATCH --account=def-rec3141
#SBATCH --job-name=plexi-saanich-mag
#SBATCH --gpus-per-node=h100:1
#SBATCH --cpus-per-task=48
#SBATCH --mem=600G
#SBATCH --time=3-00:00:00
#SBATCH --output=/home/rec3141/scratch/plexityper_saanich/logs/%x-%j.out
set -euo pipefail; source /home/rec3141/scratch/plexityper_saanich/common.sh
nvidia-smi -L
$DS/mag_analysis/run-mag-analysis.sh --container --gpu true \
  --assembly $W/ref/saanich_pooled.fasta --depths $W/ref/depths.txt --bam_dir $W/bam_dna/ \
  --outdir $W/mag --db_dir /home/rec3141/scratch/databases \
  --run_semibin true --run_comebin true --run_vamb true --run_maxbin false --run_lorbin false --run_binette false \
  --run_eukaryotic true --run_gtdbtk true \
  --binning_cpus 32 --binning_memory '400 GB' --bakta_cpus 48 --bakta_memory '400 GB' --medium_cpus 32 --medium_memory '96 GB' \
  2>&1 | tee $W/logs/mag_run_${SLURM_JOB_ID}.txt
if grep -q -i 'running on CPU\|running on the CPU image' $W/logs/mag_run_${SLURM_JOB_ID}.txt; then echo "[ERROR] GPU fallback detected"; exit 1; fi
