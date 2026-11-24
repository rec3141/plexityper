#!/bin/bash
#SBATCH --account=def-rec3141
#SBATCH --job-name=plexi-saanich-gpucheck
#SBATCH --gpus-per-node=h100:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/home/rec3141/scratch/plexityper_saanich/logs/%x-%j.out
set -uo pipefail; source /home/rec3141/scratch/plexityper_saanich/common.sh
nvidia-smi -L
ok=1
for e in binning comebin vamb; do
  d=$(apptainer exec --nv $B "$GSIF" bash -c "export PATH=$DS/mag_analysis/bin:\$PATH; torch_device.sh $e true" 2>&1 | tail -1)
  echo "env $e -> $d"; [ "$d" = cuda ] || ok=0
done
[ $ok -eq 1 ] && echo GPU_CHECK_PASS || { echo GPU_CHECK_FAIL; exit 1; }
