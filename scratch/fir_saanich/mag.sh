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
# on-node GPU guard (codex-plexi 251108215089446912): repeat the three env probes on THIS node
for e in binning comebin vamb; do
  d=$(apptainer exec --nv $B "$GSIF" bash -c "export PATH=$DS/mag_analysis/bin:\$PATH; torch_device.sh $e true" 2>&1 | tail -1)
  echo "on-node env $e -> $d"; [ "$d" = cuda ] || { echo "[ERROR] $e not cuda on this node"; exit 1; }
done
export NXF_WORK=$W/mag_work
$DS/mag_analysis/run-mag-analysis.sh --workdir $W/mag_work --container --gpu true \
  --assembly $W/ref/saanich_pooled.fasta --depths $W/ref/depths.txt --bam_dir $W/bam_dna/ \
  --outdir $W/mag --db_dir /home/rec3141/scratch/databases \
  --run_semibin true --run_comebin true --run_vamb true --run_maxbin false --run_lorbin false --run_binette false \
  --run_eukaryotic true --run_gtdbtk true \
  --binning_cpus 32 --binning_memory '400 GB' --bakta_cpus 48 --bakta_memory '400 GB' --medium_cpus 32 --medium_memory '96 GB' \
  2>&1 | tee $W/logs/mag_run_${SLURM_JOB_ID}.txt
# fallback checks: wrapper log AND every binner task's own stderr/log in the work dir
if grep -q -i 'running on CPU\|running on the CPU image\|training on CPU' $W/logs/mag_run_${SLURM_JOB_ID}.txt; then echo "[ERROR] GPU fallback in wrapper log"; exit 1; fi
if grep -r -l -i 'training on CPU' $W/mag_work --include=.command.err --include=.command.log 2>/dev/null | head -5 | grep -q .; then
  echo "[ERROR] binner task trained on CPU:"; grep -r -l -i 'training on CPU' $W/mag_work --include=.command.err --include=.command.log | head; exit 1; fi
echo GPU_GUARD_PASS
