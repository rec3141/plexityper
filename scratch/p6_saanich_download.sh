#!/bin/bash
# P6-A1 Saanich reads: one SRA run per array task (run list roster_runinfo.tsv d3f170b2f4d4).
#   sbatch --array=1-93%10 p6_saanich_download.sh
#SBATCH --account=def-rec3141
#SBATCH --job-name=plexi-saanich-dl
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=12:00:00
#SBATCH --output=/home/rec3141/scratch/plexityper_saanich/logs/%x-%A_%a.out
set -euo pipefail
W=/home/rec3141/scratch/plexityper_saanich
LINE=$(awk -F'\t' -v n=$((SLURM_ARRAY_TASK_ID+1)) 'NR==n' "$W/roster_runinfo.tsv")
KIND=$(cut -f2 <<<"$LINE" | tr A-Z a-z); RUN=$(cut -f5 <<<"$LINE"); SPOTS=$(cut -f6 <<<"$LINE"); LAYOUT=$(cut -f9 <<<"$LINE")
O=$W/reads/$KIND; mkdir -p "$O" "$W/sra_cache"
[ -s "$O/$RUN.done" ] && { echo "already done $RUN"; exit 0; }
module load sra-toolkit/3.0.9
echo "[INFO] $(date -Is) $(hostname) $RUN $KIND $LAYOUT spots=$SPOTS"
prefetch --max-size 200G -O "$W/sra_cache" "$RUN"
T=${SLURM_TMPDIR:-$W/tmp_$RUN}; mkdir -p "$T"
# prefetch may serve SRA Lite (.sralite: simplified qualities) when full .sra is unavailable
if [[ -s "$W/sra_cache/$RUN/$RUN.sra" ]]; then SRC="$W/sra_cache/$RUN/$RUN.sra"
elif [[ -s "$W/sra_cache/$RUN/$RUN.sralite" ]]; then SRC="$W/sra_cache/$RUN/$RUN.sralite"
else echo "[ERROR] no .sra/.sralite for $RUN" >&2; exit 1; fi
echo "[INFO] source $(basename "$SRC")"
fasterq-dump --split-3 --threads 4 --temp "$T" -O "$T" "$SRC"
n=0
for f in "$T"/${RUN}*.fastq; do
  c=$(( $(wc -l < "$f") / 4 )); echo "[INFO] $(basename $f) reads=$c"
  case "$f" in *_1.fastq|*/${RUN}.fastq) n=$c;; esac
  gzip -c "$f" > "$O/$(basename $f).gz.tmp" && mv "$O/$(basename $f).gz.tmp" "$O/$(basename $f).gz"
done
[ "$n" -eq "$SPOTS" ] || { echo "[ERROR] $RUN read count $n != spots $SPOTS" >&2; exit 1; }
rm -rf "$W/sra_cache/$RUN" "$W/tmp_$RUN"
echo "$RUN $KIND $LAYOUT $n $(date -Is)" > "$O/$RUN.done"
echo "[INFO] $(date -Is) done $RUN"
