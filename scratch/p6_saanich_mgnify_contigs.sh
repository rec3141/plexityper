#!/bin/bash
# Download the 39 MGnify v5 DNA assemblies (MGYS00005865) for the P6-A1 Saanich roster.
#   sbatch p6_saanich_mgnify_contigs.sh
#SBATCH --account=def-rec3141
#SBATCH --job-name=plexi-saanich-mgnify
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=06:00:00
#SBATCH --output=/home/rec3141/scratch/plexityper_saanich/logs/%x-%j.out
set -euo pipefail
W=/home/rec3141/scratch/plexityper_saanich; O=$W/mgnify_contigs; mkdir -p "$O"
tail -n +2 "$W/mgnify_contig_urls.tsv" | while IFS=$'\t' read -r bs asm run url; do
  f="$O/${asm}__${bs}__${run}.fasta.gz"
  [ -s "$f" ] && gzip -t "$f" && { echo "have $asm"; continue; }
  curl -sS --fail --retry 3 --retry-delay 20 -o "$f.tmp" "$url"
  gzip -t "$f.tmp" && mv "$f.tmp" "$f"
  echo "[INFO] $(date -Is) $asm $(stat -c %s "$f") bytes $(zcat "$f" | grep -c '^>') contigs"
done
ls "$O"/*.fasta.gz | wc -l
( cd "$O" && sha256sum *.fasta.gz ) > "$O/SHA256SUMS"
