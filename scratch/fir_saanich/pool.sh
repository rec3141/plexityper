#!/bin/bash
#SBATCH --account=def-rec3141
#SBATCH --job-name=plexi-saanich-pool
#SBATCH --cpus-per-task=16
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/home/rec3141/scratch/plexityper_saanich/logs/%x-%j.out
set -euo pipefail; source /home/rec3141/scratch/plexityper_saanich/common.sh
C=$W/mgnify_contigs; R=$W/ref; mkdir -p $R; cd $R
# reconcile exactly the 39 crosswalk files
exp=$(tail -n +2 $W/mgnify_contig_urls.tsv | awk -F'\t' '{print $2"__"$1"__"$3".fasta.gz"}' | sort)
got=$(cd $C && ls *.fasta.gz | sort)
[ "$exp" = "$got" ] || { echo "[ERROR] contig files != crosswalk"; diff <(echo "$exp") <(echo "$got"); exit 1; }
[ $(echo "$exp" | wc -l) -eq 39 ] || { echo "[ERROR] not 39"; exit 1; }
for f in $exp; do [ -s $C/$f ] && gzip -t $C/$f || { echo "[ERROR] bad $f"; exit 1; }; done
rm -f combined.fasta
for f in $exp; do a=${f%%__*}; zcat $C/$f | awk -v a="$a" '/^>/{sub(/^>/,">"a"__"); print $1; next}{print}' >> combined.fasta; done
inenv dana-mag-assembly "dedupe.sh -Xmx80g in=combined.fasta out=d100.fasta.gz sort=length uniquenames=t minidentity=100 t=16 ow=t && \
 dedupe.sh -Xmx80g in=d100.fasta.gz out=d99.fasta.gz sort=length uniquenames=t minidentity=99 t=16 ow=t && \
 dedupe.sh -Xmx80g in=d99.fasta.gz out=d98.fasta.gz sort=length uniquenames=t minidentity=98 t=16 ow=t && \
 reformat.sh in=d98.fasta.gz out=saanich_pooled.fasta minlength=500 ow=t && samtools faidx saanich_pooled.fasta"
grep -c '^>' combined.fasta saanich_pooled.fasta
sha256sum saanich_pooled.fasta > saanich_pooled.fasta.sha256; rm -f combined.fasta d100.fasta.gz d99.fasta.gz
