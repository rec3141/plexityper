W=/home/rec3141/scratch/plexityper_saanich
DS=/home/rec3141/GENICE/danaSeq
SIF=$(readlink -f $DS/mag_analysis/.danaseq-mag-analysis.sif)
GSIF=$(readlink -f $DS/mag_analysis/.danaseq-mag-analysis-gpu.sif)
module load apptainer/1.4.5 >/dev/null 2>&1 || export PATH=/opt/singularity-ce-4.4.1/bin:$PATH
B="-B /home/rec3141 -B /scratch/rec3141 -B /project/6043541"
inenv() { local env=$1; shift; apptainer exec $B "$SIF" bash -euo pipefail -c "source /opt/conda/etc/profile.d/conda.sh && conda activate /opt/conda/envs/$env && $*"; }
