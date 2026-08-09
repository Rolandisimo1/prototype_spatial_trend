#!/usr/bin/env bash
#SBATCH --job-name=svc_trend_grain3sup
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --array=1-15%10
#SBATCH --output=/home/rwkays/svc_trend_grain3sup_%A_%a.log
#SBATCH --error=/home/rwkays/svc_trend_grain3sup_%A_%a.log
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "array task: $SLURM_ARRAY_TASK_ID"; echo "node: $(hostname)"; echo "start: $(date)"
"$RS" 01d_run_grain3_supplement.R
echo "exit: $?  end: $(date)"
