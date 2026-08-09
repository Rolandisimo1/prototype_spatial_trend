#!/usr/bin/env bash
#SBATCH --job-name=eco_effort_secondary
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-30%10
#SBATCH --output=/home/rwkays/eco_effort_secondary_%A_%a.log
#SBATCH --error=/home/rwkays/eco_effort_secondary_%A_%a.log
# 01h_run_effort_secondary_sbatch.sh -- 30 reps, same model size as the
# effort sweep's 1x-camera cells (700 sites), so the same 2h budget applies.
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "array task: $SLURM_ARRAY_TASK_ID"; echo "node: $(hostname)"; echo "start: $(date)"
"$RS" 01h_run_effort_secondary.R
echo "exit: $?  end: $(date)"
