#!/usr/bin/env bash
#SBATCH --job-name=svc_trend_grain
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --array=1-45%10
#SBATCH --output=/home/rwkays/svc_trend_grain_%A_%a.log
#SBATCH --error=/home/rwkays/svc_trend_grain_%A_%a.log
# 01c_run_grain_sweep_sbatch.sh -- Phase 1 grain sweep: 3 new grain levels x
# 15 reps = 45 tasks, at the ORIGINAL design's density (210 cell50/700 sites,
# same 12G/30min-per-task sizing that worked cleanly for that design in job
# 439834 -- this sweep reuses the identical density, only the truth field's
# spatial grain differs).
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "array task: $SLURM_ARRAY_TASK_ID"; echo "node: $(hostname)"; echo "start: $(date)"
"$RS" 01c_run_grain_sweep.R
echo "exit: $?  end: $(date)"
