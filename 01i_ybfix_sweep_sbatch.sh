#!/usr/bin/env bash
#SBATCH --job-name=ybfix_sweep
#SBATCH --partition=compute
#SBATCH --qos=long
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --array=1-180
#SBATCH --output=/home/rwkays/isdm/sim_tier1/logs/ybfix_%A_%a.log
#SBATCH --error=/home/rwkays/isdm/sim_tier1/logs/ybfix_%A_%a.log
#SBATCH --nice=500
# -----------------------------------------------------------------------------
# FULL 720-row re-run of the 4-arm design (tag n30_ybfix), identical to the
# n30_nullfix run except for the metrics. sim_helpers_estimator_metrics.R now
# scores year_beta (yb_*: truth, mean, 95% bounds, bias, width, coverage,
# detection) alongside total_var_beta (tvb_*). No earlier run kept year_beta
# bounds or a separate year_beta truth, so the estimator-choice claim has
# never been validated on the reported parameter; this run does that.
# Driver, seeds and design are unchanged, so the tvb_* columns should reproduce
# n30_nullfix row for row -- that is the regression check.
# Each task runs 4 consecutive rows (720 rows / 180 tasks): QOS long caps a
# user at 1000 submitted jobs, shared with Tier 1 and the windowed refits.
# --nice=500: second priority behind Tier 1 (01k); it takes leftover slots.
# New output directory: the driver's "already done" skip would otherwise adopt
# the n30_nullfix rows, which have no yb_* columns.
# -----------------------------------------------------------------------------
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
export SIM_BASE=. N_REP=30 N_BURNIN=2000 N_ITER=8000
export SWEEP_OUTDIR="$PROJ/prototype_spatial_trend/estimator_sweep_out_n30_ybfix"
echo "node: $(hostname)  task: $SLURM_ARRAY_TASK_ID  start: $(date)"
FAIL=0
for ROW in $(seq $(( (SLURM_ARRAY_TASK_ID - 1) * 4 + 1 )) $(( SLURM_ARRAY_TASK_ID * 4 ))); do
  "$RS" 01i_run_estimator_sweep.R "$ROW"; RC=$?
  echo "row $ROW exit: $RC  $(date)"; (( RC != 0 )) && FAIL=1
done
echo "end: $(date)"; exit $FAIL
