#!/usr/bin/env bash
#SBATCH --partition=compute
#SBATCH --qos=long
#SBATCH --ntasks=1
#SBATCH --output=/home/rwkays/isdm/sim_tier1/logs/tier1_%x_%A_%a.log
# -----------------------------------------------------------------------------
# Tier 1 full run (01k_run_tier1_volume.R). Each array task runs ROWS_PER_TASK
# consecutive design rows in sequence, because QOS long caps a user at 1000
# submitted jobs and the two blocks hold 4,320 rows:
#   sbatch --job-name=t1_RW --array=<tasks> --mem=<M> --time=<T> \
#     --export=ALL,TIER1_BLOCK=RW,ROWS_PER_TASK=8,ROW_MAX=2400 01k_tier1_sbatch.sh
# Task t runs rows (t-1)*ROWS_PER_TASK+1 .. t*ROWS_PER_TASK (capped at ROW_MAX).
# Memory and walltime are set per submission from the smoke test (job
# 1389150): rows are ordered by range, window, camera count, iNat volume,
# effect, rep (build_tier1_design()), so a task range maps to contiguous design
# cells. The driver skips rows already written, so resubmitting a task resumes it.
# Output: tier1_out_<BLOCK>/row_%05d.rds; collect with 02k_collect_tier1.R.
# -----------------------------------------------------------------------------
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"; export PATH="$DST/bin:$PATH"
cd "$PROJ/prototype_spatial_trend" || exit 2
export SIM_BASE=. N_REP=30 N_BURNIN=2000 N_ITER=8000
export TIER1_OUTDIR="$PROJ/prototype_spatial_trend/tier1_out_${TIER1_BLOCK}"
K=${ROWS_PER_TASK:-1}
FIRST=$(( (SLURM_ARRAY_TASK_ID - 1) * K + 1 )); LAST=$(( SLURM_ARRAY_TASK_ID * K ))
(( LAST > ROW_MAX )) && LAST=$ROW_MAX
echo "node $(hostname) block $TIER1_BLOCK task $SLURM_ARRAY_TASK_ID rows $FIRST-$LAST start $(date)"
FAIL=0
for ROW in $(seq "$FIRST" "$LAST"); do
  Rscript 01k_run_tier1_volume.R "$ROW"; RC=$?
  echo "row $ROW exit $RC $(date)"; (( RC != 0 )) && FAIL=1
done
echo "end $(date)"; exit $FAIL
