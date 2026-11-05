#!/usr/bin/env bash
#SBATCH --job-name=tier1_smoke
#SBATCH --partition=compute
#SBATCH --qos=long
#SBATCH --ntasks=1
#SBATCH --mem=48G
#SBATCH --time=08:00:00
#SBATCH --array=1-22
#SBATCH --output=/home/rwkays/isdm/sim_tier1/logs/smoke_%A_%a.log
# -----------------------------------------------------------------------------
# Tier 1 smoke test: one short-MCMC replicate at every level of every axis
# (tier1_smoke_rows.txt: "<BLOCK> <row_id> <N_ITER>" per line). Confirms each
# axis is LIVE in the realized design and data, that the yb_* columns are
# populated and that yb_true moves with the effect level, and measures time
# and memory per design cell before the full launch. Rows with two N_ITER
# values separate fixed build/compile time from per-iteration cost.
# -----------------------------------------------------------------------------
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"; export PATH="$DST/bin:$PATH"
cd "$PROJ/prototype_spatial_trend" || exit 2
read -r BLK ROW ITER < <(sed -n "${SLURM_ARRAY_TASK_ID}p" tier1_smoke_rows.txt)
export SIM_BASE=. N_REP=30 N_BURNIN=100 N_ITER=$ITER TIER1_BLOCK=$BLK
export TIER1_OUTDIR="$PROJ/prototype_spatial_trend/tier1_smoke_it${ITER}_${BLK}"
echo "node $(hostname) task $SLURM_ARRAY_TASK_ID block $BLK row $ROW iter $ITER start $(date)"
Rscript 01k_run_tier1_volume.R "$ROW"; RC=$?; echo "exit $RC end $(date)"; exit $RC
