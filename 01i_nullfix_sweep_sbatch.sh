#!/usr/bin/env bash
#SBATCH --job-name=nullfix_sweep
#SBATCH --partition=compute
#SBATCH --qos=long
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=1-720
#SBATCH --output=/home/rwkays/nullfix_sweep_%A_%a.log
#SBATCH --error=/home/rwkays/nullfix_sweep_%A_%a.log
# -----------------------------------------------------------------------------
# FULL 720-row re-run with a GENUINELY NULL null scenario.
#
# WHAT CHANGED. 01i_run_estimator_sweep.R now zeroes truth$year_beta and
# truth$year_var under scenario == "null". Previously the null zeroed only the
# spatial deviation field (amplitude = 0), leaving the real national trend in
# place, so tvb_true = year_beta + year_var was -0.1795 in BOTH arms. The null
# arm's detect_rate was power-under-a-spatially-flat-truth, not Type I error,
# and no false-positive rate has ever been measured. This run measures it.
#
# WHY ALL 720 AND NOT JUST THE 360 NULL ROWS. The varying rows are unchanged by
# the fix and could in principle be reused, but row_id -> (estimator,
# abundance, scenario, rep) is assigned by a single expand.grid + sort over the
# whole design. Running a 361-720 subset would depend on that mapping being
# stable across the edit, which is exactly the assumption that produced the
# mislabelled-row hazard documented in 01i_abund_sweep_sbatch.sh. Recompute
# everything; the varying arm then also serves as a regression check that the
# fix changed nothing outside the null.
#
# SEPARATE OUTPUT DIRECTORY. estimator_sweep_out_n30_abund/ holds the corrected
# -ladder run (tag n30_abund, 720/720 OK) and is the current reported result.
# The driver's "already done, exiting" skip keys on files in SWEEP_OUTDIR, so
# writing there would silently adopt the old rows instead of recomputing.
# New directory, no overlap.
#
# Walltime 6h, matching the n30_abund run (which peaked well under it). Report
# if tasks hit the cap rather than raising it silently.
# -----------------------------------------------------------------------------
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
export SIM_BASE=. N_REP=30 N_BURNIN=2000 N_ITER=8000
export SWEEP_OUTDIR="$PROJ/prototype_spatial_trend/estimator_sweep_out_n30_nullfix"
echo "node: $(hostname)  task: $SLURM_ARRAY_TASK_ID  start: $(date)"
"$RS" 01i_run_estimator_sweep.R "$SLURM_ARRAY_TASK_ID"
RC=$?; echo "exit: $RC  end: $(date)"; exit $RC
