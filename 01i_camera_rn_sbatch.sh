#!/usr/bin/env bash
#SBATCH --job-name=cam_rn
#SBATCH --partition=compute
#SBATCH --qos=long
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=541-720
#SBATCH --output=/home/rwkays/cam_rn_%A_%a.log
#SBATCH --error=/home/rwkays/cam_rn_%A_%a.log
# -----------------------------------------------------------------------------
# Fourth arm only: camera-level Royle-Nichols, rows 541-720 of the EXPANDED
# 720-row design. Rows 1-540 are the already-completed camera_occ/array_occ/
# array_rn tasks -- verified unchanged by the expansion, because "camera_rn"
# sorts LAST alphabetically among the four estimator labels, so every existing
# row_id still denotes exactly the same cell. Do NOT re-run 1-540.
#   Re-verified independently: 4-arm sort gives array_occ 1-180, array_rn
#   181-360, camera_occ 361-540, camera_rn 541-720, and the (estimator,
#   abundance, scenario, rep_id) tuple for rows 1-540 is identical under the
#   3-arm and 4-arm designs.
#
# Walltime is 6h (vs the ~0.8h the other arms actually took): camera-level RN
# carries ~700 latent discrete N_cam nodes against array_rn's ~97, so mixing is
# expected to be slower. Report if tasks hit the cap; do not silently raise.
#
# Write to the SAME output directory as the n30 sweep so the collector picks up
# all 720 in one pass.
#
# PREAMBLE: copied from run_full_sweep_array.sh (job 650207, the n30 sweep)
# rather than hand-rolled. The earlier draft of this file used
# `source ~/.bashrc; conda activate nimble_env`, which this project has
# established does not work: compileNimble() shells out to `R CMD SHLIB`, so
# bare `R` must be on PATH, and the export-PATH form is what every working
# sbatch wrapper here uses. It also omitted SWEEP_OUTDIR (defaulting to the
# 3-rep PILOT directory, not the n30 directory) and requested --qos=short's
# territory without naming a QOS -- `short` caps at 02:00:00, so 6h needs
# --partition=compute --qos=long, which is what the n30 array ran under.
# -----------------------------------------------------------------------------
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
export SIM_BASE=. N_REP=30 N_BURNIN=2000 N_ITER=8000
export SWEEP_OUTDIR="$PROJ/prototype_spatial_trend/estimator_sweep_out_n30"
echo "node: $(hostname)  task: $SLURM_ARRAY_TASK_ID  start: $(date)"
"$RS" 01i_run_estimator_sweep.R "$SLURM_ARRAY_TASK_ID"
RC=$?; echo "exit: $RC  end: $(date)"; exit $RC
