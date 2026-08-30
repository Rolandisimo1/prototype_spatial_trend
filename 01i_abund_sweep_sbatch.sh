#!/usr/bin/env bash
#SBATCH --job-name=abund_sweep
#SBATCH --partition=compute
#SBATCH --qos=long
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=1-720
#SBATCH --output=/home/rwkays/abund_sweep_%A_%a.log
#SBATCH --error=/home/rwkays/abund_sweep_%A_%a.log
# -----------------------------------------------------------------------------
# FULL 720-row re-run with the abundance ladder actually applied.
#
# WHY ALL 720 AND NOT A PATCH. The design levels were renamed to
# abundance_levels_measured()$level verbatim (bobcat_like / intermediate /
# deer_like). That changes the alphabetical sort, so row_id denotes a
# different cell than it did in the inert-ladder run -- row-by-row reuse
# would silently mislabel results. Every row is recomputed.
#
# SEPARATE OUTPUT DIRECTORY. estimator_sweep_out_n30/ holds the inert-ladder
# run (tagged n30_4arm), which is still a valid n=90-per-arm single-abundance
# result and is already in use. It must not be mixed with or overwritten by
# this run -- and because row_ids now mean different things, writing here
# would ALSO trip the driver's "already done, exiting" skip and silently
# adopt 720 mislabeled files. New directory, no overlap.
#
# Walltime stays 6h. The inert-ladder run peaked at 0.83h, but every cell in
# that run was effectively bobcat abundance; deer_like carries ~4.5x the
# intensity and more latent N_cam mass, so it is genuinely untested for
# mixing speed. Report if tasks hit the cap rather than raising it silently.
# -----------------------------------------------------------------------------
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
export SIM_BASE=. N_REP=30 N_BURNIN=2000 N_ITER=8000
export SWEEP_OUTDIR="$PROJ/prototype_spatial_trend/estimator_sweep_out_n30_abund"
echo "node: $(hostname)  task: $SLURM_ARRAY_TASK_ID  start: $(date)"
"$RS" 01i_run_estimator_sweep.R "$SLURM_ARRAY_TASK_ID"
RC=$?; echo "exit: $RC  end: $(date)"; exit $RC
