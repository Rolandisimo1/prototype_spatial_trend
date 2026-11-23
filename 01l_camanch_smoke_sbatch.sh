#!/usr/bin/env bash
#SBATCH --job-name=camanch_smoke
#SBATCH --partition=compute
#SBATCH --qos=long
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --array=1-9
#SBATCH --output=/home/rwkays/isdm/sim_camanch/logs/smoke_%A_%a.log
# -----------------------------------------------------------------------------
# Smoke test for the camera-anchored regional validation (01l_run_camanch_regional.R).
# camanch_smoke_rows.txt: "<scenario> <n_site> <inat_records> <arm>" per line.
# Every level of every axis appears at least once, at the FULL MCMC budget
# (2,000 + 8,000, as in Tier 1) so per-fit runtime can size the 1,080-fit array.
# Rows 2-3 are the decisive check: a camera-only fit (inat_on = 0) of a REAL
# replicate must recover year_region, at 9,000 and at 1,000 cameras.
# 2026-11-23 (Roland): design is the FULL 900-cell grid (7 real ecoregions) with camera and iNat
# both covering the same last 5 years; rows 1-3 / 4-5 also give runtime at 1,000 / 3,000 / 9,000 cameras
# on that range for sizing the 720-fit array. Separate output directory and array from Tier 1.
# -----------------------------------------------------------------------------
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"; export PATH="$DST/bin:$PATH"
cd "$PROJ/prototype_spatial_trend" || exit 2
export SIM_BASE=. N_BURNIN=2000 N_ITER=8000 CAMANCH_SMOKE=1 RANGE_CELLS=900 WINDOW_YEARS=5
export CAMANCH_OUTDIR="$PROJ/prototype_spatial_trend/camanch_smoke_900w5"
echo "node $(hostname) task $SLURM_ARRAY_TASK_ID row $(sed -n "${SLURM_ARRAY_TASK_ID}p" camanch_smoke_rows.txt) start $(date)"
Rscript 01l_run_camanch_regional.R "$SLURM_ARRAY_TASK_ID"; RC=$?; echo "exit $RC end $(date)"; exit $RC
