#!/usr/bin/env bash
#SBATCH --job-name=svc_trend_ecoregion
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=40G
#SBATCH --time=02:30:00
#SBATCH --array=1-60%10
#SBATCH --output=/home/rwkays/svc_trend_ecoregion_%A_%a.log
#SBATCH --error=/home/rwkays/svc_trend_ecoregion_%A_%a.log
# 01e_run_ecoregion_sim_sbatch.sh -- the ecoregion-trend simulation study:
# 2 scenarios (varying, null) x 30 reps = 60 tasks, at the CAR prototype's
# "scaletest" density. Each task fits TWO models (ecoregion + national-
# scalar) sequentially in one process -- roughly 2x a single scaletest-
# density fit's cost, hence the generous mem/time vs the single-model CAR
# scaletest sbatch (which used 20G/90min). This is expected to take HOURS
# total (task's own framing) -- do not be alarmed by a long wall-clock;
# monitor via squeue/rep-file counts, not by expecting quick completion.
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "array task: $SLURM_ARRAY_TASK_ID"; echo "node: $(hostname)"; echo "start: $(date)"
"$RS" 01e_run_ecoregion_sim.R
echo "exit: $?  end: $(date)"
