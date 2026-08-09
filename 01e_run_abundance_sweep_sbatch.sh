#!/usr/bin/env bash
#SBATCH --job-name=eco_abund_sweep
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-180%10
#SBATCH --output=/home/rwkays/eco_abund_sweep_%A_%a.log
#SBATCH --error=/home/rwkays/eco_abund_sweep_%A_%a.log
# 01e_run_abundance_sweep_sbatch.sh -- 3 abundance levels x 2 trend scenarios
# x 30 reps = 180 tasks, reduced-data density (per-rep cost is a few thousand
# iterations, well under an hour -- 8G/2h per task has headroom).
#
# TWO FIXES vs the originally-drafted version of this script:
#   1. `conda activate nimble_env` inside a non-interactive sbatch script
#      relies on ~/.bashrc correctly sourcing conda's shell hook, which is
#      NOT guaranteed and is exactly the failure mode already hit and
#      debugged earlier in this project (compileNimble failed with "Failed
#      to create the shared library" because nimble_env's conda-provided
#      C++ compiler, x86_64-conda-linux-gnu-c++, wasn't on PATH -- jobs
#      435927 -> 435928). Exporting nimble_env/bin onto PATH directly and
#      invoking Rscript by full path sidesteps the whole conda-activation
#      question, matching every other sbatch script in this project.
#   2. Added a %10 concurrency throttle. The original array had none
#      (--array=1-180, all 180 tasks launched at once) -- unthrottled on a
#      shared partition isn't something to do without discussing cluster
#      impact first; %10 matches the throttle used throughout this project.
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "array task: $SLURM_ARRAY_TASK_ID"; echo "node: $(hostname)"; echo "start: $(date)"
"$RS" 01e_run_abundance_sweep.R
echo "exit: $?  end: $(date)"
