#!/usr/bin/env bash
#SBATCH --job-name=svc_trend_simval
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --array=1-50%10
#SBATCH --output=/home/rwkays/svc_trend_simval_%A_%a.log
#SBATCH --error=/home/rwkays/svc_trend_simval_%A_%a.log
# 01_run_sim_validation_sbatch.sh -- the real 50-replicate simulation
# validation (STEP 1 of the spatial-trend prototype), as a Slurm ARRAY: one
# fresh process per replicate. See 01_run_sim_validation.R header for why
# this replaced two earlier in-script-parallel attempts (mclapply fork
# collision, then PSOCK worker memory accumulation across repeated NIMBLE
# compiles). Sized from the real-scale single-process dry run (job 435980:
# 12.5 min wall, 5.5GB peak) -- 12G/task with headroom, %10 throttle means
# 50 reps / 10 concurrent =~ 5 batches * 12.5 min =~ 65 min wall to finish
# the whole array, well inside the 2h --qos=short ceiling per-task.
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "array task: $SLURM_ARRAY_TASK_ID"; echo "node: $(hostname)"; echo "start: $(date)"
"$RS" 01_run_sim_validation.R
echo "exit: $?  end: $(date)"
