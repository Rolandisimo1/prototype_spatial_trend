#!/usr/bin/env bash
#SBATCH --job-name=svc_trend_scaletest
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=01:30:00
#SBATCH --array=1-20%10
#SBATCH --output=/home/rwkays/svc_trend_scaletest_%A_%a.log
#SBATCH --error=/home/rwkays/svc_trend_scaletest_%A_%a.log
# 01b_run_sim_validation_scaletest_sbatch.sh -- confirmatory follow-up array:
# 3x the retained cell50/site counts of the original validated design (see
# 01b_run_sim_validation_scaletest.R header), 20 replicates (not 50 -- this
# is a diagnostic check of whether bias shrinks with more informed cells,
# not the final validation). Per-task mem/time bumped from the original
# array's 12G/30min since the design is ~3x larger; each array task gets
# its OWN independent Slurm allocation (unlike the earlier failed PSOCK
# attempt), so there is no shared-memory risk across concurrent tasks.
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "array task: $SLURM_ARRAY_TASK_ID"; echo "node: $(hostname)"; echo "start: $(date)"
"$RS" 01b_run_sim_validation_scaletest.R
echo "exit: $?  end: $(date)"
