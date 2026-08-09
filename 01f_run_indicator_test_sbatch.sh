#!/usr/bin/env bash
#SBATCH --job-name=eco_indicator_test
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-180%10
#SBATCH --output=/home/rwkays/eco_indicator_test_%A_%a.log
#SBATCH --error=/home/rwkays/eco_indicator_test_%A_%a.log
# 01f_run_indicator_test_sbatch.sh -- same 180-task grid as
# 01e_run_abundance_sweep_sbatch.sh (3 abundance x 2 scenario x 30 reps),
# reusing that run's stored sim_data rather than regenerating it. Each task
# fits ONE model (model_code_ecoregion_switch) but runs 2 MCMC chains
# (Gelman-Rubin needs >=2 chains -- see fit_replicate_switch() in
# sim_helpers.R) instead of the abundance sweep's 2 separate single-chain
# model fits -- time budget carried over unchanged from 01e pending the
# smoke test's actual per-task timing; adjust if the smoke test shows it's
# not enough.
#
# Same PATH-export pattern as every other sbatch script in this project
# (conda activate inside non-interactive sbatch is unreliable -- see
# 01e_run_abundance_sweep_sbatch.sh's header note), and the same %10
# concurrency throttle.
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "array task: $SLURM_ARRAY_TASK_ID"; echo "node: $(hostname)"; echo "start: $(date)"
"$RS" 01f_run_indicator_test.R
echo "exit: $?  end: $(date)"
