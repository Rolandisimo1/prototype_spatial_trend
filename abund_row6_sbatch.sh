#!/usr/bin/env bash
#SBATCH --job-name=abund_row6
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/home/rwkays/abund_row6_%j.log
#SBATCH --error=/home/rwkays/abund_row6_%j.log
# Backfill for array task 6 (bobcat_baseline/null/rep=6), which TIMEOUT'd in
# job 451229 (2h limit) while everything else completed cleanly.
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "start: $(date)"
"$DST/bin/Rscript" 01e_run_abundance_sweep.R 6
echo "exit: $?  end: $(date)"
