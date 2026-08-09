#!/usr/bin/env bash
#SBATCH --job-name=abund_single_test
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/home/rwkays/abund_single_test_%j.log
#SBATCH --error=/home/rwkays/abund_single_test_%j.log
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "start: $(date)"
"$DST/bin/Rscript" 01e_run_abundance_sweep.R 1
echo "exit: $?  end: $(date)"
