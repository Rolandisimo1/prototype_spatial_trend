#!/usr/bin/env bash
#SBATCH --job-name=eco_smoke
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=00:45:00
#SBATCH --output=/home/rwkays/eco_smoke_%j.log
#SBATCH --error=/home/rwkays/eco_smoke_%j.log
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
cd "$PROJ/prototype_spatial_trend" || exit 2
"$DST/bin/Rscript" ecoregion_smoke_test.R
echo "exit: $?"
