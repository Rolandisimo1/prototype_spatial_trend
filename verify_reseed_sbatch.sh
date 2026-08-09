#!/usr/bin/env bash
#SBATCH --job-name=svc_trend_verify_reseed
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --time=00:10:00
#SBATCH --output=/home/rwkays/svc_trend_verify_reseed_%j.log
#SBATCH --error=/home/rwkays/svc_trend_verify_reseed_%j.log
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
cd "$PROJ/prototype_spatial_trend" || exit 2
"$DST/bin/Rscript" verify_reseed.R
echo "exit: $?"
