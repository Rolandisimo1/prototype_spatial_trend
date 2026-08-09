#!/usr/bin/env bash
#SBATCH --job-name=svc_trend_psock_smoke
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=24G
#SBATCH --time=00:20:00
#SBATCH --output=/home/rwkays/svc_trend_psock_smoke_%j.log
#SBATCH --error=/home/rwkays/svc_trend_psock_smoke_%j.log
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
cd "$PROJ/prototype_spatial_trend" || exit 2
"$RS" psock_smoke_test.R
echo "exit: $?"
