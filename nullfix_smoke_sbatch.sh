#!/usr/bin/env bash
#SBATCH --job-name=nullfix_smoke
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/rwkays/nullfix_smoke_%j.log
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
cd "$PROJ/prototype_spatial_trend" || exit 2
"$DST/bin/Rscript" 00y_smoke_test.R
RC=$?; echo "smoke exit: $RC"; exit $RC
