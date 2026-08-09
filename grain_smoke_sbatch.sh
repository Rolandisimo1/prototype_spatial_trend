#!/usr/bin/env bash
#SBATCH --job-name=svc_trend_grain_smoke
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --time=00:30:00
#SBATCH --array=1,16
#SBATCH --output=/home/rwkays/svc_trend_grain_smoke_%A_%a.log
#SBATCH --error=/home/rwkays/svc_trend_grain_smoke_%A_%a.log
# grain_smoke_sbatch.sh -- 2-task smoke test (grain=0 rep=1, grain=1 rep=1)
# of the new grain-sweep pipeline before committing to the full 45-task
# array. Not a deliverable.
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
cd "$PROJ/prototype_spatial_trend" || exit 2
"$DST/bin/Rscript" 01c_run_grain_sweep.R
echo "exit: $?"
