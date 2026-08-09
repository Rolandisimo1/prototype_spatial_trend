#!/usr/bin/env bash
#SBATCH --job-name=svc_trend_dryrun
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=/home/rwkays/svc_trend_dryrun_%j.log
#SBATCH --error=/home/rwkays/svc_trend_dryrun_%j.log
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "node: $(hostname)"; echo "start: $(date)"
"$RS" dry_run.R
echo "exit: $?  end: $(date)"
