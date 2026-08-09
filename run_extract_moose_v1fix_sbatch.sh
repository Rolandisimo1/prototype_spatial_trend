#!/usr/bin/env bash
#SBATCH --job-name=moose_v1fix_posterior_extract
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=01:00:00
#SBATCH --output=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/prototype_spatial_trend/slurm_extract_moose_v1fix_%j.log
#SBATCH --error=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/prototype_spatial_trend/slurm_extract_moose_v1fix_%j.log
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "node: $(hostname)"; echo "start: $(date)"
"$RS" "$PROJ/prototype_spatial_trend/extract_moose_v1fix_posterior.R"
echo "exit: $?  end: $(date)"
