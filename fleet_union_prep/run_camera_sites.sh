#!/usr/bin/env bash
#SBATCH --job-name=camera_sites_by_mask
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=48G
#SBATCH --time=02:00:00
#SBATCH --output=/home/rwkays/isdm/fleet_v2u/slurm_camsites_%j.log
#SBATCH --error=/home/rwkays/isdm/fleet_v2u/slurm_camsites_%j.log
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/plotting_env"
export PATH="$DST/bin:$PATH" PROJ_DATA="$DST/share/proj" PROJ_LIB="$DST/share/proj"
"$DST/bin/Rscript" /home/rwkays/isdm/fleet_v2u/camera_sites_by_mask.R
RC=$?; echo "exit: $RC"; exit $RC
