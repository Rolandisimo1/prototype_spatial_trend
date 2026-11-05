#!/usr/bin/env bash
#SBATCH --job-name=win_survey
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/home/rwkays/isdm/fleet_v2u/windowed/slurm_win_survey_%j.log
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
"$PROJ/HPC/conda_envs/nimble_env/bin/Rscript" /home/rwkays/isdm/fleet_v2u/windowed/window_survey.R
