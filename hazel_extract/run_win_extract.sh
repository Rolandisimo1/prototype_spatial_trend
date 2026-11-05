#!/usr/bin/env bash
#SBATCH --job-name=win_extract
#SBATCH --partition=compute
#SBATCH --qos=long
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --output=/share/rwkays/rwkays/win_extract/logs/win_extract_%A_%a.log
# One run dir per array task. Line SLURM_ARRAY_TASK_ID of $LIST (default runs.txt) = "RUN [AFFECTED [CHAIN_GLOB]]".
# Or pass RUN (and AFFECTED / CHAIN_GLOB) via --export for a single job.
# Submit from /share/rwkays/rwkays/win_extract:  sbatch --array=1-N[%K] run_win_extract.sh
set -uo pipefail
# 2026-11-02: env moved off the full /rsstu share to usrapps; outputs default to /share scratch (chains still READ from /rsstu).
DST=/usr/local/usrapps/rwkays/conda_envs/nimble_env; export PATH="$DST/bin:$PATH"; RS="$DST/bin/Rscript"
export TMPDIR=/share/rwkays/rwkays/tmp
export OUT=${OUT:-/share/rwkays/rwkays/isdm_v2u/posteriors_20261101}; mkdir -p "$OUT"
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
W=/home/rwkays/isdm/fleet_v2u/win_extract
if [ -z "${RUN:-}" ]; then
  read -r RUN AFF GLOB <<< "$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${LIST:-$W/runs.txt})"
  export RUN; [ -n "${AFF:-}" ] && export AFFECTED=$AFF; [ -n "${GLOB:-}" ] && export CHAIN_GLOB=$GLOB
fi
echo "RUN=$RUN AFFECTED=${AFFECTED:-} CHAIN_GLOB=${CHAIN_GLOB:-} host=$(hostname) start=$(date -Is)"
"$RS" $W/win_extract.R; RC=$?
echo "exit: $RC end=$(date -Is)"; exit $RC
