#!/usr/bin/env bash
#SBATCH --job-name=eco_effort_sweep
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=1-480%10
#SBATCH --output=/home/rwkays/eco_effort_sweep_%A_%a.log
#SBATCH --error=/home/rwkays/eco_effort_sweep_%A_%a.log
# 01g_run_effort_sweep_sbatch.sh -- 4 camera levels x 4 iNat levels x 30 reps
# = 480 tasks (adjust --array upper bound to 320 if N_REPS is dialed down to
# 20 in 01g_run_effort_sweep.R -- MUST match design's nrow() exactly or rows
# beyond nrow(design) will error out with "row not found").
#
# --time=02:00:00: --qos=short caps walltime at 2h regardless (confirmed via
# a rejected 3h smoke-test submission: "QOSMaxWallDurationPerJobLimit") --
# so this is a hard ceiling, not a choice. The 4x camera level (2800 sites,
# 4x the usual 700-site baseline) fits a meaningfully larger occupancy
# submodel than anything run in this project so far; smoke-tested across the
# full cost range (cheapest to priciest cell) before committing to the full
# array specifically to confirm every cell fits inside this 2h ceiling.
#
# Same PATH-export pattern (not conda activate) and %10 throttle as every
# other sbatch script in this project -- see 01e_run_abundance_sweep_sbatch.sh
# header note for why.
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "array task: $SLURM_ARRAY_TASK_ID"; echo "node: $(hostname)"; echo "start: $(date)"
"$RS" 01g_run_effort_sweep.R
echo "exit: $?  end: $(date)"
