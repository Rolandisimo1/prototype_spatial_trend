#!/usr/bin/env bash
#SBATCH --partition=compute
#SBATCH --qos=long
#SBATCH --ntasks=1
#SBATCH --mem=110G
#SBATCH --time=96:00:00
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/HPC/$SPECIES" || exit 2
echo "node: $(hostname)"; echo "start: $(date)"
"$RS" "$PROJ/HPC/$SPECIES/HPC_run_model_single_camcount.R"
RC=$?
echo "exit: $RC  end: $(date)"
exit $RC
