#!/usr/bin/env bash
#SBATCH --job-name=svc_trend_smoke
#SBATCH --partition=compute_partners
#SBATCH --qos=short
#SBATCH --ntasks=2
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=/home/rwkays/svc_trend_smoke_%j.log
#SBATCH --error=/home/rwkays/svc_trend_smoke_%j.log
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
# nimble_env's conda-provided C++ compiler (x86_64-conda-linux-gnu-c++) is only
# on PATH via `conda activate`; since we invoke Rscript directly, prepend it
# manually or compileNimble() fails with "Failed to create the shared library"
# (confirmed via debug_compile.R, jobs 435927 -> 435928).
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "node: $(hostname)"; echo "start: $(date)"
"$RS" smoke_test.R
echo "exit: $?  end: $(date)"
