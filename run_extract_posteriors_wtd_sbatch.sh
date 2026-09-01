#!/usr/bin/env bash
#SBATCH --job-name=extract_post_wtd
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/home/rwkays/extract_post_wtd_%j.log
#SBATCH --error=/home/rwkays/extract_post_wtd_%j.log
# Posterior extraction for the four converged _single fits with no posterior
# CSVs. 64G because the moose_v1fix9 _single chains are ~977 MB compressed
# each (vs ~195 MB for the chunked v1fix files the 32G job handled); chains
# are loaded and column-subset ONE AT A TIME, so peak is one chain plus the
# accumulated subsets, but the decompressed single-chain footprint is the
# binding constraint. qos=normal caps at 4 days, so 3h is accepted (qos=short
# would cap at 02:00:00).
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
echo "node: $(hostname)  start: $(date)"
"$RS" extract_posteriors_wtd.R
RC=$?; echo "extract exit: $RC"
echo "end: $(date)"; exit $RC
