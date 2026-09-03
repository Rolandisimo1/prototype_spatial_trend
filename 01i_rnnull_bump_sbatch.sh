#!/usr/bin/env bash
#SBATCH --job-name=rnnull_bump
#SBATCH --partition=compute
#SBATCH --qos=long
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --array=1-540
#SBATCH --output=/home/rwkays/rnnull_bump_%A_%a.log
#SBATCH --error=/home/rwkays/rnnull_bump_%A_%a.log
# -----------------------------------------------------------------------------
# REPLICATE BUMP ON THE NULL ARM, BOTH ROYLE-NICHOLS ESTIMATORS.
# 540 tasks = 2 estimators x 3 abundance levels x 90 reps, null scenario only.
#
# WHY BOTH RN ARMS AND NOT JUST camera_rn. In the n30_nullfix run (716014) the
# false-positive rate rose monotonically with abundance in BOTH RN arms
# (array 0 -> 3.3 -> 6.7%; camera 0 -> 10.0 -> 13.3%) and in neither occupancy
# arm. Pooled by family at deer_like that is RN 10.0% vs occupancy 3.3%,
# Fisher exact p = 0.27 -- not established, and the pooling was post-hoc. But
# testing camera_rn alone answers a narrower question than the one that
# matters: if this is an RN-family property the report must qualify BOTH RN
# recommendations; if it is camera-level-specific the array-level
# recommendation stands clean. Those are different conclusions.
#
# WHY A 1-540 ARRAY WITH A LOOKUP, NOT THE ROW IDS DIRECTLY. Hazel's
# MaxArraySize is 1001, so a Slurm array index cannot exceed 1000 and the
# camera_rn rows (1621-2070) cannot be addressed directly -- submitting them
# fails with "Invalid job array specification". The task index is therefore
# mapped to a row_id below via six block starts. The mapping is verified
# against build_design_df() at n_rep=90 before submission, not assumed.
#
# ROW IDS ARE NOT GUESSED. build_design_df() sorts by
# (estimator, abundance, scenario, rep_id) and assigns row_id AFTER sorting,
# so the mapping changes with N_REP. The six ranges above were computed by
# rebuilding the design at n_rep=90 and reading off the null x {array_rn,
# camera_rn} rows -- not extrapolated from the n_rep=30 layout. The driver
# also stops loudly if row_id > nrow(design_df), so a missing N_REP=90 export
# fails rather than silently running the wrong cells.
#
# EVERYTHING ELSE IDENTICAL TO 716014 (commit 739f6fb): same N_ITER=8000,
# N_BURNIN=2000, same DESIGN_SEED, same unmodified 01i_run_estimator_sweep.R
# (md5 036103c34543cfbc224a83e2458554c9) carrying the
# stopifnot(tvb_true == 0) guard. Only N_REP and the row subset differ.
#
# SEED-COLLISION CAVEAT, for anyone pooling these with the 30-rep results.
# set.seed(DESIGN_SEED + row_id) runs before data generation, so identical
# row_ids across two N_REP settings give identical datasets when the design
# cell also matches. Within-arm pooling is SAFE (no old/new row_id overlap
# inside either RN arm). But new array_rn rows 541-570 reuse the exact
# datasets of OLD camera_rn/bobcat_like rows 541-570, so pooling old+new into
# a cross-arm family contrast shares 30 datasets between the two RN arms.
# Report n=90 fresh; treat n=120 as correlated across arms.
#
# SEPARATE OUTPUT DIRECTORY, per the standing rule -- the driver's
# "already done, exiting" skip keys on files in SWEEP_OUTDIR, and row_ids mean
# different things at n_rep=90 than at n_rep=30, so writing into an existing
# directory would silently adopt mislabelled rows.
# -----------------------------------------------------------------------------
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$DST/bin:$PATH"
RS="$DST/bin/Rscript"
[ -x "$RS" ] || { echo "FATAL: $RS not found"; exit 2; }
cd "$PROJ/prototype_spatial_trend" || exit 2
export SIM_BASE=. N_REP=90 N_BURNIN=2000 N_ITER=8000
export SWEEP_OUTDIR="$PROJ/prototype_spatial_trend/estimator_sweep_out_n90_rnnull"
# Task index -> row_id. Six contiguous blocks of 90, in ascending row_id
# order: array_rn {bobcat_like, deer_like, intermediate} then camera_rn ditto.
STARTS=(541 721 901 1621 1801 1981)
B=$(( (SLURM_ARRAY_TASK_ID - 1) / 90 ))
O=$(( (SLURM_ARRAY_TASK_ID - 1) % 90 ))
ROW=$(( ${STARTS[$B]} + O ))
if [ "$B" -lt 0 ] || [ "$B" -gt 5 ]; then echo "FATAL: task $SLURM_ARRAY_TASK_ID out of range"; exit 2; fi
echo "node: $(hostname)  task: $SLURM_ARRAY_TASK_ID -> row_id $ROW  start: $(date)"
"$RS" 01i_run_estimator_sweep.R "$ROW"
RC=$?; echo "exit: $RC  end: $(date)"; exit $RC
