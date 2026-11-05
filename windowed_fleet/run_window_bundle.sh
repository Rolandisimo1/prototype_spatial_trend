#!/usr/bin/env bash
#SBATCH --job-name=win_bundle
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/home/rwkays/isdm/fleet_v2u/windowed/logs/win_bundle_%x_%j.log
# Windowed national_scalar bundle (2026-10-17), both streams cut to the same years.
# Usage: sbatch --job-name=<sp>_<wtag> --export=ALL,SP=bobcat,TAG=v2u,WTAG=win10,YEAR_START=2016,YEAR_END=2025 run_window_bundle.sh
# Steps 3-4 of run_one_species_bundle.sh, unchanged except for the year window:
#  3. run_fleet_prep.R on the SAME union grid + SAME union camera cache the 18-yr bundle was built from,
#     with YEAR_START/YEAR_END/WINDOW_TAG set (pin checked on the full record, then both streams cut)
#  4. national_scalar fork, 9-cov reduction (make_reduced_input_v2.R), verify_bundle.R,
#     reorder_last_site.R if the window's last site has J < 2 (guard refuses otherwise)
# Output: HPC/<sp>_<TAG>_national_scalar_<WTAG>/input_data_<same>.RDS (pre-reduce copy in HPC/<sp>_<TAG>_<WTAG>/).
# REGRESS=1: no window env (the builder's unwindowed path) into a scratch OUT_DIR, then compare to the fitted
#   18-yr bundle; WREGRESS=1: the windowed code path with 2008-2025, same comparison. Neither writes run dirs.
set -uo pipefail
: "${SP:?set SP}" "${TAG:?set TAG}"
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
PE="$PROJ/HPC/conda_envs/plotting_env"; NE="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$PE/bin:$PATH" PROJ_DATA="$PE/share/proj" PROJ_LIB="$PE/share/proj"
W=/home/rwkays/isdm/fleet_v2u; WD=$W/windowed
export NAME_CROSSWALK=$W/species_name_crosswalk.csv
UG=$PROJ/output/inat_grids_union_fleet34; CAM=$PROJ/data/forSPO_$TAG
SRC=${SP}_${TAG}_national_scalar
cd "$PROJ" || exit 2
echo "host $(hostname) start $(date)  species $SP tag $TAG wtag ${WTAG:-} ${YEAR_START:-}-${YEAR_END:-} regress ${REGRESS:-0} wregress ${WREGRESS:-0}"
[ -f "$CAM/${SP}_4SPO.RDS" ] || { echo "FAIL: camera cache $CAM/${SP}_4SPO.RDS missing"; exit 1; }

if [ "${REGRESS:-0}" = 1 ] || [ "${WREGRESS:-0}" = 1 ]; then
  RT=$([ "${REGRESS:-0}" = 1 ] && echo regress || echo wregress)
  OD=$WD/regress/${SP}_${TAG}_$RT; mkdir -p $OD/cellmaps; rm -f $OD/input_data_${SP}_${TAG}.RDS
  WENV=(); [ "$RT" = wregress ] && WENV=(YEAR_START=2008 YEAR_END=2025 WINDOW_TAG=full18)
  env SPECIES=$SP MASK_TAG=$TAG INAT_GRID="$UG/{species}_inat_grid_union.csv" CAMERA_CACHE_DIR=$CAM \
    PINS_CSV=$W/union_pins.csv PIN_KIND=union OUT_DIR=$OD CELLMAP_DIR=$OD/cellmaps "${WENV[@]}" \
    "$PE/bin/Rscript" $W/run_fleet_prep.R || { echo "FAIL 3 prep ($RT)"; exit 1; }
  B=$OD/input_data_${SP}_${TAG}.RDS
  cp $B $OD/input_data_${SRC}.RDS
  "$NE/bin/Rscript" /home/rwkays/isdm/R/make_reduced_input_v2.R $SP $OD/input_data_${SRC}.RDS || { echo "FAIL 4 reduce ($RT)"; exit 1; }
  "$NE/bin/Rscript" $WD/compare_to_fitted.R "$SRC" "$B" "$OD/input_data_${SRC}_v2cov.RDS" "$RT"
  echo "REGRESS DONE $RT $SP end $(date)"; exit $?
fi

: "${WTAG:?set WTAG}" "${YEAR_START:?}" "${YEAR_END:?}"
TOK=${SRC}_${WTAG}; D=$PROJ/HPC/$TOK; B=$D/input_data_${TOK}.RDS; PD=$PROJ/HPC/${SP}_${TAG}_${WTAG}
[ -e "$B" ] && { echo "FAIL: $B already exists -- refusing to overwrite"; exit 1; }
[ -e "$PD/input_data_${SP}_${TAG}.RDS" ] && { echo "FAIL: $PD prep bundle exists -- refusing"; exit 1; }

echo "=== 3. bundle (window $YEAR_START-$YEAR_END) ==="
SPECIES=$SP MASK_TAG=$TAG INAT_GRID="$UG/{species}_inat_grid_union.csv" CAMERA_CACHE_DIR=$CAM \
PINS_CSV=$W/union_pins.csv PIN_KIND=union OUT_DIR=$PD \
YEAR_START=$YEAR_START YEAR_END=$YEAR_END WINDOW_TAG=$WTAG \
  "$PE/bin/Rscript" $W/run_fleet_prep.R || { echo "FAIL 3 prep"; exit 1; }

echo "=== 4. national_scalar fork, 9-cov, verify, last-site J ==="
mkdir -p $D
cp $PD/input_data_${SP}_${TAG}.RDS $B
[ -e $PD/BUNDLE_REVIEW_REQUIRED ] && cp $PD/BUNDLE_REVIEW_REQUIRED $D/
"$NE/bin/Rscript" /home/rwkays/isdm/R/make_reduced_input_v2.R $SP $B || { echo "FAIL 4 reduce"; exit 1; }
mv $B $B.bak_10cov && mv $D/input_data_${TOK}_v2cov.RDS $B
"$NE/bin/Rscript" /home/rwkays/isdm/R/verify_bundle.R $TOK || echo "WARN 4 verify_bundle FAILED for $TOK (see above)"
JL=$("$NE/bin/Rscript" -e "b <- readRDS('$B'); cat(tail(b\$constants_list\$J, 1))")
echo "J_last $JL"
if [ "$JL" -lt 2 ]; then "$NE/bin/Rscript" $W/reorder_last_site.R $TOK || { echo "FAIL 4 reorder"; exit 1; }; fi
grep -h "BUNDLE SUMMARY" $WD/logs/win_bundle_${SLURM_JOB_NAME}_${SLURM_JOB_ID}.log
ls $D/BUNDLE_REVIEW_REQUIRED 2>/dev/null && echo "REVIEW GATE SET"
echo "WINDOW BUNDLE DONE $TOK end $(date)"
