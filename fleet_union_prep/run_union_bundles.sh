#!/usr/bin/env bash
#SBATCH --job-name=union_bundles
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=96G
#SBATCH --time=12:00:00
#SBATCH --output=/home/rwkays/isdm/fleet_v2u/slurm_union_bundles_%j.log
#SBATCH --error=/home/rwkays/isdm/fleet_v2u/slurm_union_bundles_%j.log
# Union-mask (v2u) bundles with the camera filter on the SAME union mask.
# Gates, each a hard stop:
#  A. camera prep legacy mode reproduces the cached 4SPO of record (4 species)
#  B. fleet34-layer union grids byte-identical to the 09-04 union grids
#     (moose/bobcat/WTD), pins pass
#  C. union camera nsite == counts reported 9/17 (job 860979)
#  D. every bundle: 9 cov, max VIF < 10 (verify_bundle.R)
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
PE="$PROJ/HPC/conda_envs/plotting_env"; NE="$PROJ/HPC/conda_envs/nimble_env"
export PATH="$PE/bin:$PATH" PROJ_DATA="$PE/share/proj" PROJ_LIB="$PE/share/proj"
W=/home/rwkays/isdm/fleet_v2u; TAG=v2u
export NAME_CROSSWALK=$W/species_name_crosswalk.csv
LAYER=$W/fleet34_isdm/fleet34_iucn_ranges_present_day.geojson
UG=$PROJ/output/inat_grids_union_fleet34
cd "$PROJ" || exit 2
echo "host $(hostname) start $(date)"
[ -f "$NAME_CROSSWALK" ] || { echo "FAIL: crosswalk missing"; exit 1; }
declare -A NSITE=( [moose]=3044 [bobcat]=23471 [white-tailed_deer]=22238 [coyote]=24463 )
ALL="moose bobcat white-tailed_deer coyote"

echo "=== A. camera prep legacy reproduction ==="
VAL=$PROJ/data/forSPO_validate_legacy_$(date +%Y%m%d%H%M)
for sp in $ALL; do
  "$PE/bin/Rscript" -e "source('$W/prep_camera_sites.R'); invisible(prep_camera_sites('$sp', 'legacy_depid', out_dir = '$VAL'))" || { echo "FAIL A build $sp"; exit 1; }
  "$PE/bin/Rscript" $W/compare_4spo.R "$PROJ/data/forSPO/${sp}_4SPO.RDS" "$VAL/${sp}_4SPO.RDS" || { echo "FAIL A: $sp legacy camera prep does not reproduce the cache"; exit 1; }
done

echo "=== B. union grids from fleet34 layer ==="
SPECIES="$ALL" IUCN_SOURCE="$LAYER" IUCN_KEY=species \
PRES_GRID="$PROJ/output/inat_grids_v2/bobcat_inat_grid_v2.csv" \
UNMASKED_CACHE="$PROJ/data/{species}_inat_grid_unmasked_v2.csv" \
PINS_CSV=$W/union_pins.csv UNION_OUT_DIR=$UG \
  "$PE/bin/Rscript" $W/build_union_grid_fleet.R || { echo "FAIL B union grid"; exit 1; }
for sp in moose bobcat white-tailed_deer; do
  a=$(md5sum < $PROJ/output/inat_grids_union/${sp}_inat_grid_union.csv | cut -d' ' -f1)
  b=$(md5sum < $UG/${sp}_inat_grid_union.csv | cut -d' ' -f1)
  echo "  $sp 09-04 $a  fleet34 $b"; [ "$a" = "$b" ] || { echo "FAIL B: $sp union grid changed under fleet34 layer"; exit 1; }
done

echo "=== C. union camera prep + bundles ==="
for d in coyote_v2u coyote_v2u_national_scalar; do
  [ -d "$PROJ/HPC/$d" ] && mv "$PROJ/HPC/$d" "$PROJ/HPC/${d}_superseded_campoly_20260917" && echo "  moved aside HPC/$d (old IUCN camera filter)"
done
CAM=$PROJ/data/forSPO_$TAG
for sp in $ALL; do
  "$PE/bin/Rscript" -e "source('$W/prep_camera_sites.R'); d <- prep_camera_sites('$sp', 'union_cell', union_grid = '$UG/${sp}_inat_grid_union.csv', out_dir = '$CAM'); if (nrow(d\$y) != ${NSITE[$sp]}) stop('nsite ', nrow(d\$y), ' != reported ${NSITE[$sp]}')" \
    || { echo "FAIL C camera $sp"; exit 1; }
  SPECIES=$sp MASK_TAG=$TAG INAT_GRID="$UG/{species}_inat_grid_union.csv" CAMERA_CACHE_DIR=$CAM \
  PINS_CSV=$W/union_pins.csv PIN_KIND=union \
    "$PE/bin/Rscript" $W/run_fleet_prep.R || { echo "FAIL C prep $sp"; exit 1; }
done

echo "=== D. forks, 9-cov reduction, verify ==="
TOKS=""
for sp in $ALL; do
  if [ "$sp" = coyote ]; then
    ND=$PROJ/HPC/${sp}_${TAG}_national_scalar; mkdir -p $ND
    cp $PROJ/HPC/${sp}_${TAG}/input_data_${sp}_${TAG}.RDS $ND/input_data_${sp}_${TAG}_national_scalar.RDS
    FORKS="national_scalar"
  else
    "$PE/bin/Rscript" $W/join_ecoregion_fleet.R $sp $TAG || { echo "FAIL D join $sp"; exit 1; }
    FORKS="national_scalar ecoregion"
  fi
  for fk in $FORKS; do
    TOK=${sp}_${TAG}_${fk}; D=$PROJ/HPC/$TOK; B=$D/input_data_${TOK}.RDS
    [ -e $PROJ/HPC/${sp}_${TAG}/BUNDLE_REVIEW_REQUIRED ] && cp $PROJ/HPC/${sp}_${TAG}/BUNDLE_REVIEW_REQUIRED $D/
    "$NE/bin/Rscript" /home/rwkays/isdm/R/make_reduced_input_v2.R $sp $B || { echo "FAIL D reduce $TOK"; exit 1; }
    mv $B $B.bak_10cov && mv $D/input_data_${TOK}_v2cov.RDS $B
    TOKS="$TOKS $TOK"
  done
done
"$NE/bin/Rscript" /home/rwkays/isdm/R/verify_bundle.R $TOKS || { echo "FAIL D verify"; exit 1; }
grep -h "BUNDLE SUMMARY" /home/rwkays/isdm/fleet_v2u/slurm_union_bundles_${SLURM_JOB_ID}.log
echo "UNION BUNDLES DONE end $(date)"
