#!/usr/bin/env bash
#SBATCH --job-name=fleet_regress
#SBATCH --partition=compute
#SBATCH --qos=normal
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=03:00:00
#SBATCH --output=/home/rwkays/isdm/fleet_v2u/slurm_regress_%j.log
#SBATCH --error=/home/rwkays/isdm/fleet_v2u/slurm_regress_%j.log
# Regression gate for the generalised scripts. Exit 0 ONLY if
#  (1) build_union_grid_fleet.R reproduces the 2026-09-04 union grids for
#      moose/bobcat/WTD byte-for-byte (md5), with all pins passing, and
#  (2) run_fleet_prep.R at default settings rebuilds moose_v2b's base bundle
#      with every field identical (diff_bundles.R verdict EXACT REPRODUCTION).
set -uo pipefail
PROJ=/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final
DST="$PROJ/HPC/conda_envs/plotting_env"
export PATH="$DST/bin:$PATH" PROJ_DATA="$DST/share/proj" PROJ_LIB="$DST/share/proj"
W=/home/rwkays/isdm/fleet_v2u; T=$W/regress; rm -rf "$T"; mkdir -p "$T/union" "$T/moose" "$T/cellmaps"
cd "$PROJ" || exit 2
echo "host $(hostname) start $(date)"

echo "=== (1) union grids, 3 pinned species ==="
SPECIES="moose bobcat white-tailed_deer" \
IUCN_SOURCE="$PROJ/data/moose_inat_grid.csv" \
PRES_GRID="$PROJ/output/inat_grids_v2/{species}_inat_grid_v2.csv" \
UNMASKED_CACHE="$PROJ/data/{species}_inat_grid_unmasked_v2.csv" \
PINS_CSV="$W/union_pins.csv" UNION_OUT_DIR="$T/union" \
  "$DST/bin/Rscript" "$W/build_union_grid_fleet.R" || { echo "REGRESS FAIL: union script errored"; exit 1; }
for sp in moose bobcat white-tailed_deer; do
  a=$(md5sum < "$PROJ/output/inat_grids_union/${sp}_inat_grid_union.csv" | cut -d' ' -f1)
  b=$(md5sum < "$T/union/${sp}_inat_grid_union.csv" | cut -d' ' -f1)
  echo "  $sp  reference $a  rebuilt $b"
  [ "$a" = "$b" ] || { echo "REGRESS FAIL: union grid differs for $sp"; exit 1; }
  [ -e "$T/union/${sp}_REVIEW_REQUIRED" ] && { echo "REGRESS FAIL: pinned $sp flagged"; exit 1; }
done
echo "  union grids: byte-identical x3"

echo "=== (2) moose_v2b base bundle rebuild via run_fleet_prep.R defaults ==="
SPECIES=moose OUT_DIR="$T/moose" CELLMAP_DIR="$T/cellmaps" \
  "$DST/bin/Rscript" "$W/run_fleet_prep.R" || { echo "REGRESS FAIL: prep errored"; exit 1; }
[ -e "$T/moose/BUNDLE_REVIEW_REQUIRED" ] && { echo "REGRESS FAIL: pinned moose flagged"; exit 1; }
BUNDLE_A="$PROJ/HPC/moose_v2b/input_data_moose_v2b.RDS" BUNDLE_B="$T/moose/input_data_moose_v2b.RDS" \
  "$DST/bin/Rscript" /home/rwkays/isdm/windowed_v2b/diff_bundles.R | tee "$T/moose_diff.txt" | tail -8
grep -q "VERDICT: EXACT REPRODUCTION" "$T/moose_diff.txt" || { echo "REGRESS FAIL: moose bundle not exact"; exit 1; }
cmp "$PROJ/cell_maps/cell_map_50_100_moose_v2b.RDS" "$T/cellmaps/cell_map_50_100_moose_v2b.RDS" && echo "  cell map identical"
echo "REGRESS PASS  end $(date)"
