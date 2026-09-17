# Union-mask fleet prep (2026-09-17)

Scripts that build union-mask bundles for any fleet species, plus the evidence
tables behind the 9/17 decisions. Deployed copies live on Hazel in
`~/isdm/fleet_v2u/`; these are the versions of record.

## The correction
Fitted v2b bundles used two different range definitions on the two halves of one
model: the iNaturalist grid used the presence mask, while camera deployments were
filtered by the IUCN polygon in `prep_spoccupancy_data()`. The union mask is what
the reports describe, so both halves now use it, defined on the same cell50 set.

## Scripts
| file | role |
|---|---|
| `iucn_source_lib.R` | in-range cell50 from a shared-column CSV, a per-species layer, or a file template. Rasterises exactly as Arielle's `prep_inat_data_grid()` (simplifyGeom 0.01 -> project -> rasterize touches=TRUE) |
| `build_union_grid_fleet.R` | union grid per species (IUCN OR presence), pins, review flags, edge report + map |
| `prep_camera_sites.R` | camera-side prep with a parameterised range test (`legacy_depid` / `iucn_site` / `union_cell`) |
| `run_fleet_prep.R` | generalised `run_v2b_prep.R`; `MASK_TAG`, `INAT_GRID`, `CAMERA_CACHE_DIR`, pins |
| `join_ecoregion_fleet.R` | ecoregion + national-scalar forks for any mask tag |
| `build_camcount_bundle.R` | adds `y_cam` / `log_effort_cam` / `yday_mean` for the count-camera arm |
| `build_name_crosswalk.R` | builds and verifies `species_name_crosswalk.csv` |

## Validation gates (all passed 2026-09-17, jobs 856233 / 858380 / 863644)
1. `run_fleet_prep.R` at defaults rebuilds `moose_v2b`'s base bundle: **67/67 fields identical**.
2. `iucn_source_lib.R` on the old shapefiles reproduces the shared-column IUCN cell
   sets **exactly** for moose / bobcat / WTD / coyote. NOTE `grid50.tif` carries 3,364
   cells but the grids only carry the 3,322 with iNat records, so IUCN cells are
   restricted to that universe (dropped per species: bobcat 37, WTD 34, moose 3, coyote 41).
3. `prep_camera_sites.R` in `legacy_depid` mode reproduces the cached
   `forSPO/<sp>_4SPO.RDS` **field for field** for all four species.
4. fleet34-layer union grids are **md5-identical** to the 2026-09-04 grids (moose/bobcat/WTD).

## Two traps this encodes
- **The cached camera files are not what today's prep code produces.** They kept
  in-polygon (matched by `deployment_id`) & complete siteCovs & complete raster covs,
  with **no single-occasion filter**. Rerunning `prep_spoccupancy_data(redo = TRUE)`
  would silently drop ~17% more sites (bobcat 3,426, WTD 3,947, coyote 4,061) and
  confound that with any mask change.
- **Species names differ per source** (`species_name_crosswalk.csv`): armadillo is
  *Dasypus mexicanus* to iNaturalist and IUCN but *D. novemcinctus* to Wildlife
  Insights; fisher is *Martes pennanti* to IUCN only. Keying armadillo on one name
  everywhere empties one half of the model with no error.

## Evidence tables
`camera_sites_by_mask.csv` (sites per filter, all 34 species; the current-filter column
reproduces every fitted bundle's nsite), `deployment_id_fix_delta.csv` (per-site vs
deployment_id match: 0 added, 0-3 removed), `species_name_crosswalk.csv`, `union_pins.csv`.
