#!/usr/bin/env Rscript
# =============================================================================
# prep_inat_data_grid_v2.R  --  Fix 2: range-polygon mask replacement
# -----------------------------------------------------------------------------
# Forked from Arielle's `prep_inat_data_grid()` in `integration_helper.R`
# (never touches her original -- source it AFTER integration_helper.R so this
# definition of `prep_inat_data_grid` in your session shadows hers; the
# original function is left completely unmodified in her file).
#
# PROBLEM THIS FIXES (documented in team_memo_inat_pipeline_issues.md, Issue 2):
#   The original function rasterizes each species' IUCN range polygon onto
#   the 50km grid, then sets that species' count to NA (not 0) in every cell
#   outside the polygon. NA cells are dropped from the likelihood entirely --
#   "no information," not "confirmed absence." For moose this discarded 4,300
#   real observations (37% of all moose iNat records) across 149 occupied
#   grid cells, concentrated in a real, spatially coherent Colorado
#   expansion cluster (established there since the 1970s) -- exactly the
#   kind of signal a temporal-trend model should be able to detect, thrown
#   away before the model ever saw it. Quantified cross-check against the
#   model's own FITTED footprint (2026-08-09): every one of the 382 cells
#   moose_v1fix was actually fit on has in_range==TRUE under the old mask,
#   and zero excluded cells made it in -- confirming this is a data-prep-time
#   exclusion, not just a map-rendering choice.
#
# THE FIX:
#   1. DROP the range-polygon mask entirely. No NA-ing based on a historical
#      range shapefile.
#   2. ADD a captive-record filter: exclude iNat observations flagged
#      captive_cultivated == TRUE in the raw download, when that column is
#      present (it exists in the standard iNat export schema but was never
#      read by the original function).
#   3. ADD a spatial-isolation outlier filter as the mask's replacement:
#      exclude an observation only if it has ZERO same-species records
#      within ISOLATION_RADIUS_KM km and within ISOLATION_YEARS calendar
#      years. This targets true vagrants/mis-IDs without assuming the
#      historical range polygon is current or complete -- a genuine
#      expansion cluster (e.g. Colorado moose) has neighbors and survives;
#      an isolated one-off far from any other record does not.
#
# STATUS OF THE THRESHOLDS: 100km / 3yr are the starting proposal from the
# team memo, NOT a settled decision -- flagged there as needing review before
# use in a real fit. Left as named constants below so they're one edit to
# adjust after that review, not a re-derivation.
#
# USAGE:
#   source("integration_helper.R")       # defines the ORIGINAL, untouched
#   source("prep_inat_data_grid_v2.R")   # DEFINES prep_inat_data_grid_v2()
#                                          # (does NOT overwrite the original
#                                          #  name -- call v2 explicitly)
#   dat <- prep_inat_data_grid_v2(taxon_key, "moose", redo = TRUE)
#
# DEPENDENCY: requires the `sf` package (already used elsewhere in this
# pipeline, e.g. prep_data_for_spoccupancy.R, 00b_prep_ecoregion.R) for the
# spatial-isolation filter's indexed neighbor search -- see note at
# spatial_isolation_flag() below on why this replaced an earlier terra-based
# dense-matrix draft that would not have scaled past moose's small record
# count to WTD/bobcat's much larger ones.
# =============================================================================

ISOLATION_RADIUS_KM <- 100   # NOT YET REVIEWED -- team memo starting proposal
ISOLATION_YEARS      <- 3    # NOT YET REVIEWED -- team memo starting proposal

prep_inat_data_grid_v2 <- function(taxon_key, species, redo = FALSE,
                                    isolation_radius_km = ISOLATION_RADIUS_KM,
                                    isolation_years = ISOLATION_YEARS) {

  if (!exists("INAT_GRID_DIR")) INAT_GRID_DIR <- paste0(PROJ_DIR, "/output/inat_grids_v2")
  dir.create(INAT_GRID_DIR, showWarnings = FALSE, recursive = TRUE)
  outfile <- paste0(INAT_GRID_DIR, "/", species, "_inat_grid_v2.csv")

  sciname <- taxon_key$sci_name[taxon_key$common_name_clean == species]

  if (file.exists(outfile) && file.exists(GRID50_PATH) && !redo) {
    cat("  reusing precomputed v2 iNat grid:", outfile, "\n")
    return(read_csv(outfile))
  }

  raster_brick <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/covar_raster_brick.tif")
  grid100 <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/grid100.tif")
  names(grid100) <- "cell100"

  grid50 <- terra::aggregate(raster_brick[["MAP"]] *
                               raster_brick[["agriculture_pct"]] *
                               raster_brick[["soil_clay"]],
                             fact = 10, na.rm = TRUE)
  names(grid50) <- "cell50"
  terra::values(grid50)[!is.na(terra::values(grid50))] <- 1:sum(!is.na(terra::values(grid50)))

  coords50  <- as.data.frame(grid50, xy = TRUE) %>% rename(x50 = x, y50 = y)
  coords100 <- as.data.frame(grid100, xy = TRUE) %>% rename(x100 = x, y100 = y)

  inat_dat_all <- read_csv("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/inat_combo_nam_mams.csv") %>%
    filter(public_positional_accuracy < 1000) %>%
    mutate(Year = year(observed_on)) %>%
    filter(Year >= 2008 & Year <= 2025)

  # --- FIX 2, step (2): captive-record filter, only if the column exists ---
  if ("captive_cultivated" %in% colnames(inat_dat_all)) {
    n_before <- nrow(inat_dat_all)
    inat_dat_all <- inat_dat_all %>% filter(captive_cultivated == FALSE | is.na(captive_cultivated))
    cat("  captive-record filter: dropped", n_before - nrow(inat_dat_all),
        "of", n_before, "records (captive_cultivated == TRUE)\n")
  } else {
    cat("  NOTE: 'captive_cultivated' column not found in inat_combo_nam_mams.csv --",
        "skipping captive-record filter. Verify column name against the raw",
        "iNat export schema before treating this fix as complete.\n")
  }

  inat_pts <- inat_dat_all %>%
    select(longitude, latitude) %>%
    vect(geom = c("longitude", "latitude"), crs = "+proj=longlat") %>%
    project(crs(grid100))

  inat_cells <- bind_cols(
    extract(grid50, inat_pts)[, 2],
    extract(grid100, inat_pts)[, 2])
  colnames(inat_cells) <- c("cell50", "cell100")
  inat_cells$sciname <- inat_dat_all$taxon_species_name
  inat_cells$Year    <- inat_dat_all$Year
  inat_cells$lon     <- inat_dat_all$longitude
  inat_cells$lat     <- inat_dat_all$latitude

  inat_cell_summary <- inat_cells %>%
    group_by(Year) %>%
    count(cell50, cell100) %>%
    rename(effort = n) %>%
    filter(!is.na(cell50)) %>%
    left_join(coords50) %>%
    left_join(coords100)

  # --- FIX 2, step (3): spatial-isolation filter, replacing the range mask ---
  # For each species' own observations: keep a record only if >=1 other
  # same-species record exists within isolation_radius_km AND within
  # isolation_years calendar years. No range polygon is read or used here.
  #
  # SCALABILITY NOTE: an earlier draft of this function built a dense n x n
  # pairwise distance matrix (terra::distance on all points at once). That
  # is fine for moose's ~11.6k records but does not scale -- WTD and bobcat
  # have far more iNat observations, and an n x n matrix at n in the tens of
  # thousands is both a memory blowup (n^2 doubles) and a slow O(n^2) scan.
  # Rewritten to use sf::st_is_within_distance, which uses a spatial index
  # (STRtree) to return only the neighbors actually within radius_km for
  # each point -- no dense matrix is ever materialized, and the year-window
  # check only runs on that already-small candidate set per point.
  spatial_isolation_flag <- function(spec_pts_ll, years, radius_km, year_window) {
    n <- nrow(spec_pts_ll)
    if (n <= 1) return(rep(FALSE, n))  # a single record anywhere is isolated by definition

    pts_sf <- sf::st_as_sf(spec_pts_ll, coords = c("lon", "lat"), crs = 4326)
    # CONUS-centered equidistant projection so st_is_within_distance's units are meters
    pts_proj <- sf::st_transform(pts_sf, "+proj=aeqd +lat_0=45 +lon_0=-100")

    # sparse list: for each point i, indices of ALL other points within radius_km
    # (includes self; excluded below). Uses the geometry's spatial index --
    # no n x n matrix is built.
    nbr_list <- sf::st_is_within_distance(pts_proj, dist = radius_km * 1000)

    isolated <- logical(n)
    for (i in seq_len(n)) {
      candidates <- setdiff(nbr_list[[i]], i)
      if (length(candidates) == 0) {
        isolated[i] <- TRUE
        next
      }
      within_time <- abs(years[candidates] - years[i]) <= year_window
      isolated[i] <- !any(within_time)
    }
    isolated
  }

  for (i in 1:nrow(taxon_key)) {
    spec_name <- taxon_key$common_name_clean[i]
    this_spec_dat <- inat_cells %>%
      filter(sciname == taxon_key$sci_name[i]) %>%
      group_by(Year) %>%
      count(cell50, cell100)
    colnames(this_spec_dat)[4] <- spec_name

    # flag isolated (likely vagrant/mis-ID) points using the FULL per-point
    # record (not the cell-aggregated one), then drop only THOSE raw points
    # before re-aggregating to cell x year
    spec_pts <- inat_cells %>% filter(sciname == taxon_key$sci_name[i])
    if (nrow(spec_pts) > 0) {
      iso_flag <- spatial_isolation_flag(spec_pts, spec_pts$Year,
                                          isolation_radius_km, isolation_years)
      cat("  ", spec_name, ": flagged", sum(iso_flag), "of", nrow(spec_pts),
          "records as spatially isolated (excluded)\n")
      spec_pts_kept <- spec_pts[!iso_flag, ]
      this_spec_dat <- spec_pts_kept %>%
        group_by(Year) %>%
        count(cell50, cell100)
      colnames(this_spec_dat)[4] <- spec_name
    }

    inat_cell_summary <- left_join(inat_cell_summary, this_spec_dat,
                                    by = c("Year", "cell50", "cell100"))
    inat_cell_summary[[spec_name]][is.na(inat_cell_summary[[spec_name]])] <- 0

    # NOTE: NO range-polygon NA-ing step here. Every grid cell keeps a real
    # 0 or positive count -- the model gets to see the full survey area.
  }

  write_csv(inat_cell_summary, outfile)
  inat_cell_summary
}
