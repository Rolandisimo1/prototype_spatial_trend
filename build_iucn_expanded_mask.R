#!/usr/bin/env Rscript
# =============================================================================
# build_iucn_expanded_mask.R
# -----------------------------------------------------------------------------
# Third range-mask option, alongside the two already in use:
#   1. Pure IUCN polygon (the ORIGINAL mask -- discarded for moose after the
#      Colorado/Utah exclusion finding; never touched for bobcat/WTD).
#   2. Pure presence threshold (Fix 2b, prep_inat_data_grid_v2.R) -- discards
#      the IUCN polygon entirely in favor of PRESENCE_MIN_RECORDS + the
#      spatial-isolation/captive filters.
#   3. THIS SCRIPT: IUCN polygon as a base, expanded outward by a contiguous
#      chain of grid100 cells that individually clear a detection threshold --
#      per the user's prior-project convention: 100km cells, >2 detections,
#      4-connected (edge-adjacent) flood fill. A cell joins the mask only if
#      it is reachable from the IUCN polygon through an unbroken chain of
#      qualifying cells -- a genuinely isolated vagrant record far from both
#      the polygon and any other detection cluster does NOT pull in its cell,
#      even if that single cell individually has >2 detections.
#
# WHY A THIRD OPTION, NOT A REPLACEMENT: this project has not yet compared
# fleet-fit results across mask choices. Fix 2b is already validated on real
# fits (moose_v2b, converged, recovers Colorado/Utah). This script produces a
# parallel candidate so the two can be compared before either supersedes the
# other -- do not wire this into any HPC run script until that comparison has
# been made and reviewed.
#
# NOT YET RUN against real data -- the flood-fill core (below) is logic-tested
# against synthetic data only (4 cases: chain propagation, threshold block,
# isolated-eligible-cell exclusion, base-cell retention -- all passed). The
# driver below needs real per-species cell100 detection counts and the
# per-species IUCN polygon, both of which live only on Hazel:
#   - IUCN polygons: RANGE_DIR (see build_moose_unmasked_review.R for the
#     exact path convention) -- per-species, extant-only (PRESENCE==1,
#     ORIGIN==1) already extracted there.
#   - cell100 detection counts: need the SAME isolation/captive filtering
#     Fix 2b already applies (prep_inat_data_grid_v2.R) so the >2-detections
#     threshold isn't tripped by records that would already be excluded as
#     vagrants/captives elsewhere in the pipeline. Do not recompute counts
#     from raw, unfiltered records -- that would silently disagree with what
#     Fix 2b calls "presence" for the same species.
#
# USAGE (once wired to real data):
#   result <- flood_fill_expand_mask(x = cell100_df$x, y = cell100_df$y,
#                                     count = cell100_df$filtered_count,
#                                     in_range_base = cell100_df$in_iucn_polygon,
#                                     min_detections = 2)
#   cell100_df$in_range_iucn_expanded <- result
#   # then join down to cell50: a cell50 is in-range iff its parent cell100 is.
# =============================================================================

#' Flood-fill range expansion from an IUCN base, gated by a detection threshold
#'
#' @param x,y integer grid coordinates (cell100 column/row index -- NOT lon/lat)
#' @param count per-cell100 detection count, AFTER isolation/captive filtering
#'   (consistent with Fix 2b's presence-mask filters)
#' @param in_range_base logical: TRUE where the IUCN polygon covers the cell
#' @param min_detections a cell is ELIGIBLE to join the mask via expansion only
#'   if count > min_detections (strict >, matching the ">2 detections" rule)
#' @return logical vector: in_range_base OR reachable from any base cell via
#'   an unbroken 4-connected chain of eligible cells
flood_fill_expand_mask <- function(x, y, count, in_range_base, min_detections = 2) {
  n <- length(x)
  in_range <- in_range_base
  eligible <- count > min_detections
  key <- paste(x, y)

  frontier <- which(in_range)
  changed <- TRUE
  while (changed) {
    changed <- FALSE
    new_frontier <- integer(0)
    for (i in frontier) {
      neighbor_keys <- c(
        paste(x[i] + 1, y[i]), paste(x[i] - 1, y[i]),
        paste(x[i], y[i] + 1), paste(x[i], y[i] - 1)
      )
      j_vec <- match(neighbor_keys, key)  # NA for off-grid/absent neighbors
      j_vec <- j_vec[!is.na(j_vec)]
      for (j in j_vec) {
        if (!in_range[j] && eligible[j]) {
          in_range[j] <- TRUE
          new_frontier <- c(new_frontier, j)
          changed <- TRUE
        }
      }
    }
    frontier <- unique(new_frontier)
  }
  in_range
}

# =============================================================================
# Driver -- adapt paths/species before running on Hazel. Modeled on
# build_moose_unmasked_review.R's read-only, non-destructive pattern: this
# produces a labeled comparison column, it does not filter or overwrite
# anything.
# =============================================================================
build_iucn_expanded_mask_for_species <- function(sci_name, species_label,
                                                   min_detections = 2) {
  suppressPackageStartupMessages({
    library(data.table)
    library(terra)
  })

  PROJ      <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
  DATA_DIR  <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data"
  RANGE_DIR <- "/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/iucn_ranges_extantonly"

  grid100 <- rast(file.path(DATA_DIR, "grid100.tif"))

  # --- 1. IUCN polygon -> cell100 base mask ---
  # CONFIRM the exact per-species file naming convention in RANGE_DIR before
  # trusting this glob -- build_moose_unmasked_review.R references RANGE_DIR
  # but (as far as this repo shows) never actually reads a per-species file
  # from it. Verify on Hazel rather than assume the pattern below is right.
  range_files <- list.files(RANGE_DIR, pattern = sci_name, full.names = TRUE, ignore.case = TRUE)
  stopifnot("No IUCN range file found for this species -- check RANGE_DIR naming convention" =
    length(range_files) > 0)
  range_poly <- vect(range_files[1])
  range_poly <- project(range_poly, crs(grid100))

  cell100_ids_in_polygon <- unique(terra::extract(grid100, range_poly, ID = FALSE)[, 1])
  cell100_ids_in_polygon <- cell100_ids_in_polygon[!is.na(cell100_ids_in_polygon)]

  # --- 2. Real cell100 detection counts, ISOLATION/CAPTIVE-FILTERED ---
  # Must reuse Fix 2b's filtering exactly -- see prep_inat_data_grid_v2.R's
  # spatial_isolation_flag() and the captive-record check. Placeholder below
  # assumes a per-species filtered cell100 count table already exists;
  # confirm the real path/column names on Hazel before running.
  stop(paste(
    "NOT WIRED TO REAL DATA YET.",
    "Confirm on Hazel: (a) the RANGE_DIR per-species file naming convention,",
    "(b) where Fix 2b's isolation/captive-filtered cell100 counts are cached",
    "for", species_label, "(do not recompute unfiltered -- must match Fix 2b's",
    "definition of a valid detection), then replace this stop() with the real",
    "read + flood_fill_expand_mask() call + cell100->cell50 join."
  ))
}

# Example calls once wired (do not run before the stop() above is resolved):
# build_iucn_expanded_mask_for_species("Alces alces", "moose")
# build_iucn_expanded_mask_for_species("Lynx rufus", "bobcat")
# build_iucn_expanded_mask_for_species("Odocoileus virginianus", "white-tailed_deer")
