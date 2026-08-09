# sim_helpers_effort.R
# ---------------------------------------------------------------------------
# Add-on helpers for the EFFORT sweep (01g_run_effort_sweep.R). Sourced AFTER
# sim_helpers.R (reuses true_param_list_ecoregion(), simulate_replicate_data(),
# fit_replicate(), run_one_replicate(), compute_ecoregion_metrics(),
# build_reduced_constants()'s new site_replace/site_jitter_sd capability).
#
# Purpose: the abundance sweep varied the ANIMAL at fixed sampling effort.
# This sweep is the complement -- vary SAMPLING EFFORT at fixed (real
# bobcat) abundance, to ask how many cameras / how much iNat coverage a
# bobcat-density species needs for usable per-ecoregion trend resolution.
#
# Two independent levers:
#   CAMERAS: n_site_keep passed to build_reduced_constants() (already a
#     parameter there) -- detections per site are NOT set directly, they
#     emerge from the occupancy submodel at the fixed bobcat intensity, so a
#     larger camera array yields proportionally more information the same
#     way more real cameras would. Below the real reduced baseline (700
#     sites, the "1x" used throughout this project -- see build_reduced_
#     constants()'s site_replace docstring for why 700 and not the literal
#     full 20,531-site dataset) = plain subsample. Above baseline (2x/4x) =
#     resample WITH replacement (site_replace=TRUE), an explicit
#     extrapolation flagged in every design row.
#   INAT: inat_effort is take-what-you-get real data, not something more of
#     can be simulated by adding cells -- so it's scaled by a multiplier on
#     the EXISTING effort matrix instead (scale_inat_effort() below).
# ---------------------------------------------------------------------------

#' @name scale_inat_effort
#' @description Scale a (reduced) inat_effort matrix by a GLOBAL multiplier,
#'   and optionally OVERRIDE one focal ecoregion's rows with a separate
#'   multiplier instead (secondary "one low-effort region, like an
#'   under-sampled country" experiment). The focal override REPLACES the
#'   global multiplier for that region's rows -- it does not compound with
#'   it (so focal_mult=0.1 means "this region's real effort x0.1", not
#'   "global_mult x 0.1").
#' @param inat_effort Reduced effort matrix (n_cell50_kept x nyear), from
#'   build_reduced_constants()$inat_effort.
#' @param inat_cell100 constants$inat_cell100 -- cell100 id per retained
#'   cell50 row, SAME ORDER as inat_effort's rows (this is what makes the
#'   per-region row lookup below correct).
#' @param ecoregion_of_cell100 Full length-ncell100 lookup vector (cell100
#'   id -> ecoregion id), from prepped_sim_inputs.RDS.
#' @param global_mult Applied to every row (default 1 = unchanged).
#' @param focal_ecoregion_id,focal_mult Optional. If supplied, rows whose
#'   ecoregion == focal_ecoregion_id get (original real effort * focal_mult)
#'   instead of the global-multiplied value. NULL (default) = no override.
#' @return Scaled effort matrix, same dimensions as inat_effort.
scale_inat_effort <- function(inat_effort, inat_cell100, ecoregion_of_cell100,
                              global_mult = 1, focal_ecoregion_id = NULL, focal_mult = NULL) {
  out <- inat_effort * global_mult
  if (!is.null(focal_ecoregion_id)) {
    row_region <- ecoregion_of_cell100[inat_cell100]
    focal_rows <- which(row_region == focal_ecoregion_id)
    if (length(focal_rows) > 0) {
      out[focal_rows, ] <- inat_effort[focal_rows, , drop = FALSE] * focal_mult
    }
  }
  out
}

#' @name summarize_simulated_information_by_region
#' @description Per-region analogue of summarize_simulated_information()
#'   (sim_helpers_abundance.R) -- reports realized camera detections and
#'   iNat counts PER ECOREGION, so effort -> information -> resolution can
#'   actually be mapped region-by-region rather than only in aggregate. Also
#'   returns the per-region SITE and CELL50 counts actually retained (so a
#'   region with a small real camera pool is visible even before looking at
#'   detections).
#' @param sim_data List with y, y_inat from simulate_replicate_data().
#' @param constants Reduced constants_list (cell = site->cell100 lookup,
#'   inat_cell100 = cell50->cell100 lookup).
#' @param ecoregion_of_cell100 Full length-ncell100 lookup (cell100 id ->
#'   ecoregion id).
#' @param ecoregion_levels Character vector of region names, length nregion
#'   (for readable output).
#' @return data.frame, one row per ecoregion: n_site, cam_detections,
#'   n_cell50, inat_count.
summarize_simulated_information_by_region <- function(sim_data, constants,
                                                       ecoregion_of_cell100, ecoregion_levels) {
  nregion <- length(ecoregion_levels)
  site_region   <- ecoregion_of_cell100[constants$cell]
  cell50_region <- ecoregion_of_cell100[constants$inat_cell100]

  cam_det_by_site  <- rowSums(sim_data$y, na.rm = TRUE)
  inat_by_cell50   <- rowSums(sim_data$y_inat, na.rm = TRUE)

  do.call(rbind, lapply(seq_len(nregion), function(r) {
    data.frame(
      ecoregion      = ecoregion_levels[r],
      ecoregion_id   = r,
      n_site         = sum(site_region == r, na.rm = TRUE),
      cam_detections = sum(cam_det_by_site[site_region == r], na.rm = TRUE),
      n_cell50       = sum(cell50_region == r, na.rm = TRUE),
      inat_count     = sum(inat_by_cell50[cell50_region == r], na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
}

#' @name camera_levels_default
#' @description The 4 camera-axis levels for the effort sweep, expressed as
#'   n_site_keep values relative to this project's established 700-site
#'   baseline (NOT the literal full 20,531-site real dataset -- see
#'   build_reduced_constants()'s site_replace docstring). 0.5x/1x subsample
#'   without replacement (site_replace=FALSE); 2x/4x resample WITH
#'   replacement (site_replace=TRUE), an explicit extrapolation.
#' @return data.frame(level, n_site_keep, site_replace).
camera_levels_default <- function() {
  data.frame(
    level        = c("0.5x", "1x", "2x", "4x"),
    n_site_keep  = c(350, 700, 1400, 2800),
    site_replace = c(FALSE, FALSE, TRUE, TRUE),
    stringsAsFactors = FALSE
  )
}

#' @name inat_levels_default
#' @description The 4 iNat-axis levels for the effort sweep: global
#'   multipliers on the real reduced-design effort matrix.
#' @return data.frame(level, global_mult).
inat_levels_default <- function() {
  data.frame(
    level       = c("0.25x", "1x", "4x", "16x"),
    global_mult = c(0.25, 1, 4, 16),
    stringsAsFactors = FALSE
  )
}
