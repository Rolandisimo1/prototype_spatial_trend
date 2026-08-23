# sim_helpers_abundance.R
# ---------------------------------------------------------------------------
# Add-on helpers for the ABUNDANCE / OCCUPANCY sweep of the ecoregion-trend
# simulation. Sourced AFTER sim_helpers.R (reuses true_param_list(),
# simulate_replicate_data(), fit_replicate(), run_one_replicate()).
#
# Purpose: bobcat is widespread-but-rare -> a LOW-information stress case.
# Recovery of ecoregion-level trends should only get EASIER for higher-
# occupancy / higher-abundance species. This sweep places bobcat on a
# spectrum by scaling the two magnitudes that control per-cell information:
#   (1) occupancy   via link_occ_intercept  (cloglog scale, additive shift)
#   (2) iNat count  via theta0              (log scale, additive shift = log multiplier)
# Both are applied to the TRUTH parameter list BEFORE data are simulated, so
# more-abundant scenarios generate more camera detections and higher iNat
# counts through the model's own likelihoods (dOcc_v, dnbinom) -- faithful,
# not hand-rolled.
#
# NB scope: this varies abundance/detection while keeping bobcat's spatial
# arrangement, range extent, adjacency, covariates and effort fixed. It does
# NOT emulate a species with a genuinely different RANGE SHAPE (narrow-range
# endemic, single-ecoregion species). For a specific high-priority species
# with a very different distribution, re-run the whole prototype calibrated
# to that species' own fitted posterior. This sweep answers the fleet-wide
# question "how does resolvability scale with abundance", not "exactly how
# will species X behave".
# ---------------------------------------------------------------------------

#' @name scale_truth_abundance
#' @description Return a copy of a truth parameter list with occupancy and/or
#'   iNat count magnitude shifted, to emulate a more (or less) abundant/
#'   detectable species than bobcat. Shifts are on the linear-predictor scale
#'   so they compose cleanly with the existing spatial fields.
#' @param truth Named list from true_param_list() (bobcat baseline).
#' @param occ_shift numeric; ADDITIVE shift applied to every element of
#'   link_occ_intercept (cloglog scale). >0 raises occupancy. The spatial
#'   structure of the intercept field is preserved (shift is constant across
#'   cells); only its level moves. Note: large positive shifts saturate psi
#'   toward 1 (realistic for a near-ubiquitous species like deer) -- when psi
#'   saturates, camera occupancy carries little trend information, which is a
#'   useful stress in its own right and should be visible in the metrics.
#' @param count_log_mult numeric; ADDITIVE shift applied to theta0 (log scale),
#'   i.e. log() of the multiplicative factor on expected iNat count. E.g.
#'   log(3) makes the species report ~3x bobcat's iNat counts per cell-year.
#' @param label character; scenario tag carried through for the summary.
#' @return A truth list identical to `truth` except link_occ_intercept and
#'   theta0 are shifted; carries attr(,"abundance_label").
scale_truth_abundance <- function(truth, occ_shift = 0, count_log_mult = 0,
                                  label = NA_character_) {
  stopifnot(is.list(truth), "link_occ_intercept" %in% names(truth),
            "theta0" %in% names(truth))
  out <- truth
  out$link_occ_intercept <- truth$link_occ_intercept + occ_shift
  out$theta0             <- truth$theta0 + count_log_mult
  attr(out, "abundance_label") <- label
  attr(out, "occ_shift")       <- occ_shift
  attr(out, "count_log_mult")  <- count_log_mult
  out
}

#' @name abundance_levels_default
#' @description Default 3-level abundance ladder used if no deer-derived anchor
#'   is supplied. Levels are deliberately modest so the compute stays bounded
#'   (3 abundance x 2 trend scenarios). Overridden by
#'   00c_anchor_abundance_from_deer.R when the deer prep bundle is reachable.
#'   Values are placeholders on the linear-predictor scale, NOT fitted numbers.
#' @return data.frame(level, occ_shift, count_log_mult).
abundance_levels_default <- function() {
  data.frame(
    level          = c("bobcat_baseline", "moderate", "common"),
    occ_shift      = c(0.0, 0.75, 1.5),           # cloglog-scale occupancy lift
    count_log_mult = c(0.0, log(3), log(8)),      # iNat count x1, x3, x8
    stringsAsFactors = FALSE
  )
}

#' @name abundance_levels_measured
#' @description Abundance ladder anchored on MEASURED per-window independent
#'   detection counts from the raw camera data, replacing the placeholder
#'   multipliers in abundance_levels_default(). Prefer this for the estimator
#'   sweep: the estimators are expected to diverge along the abundance axis, so
#'   the axis should be real rather than notional.
#' @return data.frame(level, occ_shift, count_log_mult, target_mean_count).
#' @details Measured from combined_sequences_all.csv after (a) collapsing
#'   age/sex class rows to one row per (deployment, sequence, species),
#'   (b) a 30-minute independence filter, (c) effort bounds (0, 365] nights.
#'   Mean independent detections per occupied 10-day window:
#'
#'     bobcat              1.555   (n = 6,079 occupied windows, 28.3% > 1)
#'     moose               1.757   (n =   609 occupied windows, 34.3% > 1)
#'     white-tailed deer   6.405   (n = 60,125 occupied windows, 78.6% > 1)
#'
#'   Note bobcat and moose are nearly identical on this axis (1.56 vs 1.76) --
#'   moose is data-POOR, not low-abundance-per-window, which is why moose
#'   cannot stand in for the low rung and bobcat anchors it.
#'
#'   The intermediate rung is the GEOMETRIC mean of the two anchors (3.156),
#'   i.e. equal log-spacing, which is the right spacing because the count
#'   enters the linear predictor on the log scale.
#'
#'   Multipliers relative to the bobcat baseline: 1.00x / 2.03x / 4.12x.
#'   These are smaller than the placeholder 1x/3x/8x because they are measured
#'   rather than assumed -- the real abundance range between a widespread-rare
#'   carnivore and the most abundant ungulate is ~4x in per-window detection
#'   rate, not ~8x. Using the placeholder would have overstated the gradient
#'   the estimators are being compared across.
#'
#'   Source: rn_count_calibration.csv, and the derivation in
#'   raw_data_verification.md.
abundance_levels_measured <- function() {
  lo <- 1.555; hi <- 6.405
  mid <- sqrt(lo * hi)          # 3.156, equal log-spacing
  data.frame(
    level             = c("bobcat_like", "intermediate", "deer_like"),
    occ_shift         = c(0.0, 0.75, 1.5),
    count_log_mult    = c(0.0, log(mid / lo), log(hi / lo)),
    target_mean_count = c(lo, mid, hi),
    stringsAsFactors  = FALSE
  )
}

#' @name summarize_simulated_information
#' @description Diagnostic: given a simulated replicate's data, report the raw
#'   information content actually generated -- so we can VERIFY the abundance
#'   ladder produced the intended spread and later map it onto real species.
#'   Reports total & per-informed-cell camera detections and iNat counts.
#' @param sim_data List with y, y_inat from simulate_replicate_data().
#' @param constants Reduced constants_list (for informed-cell bookkeeping).
#' @return one-row data.frame of information summaries.
summarize_simulated_information <- function(sim_data, constants) {
  cam_det   <- sum(sim_data$y, na.rm = TRUE)
  inat_tot  <- sum(sim_data$y_inat, na.rm = TRUE)
  n_cell50  <- nrow(sim_data$y_inat)
  n_site    <- nrow(sim_data$y)
  data.frame(
    cam_detections_total = cam_det,
    cam_det_per_site     = cam_det / max(n_site, 1),
    inat_count_total     = inat_tot,
    inat_count_per_cell  = inat_tot / max(n_cell50, 1),
    n_site = n_site, n_cell50 = n_cell50,
    row.names = NULL
  )
}
