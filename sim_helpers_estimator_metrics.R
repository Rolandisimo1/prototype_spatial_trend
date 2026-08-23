# =============================================================================
# sim_helpers_estimator_metrics.R
# -----------------------------------------------------------------------------
# Metrics for the estimator-comparison sweep. Source AFTER sim_helpers.R.
#
# WHY TWO FAMILIES OF METRIC, NOT ONE
# The Brazil workshop comparison (camera-level N~1,090 sites vs array-level
# N=60 arrays, 20 species pairs) found that aggregating to arrays GENERALLY
# REDUCED statistical power but SOMETIMES INCREASED AUC. Those are different
# quantities:
#
#   power         = can we detect THAT a trend exists?      (a hypothesis test)
#   discrimination= can we say WHERE intensity is high?     (a ranking)
#
# Aggregation reduces the number of independent units (hurting power) while
# averaging out per-camera noise (potentially helping the ranking). A
# comparison that reported only one of these would give a confidently wrong
# verdict -- and reporting only power is the more tempting error, since it is
# the headline quantity for a trend paper.
#
# So every replicate returns BOTH, plus the model's own camera-corroboration
# check.
# =============================================================================

#' @name compute_estimator_metrics
#' @description Per-replicate metrics for one estimator arm: trend recovery
#'   (power family), spatial discrimination (AUC family), and the trend
#'   robustness indicator.
#' @param samples Matrix of post-burnin posterior samples (all monitored
#'   nodes), as returned by fit_replicate().
#' @param truth List of true parameter values used to simulate, from
#'   true_param_list() then scale_truth_abundance().
#' @param constants The camera-level constants list (for the cell100 mapping;
#'   note this is the CAMERA-level object even for array arms, so the
#'   discrimination metric is computed on the same spatial support in all
#'   arms and stays comparable).
#' @param cell100_geo Data frame with cell100 and region columns.
#' @return One-row data.frame of metrics.
#' @details Every quantity here is computed identically across the three
#'   estimator arms. That is the point: the arms differ only in the camera
#'   observation model, so any difference in these numbers is attributable to
#'   the estimator rather than to how it was scored.
compute_estimator_metrics <- function(samples, truth, constants, cell100_geo) {

  colm <- function(nm) if (nm %in% colnames(samples)) samples[, nm] else rep(NA_real_, nrow(samples))

  # ---- power family: the trend itself ------------------------------------
  tvb      <- colm("total_var_beta")
  yb       <- colm("year_beta")
  yv       <- colm("year_var")
  tvb_true <- truth$year_beta + truth$year_var

  tvb_mean <- mean(tvb)
  tvb_lb   <- stats::quantile(tvb, 0.025, names = FALSE)
  tvb_ub   <- stats::quantile(tvb, 0.975, names = FALSE)

  # "detected" = 95% CI excludes zero. Under the varying scenario, the rate of
  # this across replicates IS power; under null it is the false-positive rate.
  detected <- as.integer(tvb_lb > 0 | tvb_ub < 0)

  # ---- the model's own camera-corroboration check ------------------------
  # snr = year_beta / year_var; the indicator fires when the camera-anchored
  # component matches the iNat-specific one in sign and exceeds it in
  # magnitude. This is the guard against reporting an iNat effort artifact as
  # a population trend, and moose real-data currently FAILS it (0.478/0.432).
  # If an estimator raises this, that is the strongest argument for adopting
  # it -- stronger than a speed gain.
  tri <- colm("trend_robust_indicator")

  # ---- discrimination family: the spatial ranking ------------------------
  # Rank cells by fitted intensity and ask whether the ranking matches truth.
  # Computed as a rank correlation and as an AUC over the median split of
  # true intensity, so it is comparable to the Brazil AUC finding.
  disc <- .spatial_discrimination(samples, truth, cell100_geo)

  data.frame(
    # power family
    tvb_true          = tvb_true,
    tvb_mean          = tvb_mean,
    tvb_bias          = tvb_mean - tvb_true,
    tvb_ci_width      = tvb_ub - tvb_lb,
    tvb_covered       = as.integer(tvb_true >= tvb_lb & tvb_true <= tvb_ub),
    tvb_detected      = detected,
    year_beta_mean    = mean(yb),
    year_var_mean     = mean(yv),
    # camera corroboration
    tri_mean          = mean(tri),
    tri_fires         = as.integer(mean(tri) > 0.5),
    # discrimination family
    disc_spearman     = disc$spearman,
    disc_auc          = disc$auc,
    disc_rmse         = disc$rmse,
    stringsAsFactors = FALSE
  )
}

#' @name .spatial_discrimination
#' @description How well does the fitted spatial intercept field rank cells
#'   against the truth? Reported as Spearman correlation, an AUC over the
#'   median split of true intensity, and RMSE on the field itself.
#' @param samples Posterior sample matrix.
#' @param truth Truth list (must carry link_occ_intercept).
#' @param cell100_geo Data frame with cell100 column.
#' @return List with spearman, auc, rmse.
#' @details AUC is computed by the Mann-Whitney identity (the probability a
#'   randomly chosen above-median cell is ranked above a randomly chosen
#'   below-median cell), which needs no external package. Returns NA rather
#'   than a number when the required columns are absent, so a missing metric
#'   is visible as missing instead of silently defaulting.
.spatial_discrimination <- function(samples, truth, cell100_geo) {
  na_out <- list(spearman = NA_real_, auc = NA_real_, rmse = NA_real_)
  true_field <- truth$link_occ_intercept
  if (is.null(true_field)) return(na_out)

  cols <- paste0("link_occ_intercept[", seq_along(true_field), "]")
  cols <- cols[cols %in% colnames(samples)]
  if (length(cols) < 3) return(na_out)

  est <- colMeans(samples[, cols, drop = FALSE])
  idx <- as.integer(sub(".*\\[(\\d+)\\]", "\\1", cols))
  tru <- true_field[idx]

  ok <- is.finite(est) & is.finite(tru)
  if (sum(ok) < 3) return(na_out)
  est <- est[ok]; tru <- tru[ok]

  sp <- suppressWarnings(stats::cor(est, tru, method = "spearman"))

  hi <- tru > stats::median(tru)
  auc <- if (any(hi) && any(!hi)) {
    r <- rank(est)
    (sum(r[hi]) - sum(hi) * (sum(hi) + 1) / 2) / (sum(hi) * sum(!hi))
  } else NA_real_

  list(spearman = sp, auc = auc,
       rmse = sqrt(mean((est - tru)^2)))
}

#' @name summarize_estimator_sweep
#' @description Collapse per-replicate rows to per-cell summaries, with a
#'   DEGENERACY GUARD. Call this in the collector, not in the driver.
#' @param rows data.frame of stacked per-replicate results.
#' @return data.frame, one row per (estimator x abundance x scenario).
#' @details The guard exists because this project has twice shipped a sweep in
#'   which every replicate within a cell was byte-identical (a missing
#'   per-replicate reseed), making an n=1 result look like n=30. Rates and
#'   coverages from such a run are meaningless. This function STOPS rather than
#'   producing a plausible-looking summary from degenerate input.
summarize_estimator_sweep <- function(rows) {
  ok <- rows[rows$status == "OK", , drop = FALSE]
  if (nrow(ok) == 0) stop("no successful replicates to summarize")

  key <- paste(ok$estimator, ok$abundance, ok$scenario, sep = "|")
  for (k in unique(key)) {
    sub <- ok[key == k, , drop = FALSE]
    if (nrow(sub) > 2 && length(unique(round(sub$tvb_mean, 12))) == 1L) {
      stop(sprintf(
        "DEGENERATE: all %d replicates in cell '%s' have identical tvb_mean. ",
        nrow(sub), k),
        "This is the missing-reseed bug (see 01i header). Rates from this run ",
        "would be n=1 masquerading as n=", nrow(sub), ". Fix the reseed and ",
        "re-run; do not summarize.")
    }
  }

  agg <- function(f, col) tapply(ok[[col]], key, f, simplify = TRUE)
  parts <- strsplit(names(agg(mean, "tvb_mean")), "|", fixed = TRUE)

  data.frame(
    estimator   = vapply(parts, `[`, character(1), 1),
    abundance   = vapply(parts, `[`, character(1), 2),
    scenario    = vapply(parts, `[`, character(1), 3),
    n_rep       = as.integer(agg(length, "tvb_mean")),
    # power family. Under scenario=="varying" detect_rate IS power;
    # under "null" the SAME column is the false-positive rate.
    detect_rate = as.numeric(agg(mean, "tvb_detected")),
    bias        = as.numeric(agg(mean, "tvb_bias")),
    ci_width    = as.numeric(agg(mean, "tvb_ci_width")),
    coverage    = as.numeric(agg(mean, "tvb_covered")),
    # camera corroboration
    tri_rate    = as.numeric(agg(mean, "tri_fires")),
    # discrimination family
    auc         = as.numeric(agg(mean, "disc_auc")),
    spearman    = as.numeric(agg(mean, "disc_spearman")),
    field_rmse  = as.numeric(agg(mean, "disc_rmse")),
    elapsed_med = as.numeric(agg(stats::median, "elapsed_sec")),
    stringsAsFactors = FALSE, row.names = NULL
  )
}
