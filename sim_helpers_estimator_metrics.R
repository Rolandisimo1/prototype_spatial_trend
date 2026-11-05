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
  # yb_true is recorded on its OWN, not only inside the sum (2026-10-17). The
  # reportable ecological trend is year_beta, the camera-anchored slope;
  # year_var is the iNaturalist-only sink. Every run before this date scored
  # tvb against tvb_true and kept only year_beta_mean -- no interval and no
  # separate truth -- so year_beta bias/coverage/power cannot be recovered
  # from those rows and needs a re-run. The tvb_* family is kept alongside so
  # the two can be compared: the gap between them is what the sink absorbs.
  yb_true  <- truth$year_beta

  tvb_mean <- mean(tvb)
  tvb_lb   <- stats::quantile(tvb, 0.025, names = FALSE)
  tvb_ub   <- stats::quantile(tvb, 0.975, names = FALSE)

  # "detected" = 95% CI excludes zero. Under the varying scenario, the rate of
  # this across replicates IS power; under null it is the false-positive rate.
  detected <- as.integer(tvb_lb > 0 | tvb_ub < 0)

  # same four quantities for year_beta, the reported parameter
  yb_mean <- mean(yb)
  yb_lb   <- stats::quantile(yb, 0.025, names = FALSE)
  yb_ub   <- stats::quantile(yb, 0.975, names = FALSE)
  # split-R-hat on year_beta (2026-10-17). Each replicate is ONE chain, so the
  # chain is halved and the two halves treated as two chains (BDA3 split-R-hat).
  # It detects drift / non-stationarity within the chain; it cannot see a chain
  # that sits stably in the wrong mode. Converged = yb_rhat_split < 1.1.
  yb_rhat_split <- .split_rhat(yb)
  # Effective sample size alongside it: on one chain split-R-hat is blind to
  # slow but stationary mixing (an AR(0.99) chain of 8,000 scores ~1.004), so
  # low ESS is the signal for that failure. Recorded, not used as a filter.
  yb_ess <- if (requireNamespace("coda", quietly = TRUE) && all(is.finite(yb)))
    unname(coda::effectiveSize(yb)) else NA_real_

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
    # power family for year_beta (the reported, camera-anchored trend)
    yb_true           = yb_true,
    yb_mean           = yb_mean,
    yb_ci_lb          = yb_lb,
    yb_ci_ub          = yb_ub,
    yb_bias           = yb_mean - yb_true,
    yb_ci_width       = yb_ub - yb_lb,
    yb_covered        = as.integer(yb_true >= yb_lb & yb_true <= yb_ub),
    yb_detected       = as.integer(yb_lb > 0 | yb_ub < 0),
    yb_rhat_split     = yb_rhat_split,
    yb_ess            = yb_ess,
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

#' @name .split_rhat
#' @description Split-R-hat (Gelman et al. BDA3 eq. 11.4) of one chain: the
#'   chain is cut into two equal halves treated as separate chains.
#' @param x Numeric vector of post-burnin draws.
#' @return R-hat, or NA if the chain is constant, too short or not finite.
.split_rhat <- function(x) {
  x <- x[is.finite(x)]; n <- floor(length(x) / 2)
  if (n < 2) return(NA_real_)
  h <- cbind(x[1:n], x[(n + 1):(2 * n)])
  W <- mean(apply(h, 2, stats::var))
  if (!is.finite(W) || W <= 0) return(NA_real_)
  B <- n * stats::var(colMeans(h))
  sqrt(((n - 1) / n * W + B / n) / W)
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
summarize_estimator_sweep <- function(rows, rhat_max = 1.1) {
  ok <- rows[rows$status == "OK", , drop = FALSE]
  if (nrow(ok) == 0) stop("no successful replicates to summarize")
  # Every metric is reported twice (2026-10-17): over all replicates, and over
  # replicates with split-R-hat(year_beta) < rhat_max. If the two agree,
  # under-convergence does not explain the power numbers; if they diverge,
  # that is the result. Runs scored before 2026-10-17 have no yb_rhat_split
  # and get the all-replicate set only.
  all_s <- .summarize_set(ok)
  all_s$replicate_set <- "all"
  if (!"yb_rhat_split" %in% names(ok)) return(all_s)
  conv <- ok[!is.na(ok$yb_rhat_split) & ok$yb_rhat_split < rhat_max, , drop = FALSE]
  frac <- tapply(!is.na(ok$yb_rhat_split) & ok$yb_rhat_split < rhat_max,
                 paste(ok$estimator, ok$abundance, ok$scenario, sep = "|"), mean)
  all_s$frac_converged <- as.numeric(frac[paste(all_s$estimator, all_s$abundance, all_s$scenario, sep = "|")])
  if (nrow(conv) == 0) return(all_s)
  conv_s <- .summarize_set(conv, guard = FALSE)
  conv_s$replicate_set <- paste0("rhat_split<", rhat_max)
  conv_s$frac_converged <- NA_real_
  rbind(all_s, conv_s)
}

.summarize_set <- function(ok, guard = TRUE) {

  key <- paste(ok$estimator, ok$abundance, ok$scenario, sep = "|")
  for (k in if (guard) unique(key) else character(0)) {
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
    # year_beta family; NA for runs scored before 2026-10-17 (no yb_* columns)
    yb_detect_rate = if ("yb_detected" %in% names(ok)) as.numeric(agg(mean, "yb_detected")) else NA_real_,
    yb_bias        = if ("yb_bias"     %in% names(ok)) as.numeric(agg(mean, "yb_bias"))     else NA_real_,
    yb_ci_width    = if ("yb_ci_width" %in% names(ok)) as.numeric(agg(mean, "yb_ci_width")) else NA_real_,
    yb_coverage    = if ("yb_covered"  %in% names(ok)) as.numeric(agg(mean, "yb_covered"))  else NA_real_,
    yb_ess_median  = if ("yb_ess" %in% names(ok)) as.numeric(agg(function(v) stats::median(v, na.rm = TRUE), "yb_ess")) else NA_real_,
    yb_rhat_median = if ("yb_rhat_split" %in% names(ok)) as.numeric(agg(function(v) stats::median(v, na.rm = TRUE), "yb_rhat_split")) else NA_real_,
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
