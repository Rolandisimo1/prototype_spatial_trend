#!/usr/bin/env Rscript
# =============================================================================
# 02k_collect_tier1.R -- stack and summarize the Tier 1 data-volume sweep
# (01k_run_tier1_volume.R). One summary row per design cell x effect level.
#
# For the reported trend (year_beta, camera-anchored) and for total_var_beta:
#   detect_rate  share of replicates whose 95% interval excludes zero. Under an
#                effect this is POWER; under the null it is the FALSE-POSITIVE
#                rate (same column, read by the effect level).
#   bias_pp      mean (posterior mean - truth), in %/yr
#   ci_width_pp  mean 95% interval width, in %/yr
#   coverage     share of intervals containing the truth
# %/yr conversion is the fleet's own: 100*(exp(x * s_occ) - 1) for year_beta,
# 100*(exp(x * s_vals) - 1) for total_var_beta (trend_table_fleet.R).
# Every metric is reported TWICE (2026-10-17): replicate_set "all", and
# "rhat_split<1.1" = replicates whose split-R-hat on year_beta is below 1.1
# (one chain per replicate, halved). frac_converged and median ESS are given
# per cell. Agreement between the two sets does NOT establish convergence (single-chain split-R-hat misses slow
# mixing); the per-cell ESS distribution (min / p10 / median, share < 100) and the ess_flag column are the check,
# and the low-ESS concentration tables say which design factors drive it.
# Stops on a degenerate cell (identical estimates across replicates = the
# missing-reseed bug) rather than summarizing it.
# =============================================================================
args <- commandArgs(trailingOnly = TRUE)
dirs <- if (length(args)) args else c("tier1_out_RW", "tier1_out_VOL")
rows <- do.call(rbind, lapply(dirs, function(d) {
  fs <- list.files(d, "^row_\\d+\\.rds$", full.names = TRUE)
  do.call(rbind, lapply(fs, function(f) { x <- readRDS(f); if (is.null(x$error_msg)) x$error_msg <- NA_character_; x }))
}))
cat("rows:", nrow(rows), " OK:", sum(rows$status == "OK"), " errored:", sum(rows$status != "OK"), "\n")
if (any(rows$status != "OK")) print(table(rows$error_msg[rows$status != "OK"]))
ok <- rows[rows$status == "OK", ]

pct <- function(x, s) 100 * (exp(x * s) - 1)
ok$yb_pct_mean <- pct(ok$yb_mean, ok$d_s_occ);  ok$yb_pct_lb <- pct(ok$yb_ci_lb, ok$d_s_occ)
ok$yb_pct_ub   <- pct(ok$yb_ci_ub, ok$d_s_occ); ok$yb_pct_true <- pct(ok$yb_true, ok$d_s_occ)

cellvars <- c("block", "range_cells", "window_frac", "n_site", "inat_records", "effect")
key <- do.call(paste, c(ok[cellvars], sep = "|"))
for (k in unique(key)) {
  s <- ok[key == k, ]
  if (nrow(s) > 2 && length(unique(round(s$yb_mean, 12))) == 1L)
    stop("DEGENERATE cell ", k, ": identical yb_mean across ", nrow(s), " replicates (missing reseed)")
}
RHAT_MAX <- 1.1
# ESS is the convergence check that matters here (2026-10-19): split-R-hat on one halved chain is blind to slow
# stationary mixing. A cell whose ESS on year_beta is systematically low has coverage / detection numbers that
# measure the sampler, not the estimator, and needs more iterations before its power number is reportable.
ESS_MIN <- 100          # per-replicate floor for a 95% interval
ESS_CELL_MED <- 400     # cell median below this -> flag the cell
stopifnot("yb_rhat_split" %in% names(ok))
ok$converged <- !is.na(ok$yb_rhat_split) & ok$yb_rhat_split < RHAT_MAX
summ_set <- function(ok, set) do.call(rbind, lapply(split(ok, do.call(paste, c(ok[cellvars], sep = "|"))), function(s) {
  data.frame(s[1, cellvars], replicate_set = set, n_rep = nrow(s),
    yb_detect_rate = mean(s$yb_detected), yb_bias_pp = mean(s$yb_pct_mean - s$yb_pct_true),
    yb_ci_width_pp = mean(s$yb_pct_ub - s$yb_pct_lb), yb_coverage = mean(s$yb_covered),
    tvb_detect_rate = mean(s$tvb_detected), tvb_coverage = mean(s$tvb_covered),
    sink_mean = mean(s$year_var_mean), sink_true = s$yv_true[1],
    det_sites = median(s$r_cam_det_sites), inat_records = median(s$r_inat_records),
    ncell50 = s$d_ncell50[1], site_pool = s$d_site_pool[1],
    yb_rhat_med = median(s$yb_rhat_split, na.rm = TRUE),
    yb_ess_min = min(s$yb_ess, na.rm = TRUE), yb_ess_p10 = unname(quantile(s$yb_ess, .1, na.rm = TRUE)),
    yb_ess_med = median(s$yb_ess, na.rm = TRUE), frac_ess_lt100 = mean(s$yb_ess < ESS_MIN, na.rm = TRUE),
    elapsed_h_med = median(s$elapsed_sec) / 3600, row.names = NULL)
}))
sa <- summ_set(ok, "all")
sa$frac_converged <- as.numeric(tapply(ok$converged, key, mean)[do.call(paste, c(sa[cellvars], sep = "|"))])
sc <- summ_set(ok[ok$converged, ], paste0("rhat_split<", RHAT_MAX)); sc$frac_converged <- NA_real_
sa$ess_flag <- ifelse(sa$yb_ess_med < ESS_CELL_MED | sa$frac_ess_lt100 > 0.10,
                      "LOW ESS: more iterations before reporting", "")
sc$ess_flag <- NA_character_
summ <- rbind(sa, sc)
summ <- summ[do.call(order, summ[c(cellvars, "replicate_set")]), ]
options(width = 220); print(summ, row.names = FALSE, digits = 3)
# Where do low-ESS replicates concentrate? Share of all replicates vs share of the ESS < ESS_MIN replicates, by
# each design factor separately and by cell (expectation: range 900 and the 20% window, if anywhere).
ok$low_ess <- !is.na(ok$yb_ess) & ok$yb_ess < ESS_MIN
cat(sprintf("\nESS on year_beta: %d of %d replicates below %d (%.1f%%); overall quantiles:\n",
            sum(ok$low_ess), nrow(ok), ESS_MIN, 100 * mean(ok$low_ess)))
print(round(quantile(ok$yb_ess, c(0, .01, .05, .1, .25, .5), na.rm = TRUE), 1))
for (f in c("block", "range_cells", "window_frac", "n_site", "inat_records", "effect")) {
  tb <- do.call(rbind, lapply(split(ok, ok[[f]]), function(s) data.frame(level = s[[f]][1], n_rep = nrow(s),
    share_of_reps = nrow(s) / nrow(ok), share_of_low = sum(s$low_ess) / max(1, sum(ok$low_ess)),
    frac_low = mean(s$low_ess), ess_med = median(s$yb_ess, na.rm = TRUE), ess_p10 = unname(quantile(s$yb_ess, .1, na.rm = TRUE)))))
  cat("\n-- by", f, "\n"); print(tb, row.names = FALSE, digits = 3)
}
cat("\n-- flagged cells (median ESS <", ESS_CELL_MED, "or >10% of replicates below", ESS_MIN, ")\n")
fl <- sa[sa$ess_flag != "", c(cellvars, "n_rep", "yb_ess_min", "yb_ess_p10", "yb_ess_med", "frac_ess_lt100",
                              "yb_coverage", "yb_detect_rate", "frac_converged")]
if (nrow(fl)) print(fl, row.names = FALSE, digits = 3) else cat("none\n")
write.csv(ok, "tier1_rows.csv", row.names = FALSE)
write.csv(summ, "tier1_summary.csv", row.names = FALSE)
cat("wrote tier1_rows.csv, tier1_summary.csv\n")
