#!/usr/bin/env Rscript
# =============================================================================
# 02f_collect_indicator_test.R
# Walks sim_results_indicator/<abundance>_<scenario>/rep_*.RDS (from
# 01f_run_indicator_test.R) into one tidy table, then FP/FN rates of the
# RJMCMC switch indicator at support thresholds 0.5 and 0.9 -- the direct
# replacement for the WAIC gate shown unusable by the corrected abundance
# sweep (see abundance_sweep_evaluation.md, rjmcmc_indicator_test_plan.md).
#   null scenario:    false-positive rate = frac(reps with P(gamma=1) >= threshold)
#   varying scenario: false-negative rate = frac(reps with P(gamma=1) <  threshold)
# Compared against the Goldstein et al. published single-effect benchmark:
# FP 9.1% / FN 14.9% at 0.5; FP 0.7% / FN 38.1% at 0.9.
#
# RUN (needs ggplot2 -- plotting_env, not nimble_env):
#   $PROJ/HPC/conda_envs/plotting_env/bin/Rscript 02f_collect_indicator_test.R
# Writes: indicator_test_summary.csv, indicator_test_recovery.png/.tif
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
setwd(PROTO_DIR)

RESULTS_ROOT <- paste0(PROTO_DIR, "/sim_results_indicator")
dirs <- list.dirs(RESULTS_ROOT, recursive = FALSE)
if (length(dirs) == 0) stop("No sim_results_indicator/<abundance>_<scenario>/ dirs found -- has the array run yet?")

# ------------------------------ walk every replicate file ---------------------
read_one <- function(f) {
  x <- readRDS(f)
  base <- data.frame(
    rep = x$rep, status = x$status,
    abundance = if (!is.null(x$abundance)) x$abundance else NA_character_,
    scenario  = if (!is.null(x$scenario))  x$scenario  else NA_character_,
    stringsAsFactors = FALSE
  )
  if (!identical(x$status, "ok")) {
    base$detail <- if (!is.null(x$detail)) x$detail else NA_character_
    return(base)
  }
  cbind(base,
       gamma_mean  = x$gamma_mean,
       rhat_logLik = x$rhat_logLik,
       ess_logLik  = x$ess_logLik,
       gamma_by_chain = x$gamma_by_chain,
       elapsed_sec = x$elapsed_sec)
}

files <- unlist(lapply(dirs, list.files, pattern = "^rep_.*\\.RDS$", full.names = TRUE))
cat("found", length(files), "replicate result files across", length(dirs), "abundance x scenario dirs\n")

rows <- lapply(files, read_one)
all_cols <- unique(unlist(lapply(rows, names)))
combined <- do.call(rbind, lapply(rows, function(r) {
  missing <- setdiff(all_cols, names(r))
  for (m in missing) r[[m]] <- NA
  r[all_cols]
}))

# ------------------------------ guard: refuse degenerate (unreseeded) replicates ---
# Same anti-regression check added to 02e_collect_abundance_sweep.R after the
# abundance sweep's seeding regression (see abundance_sweep_seed_diagnosis.md):
# if every replicate in an (abundance, scenario) cell reports the exact same
# gamma_mean, either the underlying sim_data was byte-identical across reps
# (a seeding regression) or the RJMCMC sampler is stuck (e.g. never actually
# proposing a jump) -- either way, refuse to produce a summary/plot through it
# silently.
ok_check <- combined[combined$status == "ok", ]
degenerate_check <- aggregate(
  cbind(n_unique_gamma = gamma_mean) ~ abundance + scenario,
  data = ok_check, FUN = function(x) length(unique(x)))
n_reps_check <- aggregate(rep ~ abundance + scenario, data = ok_check, FUN = length)
names(n_reps_check)[3] <- "n_reps"
degenerate_check <- merge(degenerate_check, n_reps_check, by = c("abundance", "scenario"))
degenerate_cells <- degenerate_check[degenerate_check$n_reps > 1 & degenerate_check$n_unique_gamma == 1, ]
if (nrow(degenerate_cells) > 0) {
  cell_desc <- apply(degenerate_cells, 1, function(r)
    sprintf("  abundance=%s scenario=%s (n_reps=%s, unique gamma_mean=%s)",
           r["abundance"], r["scenario"], r["n_reps"], r["n_unique_gamma"]))
  stop("DEGENERATE REPLICATES DETECTED -- refusing to produce a summary/plot.\n",
      "The following (abundance, scenario) cell(s) have all-identical gamma_mean\n",
      "across >1 replicate -- either the underlying simulated data was identical\n",
      "across reps (a seeding regression -- see abundance_sweep_seed_diagnosis.md\n",
      "for the exact failure mode this guards against) or the RJMCMC sampler is\n",
      "stuck. Investigate before collecting.\n",
      paste(cell_desc, collapse = "\n"))
}
cat("guard passed: no degenerate (all-identical) replicate cells detected\n\n")

write.csv(combined, "indicator_test_raw.csv", row.names = FALSE)
cat("wrote indicator_test_raw.csv (", nrow(combined), "rows )\n\n")

cat("status counts:\n"); print(table(combined$abundance, combined$scenario, combined$status))

# ------------------------------ convergence summary (data log-lik, NOT gamma) --
# rjmcmc_indicator_test_plan.md gotcha #1: R-hat on gamma is meaningless (a
# Bernoulli indicator can sit at 0/1 for long stretches while mixing well).
# Convergence here is judged on rhat_logLik/ess_logLik (Gelman-Rubin + ESS on
# total_data_logLik, computed in fit_replicate_switch()) -- gamma_mean is
# reported separately, purely as the posterior inclusion probability.
ok <- combined[combined$status == "ok", ]
cat("==================== CONVERGENCE (judged on total_data_logLik, not gamma) ====================\n")
cat("rhat_logLik summary (Gelman-Rubin on the data log-likelihood, across all", nrow(ok), "ok replicates):\n")
print(summary(ok$rhat_logLik))
n_bad_rhat <- sum(ok$rhat_logLik > 1.1, na.rm = TRUE)
cat("replicates with rhat_logLik > 1.1 (conventional non-convergence flag):", n_bad_rhat, "/", nrow(ok), "\n")
cat("\ness_logLik summary (effective sample size on the data log-likelihood, pooled across chains):\n")
print(summary(ok$ess_logLik))
cat("\n(gamma itself is NOT used for this convergence check -- see gamma_by_chain column in\n",
   "indicator_test_raw.csv for a per-chain sanity check that both starting points, gamma=1 and\n",
   "gamma=0, ended up agreeing on gamma_mean.)\n\n")

# CONVERGENCE CAVEAT, surfaced deliberately rather than smoothed over: a
# nontrivial fraction of replicates fail the conventional rhat_logLik<=1.1
# bar (inherited MCMC budget -- 1000 burnin + 3000 iter per chain -- was
# never certified for per-replicate convergence even in the original
# ecoregion model, which only ever ran 1 chain; this indicator test is the
# first time that's been checked at all). Rather than silently drop or
# silently trust those replicates, the FP/FN table below is computed BOTH
# ways -- all "ok" replicates, and restricted to rhat_logLik<=1.1 only -- so
# the reader can see whether/how much the convergence issue moves the
# headline numbers.
ok$converged <- !is.na(ok$rhat_logLik) & ok$rhat_logLik <= 1.1
n_nonconverged <- sum(!ok$converged)
cat("replicates flagged non-converged (rhat_logLik > 1.1 or NA):", n_nonconverged, "/", nrow(ok), "\n\n")

# ------------------------------ FP/FN at thresholds 0.5 and 0.9 ----------------
ABUND_LEVELS <- c("bobcat_baseline", "moderate", "common_deerlike")
ok$abundance <- factor(ok$abundance, levels = intersect(ABUND_LEVELS, unique(ok$abundance)))
THRESHOLDS <- c(0.5, 0.9)

make_summary <- function(dat) {
  do.call(rbind, lapply(THRESHOLDS, function(thr) {
    do.call(rbind, lapply(levels(dat$abundance), function(ab) {
      null_cell <- dat[dat$abundance == ab & dat$scenario == "null", ]
      varying_cell <- dat[dat$abundance == ab & dat$scenario == "varying", ]
      data.frame(
        abundance = ab, threshold = thr,
        n_null = nrow(null_cell),
        fp_rate = if (nrow(null_cell) > 0) mean(null_cell$gamma_mean >= thr) else NA_real_,
        fp_count = if (nrow(null_cell) > 0) sprintf("%d/%d", sum(null_cell$gamma_mean >= thr), nrow(null_cell)) else NA_character_,
        n_varying = nrow(varying_cell),
        fn_rate = if (nrow(varying_cell) > 0) mean(varying_cell$gamma_mean < thr) else NA_real_,
        fn_count = if (nrow(varying_cell) > 0) sprintf("%d/%d", sum(varying_cell$gamma_mean < thr), nrow(varying_cell)) else NA_character_,
        mean_gamma_null = if (nrow(null_cell) > 0) mean(null_cell$gamma_mean) else NA_real_,
        mean_gamma_varying = if (nrow(varying_cell) > 0) mean(varying_cell$gamma_mean) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  }))
}

summary_rows <- make_summary(ok)
summary_rows <- summary_rows[order(summary_rows$threshold, summary_rows$abundance), ]
summary_rows$convergence_filter <- "all_ok_replicates"

summary_rows_converged <- make_summary(ok[ok$converged, ])
summary_rows_converged <- summary_rows_converged[order(summary_rows_converged$threshold, summary_rows_converged$abundance), ]
summary_rows_converged$convergence_filter <- "rhat_logLik_le_1.1_only"

summary_both <- rbind(summary_rows, summary_rows_converged)
write.csv(summary_both, "indicator_test_summary.csv", row.names = FALSE)
cat("wrote indicator_test_summary.csv (", nrow(summary_both), "rows -- both convergence filters )\n\n")

cat("==================== INDICATOR TEST: FP/FN SUMMARY, ALL OK REPLICATES ====================\n")
print(summary_rows[, c("abundance","threshold","fp_count","fp_rate","fn_count","fn_rate")],
     row.names = FALSE, digits = 3)
cat("\n==================== INDICATOR TEST: FP/FN SUMMARY, CONVERGED ONLY (rhat_logLik<=1.1) ====================\n")
print(summary_rows_converged[, c("abundance","threshold","fp_count","fp_rate","fn_count","fn_rate")],
     row.names = FALSE, digits = 3)

# ------------------------------ figure -----------------------------------------
# Published single-effect benchmark (Goldstein et al.): FP 9.1% / FN 14.9% at
# threshold 0.5; FP 0.7% / FN 38.1% at threshold 0.9.
benchmark <- data.frame(
  threshold = c(0.5, 0.5, 0.9, 0.9),
  metric    = c("FP", "FN", "FP", "FN"),
  rate      = c(0.091, 0.149, 0.007, 0.381)
)

plot_df <- rbind(
  data.frame(abundance = summary_rows$abundance, threshold = summary_rows$threshold,
            metric = "FP", rate = summary_rows$fp_rate),
  data.frame(abundance = summary_rows$abundance, threshold = summary_rows$threshold,
            metric = "FN", rate = summary_rows$fn_rate)
)
plot_df$abundance <- factor(plot_df$abundance, levels = ABUND_LEVELS)
plot_df$threshold_lab <- factor(paste0("threshold = ", plot_df$threshold),
                                levels = c("threshold = 0.5", "threshold = 0.9"))
benchmark$threshold_lab <- factor(paste0("threshold = ", benchmark$threshold),
                                  levels = c("threshold = 0.5", "threshold = 0.9"))

p <- ggplot(plot_df, aes(x = abundance, y = rate, color = metric, group = metric)) +
  geom_hline(data = benchmark, aes(yintercept = rate, color = metric),
            linetype = "dashed", alpha = 0.6) +
  geom_point(size = 2.5) +
  geom_line() +
  facet_wrap(~ threshold_lab) +
  scale_y_continuous(labels = scales::percent, limits = c(0, NA)) +
  scale_color_manual(values = c(FP = "firebrick", FN = "steelblue")) +
  labs(x = "abundance level (bobcat_baseline = actual bobcat magnitude)",
      y = "rate", color = "metric",
      title = "RJMCMC switch indicator: false-positive / false-negative rates vs abundance",
      subtitle = sprintf("dashed lines = Goldstein et al. published benchmark (FP 9.1%%/0.7%%, FN 14.9%%/38.1%%). All %d ok replicates shown (%d/%d flagged rhat_logLik>1.1 -- see indicator_test_summary.csv for converged-only rates).",
                         nrow(ok), n_nonconverged, nrow(ok))) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

ggsave("indicator_test_recovery.png", p, width = 10, height = 6, dpi = 150, bg = "white")
cat("\nwrote indicator_test_recovery.png\n")
