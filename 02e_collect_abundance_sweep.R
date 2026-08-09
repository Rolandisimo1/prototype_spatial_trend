#!/usr/bin/env Rscript
# =============================================================================
# 02e_collect_abundance_sweep.R
# Walks sim_results_abundance/<abundance>_<scenario>/rep_*.RDS (from
# 01e_run_abundance_sweep.R) into one tidy table: per-region recovery
# metrics stratified by (a) abundance level and (b) realized information
# (cam_detections_total / inat_count_total, actually generated per
# replicate -- the abundance ladder's intended effect, VERIFIED rather than
# assumed), plus the ecoregion-vs-national-scalar WAIC delta. Writes
# abundance_sweep_summary.csv, then a figure: sign-recovery rate and bias
# against abundance level / realized information, faceted by trend scenario,
# bobcat's baseline rung marked -- "how abundant does a species need to be
# to resolve ecoregion-level trends."
#
# bobcat_baseline (occ_shift=0, count_log_mult=0) IS the pure ecoregion-
# trend validation (identical truth to a standalone ecoregion run) -- by
# design, no separate ecoregion-only array was run; this collector's
# bobcat_baseline rows serve double duty as that result too.
#
# RUN (needs ggplot2 -- plotting_env, not nimble_env):
#   $PROJ/HPC/conda_envs/plotting_env/bin/Rscript 02e_collect_abundance_sweep.R
# Writes: abundance_sweep_summary.csv, abundance_sweep_recovery.png/.tif
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
setwd(PROTO_DIR)

RESULTS_ROOT <- paste0(PROTO_DIR, "/sim_results_abundance")
dirs <- list.dirs(RESULTS_ROOT, recursive = FALSE)
if (length(dirs) == 0) stop("No sim_results_abundance/<abundance>_<scenario>/ dirs found -- has the array run yet?")

# ------------------------------ walk every replicate file ---------------------
read_one <- function(f) {
  x <- readRDS(f)
  base <- data.frame(
    rep = x$rep, status = x$status,
    abundance = if (!is.null(x$abundance)) x$abundance else NA_character_,
    scenario  = if (!is.null(x$scenario))  x$scenario  else NA_character_,
    occ_shift = if (!is.null(x$occ_shift)) x$occ_shift else NA_real_,
    count_log_mult = if (!is.null(x$count_log_mult)) x$count_log_mult else NA_real_,
    stringsAsFactors = FALSE
  )
  if (!identical(x$status, "ok")) {
    base$detail <- if (!is.null(x$detail)) x$detail else NA_character_
    return(base)
  }
  # metrics from compute_ecoregion_metrics() -- take the scalar columns only
  # (per-region est_*/err_*/cov_* columns kept too, just appended as-is)
  scalar_cols <- c("bias_all", "rmse_all", "coverage_all", "frac_sign_correct",
                   "sigma_region_mean", "sigma_region_lb", "sigma_region_ub",
                   "waic_primary", "waic_national", "waic_diff", "elapsed_sec")
  metrics <- x[intersect(scalar_cols, names(x))]
  region_cols <- grep("^(est_|err_|cov_)", names(x), value = TRUE)
  region_metrics <- x[region_cols]

  # information is a nested data.frame column -- flatten to top-level scalars
  info <- if (!is.null(x$information)) x$information else NULL
  info_flat <- if (!is.null(info)) {
    setNames(as.data.frame(info), paste0("info_", names(info)))
  } else {
    data.frame()
  }

  cbind(base, as.data.frame(metrics), as.data.frame(region_metrics), info_flat)
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
# Catches the exact failure mode documented in abundance_sweep_seed_diagnosis.md:
# build_reduced_constants() calls set.seed() internally with a FIXED seed (by
# design, so every replicate shares the same retained-cell/site subsample),
# and if the driver doesn't reseed per-replicate afterward, R's global RNG is
# left identical across every array task -- every "replicate" within a cell
# then simulates and fits the exact same dataset (n=1 dressed up as n=30).
# If this ever regresses again (e.g. a future driver forgets the reseed the
# same way 01e_run_abundance_sweep.R once did), refuse to produce a summary
# or plot through it silently -- stop() naming the exact offending cell.
ok_check <- combined[combined$status == "ok", ]
degenerate_check <- aggregate(
  cbind(n_unique_rmse = rmse_all, n_unique_waic = waic_primary) ~ abundance + scenario,
  data = ok_check, FUN = function(x) length(unique(x)))
n_reps_check <- aggregate(rep ~ abundance + scenario, data = ok_check, FUN = length)
names(n_reps_check)[3] <- "n_reps"
degenerate_check <- merge(degenerate_check, n_reps_check, by = c("abundance", "scenario"))
degenerate_cells <- degenerate_check[
  degenerate_check$n_reps > 1 &
  (degenerate_check$n_unique_rmse == 1 | degenerate_check$n_unique_waic == 1), ]
if (nrow(degenerate_cells) > 0) {
  cell_desc <- apply(degenerate_cells, 1, function(r)
    sprintf("  abundance=%s scenario=%s (n_reps=%s, unique rmse_all=%s, unique waic_primary=%s)",
           r["abundance"], r["scenario"], r["n_reps"], r["n_unique_rmse"], r["n_unique_waic"]))
  stop("DEGENERATE REPLICATES DETECTED -- refusing to produce a summary/plot.\n",
      "The following (abundance, scenario) cell(s) have all-identical rmse_all and/or\n",
      "waic_primary across >1 replicate -- every 'replicate' simulated and fit the\n",
      "SAME dataset (a per-replicate set.seed() is likely missing after\n",
      "build_reduced_constants(), which resets the RNG internally -- see\n",
      "abundance_sweep_seed_diagnosis.md). Fix the driver and re-run before collecting.\n",
      paste(cell_desc, collapse = "\n"))
}
cat("guard passed: no degenerate (all-identical) replicate cells detected\n\n")

write.csv(combined, "abundance_sweep_summary.csv", row.names = FALSE)
cat("wrote abundance_sweep_summary.csv (", nrow(combined), "rows )\n\n")

cat("status counts:\n"); print(table(combined$abundance, combined$scenario, combined$status))

# ------------------------------ aggregate summary table ------------------------
ok <- combined[combined$status == "ok", ]
ABUND_LEVELS <- c("bobcat_baseline", "moderate", "common_deerlike")
ok$abundance <- factor(ok$abundance, levels = intersect(ABUND_LEVELS, unique(ok$abundance)))

agg <- aggregate(
  cbind(bias_all, rmse_all, coverage_all, frac_sign_correct, waic_diff,
       info_cam_detections_total, info_inat_count_total) ~ abundance + scenario,
  data = ok, FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass)
n_reps <- aggregate(rep ~ abundance + scenario, data = ok, FUN = length)
names(n_reps)[3] <- "n_reps"
agg <- merge(agg, n_reps, by = c("abundance", "scenario"))
agg <- agg[order(agg$scenario, agg$abundance), ]

cat("\n==================== ABUNDANCE SWEEP SUMMARY ====================\n")
print(agg, row.names = FALSE, digits = 3)

# WAIC preference rate: fraction of replicates where the ecoregion model was
# actually preferred (waic_diff > 0), per (abundance, scenario) -- a cleaner
# "was the complexity warranted" readout than the mean diff alone, since
# WAIC differences aren't on a bounded/comparable scale across replicates.
waic_pref <- aggregate(waic_diff ~ abundance + scenario, data = ok,
                       FUN = function(x) mean(x > 0, na.rm = TRUE))
names(waic_pref)[3] <- "frac_ecoregion_preferred"
cat("\nfraction of replicates where the ecoregion model beat national-scalar on WAIC:\n")
print(merge(agg[, c("abundance","scenario","n_reps")], waic_pref), row.names = FALSE)

# ------------------------------ figure -----------------------------------------
# panel 1: recovery vs abundance level (categorical), faceted by scenario
p_abund <- ggplot(agg, aes(x = abundance, group = scenario)) +
  geom_col(aes(y = frac_sign_correct), fill = "steelblue", alpha = 0.6, na.rm = TRUE) +
  geom_point(aes(y = abs(bias_all) * 2), color = "firebrick", size = 2, na.rm = TRUE) +
  geom_line(aes(y = abs(bias_all) * 2), color = "firebrick", group = 1, na.rm = TRUE) +
  scale_y_continuous(
    name = "sign-recovery rate (bars)",
    sec.axis = sec_axis(~ . / 2, name = "|bias| (points/line)")) +
  facet_wrap(~ scenario, labeller = labeller(scenario = c(varying = "VARYING (real regional trend)",
                                                          null = "NULL (no regional trend)"))) +
  labs(x = "abundance level", title = "Ecoregion trend recovery vs. species abundance",
      subtitle = "bobcat_baseline = actual bobcat magnitude; moderate/common_deerlike = deer-anchored scaling") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1))

# panel 2: recovery vs REALIZED information (continuous -- verifies the
# ladder actually produced more data, not just assumed it), faceted by scenario
p_info <- ggplot(ok, aes(x = info_inat_count_total, y = frac_sign_correct, color = abundance)) +
  geom_point(alpha = 0.6, na.rm = TRUE) +
  geom_smooth(aes(group = scenario), method = "loess", se = FALSE, color = "black",
             linewidth = 0.6, na.rm = TRUE) +
  facet_wrap(~ scenario) +
  labs(x = "realized iNat count (total, this replicate)", y = "sign-recovery rate (per-replicate)",
      color = "abundance level",
      title = "Recovery vs. realized information (verifies the abundance ladder, not just assumed)") +
  theme_bw()

combined_plot <- cowplot::plot_grid(p_abund, p_info, ncol = 1, align = "v")
ggsave("abundance_sweep_recovery.png", combined_plot, width = 10, height = 10, dpi = 150, bg = "white")
cat("\nwrote abundance_sweep_recovery.png\n")
