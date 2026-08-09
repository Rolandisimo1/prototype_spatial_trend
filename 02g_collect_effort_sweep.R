#!/usr/bin/env Rscript
# =============================================================================
# 02g_collect_effort_sweep.R
# Walks sim_results_effort/<camera_level>_<inat_level>/rep_*.RDS (from
# 01g_run_effort_sweep.R) into a per-region long table, then per-region RMSE
# / coverage / sign-recovery as a function of realized effort (cameras and
# iNat counts) -- "how much sampling does a bobcat-density species need for
# usable per-ecoregion trend resolution."
#
# RUN (needs ggplot2 -- plotting_env, not nimble_env):
#   $PROJ/HPC/conda_envs/plotting_env/bin/Rscript 02g_collect_effort_sweep.R
# Writes: effort_sweep_summary.csv, effort_sweep_by_region.csv,
#         effort_sweep_target_rmse_table.csv,
#         effort_sweep_rmse_vs_cameras.png, effort_sweep_rmse_vs_inat.png
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2) })

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
setwd(PROTO_DIR)

prepped <- readRDS("prepped_sim_inputs.RDS")
ecoregion_levels <- prepped$ecoregion_levels
region_names_safe <- make.names(ecoregion_levels)

RESULTS_ROOT <- paste0(PROTO_DIR, "/sim_results_effort")
dirs <- list.dirs(RESULTS_ROOT, recursive = FALSE)
if (length(dirs) == 0) stop("No sim_results_effort/<camera_level>_<inat_level>/ dirs found -- has the array run yet?")

CAMERA_LEVELS <- c("0.5x", "1x", "2x", "4x")
INAT_LEVELS   <- c("0.25x", "1x", "4x", "16x")

# ------------------------------ walk every replicate file: WIDE (status/guard) ---
read_one_wide <- function(f) {
  x <- readRDS(f)
  base <- data.frame(
    rep = x$rep, status = x$status,
    camera_level = if (!is.null(x$camera_level)) x$camera_level else NA_character_,
    n_site_keep  = if (!is.null(x$n_site_keep))  x$n_site_keep  else NA_integer_,
    site_replace = if (!is.null(x$site_replace)) x$site_replace else NA,
    inat_level   = if (!is.null(x$inat_level))   x$inat_level   else NA_character_,
    global_mult  = if (!is.null(x$global_mult))  x$global_mult  else NA_real_,
    stringsAsFactors = FALSE
  )
  if (!identical(x$status, "ok")) {
    base$detail <- if (!is.null(x$detail)) x$detail else NA_character_
    return(base)
  }
  cbind(base,
       bias_all = x$bias_all, rmse_all = x$rmse_all, coverage_all = x$coverage_all,
       frac_sign_correct = x$frac_sign_correct,
       sigma_region_mean = x$sigma_region_mean, elapsed_sec = x$elapsed_sec)
}

# ------------------------------ walk every replicate file: LONG (per-region) -----
read_one_long <- function(f) {
  x <- readRDS(f)
  if (!identical(x$status, "ok")) return(NULL)
  info <- x$info_by_region[[1]]
  info$region_safe <- make.names(info$ecoregion)
  do.call(rbind, lapply(seq_along(region_names_safe), function(r) {
    rn <- region_names_safe[r]
    info_row <- info[info$region_safe == rn, ]
    data.frame(
      rep = x$rep,
      camera_level = x$camera_level, n_site_keep = x$n_site_keep, site_replace = x$site_replace,
      inat_level = x$inat_level, global_mult = x$global_mult,
      sigma_region_mean = x$sigma_region_mean,
      ecoregion = ecoregion_levels[r],
      est = x[[paste0("est_", rn)]],
      err = x[[paste0("err_", rn)]],
      covered = x[[paste0("cov_", rn)]],
      n_site_region = if (nrow(info_row) > 0) info_row$n_site else NA_integer_,
      cam_detections_region = if (nrow(info_row) > 0) info_row$cam_detections else NA_real_,
      n_cell50_region = if (nrow(info_row) > 0) info_row$n_cell50 else NA_integer_,
      inat_count_region = if (nrow(info_row) > 0) info_row$inat_count else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
}

files <- unlist(lapply(dirs, list.files, pattern = "^rep_.*\\.RDS$", full.names = TRUE))
cat("found", length(files), "replicate result files across", length(dirs), "camera x inat dirs\n")

rows_wide <- lapply(files, read_one_wide)
all_cols <- unique(unlist(lapply(rows_wide, names)))
combined <- do.call(rbind, lapply(rows_wide, function(r) {
  missing <- setdiff(all_cols, names(r))
  for (m in missing) r[[m]] <- NA
  r[all_cols]
}))

# ------------------------------ guard: refuse degenerate (unreseeded) replicates ---
ok_check <- combined[combined$status == "ok", ]
degenerate_check <- aggregate(
  cbind(n_unique_rmse = rmse_all, n_unique_bias = bias_all) ~ camera_level + inat_level,
  data = ok_check, FUN = function(x) length(unique(x)))
n_reps_check <- aggregate(rep ~ camera_level + inat_level, data = ok_check, FUN = length)
names(n_reps_check)[3] <- "n_reps"
degenerate_check <- merge(degenerate_check, n_reps_check, by = c("camera_level", "inat_level"))
degenerate_cells <- degenerate_check[
  degenerate_check$n_reps > 1 &
  (degenerate_check$n_unique_rmse == 1 | degenerate_check$n_unique_bias == 1), ]
if (nrow(degenerate_cells) > 0) {
  cell_desc <- apply(degenerate_cells, 1, function(r)
    sprintf("  camera=%s inat=%s (n_reps=%s, unique rmse_all=%s, unique bias_all=%s)",
           r["camera_level"], r["inat_level"], r["n_reps"], r["n_unique_rmse"], r["n_unique_bias"]))
  stop("DEGENERATE REPLICATES DETECTED -- refusing to produce a summary/plot.\n",
      "The following (camera_level, inat_level) cell(s) have all-identical rmse_all\n",
      "and/or bias_all across >1 replicate -- every 'replicate' simulated and fit the\n",
      "SAME dataset (a per-replicate set.seed() is likely missing after\n",
      "build_reduced_constants(), which resets the RNG internally -- the same bug\n",
      "class documented in abundance_sweep_seed_diagnosis.md). Fix the driver and\n",
      "re-run before collecting.\n",
      paste(cell_desc, collapse = "\n"))
}
cat("guard passed: no degenerate (all-identical) replicate cells detected\n\n")

write.csv(combined, "effort_sweep_summary.csv", row.names = FALSE)
cat("wrote effort_sweep_summary.csv (", nrow(combined), "rows )\n\n")
cat("status counts:\n"); print(table(combined$camera_level, combined$inat_level, combined$status))

# ------------------------------ per-region long table + aggregation --------------
rows_long <- lapply(files, read_one_long)
long_all <- do.call(rbind, rows_long)
write.csv(long_all, "effort_sweep_by_region_raw.csv", row.names = FALSE)

long_all$camera_level <- factor(long_all$camera_level, levels = CAMERA_LEVELS)
long_all$inat_level   <- factor(long_all$inat_level, levels = INAT_LEVELS)

# ------------------------------ divergence flag (post-hoc, no re-run needed) -----
# fit_replicate() runs a SINGLE chain with no convergence diagnostic at all (the
# indicator test's fit_replicate_switch() added Gelman-Rubin on total_data_logLik
# specifically because this project had never checked convergence before -- this
# sweep's driver, 01g_run_effort_sweep.R, reuses the older single-chain
# fit_replicate() and inherited the same gap). sigma_region_mean is a cheap,
# already-saved proxy: true target_sd=0.2, prior is dexp(1) (mean 1) -- a posterior
# mean far above 1 is a strong sign of a runaway/diverged chain, not a real
# estimate. Flagging (not silently dropping) so the report can show both.
long_all$divergent <- long_all$sigma_region_mean > 1.0
n_divergent <- sum(long_all$divergent) / length(unique(long_all$ecoregion))  # per-replicate, not per-region-row
cat("replicates flagged divergent (sigma_region_mean > 1.0):", n_divergent, "/",
   length(unique(paste(long_all$camera_level, long_all$inat_level, long_all$rep))), "\n\n")

make_agg_region <- function(dat) {
  out <- do.call(rbind, lapply(split(dat, list(dat$camera_level, dat$inat_level, dat$ecoregion), drop = TRUE), function(d) {
    data.frame(
      camera_level = d$camera_level[1], n_site_keep = d$n_site_keep[1],
      inat_level = d$inat_level[1], global_mult = d$global_mult[1],
      ecoregion = d$ecoregion[1],
      n_reps = nrow(d),
      rmse = sqrt(mean(d$err^2, na.rm = TRUE)),
      bias = mean(d$err, na.rm = TRUE),
      coverage = mean(d$covered, na.rm = TRUE),
      mean_n_site_region = mean(d$n_site_region, na.rm = TRUE),
      mean_cam_detections_region = mean(d$cam_detections_region, na.rm = TRUE),
      mean_inat_count_region = mean(d$inat_count_region, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  out[order(out$camera_level, out$inat_level, out$ecoregion), ]
}

agg_region_all       <- make_agg_region(long_all)
agg_region_converged <- make_agg_region(long_all[!long_all$divergent, ])
agg_region_all$convergence_filter       <- "all_reps"
agg_region_converged$convergence_filter <- "sigma_region_mean_le_1_only"
write.csv(rbind(agg_region_all, agg_region_converged), "effort_sweep_by_region.csv", row.names = FALSE)
cat("wrote effort_sweep_by_region.csv (", nrow(agg_region_all) + nrow(agg_region_converged),
   "rows -- one per camera x inat x ecoregion cell, both convergence filters )\n\n")

# per-cell divergence rate, for the report -- concentrated at the grid's extremes
div_by_cell <- aggregate(divergent ~ camera_level + inat_level, data = long_all[!duplicated(long_all[,c("camera_level","inat_level","rep")]),],
                         FUN = function(x) sprintf("%d/%d (%.0f%%)", sum(x), length(x), 100*mean(x)))
cat("==================== DIVERGENCE RATE BY CELL (sigma_region_mean > 1.0) ====================\n")
print(div_by_cell[order(div_by_cell$camera_level, div_by_cell$inat_level), ], row.names = FALSE)

cat("\n==================== EFFORT SWEEP: PER-REGION RMSE, ALL REPS (rounded) ====================\n")
print(agg_region_all[, c("camera_level","inat_level","ecoregion","rmse","coverage","mean_inat_count_region")],
     row.names = FALSE, digits = 3)
cat("\n==================== EFFORT SWEEP: PER-REGION RMSE, DIVERGENCE-FILTERED (rounded) ====================\n")
print(agg_region_converged[, c("camera_level","inat_level","ecoregion","rmse","coverage","mean_inat_count_region")],
     row.names = FALSE, digits = 3)

# Use the DIVERGENCE-FILTERED table as primary for the target-RMSE table and
# figures below -- the raw/all-reps numbers are still in effort_sweep_by_region.csv
# for transparency, but reporting the outlier-contaminated version as the headline
# would be misleading (a handful of runaway chains dominate a squared-error metric).
agg_region <- agg_region_converged

# ------------------------------ inverted "effort needed for target RMSE" table ---
TARGET_RMSE <- c(0.10, 0.15, 0.20)
target_table <- do.call(rbind, lapply(TARGET_RMSE, function(target) {
  do.call(rbind, lapply(ecoregion_levels, function(reg) {
    d <- agg_region[agg_region$ecoregion == reg & agg_region$rmse <= target, ]
    if (nrow(d) == 0) {
      return(data.frame(ecoregion = reg, target_rmse = target,
                        camera_level = NA_character_, inat_level = NA_character_,
                        achieved_rmse = NA_real_, note = "not achieved anywhere in grid",
                        stringsAsFactors = FALSE))
    }
    # cheapest qualifying cell, cost = n_site_keep * global_mult (simple combined proxy)
    d$cost <- d$n_site_keep * d$global_mult
    best <- d[which.min(d$cost), ]
    data.frame(ecoregion = reg, target_rmse = target,
              camera_level = best$camera_level, inat_level = best$inat_level,
              achieved_rmse = best$rmse, note = "", stringsAsFactors = FALSE)
  }))
}))
write.csv(target_table, "effort_sweep_target_rmse_table.csv", row.names = FALSE)
cat("\n==================== EFFORT NEEDED FOR TARGET RMSE (min-cost qualifying cell) ====================\n")
print(target_table, row.names = FALSE, digits = 3)

# ------------------------------ figures -----------------------------------------
p_cam <- ggplot(agg_region, aes(x = camera_level, y = rmse, color = ecoregion, group = ecoregion)) +
  geom_point(size = 2) + geom_line() +
  facet_wrap(~ inat_level, labeller = label_both) +
  labs(x = "camera level (n_site_keep multiplier of the 700-site baseline)",
      y = "per-region RMSE (trend deviation)",
      title = "Per-region trend RMSE vs camera effort, faceted by iNat effort level",
      subtitle = sprintf("bobcat abundance fixed at real magnitude; varying-trend scenario. Divergence-filtered (%d/480 replicates excluded, sigma_region_mean>1.0 -- see console output for per-cell rates).", n_divergent)) +
  theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7)) +
  guides(color = guide_legend(nrow = 3))
ggsave("effort_sweep_rmse_vs_cameras.png", p_cam, width = 11, height = 8, dpi = 150, bg = "white")

p_inat <- ggplot(agg_region, aes(x = mean_inat_count_region, y = rmse, color = ecoregion, group = ecoregion)) +
  geom_point(size = 2) + geom_line() +
  facet_wrap(~ camera_level, labeller = label_both, scales = "free_x") +
  scale_x_log10() +
  labs(x = "realized iNat count (this region, mean over reps, log scale)",
      y = "per-region RMSE (trend deviation)",
      title = "Per-region trend RMSE vs realized iNat effort, faceted by camera level",
      subtitle = sprintf("bobcat abundance fixed at real magnitude; varying-trend scenario. Divergence-filtered (%d/480 replicates excluded, sigma_region_mean>1.0 -- see console output for per-cell rates).", n_divergent)) +
  theme_bw() + theme(legend.position = "bottom", legend.text = element_text(size = 7)) +
  guides(color = guide_legend(nrow = 3))
ggsave("effort_sweep_rmse_vs_inat.png", p_inat, width = 11, height = 8, dpi = 150, bg = "white")

cat("\nwrote effort_sweep_rmse_vs_cameras.png, effort_sweep_rmse_vs_inat.png\n")
