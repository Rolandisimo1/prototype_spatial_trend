#!/usr/bin/env Rscript
# =============================================================================
# 02_collect_sim_summary.R
# Run AFTER the 01_run_sim_validation.R array job completes. Aggregates
# sim_results/rep_*.RDS into sim_summary.csv and prints the bias/RMSE/
# coverage/sign-recovery verdict that determines whether Step 2 (real
# bobcat refit) is warranted.
#
# RUN: $PROJ/HPC/conda_envs/nimble_env/bin/Rscript 02_collect_sim_summary.R
# =============================================================================

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
RESULTS_DIR <- paste0(PROTO_DIR, "/sim_results")

files <- list.files(RESULTS_DIR, pattern = "^rep_.*\\.RDS$", full.names = TRUE)
if (length(files) == 0) stop("No replicate result files found under ", RESULTS_DIR)
cat("found", length(files), "replicate result files\n")

results_list <- lapply(files, readRDS)
all_cols <- unique(unlist(lapply(results_list, names)))
results_df <- do.call(rbind, lapply(results_list, function(x) {
  missing <- setdiff(all_cols, names(x))
  for (m in missing) x[[m]] <- NA
  x[all_cols]
}))
results_df <- results_df[order(results_df$rep), ]

write.csv(results_df, paste0(PROTO_DIR, "/sim_summary.csv"), row.names = FALSE)

cat("\n==================== SIM VALIDATION SUMMARY ====================\n")
cat("status counts:\n"); print(table(results_df$status))

ok <- results_df[results_df$status == "ok", ]
if (nrow(ok) > 0) {
  cat("\nmean over", nrow(ok), "successful replicates (of", nrow(results_df), "total):\n")
  cat("  bias (all informed cells):    ", round(mean(ok$bias_all, na.rm=TRUE), 4), "\n")
  cat("  RMSE (all informed cells):    ", round(mean(ok$rmse_all, na.rm=TRUE), 4), "\n")
  cat("  coverage (all informed, 95%): ", round(mean(ok$coverage_all, na.rm=TRUE), 3), "\n")
  cat("  bias / RMSE / coverage, NE:   ", round(mean(ok$bias_NE, na.rm=TRUE), 4), "/",
      round(mean(ok$rmse_NE, na.rm=TRUE), 4), "/", round(mean(ok$coverage_NE, na.rm=TRUE), 3), "\n")
  cat("  bias / RMSE / coverage, West: ", round(mean(ok$bias_West, na.rm=TRUE), 4), "/",
      round(mean(ok$rmse_West, na.rm=TRUE), 4), "/", round(mean(ok$coverage_West, na.rm=TRUE), 3), "\n")
  cat("  bias / RMSE / coverage, other:", round(mean(ok$bias_other, na.rm=TRUE), 4), "/",
      round(mean(ok$rmse_other, na.rm=TRUE), 4), "/", round(mean(ok$coverage_other, na.rm=TRUE), 3), "\n")
  cat("  NE-vs-West sign recovered:    ", sum(ok$sign_recovered), "/", nrow(ok),
      " replicates (", round(100 * mean(ok$sign_recovered), 1), "%)\n")
} else {
  cat("\nNO successful replicates -- see sim_summary.csv 'detail' column for errors.\n")
}

not_ok <- results_df[results_df$status != "ok", ]
if (nrow(not_ok) > 0) {
  cat("\n--- non-ok replicates ---\n")
  print(not_ok[, intersect(c("rep","status","detail"), names(not_ok))])
}

cat("\nwrote:", paste0(PROTO_DIR, "/sim_summary.csv"), "\n")
