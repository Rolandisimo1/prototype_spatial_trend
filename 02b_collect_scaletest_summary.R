#!/usr/bin/env Rscript
# =============================================================================
# 02b_collect_scaletest_summary.R
# Aggregates sim_results_scaletest/rep_*.RDS (the 3x-larger confirmatory
# design) and prints the same bias/RMSE/coverage/sign-recovery summary as
# 02_collect_sim_summary.R, side by side against the original 50-rep result,
# to check whether bias shrank as the data-thinning-artifact hypothesis predicts.
# =============================================================================

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")

summarize_dir <- function(results_dir, label) {
  files <- list.files(results_dir, pattern = "^rep_.*\\.RDS$", full.names = TRUE)
  if (length(files) == 0) { cat(label, ": no result files found\n"); return(invisible(NULL)) }
  results_list <- lapply(files, readRDS)
  all_cols <- unique(unlist(lapply(results_list, names)))
  df <- do.call(rbind, lapply(results_list, function(x) {
    missing <- setdiff(all_cols, names(x)); for (m in missing) x[[m]] <- NA; x[all_cols]
  }))
  ok <- df[df$status == "ok", ]
  cat("\n===", label, "(", nrow(ok), "/", nrow(df), "ok ) ===\n")
  if (nrow(ok) == 0) return(invisible(df))
  cat("  bias/RMSE/coverage, all:  ", round(mean(ok$bias_all,na.rm=TRUE),4), "/",
      round(mean(ok$rmse_all,na.rm=TRUE),4), "/", round(mean(ok$coverage_all,na.rm=TRUE),3), "\n")
  cat("  bias/RMSE/coverage, NE:   ", round(mean(ok$bias_NE,na.rm=TRUE),4), "/",
      round(mean(ok$rmse_NE,na.rm=TRUE),4), "/", round(mean(ok$coverage_NE,na.rm=TRUE),3), "\n")
  cat("  bias/RMSE/coverage, West: ", round(mean(ok$bias_West,na.rm=TRUE),4), "/",
      round(mean(ok$rmse_West,na.rm=TRUE),4), "/", round(mean(ok$coverage_West,na.rm=TRUE),3), "\n")
  cat("  sign recovered:           ", sum(ok$sign_recovered), "/", nrow(ok),
      "(", round(100*mean(ok$sign_recovered),1), "%)\n")
  cat("  n_informed cells:         ", ok$n_informed[1], "\n")
  df
}

df_orig  <- summarize_dir(paste0(PROTO_DIR, "/sim_results"), "ORIGINAL (210 cell50 / 700 sites, 178 informed cells, 50 reps)")
df_scale <- summarize_dir(paste0(PROTO_DIR, "/sim_results_scaletest"), "SCALE-TEST (630 cell50 / 2000 sites, 20 reps)")

if (!is.null(df_scale)) write.csv(df_scale, paste0(PROTO_DIR, "/sim_summary_scaletest.csv"), row.names = FALSE)
cat("\nwrote:", paste0(PROTO_DIR, "/sim_summary_scaletest.csv"), "\n")
