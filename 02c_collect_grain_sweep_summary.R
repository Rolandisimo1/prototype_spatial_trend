#!/usr/bin/env Rscript
# =============================================================================
# 02c_collect_grain_sweep_summary.R
# Combines the 3 NEW grain-sweep levels (sim_results_grain/grain_{0,1,8}/)
# with the EXISTING original 50-replicate run (sim_results/, which used the
# exact same smoothing recipe as n_smooth_iter=3 -- see
# make_true_year_effect_grain()'s docstring) into a single grain-vs-recovery
# table, without re-simulating the grain=3 point.
#
# NOTE: the original run's saved rep_*.RDS files only contain computed
# metrics (bias_all/rmse_all/coverage_all), not the raw posterior samples --
# those were never written to disk (only held in-memory during that job) --
# so spatial_cor (which needs the raw per-cell posterior means) cannot be
# retroactively computed for the grain=3 point. bias_all/rmse_all/coverage_all
# ARE directly comparable across both metric functions (identical formulas);
# only spatial_cor is missing for that one point, marked NA below.
# =============================================================================

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")

read_results <- function(dir, grain_value = NA, label = NA) {
  files <- list.files(dir, pattern = "^rep_.*\\.RDS$", full.names = TRUE)
  if (length(files) == 0) return(NULL)
  rows <- lapply(files, function(f) {
    x <- readRDS(f)
    data.frame(
      rep = x$rep, status = x$status,
      grain = if (!is.null(x$grain)) x$grain else grain_value,
      label = label,
      bias_all = if (!is.null(x$bias_all)) x$bias_all else NA,
      rmse_all = if (!is.null(x$rmse_all)) x$rmse_all else NA,
      coverage_all = if (!is.null(x$coverage_all)) x$coverage_all else NA,
      n_informed = if (!is.null(x$n_informed)) x$n_informed else NA,
      spatial_cor = if (!is.null(x$spatial_cor)) x$spatial_cor else NA_real_
    )
  })
  do.call(rbind, rows)
}

new_grain <- do.call(rbind, lapply(c(0, 1, 8), function(g) {
  read_results(paste0(PROTO_DIR, "/sim_results_grain/grain_", g), grain_value = g,
              label = paste0("grain=", g, " (noise-seeded)"))
}))

# TWO distinct grain=3 data points -- NOT the same truth field, both kept
# separate rather than merged:
#   original: make_true_year_effect() -- regionally-structured NE+/West- seed
#     (the already-validated 50-rep design, job 439834). No spatial_cor
#     available (posterior samples weren't saved).
#   supplement: make_true_year_effect_grain(n_smooth_iter=3) -- white-noise
#     seed, same smoothing recipe. Run specifically to get a spatial_cor
#     value directly comparable to grain 0/1/8 (job 443144).
original <- read_results(paste0(PROTO_DIR, "/sim_results"), grain_value = 3,
                         label = "grain=3 (original, regional-structure truth)")
supplement <- read_results(paste0(PROTO_DIR, "/sim_results_grain/grain_3"), grain_value = 3,
                           label = "grain=3 (noise-seeded, supplement)")

combined <- rbind(new_grain, original, supplement)
combined <- combined[combined$status == "ok", ]
write.csv(combined, paste0(PROTO_DIR, "/sim_summary_grain_sweep.csv"), row.names = FALSE)

cat("==================== GRAIN SWEEP SUMMARY ====================\n")
cat("(grain = smoothing passes on the truth field; higher = broader/smoother,\n")
cat(" 0 = finest/patchiest. Two grain=3 rows are kept separate -- different\n")
cat(" specific truth fields, see file header.)\n\n")

# na.action = na.pass: aggregate()'s formula interface defaults to na.omit,
# which drops the ENTIRE ROW if ANY response column is NA -- since spatial_cor
# is NA for every original-grain=3 row (no saved posterior samples to compute
# it retroactively, see file header), that default silently dropped that
# whole group before aggregation ran, even though its bias/RMSE/coverage
# were all valid. na.pass keeps every row; FUN's na.rm=TRUE handles the NAs
# per-column instead.
agg <- aggregate(cbind(bias_all, rmse_all, coverage_all, spatial_cor) ~ label + grain,
                 data = combined, FUN = function(x) mean(x, na.rm = TRUE),
                 na.action = na.pass)
n_reps <- aggregate(rep ~ label, data = combined, FUN = length)
names(n_reps)[2] <- "n_reps"
agg <- merge(agg, n_reps, by = "label")
agg <- agg[order(agg$grain), ]
agg[, c("bias_all","rmse_all","coverage_all","spatial_cor")] <-
  round(agg[, c("bias_all","rmse_all","coverage_all","spatial_cor")], 4)

print(agg[, c("label","bias_all","rmse_all","coverage_all","spatial_cor","n_reps")], row.names = FALSE)

cat("\nwrote:", paste0(PROTO_DIR, "/sim_summary_grain_sweep.csv"), "\n")
