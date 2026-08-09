#!/usr/bin/env Rscript
# =============================================================================
# 02h_collect_effort_secondary.R
# Collects sim_results_effort_secondary/rep_*.RDS (from
# 01h_run_effort_secondary.R): one ecoregion's iNat effort held at 0.1x
# while every other region (including the comparison region) stays at 1x,
# cameras fixed at 1x throughout. Compares FOCAL vs COMPARISON region RMSE
# from the SAME 30 replicates (same-run contrast, not cross-run).
#
# RUN (plotting_env):
#   $PROJ/HPC/conda_envs/plotting_env/bin/Rscript 02h_collect_effort_secondary.R
# Writes: effort_secondary_summary.csv
# =============================================================================

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
setwd(PROTO_DIR)

prepped <- readRDS("prepped_sim_inputs.RDS")
ecoregion_levels <- prepped$ecoregion_levels
region_names_safe <- make.names(ecoregion_levels)

RESULTS_DIR <- paste0(PROTO_DIR, "/sim_results_effort_secondary")
files <- list.files(RESULTS_DIR, pattern = "^rep_.*\\.RDS$", full.names = TRUE)
if (length(files) == 0) stop("No sim_results_effort_secondary/rep_*.RDS found -- has it run yet?")
cat("found", length(files), "replicate result files\n")

xs <- lapply(files, readRDS)
statuses <- sapply(xs, function(x) x$status)
cat("status counts:\n"); print(table(statuses))
ok <- xs[statuses == "ok"]
if (length(ok) == 0) stop("No ok replicates -- nothing to collect.")

FOCAL_ECOREGION_NAME      <- ok[[1]]$focal_ecoregion
COMPARISON_ECOREGION_NAME <- ok[[1]]$comparison_ecoregion
FOCAL_MULT                <- ok[[1]]$focal_mult
cat(sprintf("FOCAL: %s (x%.2f)  COMPARISON: %s (x1)\n",
            FOCAL_ECOREGION_NAME, FOCAL_MULT, COMPARISON_ECOREGION_NAME))

# ------------------------------ guard: refuse degenerate replicates --------------
rmse_all_vec <- sapply(ok, function(x) x$rmse_all)
if (length(ok) > 1 && length(unique(rmse_all_vec)) == 1) {
  stop("DEGENERATE REPLICATES DETECTED in sim_results_effort_secondary -- all ",
      length(ok), " ok replicates have identical rmse_all (", rmse_all_vec[1], "). ",
      "Every 'replicate' likely simulated and fit the SAME dataset (missing per-rep ",
      "reseed after build_reduced_constants() -- see abundance_sweep_seed_diagnosis.md). ",
      "Fix 01h_run_effort_secondary.R and re-run before collecting.")
}
cat("guard passed: rmse_all not degenerate across", length(ok), "ok replicates\n\n")

extract_region <- function(x, region_name) {
  rn <- make.names(region_name)
  data.frame(
    rep = x$rep,
    ecoregion = region_name,
    est = x[[paste0("est_", rn)]],
    err = x[[paste0("err_", rn)]],
    covered = x[[paste0("cov_", rn)]],
    stringsAsFactors = FALSE
  )
}

focal_rows      <- do.call(rbind, lapply(ok, extract_region, region_name = FOCAL_ECOREGION_NAME))
comparison_rows <- do.call(rbind, lapply(ok, extract_region, region_name = COMPARISON_ECOREGION_NAME))

info_rows <- do.call(rbind, lapply(ok, function(x) {
  info <- x$info_by_region[[1]]
  rbind(
    data.frame(rep = x$rep, ecoregion = FOCAL_ECOREGION_NAME,
              inat_count = info$inat_count[info$ecoregion == FOCAL_ECOREGION_NAME],
              n_cell50 = info$n_cell50[info$ecoregion == FOCAL_ECOREGION_NAME]),
    data.frame(rep = x$rep, ecoregion = COMPARISON_ECOREGION_NAME,
              inat_count = info$inat_count[info$ecoregion == COMPARISON_ECOREGION_NAME],
              n_cell50 = info$n_cell50[info$ecoregion == COMPARISON_ECOREGION_NAME])
  )
}))

summary_tbl <- do.call(rbind, lapply(list(
  list(name = FOCAL_ECOREGION_NAME, rows = focal_rows, role = sprintf("FOCAL (iNat x%.2f)", FOCAL_MULT)),
  list(name = COMPARISON_ECOREGION_NAME, rows = comparison_rows, role = "COMPARISON (iNat x1)")
), function(g) {
  info_g <- info_rows[info_rows$ecoregion == g$name, ]
  data.frame(
    ecoregion = g$name, role = g$role, n_reps = nrow(g$rows),
    rmse = sqrt(mean(g$rows$err^2, na.rm = TRUE)),
    bias = mean(g$rows$err, na.rm = TRUE),
    coverage = mean(g$rows$covered, na.rm = TRUE),
    mean_inat_count = mean(info_g$inat_count, na.rm = TRUE),
    mean_n_cell50 = mean(info_g$n_cell50, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}))

write.csv(summary_tbl, "effort_secondary_summary.csv", row.names = FALSE)
cat("\n==================== SECONDARY EFFORT SCENARIO: FOCAL vs COMPARISON ====================\n")
print(summary_tbl, row.names = FALSE, digits = 3)
cat("\nwrote effort_secondary_summary.csv\n")
