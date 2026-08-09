#!/usr/bin/env Rscript
# ecoregion_smoke_test.R -- GATE 1 (per task spec): confirm the ecoregion
# model compiles and samples year_region/sigma_region without dimension
# errors, exercising run_one_replicate(fit_null_scalar=TRUE) (both the
# ecoregion model and model_code_national_scalar) on tiny fake-scale data.
# Not a deliverable.
suppressPackageStartupMessages({
  library(nimble); library(nimbleEcology); library(coda)
})

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
setwd(PROTO_DIR)

source(paste0(PROJ, "/HPC/bobcat/integration_helper.R"))
source("model_code_ecoregion_trend.R")
source("model_code_national_scalar.R")
source("sim_helpers.R")

prepped <- readRDS("prepped_sim_inputs.RDS")
cl <- prepped$constants_list

cat("=== building TINY reduced design (ecoregion-stratified) ===\n")
reduced <- build_reduced_constants(
  cl = cl, inat_effort_real = prepped$inat_effort_real, real_y_template = prepped$real_y_template,
  cell100_geo = prepped$cell100_geo,
  n_cell50_keep = 30, max_subcell_per_cell = 5, n_site_keep = 60, seed = 1,
  stratify_by = "ecoregion",
  ecoregion_of_cell100_full = prepped$ecoregion_of_cell100, nregion = prepped$nregion)
cat("ncell50=", reduced$constants$ncell50, " nsite=", reduced$constants$nsite, "\n")

year_region_true <- make_true_year_region(prepped$cell100_geo, scenario = "varying")
truth <- true_param_list_ecoregion(prepped$real_post_means, year_region_true)

cat("\n=== rep 1 (varying scenario, fit_null_scalar=TRUE) ===\n")
t0 <- Sys.time()
r1 <- run_one_replicate(
  rep_id = 1, model_code = model_code_ecoregion,
  constants = reduced$constants, inat_effort = reduced$inat_effort, y_ncol = reduced$y_ncol,
  truth = truth, base_inits = prepped$base_inits, n_burnin = 100, n_iter = 200,
  metrics_fn = compute_ecoregion_metrics,
  trend_inits = list(year_beta = 0, year_var = 0, sigma_region = 1,
                    year_region = rep(0, prepped$nregion)),
  extra_monitors = c("year_region", "sigma_region"),
  fit_null_scalar = TRUE,
  year_region_true = year_region_true, ecoregion_levels = prepped$ecoregion_levels)
cat("elapsed:", as.numeric(difftime(Sys.time(), t0, units="secs")), "sec\n")
print(r1[, setdiff(names(r1), "sim_data")])

cat("\n=== rep 1 (null scenario) ===\n")
year_region_null <- make_true_year_region(prepped$cell100_geo, scenario = "null")
truth_null <- true_param_list_ecoregion(prepped$real_post_means, year_region_null)
r2 <- run_one_replicate(
  rep_id = 1, model_code = model_code_ecoregion,
  constants = reduced$constants, inat_effort = reduced$inat_effort, y_ncol = reduced$y_ncol,
  truth = truth_null, base_inits = prepped$base_inits, n_burnin = 100, n_iter = 200,
  metrics_fn = compute_ecoregion_metrics,
  trend_inits = list(year_beta = 0, year_var = 0, sigma_region = 1,
                    year_region = rep(0, prepped$nregion)),
  extra_monitors = c("year_region", "sigma_region"),
  fit_null_scalar = TRUE,
  year_region_true = year_region_null, ecoregion_levels = prepped$ecoregion_levels)
print(r2[, setdiff(names(r2), "sim_data")])

cat("\n=== abundance scaling check (via run_one_replicate, moderate abundance) ===\n")
source("sim_helpers_abundance.R")
truth_abund <- scale_truth_abundance(truth, occ_shift = 0.75, count_log_mult = log(3), label = "moderate")
r3 <- run_one_replicate(
  rep_id = 1, model_code = model_code_ecoregion,
  constants = reduced$constants, inat_effort = reduced$inat_effort, y_ncol = reduced$y_ncol,
  truth = truth_abund, base_inits = prepped$base_inits, n_burnin = 100, n_iter = 200,
  metrics_fn = compute_ecoregion_metrics,
  trend_inits = list(year_beta = 0, year_var = 0, sigma_region = 1,
                    year_region = rep(0, prepped$nregion)),
  extra_monitors = c("year_region", "sigma_region"),
  fit_null_scalar = TRUE,
  year_region_true = year_region_true, ecoregion_levels = prepped$ecoregion_levels)
info <- summarize_simulated_information(r3$sim_data[[1]], reduced$constants)
cat("status:", r3$status, "\n")
print(info)

cat("\nSMOKE TEST DONE\n")
