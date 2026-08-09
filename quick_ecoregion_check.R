#!/usr/bin/env Rscript
# quick_ecoregion_check.R -- cheap, nimble-free sanity check of the
# reconciled ecoregion interfaces before spending a real MCMC smoke test.
# Not a deliverable.
PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
setwd(PROTO_DIR)
source("sim_helpers.R")

prepped <- readRDS("prepped_sim_inputs.RDS")
cat("names(prepped):", paste(names(prepped), collapse=", "), "\n")
stopifnot(all(c("ecoregion_of_cell100","nregion","ecoregion_levels",
               "inat_effort_real","real_y_template","base_inits") %in% names(prepped)))
stopifnot("ecoregion_id" %in% names(prepped$cell100_geo))
cat("nregion:", prepped$nregion, " ecoregion_levels:", paste(prepped$ecoregion_levels, collapse=", "), "\n\n")

cat("=== make_true_year_region() ===\n")
yr_varying <- make_true_year_region(prepped$cell100_geo, scenario = "varying")
yr_null    <- make_true_year_region(prepped$cell100_geo, scenario = "null")
cat("varying:", paste(round(yr_varying, 3), collapse=", "), " (expect spread, mean~0)\n")
cat("null:   ", paste(yr_null, collapse=", "), " (expect all 0)\n")
stopifnot(length(yr_varying) == prepped$nregion, length(yr_null) == prepped$nregion,
         all(yr_null == 0), abs(mean(yr_varying)) < 1e-10)

cat("\n=== true_param_list_ecoregion() ===\n")
truth <- true_param_list_ecoregion(prepped$real_post_means, yr_varying)
cat("names(truth):", paste(names(truth), collapse=", "), "\n")
stopifnot("year_region" %in% names(truth), "sigma_region" %in% names(truth),
         "link_occ_intercept" %in% names(truth), "theta0" %in% names(truth))
cat("year_region matches input:", isTRUE(all.equal(truth$year_region, yr_varying)), "\n")

cat("\n=== build_reduced_constants(stratify_by='ecoregion', ...) ===\n")
cl <- prepped$constants_list
reduced <- build_reduced_constants(
  cl = cl, inat_effort_real = prepped$inat_effort_real, real_y_template = prepped$real_y_template,
  cell100_geo = prepped$cell100_geo,
  n_cell50_keep = 30, max_subcell_per_cell = 5, n_site_keep = 60, seed = 1,
  stratify_by = "ecoregion",
  ecoregion_of_cell100_full = prepped$ecoregion_of_cell100, nregion = prepped$nregion)
cat("ncell50=", reduced$constants$ncell50, " nsite=", reduced$constants$nsite,
    " has ecoregion_of_cell100:", !is.null(reduced$constants$ecoregion_of_cell100),
    " nregion in constants:", reduced$constants$nregion, "\n")
stopifnot(!is.null(reduced$constants$ecoregion_of_cell100), reduced$constants$nregion == prepped$nregion,
         length(reduced$constants$ecoregion_of_cell100) == cl$ncell100)

cat("\n=== scale_truth_abundance() ===\n")
source("sim_helpers_abundance.R")
truth_scaled <- scale_truth_abundance(truth, occ_shift = 0.75, count_log_mult = log(3), label = "moderate")
cat("link_occ_intercept shift check:",
   isTRUE(all.equal(truth_scaled$link_occ_intercept, truth$link_occ_intercept + 0.75)), "\n")
cat("theta0 shift check:", isTRUE(all.equal(truth_scaled$theta0, truth$theta0 + log(3))), "\n")

cat("\nQUICK ECOREGION CHECK DONE\n")
