#!/usr/bin/env Rscript
# dry_run.R -- ONE replicate at the real planned design size/budget (not the
# tiny smoke-test scale), to get an actual per-replicate runtime for capacity
# planning of the 50-replicate array, and a more meaningful (though still
# N=1) bias/coverage/sign sanity check. Not a deliverable itself.
suppressPackageStartupMessages({
  library(nimble); library(nimbleEcology); library(coda)
})

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
setwd(PROTO_DIR)

source(paste0(PROJ, "/HPC/bobcat/integration_helper.R"))
source("model_code_spatial_trend.R")
source("sim_helpers.R")

prepped <- readRDS("prepped_sim_inputs.RDS")
inp <- readRDS(paste0(PROJ, "/HPC/bobcat/input_data_bobcat.RDS"))
cl  <- prepped$constants_list

cat("=== building REAL-SCALE reduced design ===\n")
t_build <- Sys.time()
reduced <- build_reduced_constants(
  cl = cl, inat_effort_real = inp$inat_effort, real_y_template = inp$real_data$y,
  cell100_geo = prepped$cell100_geo,
  n_cell50_keep = 210, max_subcell_per_cell = 15, n_site_keep = 700, seed = 20260712)
cat("ncell50=", reduced$constants$ncell50, " nsite=", reduced$constants$nsite,
    " subcell rows=", dim(reduced$constants$xdat_inat)[1],
    " (build took", round(as.numeric(difftime(Sys.time(), t_build, units="secs")),1), "sec)\n")

year_effect_true <- make_true_year_effect(prepped$cell100_geo, cl$adj, cl$num, amplitude = 0.3)
truth <- true_param_list(prepped$real_post_means, year_effect_true)
informed_cell100 <- unique(reduced$constants$inat_cell100)
cat("informed cell100 cells:", length(informed_cell100), "/", cl$ncell100, "\n")

cat("\n=== dry-run replicate (n_burnin=1000, n_iter=3000) ===\n")
t0 <- Sys.time()
r1 <- run_one_replicate(1, model_code_spatial_trend, reduced$constants, reduced$inat_effort,
                        reduced$y_ncol, truth, year_effect_true, prepped$cell100_geo,
                        informed_cell100, inp$inits_list, n_burnin = 1000, n_iter = 3000)
cat("\ntotal elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "sec\n")
print(r1)
cat("\nDRY RUN DONE\n")
