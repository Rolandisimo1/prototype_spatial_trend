#!/usr/bin/env Rscript
# smoke_test.R -- tiny end-to-end run (2 reps, small reduced design, short
# MCMC) to catch API/dimension errors before committing to the real
# 50-replicate array via 01_run_sim_validation.R. Not a deliverable, just a
# pre-flight check.
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

cat("=== building TINY reduced design ===\n")
reduced <- build_reduced_constants(
  cl = cl, inat_effort_real = inp$inat_effort, real_y_template = inp$real_data$y,
  cell100_geo = prepped$cell100_geo,
  n_cell50_keep = 30, max_subcell_per_cell = 5, n_site_keep = 60,
  seed = 1)
cat("ncell50=", reduced$constants$ncell50, " nsite=", reduced$constants$nsite,
    " subcell rows=", dim(reduced$constants$xdat_inat)[1], "\n")

year_effect_true <- make_true_year_effect(prepped$cell100_geo, cl$adj, cl$num, amplitude = 0.3)
truth <- true_param_list(prepped$real_post_means, year_effect_true)
informed_cell100 <- unique(reduced$constants$inat_cell100)

cat("\n=== rep 1 ===\n")
t0 <- Sys.time()
r1 <- run_one_replicate(1, model_code_spatial_trend, reduced$constants, reduced$inat_effort,
                        reduced$y_ncol, truth, year_effect_true, prepped$cell100_geo,
                        informed_cell100, inp$inits_list, n_burnin = 100, n_iter = 200)
cat("elapsed:", as.numeric(difftime(Sys.time(), t0, units="secs")), "sec\n")
print(r1)

cat("\n=== rep 2 ===\n")
r2 <- run_one_replicate(2, model_code_spatial_trend, reduced$constants, reduced$inat_effort,
                        reduced$y_ncol, truth, year_effect_true, prepped$cell100_geo,
                        informed_cell100, inp$inits_list, n_burnin = 100, n_iter = 200)
print(r2)

cat("\nSMOKE TEST DONE\n")
