#!/usr/bin/env Rscript
# verify_reseed.R -- cheap check (no MCMC/compile) that rep-specific
# reseeding actually produces different simulated data across replicates.
# Not a deliverable.
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
SIM_SEED <- 20260712

sim_for_rep <- function(rep_id) {
  reduced <- build_reduced_constants(
    cl = cl, inat_effort_real = inp$inat_effort, real_y_template = inp$real_data$y,
    cell100_geo = prepped$cell100_geo,
    n_cell50_keep = 210, max_subcell_per_cell = 15, n_site_keep = 700, seed = SIM_SEED)
  year_effect_true <- make_true_year_effect(prepped$cell100_geo, cl$adj, cl$num, amplitude = 0.3)
  truth <- true_param_list(prepped$real_post_means, year_effect_true)
  set.seed(SIM_SEED + rep_id)
  simulate_replicate_data(model_code_spatial_trend, reduced$constants, reduced$inat_effort,
                          reduced$y_ncol, truth)
}

s1 <- sim_for_rep(1)
s2 <- sim_for_rep(2)
s1b <- sim_for_rep(1)  # same rep_id should reproduce identically

cat("sum(y_inat) rep1:", sum(s1$y_inat, na.rm=TRUE), "\n")
cat("sum(y_inat) rep2:", sum(s2$y_inat, na.rm=TRUE), "\n")
cat("sum(y_inat) rep1 (rerun):", sum(s1b$y_inat, na.rm=TRUE), "\n")
cat("rep1 == rep2 (should be FALSE):", isTRUE(all.equal(s1$y_inat, s2$y_inat)), "\n")
cat("rep1 == rep1-rerun (should be TRUE, reproducible):", isTRUE(all.equal(s1$y_inat, s1b$y_inat)), "\n")
