#!/usr/bin/env Rscript
# debug_compile.R -- isolate the compileNimble failure seen in smoke_test.R
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

reduced <- build_reduced_constants(
  cl = cl, inat_effort_real = inp$inat_effort, real_y_template = inp$real_data$y,
  cell100_geo = prepped$cell100_geo,
  n_cell50_keep = 30, max_subcell_per_cell = 5, n_site_keep = 60, seed = 1)

year_effect_true <- make_true_year_effect(prepped$cell100_geo, cl$adj, cl$num, amplitude = 0.3)
truth <- true_param_list(prepped$real_post_means, year_effect_true)

cat("=== building truth model + simulating ===\n")
sim_data <- simulate_replicate_data(model_code_spatial_trend, reduced$constants,
                                    reduced$inat_effort, reduced$y_ncol, truth)
cat("simulated y dim:", paste(dim(sim_data$y), collapse="x"),
    " y_inat dim:", paste(dim(sim_data$y_inat), collapse="x"), "\n")
cat("y_inat range:", paste(range(sim_data$y_inat, na.rm=TRUE), collapse=" to "), "\n")
cat("y range:", paste(range(sim_data$y, na.rm=TRUE), collapse=" to "), "\n")

cat("\n=== building fit_model ===\n")
fit_model <- nimbleModel(
  model_code_spatial_trend, constants = reduced$constants,
  data = list(y = sim_data$y, y_inat = sim_data$y_inat, inat_effort = reduced$inat_effort),
  inits = list(year_beta = 0, year_var = 0, tau_year = 1,
              year_effect = rep(0, reduced$constants$ncell100)),
  calculate = FALSE)

cat("\n=== compiling (catching error) ===\n")
Cmodel <- tryCatch(compileNimble(fit_model), error = function(e) {
  cat("COMPILE ERROR:", conditionMessage(e), "\n")
  NULL
})

if (is.null(Cmodel)) {
  cat("\n=== printErrors() ===\n")
  print(printErrors())
} else {
  cat("\nCOMPILE SUCCEEDED\n")
}
