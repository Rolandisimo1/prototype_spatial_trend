#!/usr/bin/env Rscript
# psock_smoke_test.R -- verify the PSOCK cluster fix actually avoids the
# mclapply fork-collision seen in job 436025, before committing to the full
# ~65-min 50-replicate run. 4 reps on 4 workers, short MCMC budget -- only
# checking that parallel compilation succeeds, not inference quality.
# run_rep_from_scratch here must stay IDENTICAL to the one in
# 01_run_sim_validation.R. Not a deliverable, just a pre-flight check.
suppressPackageStartupMessages(library(parallel))

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
setwd(PROTO_DIR)

run_rep_from_scratch <- function(rep_id, proj, proto_dir, n_burnin, n_iter, sim_seed) {
  suppressPackageStartupMessages({
    library(nimble); library(nimbleEcology); library(coda)
  })
  setwd(proto_dir)
  source(paste0(proj, "/HPC/bobcat/integration_helper.R"))
  source("model_code_spatial_trend.R")
  source("sim_helpers.R")

  prepped <- readRDS("prepped_sim_inputs.RDS")
  inp <- readRDS(paste0(proj, "/HPC/bobcat/input_data_bobcat.RDS"))
  cl  <- prepped$constants_list

  reduced <- build_reduced_constants(
    cl = cl, inat_effort_real = inp$inat_effort, real_y_template = inp$real_data$y,
    cell100_geo = prepped$cell100_geo,
    n_cell50_keep = 30, max_subcell_per_cell = 5, n_site_keep = 60,  # tiny -- speed only
    seed = sim_seed)

  year_effect_true <- make_true_year_effect(prepped$cell100_geo, cl$adj, cl$num, amplitude = 0.3)
  truth <- true_param_list(prepped$real_post_means, year_effect_true)
  informed_cell100 <- unique(reduced$constants$inat_cell100)

  res <- run_one_replicate(
    rep_id = rep_id, model_code = model_code_spatial_trend,
    constants = reduced$constants, inat_effort = reduced$inat_effort,
    y_ncol = reduced$y_ncol, truth = truth, year_effect_true = year_effect_true,
    cell100_geo = prepped$cell100_geo, informed_cell100 = informed_cell100,
    base_inits = inp$inits_list, n_burnin = n_burnin, n_iter = n_iter)
  res
}

cat("=== PSOCK smoke test: 4 reps, 4 workers, tiny budget ===\n")
cl_psock <- makeCluster(4, type = "PSOCK")
clusterSetRNGStream(cl_psock, 1)
t0 <- Sys.time()
results <- parLapply(cl_psock, 1:4, run_rep_from_scratch,
                     proj = PROJ, proto_dir = PROTO_DIR,
                     n_burnin = 100, n_iter = 200, sim_seed = 1)
stopCluster(cl_psock)
cat("elapsed:", round(as.numeric(difftime(Sys.time(), t0, units="secs")), 1), "sec\n\n")

for (r in results) print(r[, c("rep","status")])
statuses <- sapply(results, function(x) x$status)
cat("\nstatuses:", paste(statuses, collapse=", "), "\n")
if (all(statuses == "ok")) {
  cat("PSOCK SMOKE TEST: ALL OK\n")
} else {
  cat("PSOCK SMOKE TEST: FAILURES PRESENT\n")
  for (r in results) if (r$status != "ok") { cat("---\n"); print(r) }
}
