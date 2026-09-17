#!/usr/bin/env Rscript
# =============================================================================
# HPC_run_model_single_camcount.R -- single-shot chain for the COUNT camera arm.
#
# Fork of HPC/white-tailed_deer_v2b_national_scalar/HPC_run_model_single.R
# (the occupancy arm it is compared against). Identical: nimbleModel ->
# compile -> configureMCMC(enableWAIC = TRUE) -> ONE Cmcmc$run(niter + nburnin,
# nburnin), seed set.seed(1000 + chain_id), NITER 50000 / NBURNIN 5000, same
# monitor list. Different, and only these:
#   - model: model_code_national_scalar_camcount, sourced from the copy of
#     model_code_national_scalar_camcount.R in this directory (md5 logged)
#   - data: y_cam replaces y; constants add log_effort_cam / yday_mean
#   - monitors: + overdisp_cam, overdisp_inat (both named explicitly)
#   - data-shape checks run BEFORE nimbleModel(), so a wiring error fails in
#     seconds instead of after the multi-hour graph build
# ENV: SPECIES (run dir token), CHAIN_ID, NITER, NBURNIN.
# =============================================================================
library(nimbleEcology)
library(tidyverse)
library(coda)

PROJ     <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
species  <- Sys.getenv("SPECIES")
chain_id <- as.integer(Sys.getenv("CHAIN_ID"))
NITER    <- as.integer(Sys.getenv("NITER",   "50000"))
NBURNIN  <- as.integer(Sys.getenv("NBURNIN", "5000"))
stopifnot(nzchar(species), !is.na(chain_id), chain_id %in% 1:3)

setwd(file.path(PROJ, "HPC", species))
project_dir <- getwd()
if (file.exists("BUNDLE_REVIEW_REQUIRED")) stop("BUNDLE_REVIEW_REQUIRED present -- bundle not signed off")
source("pipeline_helper.R")
source("integration_helper_v2.R")

cat("=== SINGLE-SHOT run (count camera arm) ===\n")
cat("species:", species, " chain:", chain_id, " nburnin:", NBURNIN, " niter(recorded):", NITER, "\n")
cat("start:", format(Sys.time()), "\n")

model_file <- "model_code_national_scalar_camcount.R"
cat("model file md5:", unname(tools::md5sum(model_file)), "\n")
source(model_file)
model_code <- model_code_national_scalar_camcount

mons <- c("occ_beta", "link_occ_intercept", "year_beta", "year_var", "total_var_beta",
          "MWMT_effect", "MCMT_effect", "trend_robust_indicator",
          "overdisp_cam", "overdisp_inat")
cat("monitors:", paste(mons, collapse = ", "), "\n")

input_data <- readRDS(paste0("input_data_", species, ".RDS"))
cl <- input_data$constants_list
n  <- cl$nsite
stopifnot(
  "y_cam length"           = length(input_data$y_cam) == n,
  "y_cam integer >= 0"     = all(input_data$y_cam >= 0) && all(input_data$y_cam == round(input_data$y_cam)),
  "log_effort_cam length"  = length(cl$log_effort_cam) == n && all(is.finite(cl$log_effort_cam)),
  "yday_mean length"       = length(cl$yday_mean) == n && !anyNA(cl$yday_mean),
  "canopy/road length"     = length(cl$canopy_height) == n && length(cl$log_roaddist) == n,
  "occ_covars rows"        = nrow(cl$occ_covars) == n && ncol(cl$occ_covars) == cl$numOccCovars,
  "no det_intercept init"  = is.null(input_data$inits_list$det_intercept),
  "overdisp_cam init"      = !is.null(input_data$inits_list$overdisp_cam)
)
cat("pre-build data checks passed: nsite", n, " sum(y_cam)", sum(input_data$y_cam), "\n")

out_file <- file.path(project_dir, paste0("chain_", species, "_", chain_id, "_single.RDS"))
if (file.exists(out_file)) stop("output already exists, refusing to overwrite: ", out_file)

t0 <- Sys.time()
model <- nimbleModel(model_code,
                     constants = cl,
                     data = list(y_cam = input_data$y_cam,
                                 y_inat = input_data$inat_y,
                                 inat_effort = input_data$inat_effort),
                     inits = input_data$inits_list,
                     calculate = FALSE)
cat("model defined, seconds:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "\n")
Cmodel <- compileNimble(model)
conf   <- configureMCMC(model, enableWAIC = TRUE)
do.call(conf$addMonitors, as.list(mons))
mcmc   <- buildMCMC(conf)
Cmcmc  <- compileNimble(mcmc, project = model)
t_build <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
cat("build+compile seconds:", round(t_build, 1), "\n")
cat("initial logProb:", Cmodel$calculate(), "\n")

set.seed(1000 + chain_id)
t1 <- Sys.time()
cat("running", NBURNIN + NITER, "iterations in ONE call (nburnin =", NBURNIN, ")\n")
Cmcmc$run(niter = NBURNIN + NITER, nburnin = NBURNIN)
t_mcmc <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
cat("mcmc seconds:", round(t_mcmc, 1), " sec/iter:", round(t_mcmc / (NBURNIN + NITER), 4), "\n")

samples <- as.matrix(Cmcmc$mvSamples)
cat("samples:", nrow(samples), "x", ncol(samples), "\n")
stopifnot(nrow(samples) == NITER)
waic <- tryCatch(Cmcmc$getWAIC(), error = function(e) { cat("WAIC extraction failed:", conditionMessage(e), "\n"); NULL })

saveRDS(list(samples = samples, iter_total = NITER, nburnin = NBURNIN, waic = waic,
             chain_id = chain_id, species = species, single_shot = TRUE,
             model_file_md5 = unname(tools::md5sum(model_file)),
             build_sec = t_build, mcmc_sec = t_mcmc, finished_at = format(Sys.time())),
        file = out_file)
cat("wrote", out_file, "\nend:", format(Sys.time()), "\nDONE\n")
