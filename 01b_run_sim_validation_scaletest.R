#!/usr/bin/env Rscript
# =============================================================================
# 01b_run_sim_validation_scaletest.R
# CONFIRMATORY follow-up to 01_run_sim_validation.R. The first valid
# 50-replicate run (job 439834) showed good coverage (~99.5%) but real
# shrinkage bias toward zero in the informative regions (NE bias -0.27,
# West bias +0.13 against a true +-0.3 signal) and only 80% NE-vs-West
# sign recovery. Hypothesis: this is a reduced-simulation-DESIGN artifact,
# not a model-formulation problem -- the original design only directly
# informs 178/908 CAR cells (~20%), far sparser than the real deployment
# (all 908 cells have camera and/or iNat data), so the CAR prior's
# zero-mean shrinkage is expected to bite harder here than it will on the
# real refit.
#
# This script triples the retained cell50 count and camera sites (same
# per-cell subcell cap) to roughly triple the informed-cell count, and
# checks whether bias/sign-recovery improve as predicted. Fewer replicates
# (20, not 50) since this is a diagnostic check, not the final validation
# -- if bias shrinks as expected here, that's strong evidence the full
# real-data refit (Step 2) will do even better, without paying for another
# full 50-replicate array at 3x the per-replicate cost.
#
# RUN: as a Slurm array (see 01b_run_sim_validation_scaletest_sbatch.sh).
# Writes: sim_results_scaletest/rep_<id>.RDS
# =============================================================================

suppressPackageStartupMessages({
  library(nimble); library(nimbleEcology); library(coda)
})

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
setwd(PROTO_DIR)

source(paste0(PROJ, "/HPC/bobcat/integration_helper.R"))
source("model_code_spatial_trend.R")
source("sim_helpers.R")

N_CELL50_KEEP <- 630   # 3x the original 210
MAX_SUBCELL   <- 15    # unchanged
N_SITE_KEEP   <- 2000  # ~2.9x the original 700
N_BURNIN <- 1000
N_ITER   <- 3000
SIM_SEED <- 20260713   # different from the original run's design seed (20260712) -- a different, larger subsample, not a superset of the same one

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
args <- commandArgs(trailingOnly = TRUE)
rep_id <- if (!is.na(task_id)) as.integer(task_id) else if (length(args) >= 1) as.integer(args[1]) else stop(
  "No SLURM_ARRAY_TASK_ID and no CLI rep id given -- this script runs ONE replicate per invocation.")

RESULTS_DIR <- paste0(PROTO_DIR, "/sim_results_scaletest")
dir.create(RESULTS_DIR, showWarnings = FALSE)

cat("=== scale-test replicate", rep_id, "===\n")

prepped <- readRDS("prepped_sim_inputs.RDS")
inp <- readRDS(paste0(PROJ, "/HPC/bobcat/input_data_bobcat.RDS"))
cl  <- prepped$constants_list

reduced <- build_reduced_constants(
  cl = cl, inat_effort_real = inp$inat_effort, real_y_template = inp$real_data$y,
  cell100_geo = prepped$cell100_geo,
  n_cell50_keep = N_CELL50_KEEP, max_subcell_per_cell = MAX_SUBCELL, n_site_keep = N_SITE_KEEP,
  seed = SIM_SEED)
cat("ncell50=", reduced$constants$ncell50, " nsite=", reduced$constants$nsite,
    " subcell rows=", dim(reduced$constants$xdat_inat)[1], "\n")

year_effect_true <- make_true_year_effect(prepped$cell100_geo, cl$adj, cl$num, amplitude = 0.3)
truth <- true_param_list(prepped$real_post_means, year_effect_true)
informed_cell100 <- unique(reduced$constants$inat_cell100)
cat("informed cell100 cells:", length(informed_cell100), "/", cl$ncell100,
    "(original design had 178)\n")

set.seed(SIM_SEED + rep_id)   # see 01_run_sim_validation.R for why this matters

res <- run_one_replicate(
  rep_id = rep_id, model_code = model_code_spatial_trend,
  constants = reduced$constants, inat_effort = reduced$inat_effort,
  y_ncol = reduced$y_ncol, truth = truth, year_effect_true = year_effect_true,
  cell100_geo = prepped$cell100_geo, informed_cell100 = informed_cell100,
  base_inits = inp$inits_list, n_burnin = N_BURNIN, n_iter = N_ITER)

saveRDS(res, paste0(RESULTS_DIR, "/rep_", rep_id, ".RDS"))
cat("\nreplicate", rep_id, "status:", res$status, " elapsed_sec:", round(res$elapsed_sec, 1), "\n")
print(res)
