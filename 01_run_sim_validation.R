#!/usr/bin/env Rscript
# =============================================================================
# 01_run_sim_validation.R
# STEP 1 of the spatially-varying-trend prototype: ONE replicate of the
# simulation-based validation. Run as a Slurm ARRAY job (see
# 01_run_sim_validation_sbatch.sh, --array=1-50) so every replicate gets a
# genuinely fresh OS process. After the array completes, run
# 02_collect_sim_summary.R to aggregate bias/RMSE/coverage/sign-recovery.
#
# WHY AN ARRAY, NOT AN IN-SCRIPT PARALLEL POOL: two earlier attempts failed
# for different reasons, both instructive:
#   - job 436025 used mclapply (fork-based): all forked workers inherit the
#     SAME tempdir() and the SAME nimble model-ID counter from the parent,
#     so concurrent compileNimble() calls raced to write the identical
#     generated-code path and corrupted each other's shared library (every
#     one of 45 completed replicates errored: "no such symbol
#     new_model_code_MID_2").
#   - job 436626 switched to a PSOCK cluster (genuinely separate R sessions,
#     confirmed correct: first 10 replicates all completed "ok"), but each
#     PSOCK worker is PERSISTENT and processes many replicates sequentially
#     within the same R session. NIMBLE compiles a fresh model from scratch
#     every replicate, and its generated shared libraries stay loaded for
#     the life of the process (never freed) -- so per-worker memory grows
#     with every additional replicate it handles, and the job was OOM-killed
#     partway through the second batch (80G requested, all consumed).
# A Slurm array sidesteps both: every replicate is a brand-new process that
# compiles exactly once and exits, so there is no shared tempdir/UID state to
# race on and no accumulation to OOM on.
#
# Real-scale dry run (job 435980, single fresh process): 12.5 min wall,
# 5.5GB peak, bias=-0.013, coverage=100%, sign_recovered=TRUE. Reduced
# design (build_reduced_constants(): real 908-cell CAR graph kept exactly,
# stratified iNat/camera subsample -- see sim_helpers.R) confirmed with the
# user before writing this script.
#
# RUN (as a Slurm array task -- SLURM_ARRAY_TASK_ID selects the replicate id;
# falls back to a CLI arg for a one-off manual rerun):
#   $PROJ/HPC/conda_envs/nimble_env/bin/Rscript 01_run_sim_validation.R
# Writes: sim_results/rep_<id>.RDS
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

N_CELL50_KEEP <- 210
MAX_SUBCELL   <- 15
N_SITE_KEEP   <- 700
N_BURNIN <- 1000
N_ITER   <- 3000
SIM_SEED <- 20260712   # fixed design (which cells/sites retained) -- same for every replicate

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
args <- commandArgs(trailingOnly = TRUE)
rep_id <- if (!is.na(task_id)) as.integer(task_id) else if (length(args) >= 1) as.integer(args[1]) else stop(
  "No SLURM_ARRAY_TASK_ID and no CLI rep id given -- this script runs ONE replicate per invocation.")

RESULTS_DIR <- paste0(PROTO_DIR, "/sim_results")
dir.create(RESULTS_DIR, showWarnings = FALSE)

cat("=== replicate", rep_id, "===\n")

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

# CRITICAL: build_reduced_constants() calls set.seed(SIM_SEED) internally so
# every replicate shares the identical retained-cell/site DESIGN (intended).
# But that leaves R's global RNG state deterministic and IDENTICAL across
# every array task at this point -- if the data simulation below ran with
# that leftover state, every "replicate" would draw the exact same y/y_inat
# (caught in job 436723: rep 1 and rep 2 had bit-identical bias_all). Reseed
# per-replicate now, after the shared design is fixed, so each replicate's
# simulated data is genuinely independent.
set.seed(SIM_SEED + rep_id)

res <- run_one_replicate(
  rep_id = rep_id, model_code = model_code_spatial_trend,
  constants = reduced$constants, inat_effort = reduced$inat_effort,
  y_ncol = reduced$y_ncol, truth = truth, year_effect_true = year_effect_true,
  cell100_geo = prepped$cell100_geo, informed_cell100 = informed_cell100,
  base_inits = inp$inits_list, n_burnin = N_BURNIN, n_iter = N_ITER)

saveRDS(res, paste0(RESULTS_DIR, "/rep_", rep_id, ".RDS"))
cat("\nreplicate", rep_id, "status:", res$status, " elapsed_sec:", round(res$elapsed_sec, 1), "\n")
print(res)
