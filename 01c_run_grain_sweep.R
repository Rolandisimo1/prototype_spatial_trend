#!/usr/bin/env Rscript
# =============================================================================
# 01c_run_grain_sweep.R
# PHASE 1 of the follow-up study: what spatial grain of trend heterogeneity
# can actually be detected given current data density? Simulates truth fields
# at varying spatial grain (see make_true_year_effect_grain() in
# sim_helpers.R: white noise diffused over the real CAR graph, n_smooth_iter
# passes controls grain -- fewer passes = patchier/finer, more = smoother/
# broader) and checks recovery.
#
# n_smooth_iter = 3 exactly reproduces the grain of the already-validated
# continental-scale truth field (make_true_year_effect(), used by the
# original 50-replicate run, job 439834). That existing result is REUSED as
# one point on this grain curve (see 02c_collect_grain_sweep_summary.R) --
# this script only runs the 3 NEW grain levels, saving ~25% of the compute.
#
# Density held FIXED at the original design's level (210 cell50 / 700 sites,
# same SIM_SEED so it's the identical retained-cell/site subsample as the
# original run) so grain is the only thing varying across this sweep --
# density is phases 2-3's job, not this one.
#
# ONE (grain_level, rep) combination per Slurm array task -- see the
# combined design_df built below, matching the established one-fresh-process-
# per-replicate pattern (see 01_run_sim_validation.R header for why: avoids
# both the mclapply fork-collision and the PSOCK memory-accumulation bugs
# hit earlier in this project).
#
# RUN: as a Slurm array (see 01c_run_grain_sweep_sbatch.sh).
# Writes: sim_results_grain/grain_<n_smooth_iter>/rep_<id>.RDS
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

# ------------------------------ config ---------------------------------------
GRAIN_LEVELS  <- c(0, 1, 8)     # smoothing passes; 3 already covered by job 439834
N_REPS        <- 15
N_CELL50_KEEP <- 210            # matches the original design exactly (job 439834)
MAX_SUBCELL   <- 15
N_SITE_KEEP   <- 700
N_BURNIN <- 1000
N_ITER   <- 3000
DESIGN_SEED <- 20260712         # SAME as the original run -- identical retained cell/site subsample

grain_design_df <- expand.grid(grain = GRAIN_LEVELS, rep = seq_len(N_REPS))
grain_design_df <- grain_design_df[order(grain_design_df$grain, grain_design_df$rep), ]
grain_design_df$row <- seq_len(nrow(grain_design_df))

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
args <- commandArgs(trailingOnly = TRUE)
row_id <- if (!is.na(task_id)) as.integer(task_id) else if (length(args) >= 1) as.integer(args[1]) else stop(
  "No SLURM_ARRAY_TASK_ID and no CLI row id given -- this script runs ONE (grain, rep) combination per invocation.")

this_row <- grain_design_df[grain_design_df$row == row_id, ]
n_smooth_iter <- this_row$grain
rep_id <- this_row$rep

RESULTS_DIR <- paste0(PROTO_DIR, "/sim_results_grain/grain_", n_smooth_iter)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== grain sweep: n_smooth_iter =", n_smooth_iter, " rep =", rep_id, "===\n")

prepped <- readRDS("prepped_sim_inputs.RDS")
inp <- readRDS(paste0(PROJ, "/HPC/bobcat/input_data_bobcat.RDS"))
cl  <- prepped$constants_list

reduced <- build_reduced_constants(
  cl = cl, inat_effort_real = inp$inat_effort, real_y_template = inp$real_data$y,
  cell100_geo = prepped$cell100_geo,
  n_cell50_keep = N_CELL50_KEEP, max_subcell_per_cell = MAX_SUBCELL, n_site_keep = N_SITE_KEEP,
  seed = DESIGN_SEED)
cat("ncell50=", reduced$constants$ncell50, " nsite=", reduced$constants$nsite, "\n")

year_effect_true <- make_true_year_effect_grain(
  prepped$cell100_geo, cl$adj, cl$num, n_smooth_iter = n_smooth_iter, seed = 20260713)
cat("true year_effect (grain=", n_smooth_iter, "): sd=", round(sd(year_effect_true), 4),
    " range=", paste(round(range(year_effect_true), 3), collapse=" to "),
    " mean=", round(mean(year_effect_true), 6), "\n")

truth <- true_param_list(prepped$real_post_means, year_effect_true)
informed_cell100 <- unique(reduced$constants$inat_cell100)

set.seed(DESIGN_SEED + n_smooth_iter * 100 + rep_id)   # see 01_run_sim_validation.R for why this matters

res <- run_one_replicate(
  rep_id = rep_id, model_code = model_code_spatial_trend,
  constants = reduced$constants, inat_effort = reduced$inat_effort,
  y_ncol = reduced$y_ncol, truth = truth, year_effect_true = year_effect_true,
  cell100_geo = prepped$cell100_geo, informed_cell100 = informed_cell100,
  base_inits = inp$inits_list, n_burnin = N_BURNIN, n_iter = N_ITER,
  metrics_fn = compute_grain_metrics)
res$grain <- n_smooth_iter

saveRDS(res, paste0(RESULTS_DIR, "/rep_", rep_id, ".RDS"))
cat("\ngrain=", n_smooth_iter, " rep=", rep_id, " status:", res$status,
    " elapsed_sec:", round(res$elapsed_sec, 1), "\n")
print(res)
