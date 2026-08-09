#!/usr/bin/env Rscript
# =============================================================================
# 01d_run_grain3_supplement.R
# Supplementary run: 15 replicates at n_smooth_iter=3 using the SAME
# white-noise-seeded truth-field generator as the grain=0/1/8 sweep
# (make_true_year_effect_grain()), so spatial_cor becomes available at that
# grain level too. NOT a duplicate of the original 50-rep grain=3 result --
# that one used make_true_year_effect() (a regionally-structured NE+/West-
# seed, not white noise), so while both share "3 smoothing passes," they are
# different specific truth fields. This run exists specifically to complete
# the grain curve on a consistent (noise-seeded) basis, since the surprising
# finding that spatial_cor stayed near-zero across grain 0/1/8 raised the
# question of whether the original grain=3 result (strong NE/West sign
# recovery) reflects better SPATIAL pattern recovery, or just a coarser,
# more forgiving two-group-average test -- and the original run didn't save
# posterior samples, so spatial_cor can't be computed retroactively for it.
#
# rep_id = SLURM_ARRAY_TASK_ID directly (1-15), no combined design_df needed
# since this is a single grain level.
#
# RUN: as a Slurm array (see 01d_run_grain3_supplement_sbatch.sh).
# Writes: sim_results_grain/grain_3/rep_<id>.RDS
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

N_SMOOTH_ITER <- 3
N_CELL50_KEEP <- 210
MAX_SUBCELL   <- 15
N_SITE_KEEP   <- 700
N_BURNIN <- 1000
N_ITER   <- 3000
DESIGN_SEED <- 20260712   # same retained cell/site subsample as the rest of the sweep

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
args <- commandArgs(trailingOnly = TRUE)
rep_id <- if (!is.na(task_id)) as.integer(task_id) else if (length(args) >= 1) as.integer(args[1]) else stop(
  "No SLURM_ARRAY_TASK_ID and no CLI rep id given.")

RESULTS_DIR <- paste0(PROTO_DIR, "/sim_results_grain/grain_", N_SMOOTH_ITER)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== grain=3 supplement (noise-seeded): rep =", rep_id, "===\n")

prepped <- readRDS("prepped_sim_inputs.RDS")
inp <- readRDS(paste0(PROJ, "/HPC/bobcat/input_data_bobcat.RDS"))
cl  <- prepped$constants_list

reduced <- build_reduced_constants(
  cl = cl, inat_effort_real = inp$inat_effort, real_y_template = inp$real_data$y,
  cell100_geo = prepped$cell100_geo,
  n_cell50_keep = N_CELL50_KEEP, max_subcell_per_cell = MAX_SUBCELL, n_site_keep = N_SITE_KEEP,
  seed = DESIGN_SEED)

year_effect_true <- make_true_year_effect_grain(
  prepped$cell100_geo, cl$adj, cl$num, n_smooth_iter = N_SMOOTH_ITER, seed = 20260713)
cat("true year_effect (grain=3, noise-seeded): sd=", round(sd(year_effect_true), 4),
    " range=", paste(round(range(year_effect_true), 3), collapse=" to "), "\n")

truth <- true_param_list(prepped$real_post_means, year_effect_true)
informed_cell100 <- unique(reduced$constants$inat_cell100)

set.seed(DESIGN_SEED + N_SMOOTH_ITER * 100 + rep_id)

res <- run_one_replicate(
  rep_id = rep_id, model_code = model_code_spatial_trend,
  constants = reduced$constants, inat_effort = reduced$inat_effort,
  y_ncol = reduced$y_ncol, truth = truth, year_effect_true = year_effect_true,
  cell100_geo = prepped$cell100_geo, informed_cell100 = informed_cell100,
  base_inits = inp$inits_list, n_burnin = N_BURNIN, n_iter = N_ITER,
  metrics_fn = compute_grain_metrics)
res$grain <- N_SMOOTH_ITER

saveRDS(res, paste0(RESULTS_DIR, "/rep_", rep_id, ".RDS"))
cat("\nrep=", rep_id, " status:", res$status, " elapsed_sec:", round(res$elapsed_sec, 1), "\n")
print(res)
