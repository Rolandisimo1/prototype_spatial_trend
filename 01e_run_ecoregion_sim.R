#!/usr/bin/env Rscript
# =============================================================================
# 01e_run_ecoregion_sim.R
# The ecoregion-trend simulation study's core deliverable: ONE (scenario, rep)
# combination per Slurm array task. For each replicate, simulates data ONCE
# from the ecoregion model's own likelihood (model_code_ecoregion, defined in
# model_code_ecoregion_trend.R) at a known truth, then fits BOTH candidate
# models to that same data -- the national-scalar production structure
# (model_code_national_scalar.R) and the ecoregion extension -- and compares
# via WAIC (run_one_replicate(..., fit_null_scalar=TRUE)). Two scenarios:
#   "varying": each of the K=8 ecoregions gets its own true trend deviation
#     (deterministic spread from clearly-negative through ~0 to
#     clearly-positive -- see make_true_year_region()) -- tests recovery of
#     real regional structure.
#   "null": every ecoregion gets the SAME trend (all year_region=0) -- tests
#     whether the model invents regional structure that isn't there.
#
# This is the FIXED-abundance (bobcat baseline) two-scenario sweep; see
# 01e_run_abundance_sweep.R for the abundance x scenario extension layered
# on top (same underlying model/truth-builder, varying occupancy/detection
# magnitude in addition to the trend scenario).
#
# Density fixed at the CAR prototype's higher "scaletest" density (630
# cell50 / 2000 sites, ~3x the original design) per the task spec, since
# that's closer to the real bobcat fit's information content than the
# original design used for most of the CAR-field sweep.
#
# RUN: as a Slurm array (see 01e_run_ecoregion_sim_sbatch.sh).
# Writes: sim_results_ecoregion/<scenario>/rep_<id>.RDS
# =============================================================================

suppressPackageStartupMessages({
  library(nimble); library(nimbleEcology); library(coda)
})

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
setwd(PROTO_DIR)

source(paste0(PROJ, "/HPC/bobcat/integration_helper.R"))
source("model_code_ecoregion_trend.R")   # defines model_code_ecoregion
source("model_code_national_scalar.R")   # defines model_code_national_scalar (needed by fit_null_scalar=TRUE)
source("sim_helpers.R")

# ------------------------------ config ---------------------------------------
SCENARIOS <- c("varying", "null")
N_REPS <- 30
N_CELL50_KEEP <- 630   # CAR prototype's "scaletest" (higher) density
MAX_SUBCELL   <- 15
N_SITE_KEEP   <- 2000
N_BURNIN <- 1000
N_ITER   <- 3000
DESIGN_SEED <- 20260712   # same design seed as the rest of the prototype

design_df <- expand.grid(scenario = SCENARIOS, rep = seq_len(N_REPS), stringsAsFactors = FALSE)
design_df <- design_df[order(design_df$scenario, design_df$rep), ]
design_df$row <- seq_len(nrow(design_df))

task_id <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = NA)
args <- commandArgs(trailingOnly = TRUE)
row_id <- if (!is.na(task_id)) as.integer(task_id) else if (length(args) >= 1) as.integer(args[1]) else stop(
  "No SLURM_ARRAY_TASK_ID and no CLI row id given -- this script runs ONE (scenario, rep) combination per invocation.")

this_row <- design_df[design_df$row == row_id, ]
scenario <- this_row$scenario
rep_id <- this_row$rep

RESULTS_DIR <- paste0(PROTO_DIR, "/sim_results_ecoregion/", scenario)
dir.create(RESULTS_DIR, showWarnings = FALSE, recursive = TRUE)

cat("=== ecoregion sim: scenario =", scenario, " rep =", rep_id, "===\n")

prepped <- readRDS("prepped_sim_inputs.RDS")   # ecoregion fields baked in by 00b_prep_ecoregion.R
stopifnot("ecoregion_id" %in% names(prepped$cell100_geo),
         "ecoregion_of_cell100" %in% names(prepped))
cl  <- prepped$constants_list

reduced <- build_reduced_constants(
  cl = cl, inat_effort_real = prepped$inat_effort_real, real_y_template = prepped$real_y_template,
  cell100_geo = prepped$cell100_geo,
  n_cell50_keep = N_CELL50_KEEP, max_subcell_per_cell = MAX_SUBCELL, n_site_keep = N_SITE_KEEP,
  seed = DESIGN_SEED,
  stratify_by = "ecoregion",
  ecoregion_of_cell100_full = prepped$ecoregion_of_cell100, nregion = prepped$nregion)
cat("ncell50=", reduced$constants$ncell50, " nsite=", reduced$constants$nsite,
    " nregion=", reduced$constants$nregion, "\n")

year_region_true <- make_true_year_region(prepped$cell100_geo, scenario = scenario)
cat("true year_region (scenario=", scenario, "):",
    paste(round(year_region_true, 3), collapse=", "), "\n")

truth <- true_param_list_ecoregion(prepped$real_post_means, year_region_true)

set.seed(DESIGN_SEED + (scenario == "null") * 1000 + rep_id)   # see 01_run_sim_validation.R for why this matters

res <- run_one_replicate(
  rep_id = rep_id, model_code = model_code_ecoregion,
  constants = reduced$constants, inat_effort = reduced$inat_effort, y_ncol = reduced$y_ncol,
  truth = truth, base_inits = prepped$base_inits, n_burnin = N_BURNIN, n_iter = N_ITER,
  metrics_fn = compute_ecoregion_metrics,
  trend_inits = list(year_beta = 0, year_var = 0, sigma_region = 1,
                    year_region = rep(0, prepped$nregion)),
  extra_monitors = c("year_region", "sigma_region"),
  fit_null_scalar = TRUE,
  year_region_true = year_region_true, ecoregion_levels = prepped$ecoregion_levels)
res$scenario <- scenario

saveRDS(res, paste0(RESULTS_DIR, "/rep_", rep_id, ".RDS"))
cat("\nscenario=", scenario, " rep=", rep_id, " status:", res$status,
    " elapsed_sec:", round(res$elapsed_sec, 1), "\n")
print(res[, setdiff(names(res), "sim_data")])
