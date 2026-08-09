# 01f_run_indicator_test.R
# ---------------------------------------------------------------------------
# RJMCMC indicator test: replaces the WAIC gate (shown by the corrected
# abundance sweep to be unusable -- null/varying WAIC distributions overlap
# almost completely, p = 0.80/0.91/0.61, see abundance_sweep_evaluation.md)
# with the Goldstein et al. (bioRxiv 2025.01.17.633640) approach: a single
# Bernoulli "switch" indicator gamma on the ecoregion trend deviation, fit
# with reversible-jump MCMC, read via posterior inclusion probability
# P(gamma=1). See rjmcmc_indicator_test_plan.md.
#
# Same 180-task array structure as 01e_run_abundance_sweep.R (abundance x
# scenario x rep), same design grid, same row_id -> design row mapping --
# BUT DOES NOT REGENERATE DATA. Reuses the exact simulated y/y_inat already
# stored in sim_results_abundance/<abundance>_<scenario>/rep_<id>.RDS's
# res$sim_data[[1]] list-column, from the corrected (reseeded, guard-passed)
# 180-replicate run. Only falls back to re-simulating (with the identical
# per-replicate seed convention, set.seed(20260712 + row_id) after
# build_reduced_constants()) if that file or its sim_data is missing.
#
# Depends on (source in this order):
#   model_code_ecoregion_switch.R  (defines model_code_ecoregion_switch --
#                                   gamma indicator, deterministic zeroing
#                                   gate, total_data_logLik convergence node)
#   model_code_ecoregion_trend.R   (defines model_code_ecoregion -- needed
#                                   ONLY by the resimulation fallback path,
#                                   same truth-generation model as 01e)
#   sim_helpers.R                  (fit_replicate_switch, run_one_indicator_replicate,
#                                   build_reduced_constants, make_true_year_region,
#                                   true_param_list_ecoregion -- fallback path)
#   sim_helpers_abundance.R        (scale_truth_abundance, abundance ladder --
#                                   fallback path only)
#   prepped_sim_inputs.RDS         (+ ecoregion fields, baked in by 00b)
#   abundance_anchor.RDS           (optional, from 00c; else default ladder)
#
# pi_gamma (prior inclusion probability) and slab_sd (FIXED sd for
# year_region_raw's prior -- RJMCMC requires a constant target-node
# hyperparameter, see model_code_ecoregion_switch.R header) are set as
# constants here, not estimated: pi_gamma = 0.5 (plan default), slab_sd =
# 0.227 (mean posterior sigma_region_mean from the corrected abundance
# sweep's varying scenario -- close to the 0.2 truth-generation target).
#
# Writes: sim_results_indicator/<abundance>_<scenario>/rep_<id>.RDS
# ---------------------------------------------------------------------------

suppressMessages({ library(nimble); library(nimbleEcology); library(dplyr); library(coda) })

PROTO <- getwd()
PROJ  <- dirname(PROTO)
source(file.path(PROJ, "HPC", "bobcat", "integration_helper.R"))
source(file.path(PROTO, "sim_helpers.R"))
source(file.path(PROTO, "sim_helpers_abundance.R"))
source(file.path(PROTO, "model_code_ecoregion_trend.R"))    # defines model_code_ecoregion (fallback path only)
source(file.path(PROTO, "model_code_ecoregion_switch.R"))   # defines model_code_ecoregion_switch

PI_GAMMA <- 0.5
SLAB_SD  <- 0.227

prepped <- readRDS(file.path(PROTO, "prepped_sim_inputs.RDS"))
stopifnot("ecoregion_id" %in% names(prepped$cell100_geo),
          "ecoregion_of_cell100" %in% names(prepped))

# ---- design (IDENTICAL to 01e_run_abundance_sweep.R's grid/ordering) ------
N_REPS   <- 30
abund_path <- file.path(PROTO, "abundance_anchor.RDS")
abund_levels <- if (file.exists(abund_path)) readRDS(abund_path) else abundance_levels_default()

trend_scenarios <- c("varying", "null")

design <- expand.grid(abundance = abund_levels$level,
                      scenario  = trend_scenarios,
                      rep       = seq_len(N_REPS),
                      stringsAsFactors = FALSE)
design <- design[order(design$abundance, design$scenario, design$rep), ]
design$row <- seq_len(nrow(design))

# ---- pick this task's row --------------------------------------------------
arg <- Sys.getenv("SLURM_ARRAY_TASK_ID", NA)
if (is.na(arg)) { a <- commandArgs(trailingOnly = TRUE); arg <- if (length(a)) a[1] else NA }
stopifnot("Need SLURM_ARRAY_TASK_ID or CLI row id (1 combo per task)" = !is.na(arg))
row_id <- as.integer(arg)
this <- design[design$row == row_id, ]
abn  <- abund_levels[abund_levels$level == this$abundance, ]
cat(sprintf("=== abundance=%s  scenario=%s  rep=%d (row %d/%d) ===\n",
            this$abundance, this$scenario, this$rep, row_id, nrow(design)))

# ---- reduced design (ecoregion-stratified) -- IDENTICAL fixed-seed subsample
# as 01e (build_reduced_constants() seeds internally with a fixed seed by
# design), so dimensions line up exactly with the stored sim_data.
cl <- build_reduced_constants(prepped$constants_list, prepped$inat_effort_real,
                              prepped$real_y_template, prepped$cell100_geo,
                              stratify_by = "ecoregion",
                              ecoregion_of_cell100_full = prepped$ecoregion_of_cell100,
                              nregion = prepped$nregion)
cl$constants$pi_gamma <- PI_GAMMA
cl$constants$slab_sd  <- SLAB_SD

# ---- reuse stored sim_data from the corrected abundance sweep -------------
stored_path <- file.path(PROTO, "sim_results_abundance",
                         paste0(this$abundance, "_", this$scenario),
                         sprintf("rep_%03d.RDS", this$rep))
sim_data <- NULL
truth <- NULL
if (file.exists(stored_path)) {
  stored <- readRDS(stored_path)
  if (identical(stored$status, "ok") && !is.null(stored$sim_data) && length(stored$sim_data) >= 1) {
    sim_data <- stored$sim_data[[1]]
    cat("Reusing stored sim_data from", stored_path, "\n")
  }
}
if (is.null(sim_data)) {
  cat("WARNING: stored sim_data not usable at", stored_path, "-- falling back to re-simulation ",
      "with the identical per-replicate seed convention used by 01e_run_abundance_sweep.R.\n")
  year_region_true <- make_true_year_region(prepped$cell100_geo, scenario = this$scenario)
  truth <- true_param_list_ecoregion(prepped$real_post_means, year_region_true)
  truth <- scale_truth_abundance(truth,
                                 occ_shift      = abn$occ_shift,
                                 count_log_mult = abn$count_log_mult,
                                 label          = this$abundance)
  set.seed(20260712 + row_id)   # same convention as 01e -- reproduces byte-identical data
}

# ---- fit the switch model with explicit RJMCMC -----------------------------
res <- run_one_indicator_replicate(
  rep_id      = this$rep,
  constants   = cl$constants,
  inat_effort = cl$inat_effort,
  y_ncol      = cl$y_ncol,
  truth       = truth,
  base_inits  = prepped$base_inits,
  sim_data    = sim_data,
  n_chains    = 2)

res$abundance      <- this$abundance
res$scenario       <- this$scenario
res$occ_shift      <- abn$occ_shift
res$count_log_mult <- abn$count_log_mult

OUT <- file.path(PROTO, "sim_results_indicator",
                 paste0(this$abundance, "_", this$scenario))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(OUT, sprintf("rep_%03d.RDS", this$rep)))
cat(sprintf("status: %s | P(gamma=1)=%s | rhat_logLik=%s | ess_logLik=%s\n",
            res$status, round(res$gamma_mean, 4), round(res$rhat_logLik, 4), round(res$ess_logLik, 1)))
