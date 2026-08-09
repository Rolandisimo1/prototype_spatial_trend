# 01e_run_abundance_sweep.R
# ---------------------------------------------------------------------------
# Ecoregion-trend simulation, ABUNDANCE x TREND-SCENARIO sweep.
# One (abundance_level, trend_scenario, rep) combination per Slurm array task.
#
# Crosses:
#   abundance_level : bobcat_baseline / moderate / common(_deerlike)
#   trend_scenario  : "varying" (each ecoregion its own trend)  vs
#                     "null"    (all ecoregions share national trend; year_region = 0)
#
# Depends on (source in this order):
#   model_code_ecoregion_trend.R   (forked model: defines model_code_ecoregion,
#                                   year_region[r] ~ dnorm(0, sigma_region) random
#                                   effect, plus a deterministic year_effect[c] <-
#                                   year_region[ecoregion_of_cell100[c]] pass-through
#                                   so the calcIntensity_SVC call site matches the
#                                   CAR model byte-for-byte)
#   model_code_national_scalar.R   (defines model_code_national_scalar -- needed
#                                   by run_one_replicate(fit_null_scalar=TRUE))
#   sim_helpers.R                  (build_reduced_constants, simulate_replicate_data,
#                                   fit_replicate, run_one_replicate, metrics)
#   sim_helpers_abundance.R        (scale_truth_abundance, abundance ladder,
#                                   summarize_simulated_information)
#   prepped_sim_inputs.RDS         (+ ecoregion fields + inat_effort_real /
#                                   real_y_template / base_inits, all baked in
#                                   by 00b_prep_ecoregion.R)
#   abundance_anchor.RDS           (optional, from 00c; else default ladder)
#
# NOTE on the truth-building step: year_region[r] is the model's actual
# stochastic/estimated node -- year_effect[c] is a DETERMINISTIC pass-through
# from it, not something you set an init on directly. So the truth list has
# to be built with true_param_list_ecoregion(real_post_means, year_region_true)
# (sets year_region), not the CAR model's true_param_list() with a
# per-cell-mapped vector (which would try to init a node -- year_effect --
# that no longer accepts one in this model).
#
# Writes: sim_results_abundance/<abundance>_<scenario>/rep_<id>.RDS
# ---------------------------------------------------------------------------

suppressMessages({ library(nimble); library(nimbleEcology); library(dplyr) })

PROTO <- getwd()
PROJ  <- dirname(PROTO)
# calcIntensity_SVC/calcIntensity_noSVC (referenced by both model_code_ecoregion
# and model_code_national_scalar) are defined here -- missing this source line
# was caught the hard way: array job 444455 failed EVERY task with "R function
# 'calcIntensity_SVC' ... does not exist" (fast failures, ~1.5 min each, before
# any real MCMC ran). Every other script in this project already sources this;
# this file was the one gap.
source(file.path(PROJ, "HPC", "bobcat", "integration_helper.R"))
source(file.path(PROTO, "sim_helpers.R"))
source(file.path(PROTO, "sim_helpers_abundance.R"))
source(file.path(PROTO, "model_code_ecoregion_trend.R"))   # defines model_code_ecoregion
source(file.path(PROTO, "model_code_national_scalar.R"))   # defines model_code_national_scalar

prepped <- readRDS(file.path(PROTO, "prepped_sim_inputs.RDS"))
stopifnot("ecoregion_id" %in% names(prepped$cell100_geo),
          "ecoregion_of_cell100" %in% names(prepped))

# ---- design ---------------------------------------------------------------
N_REPS   <- 30
abund_path <- file.path(PROTO, "abundance_anchor.RDS")
abund_levels <- if (file.exists(abund_path)) readRDS(abund_path) else abundance_levels_default()
cat("Abundance ladder in use:\n"); print(abund_levels)

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

# ---- build the known-truth ecoregion trend field --------------------------
# make_true_year_region(): "varying" -> each ecoregion_id gets its own trend
# (deterministic spread of +/-/~0); "null" -> all ecoregion deviations = 0
# (national trend only).
year_region_true <- make_true_year_region(prepped$cell100_geo, scenario = this$scenario)
cat("true year_region by ecoregion:\n"); print(round(year_region_true, 3))

# ---- reduced design (ecoregion-stratified) ---------------------------------
cl <- build_reduced_constants(prepped$constants_list, prepped$inat_effort_real,
                              prepped$real_y_template, prepped$cell100_geo,
                              stratify_by = "ecoregion",
                              ecoregion_of_cell100_full = prepped$ecoregion_of_cell100,
                              nregion = prepped$nregion)

# ---- baseline (bobcat) truth, then APPLY ABUNDANCE SCALING ----------------
# year_region (NOT a per-cell year_effect vector) is the real init here --
# see file header note.
truth <- true_param_list_ecoregion(prepped$real_post_means, year_region_true)
truth <- scale_truth_abundance(truth,
                               occ_shift      = abn$occ_shift,
                               count_log_mult = abn$count_log_mult,
                               label          = this$abundance)

set.seed(20260712 + row_id)   # unique per array task -> independent simulated data per replicate

# ---- simulate data, record realized information, fit BOTH models ----------
res <- run_one_replicate(
  rep_id      = this$rep,
  model_code  = model_code_ecoregion,
  constants   = cl$constants,
  inat_effort = cl$inat_effort,
  y_ncol      = cl$y_ncol,
  truth       = truth,
  base_inits  = prepped$base_inits,
  metrics_fn  = compute_ecoregion_metrics,   # per-region bias/cov/frac_sign_correct + WAIC
  trend_inits = list(year_beta = 0, year_var = 0, sigma_region = 1,
                    year_region = rep(0, prepped$nregion)),
  extra_monitors = c("year_region", "sigma_region"),
  fit_null_scalar = TRUE,                     # also fit the no-year_region model for WAIC compare
  year_region_true = year_region_true, ecoregion_levels = prepped$ecoregion_levels)

# attach realized information + design tags for the collector
# (res$sim_data is a list-column -- [[1]] unwraps the one row's value)
info <- summarize_simulated_information(res$sim_data[[1]], cl$constants)
res$abundance      <- this$abundance
res$scenario       <- this$scenario
res$occ_shift      <- abn$occ_shift
res$count_log_mult <- abn$count_log_mult
res$information    <- info

OUT <- file.path(PROTO, "sim_results_abundance",
                 paste0(this$abundance, "_", this$scenario))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(OUT, sprintf("rep_%03d.RDS", this$rep)))
cat(sprintf("status: %s  | cam_det=%s inat=%s\n", res$status,
            round(info$cam_detections_total), round(info$inat_count_total)))
