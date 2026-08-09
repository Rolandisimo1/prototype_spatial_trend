# 01h_run_effort_secondary.R
# ---------------------------------------------------------------------------
# Secondary per-region effort scenario (01g_run_effort_sweep.R's companion):
# holds cameras at 1x (700 sites, no resampling) and GLOBAL iNat effort at
# 1x, but overrides ONE focal ecoregion's iNat effort down to 0.1x -- a
# stand-in for "an under-sampled region/country" at the same (bobcat) animal
# abundance as everywhere else in the design. Every OTHER ecoregion,
# including the comparison region below, is untouched (stays at the real
# 1x effort) -- so FOCAL vs COMPARISON is a same-replicate, same-run
# contrast, not a cross-run one.
#
# FOCAL_ECOREGION = "NORTH AMERICAN DESERTS" (LOWEST real mean iNat effort
#   of all 8 ecoregions, ~11.1 in the reduced design -- verified directly by
#   inspecting build_reduced_constants()$inat_effort per region before
#   picking, not guessed. NOTE: an earlier choice of "MARINE WEST COAST
#   FOREST" -- picked because it has the smallest real CAMERA pool -- turned
#   out to have the HIGHEST real iNat effort (~148, likely reflecting dense
#   citizen-scientist populations near Pacific coast cities), the opposite
#   of a low-effort stand-in. Camera density and iNat observer density are
#   NOT the same thing in this dataset -- caught via the smoke test's
#   single-replicate result looking backwards (focal region showed MORE
#   simulated iNat counts than the comparison region despite the 0.1x cut),
#   traced to the real data before assuming it was a bug, and fixed here by
#   picking the actual lowest-real-effort region instead.
# COMPARISON_ECOREGION = "MARINE WEST COAST FOREST" (HIGHEST real mean iNat
#   effort, ~148 -- kept at full 1x, now used as the well-sampled reference
#   it was originally intended to be an alternative to).
#
# Depends on the same sources as 01g_run_effort_sweep.R.
# Writes: sim_results_effort_secondary/rep_<id>.RDS
# ---------------------------------------------------------------------------

suppressMessages({ library(nimble); library(nimbleEcology); library(dplyr) })

PROTO <- getwd()
PROJ  <- dirname(PROTO)
source(file.path(PROJ, "HPC", "bobcat", "integration_helper.R"))
source(file.path(PROTO, "sim_helpers.R"))
source(file.path(PROTO, "sim_helpers_effort.R"))
source(file.path(PROTO, "model_code_ecoregion_trend.R"))

SECONDARY_SEED <- 20260717
FOCAL_ECOREGION_NAME      <- "NORTH AMERICAN DESERTS"
COMPARISON_ECOREGION_NAME <- "MARINE WEST COAST FOREST"
FOCAL_MULT <- 0.1

prepped <- readRDS(file.path(PROTO, "prepped_sim_inputs.RDS"))
stopifnot("ecoregion_id" %in% names(prepped$cell100_geo),
          "ecoregion_of_cell100" %in% names(prepped))

focal_ecoregion_id <- which(prepped$ecoregion_levels == FOCAL_ECOREGION_NAME)
stopifnot(length(focal_ecoregion_id) == 1)
cat(sprintf("FOCAL ecoregion: %s (id=%d), effort x%.2f | COMPARISON ecoregion: %s (unchanged, x1)\n",
            FOCAL_ECOREGION_NAME, focal_ecoregion_id, FOCAL_MULT, COMPARISON_ECOREGION_NAME))

N_REPS <- 30
arg <- Sys.getenv("SLURM_ARRAY_TASK_ID", NA)
if (is.na(arg)) { a <- commandArgs(trailingOnly = TRUE); arg <- if (length(a)) a[1] else NA }
stopifnot("Need SLURM_ARRAY_TASK_ID or CLI row id (1 rep per task)" = !is.na(arg))
rep_id <- as.integer(arg)
stopifnot(rep_id >= 1, rep_id <= N_REPS)
cat(sprintf("=== secondary effort scenario: rep=%d/%d ===\n", rep_id, N_REPS))

# ---- fixed bobcat abundance, varying trend scenario, 1x cameras -----------
year_region_true <- make_true_year_region(prepped$cell100_geo, scenario = "varying")
truth <- true_param_list_ecoregion(prepped$real_post_means, year_region_true)

cl <- build_reduced_constants(prepped$constants_list, prepped$inat_effort_real,
                              prepped$real_y_template, prepped$cell100_geo,
                              n_site_keep = 700,
                              stratify_by = "ecoregion",
                              ecoregion_of_cell100_full = prepped$ecoregion_of_cell100,
                              nregion = prepped$nregion,
                              site_replace = FALSE)

# ---- iNat effort: global 1x, FOCAL ecoregion overridden to 0.1x -----------
inat_effort_scaled <- scale_inat_effort(cl$inat_effort, cl$constants$inat_cell100,
                                        prepped$ecoregion_of_cell100,
                                        global_mult = 1,
                                        focal_ecoregion_id = focal_ecoregion_id,
                                        focal_mult = FOCAL_MULT)

set.seed(SECONDARY_SEED + rep_id)   # per-rep reseed AFTER build_reduced_constants, BEFORE run_one_replicate

res <- run_one_replicate(
  rep_id      = rep_id,
  model_code  = model_code_ecoregion,
  constants   = cl$constants,
  inat_effort = inat_effort_scaled,
  y_ncol      = cl$y_ncol,
  truth       = truth,
  base_inits  = prepped$base_inits,
  metrics_fn  = compute_ecoregion_metrics,
  trend_inits = list(year_beta = 0, year_var = 0, sigma_region = 1,
                    year_region = rep(0, prepped$nregion)),
  extra_monitors = c("year_region", "sigma_region"),
  year_region_true = year_region_true, ecoregion_levels = prepped$ecoregion_levels)

info_by_region <- summarize_simulated_information_by_region(
  res$sim_data[[1]], cl$constants, prepped$ecoregion_of_cell100, prepped$ecoregion_levels)

res$focal_ecoregion      <- FOCAL_ECOREGION_NAME
res$comparison_ecoregion <- COMPARISON_ECOREGION_NAME
res$focal_mult           <- FOCAL_MULT
res$info_by_region       <- list(info_by_region)

OUT <- file.path(PROTO, "sim_results_effort_secondary")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(OUT, sprintf("rep_%03d.RDS", rep_id)))
cat(sprintf("status: %s | bias_all=%s rmse_all=%s\n",
            res$status, round(res$bias_all, 4), round(res$rmse_all, 4)))
