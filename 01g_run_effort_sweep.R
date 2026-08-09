# 01g_run_effort_sweep.R
# ---------------------------------------------------------------------------
# Effort sweep: how much SAMPLING (not animal abundance) does a bobcat-
# density species need for usable per-ecoregion trend resolution? Complement
# of 01e_run_abundance_sweep.R, which varied the animal at fixed sampling --
# this fixes the animal at real bobcat magnitude (occ_shift=0,
# count_log_mult=0) and varies:
#   CAMERAS (n_site_keep): 0.5x/1x/2x/4x of this project's established
#     700-site reduced baseline (NOT the literal 20,531-site full real
#     dataset -- confirmed with the user; see build_reduced_constants()'s
#     site_replace docstring). 0.5x/1x subsample without replacement (same
#     as every prior script); 2x/4x resample WITH replacement
#     (site_replace=TRUE), an explicit, flagged extrapolation beyond the
#     real camera density.
#   INAT (global effort multiplier): 0.25x/1x/4x/16x applied to the real
#     reduced-design effort matrix via scale_inat_effort().
# Full 4x4 cross x N_REPS reps, scenario FIXED at "varying" (real regional
# trend -- makes per-region RMSE meaningful; a "null" companion run to
# re-check false-positive behavior under low effort was discussed as
# optional and is NOT part of this core deliverable).
#
# Depends on (source in this order):
#   model_code_ecoregion_trend.R   (defines model_code_ecoregion)
#   sim_helpers.R                  (build_reduced_constants -- now with
#                                   site_replace/site_jitter_sd,
#                                   true_param_list_ecoregion,
#                                   make_true_year_region, run_one_replicate,
#                                   compute_ecoregion_metrics)
#   sim_helpers_effort.R           (scale_inat_effort,
#                                   summarize_simulated_information_by_region,
#                                   camera_levels_default, inat_levels_default)
#   prepped_sim_inputs.RDS
#
# GUARDRAILS carried from prior runs (both bugs this project has hit twice):
#   - per-rep reseed set.seed(EFFORT_SEED + row_id) AFTER
#     build_reduced_constants() (which reseeds internally with a FIXED seed
#     by design) and BEFORE run_one_replicate() -- the "Fixed-seed RNG leak"
#     bug class, see README bug #4 and abundance_sweep_seed_diagnosis.md.
#   - collector (02g) carries the same anti-degeneracy guard added to 02e/02f.
#
# Writes: sim_results_effort/<camera_level>_<inat_level>/rep_<id>.RDS
# ---------------------------------------------------------------------------

suppressMessages({ library(nimble); library(nimbleEcology); library(dplyr) })

PROTO <- getwd()
PROJ  <- dirname(PROTO)
source(file.path(PROJ, "HPC", "bobcat", "integration_helper.R"))
source(file.path(PROTO, "sim_helpers.R"))
source(file.path(PROTO, "sim_helpers_effort.R"))
source(file.path(PROTO, "model_code_ecoregion_trend.R"))   # defines model_code_ecoregion

EFFORT_SEED <- 20260716
SITE_JITTER_SD <- 0.05   # only applied when site_replace=TRUE (2x/4x camera levels)

prepped <- readRDS(file.path(PROTO, "prepped_sim_inputs.RDS"))
stopifnot("ecoregion_id" %in% names(prepped$cell100_geo),
          "ecoregion_of_cell100" %in% names(prepped))

# ---- design ---------------------------------------------------------------
N_REPS <- 30
camera_levels <- camera_levels_default()
inat_levels   <- inat_levels_default()

design <- expand.grid(camera_level = camera_levels$level,
                      inat_level   = inat_levels$level,
                      rep          = seq_len(N_REPS),
                      stringsAsFactors = FALSE)
design <- design[order(design$camera_level, design$inat_level, design$rep), ]
design$row <- seq_len(nrow(design))

# ---- pick this task's row --------------------------------------------------
arg <- Sys.getenv("SLURM_ARRAY_TASK_ID", NA)
if (is.na(arg)) { a <- commandArgs(trailingOnly = TRUE); arg <- if (length(a)) a[1] else NA }
stopifnot("Need SLURM_ARRAY_TASK_ID or CLI row id (1 combo per task)" = !is.na(arg))
row_id <- as.integer(arg)
this <- design[design$row == row_id, ]
cam  <- camera_levels[camera_levels$level == this$camera_level, ]
inat <- inat_levels[inat_levels$level == this$inat_level, ]
cat(sprintf("=== camera=%s (n_site_keep=%d, replace=%s)  inat=%s (mult=%.2f)  rep=%d (row %d/%d) ===\n",
            this$camera_level, cam$n_site_keep, cam$site_replace,
            this$inat_level, inat$global_mult, this$rep, row_id, nrow(design)))

# ---- fixed bobcat abundance, varying trend scenario ------------------------
year_region_true <- make_true_year_region(prepped$cell100_geo, scenario = "varying")
truth <- true_param_list_ecoregion(prepped$real_post_means, year_region_true)

# ---- reduced design at THIS row's camera density ---------------------------
cl <- build_reduced_constants(prepped$constants_list, prepped$inat_effort_real,
                              prepped$real_y_template, prepped$cell100_geo,
                              n_site_keep = cam$n_site_keep,
                              stratify_by = "ecoregion",
                              ecoregion_of_cell100_full = prepped$ecoregion_of_cell100,
                              nregion = prepped$nregion,
                              site_replace = cam$site_replace,
                              site_jitter_sd = SITE_JITTER_SD)

# ---- scale iNat effort at THIS row's inat level -----------------------------
inat_effort_scaled <- scale_inat_effort(cl$inat_effort, cl$constants$inat_cell100,
                                        prepped$ecoregion_of_cell100,
                                        global_mult = inat$global_mult)

set.seed(EFFORT_SEED + row_id)   # unique per array task -> independent simulated data per replicate

# ---- fit the ecoregion model, score per-region recovery --------------------
res <- run_one_replicate(
  rep_id      = this$rep,
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

# per-region realized information (res$sim_data is a list-column -- [[1]] unwraps it)
info_by_region <- summarize_simulated_information_by_region(
  res$sim_data[[1]], cl$constants, prepped$ecoregion_of_cell100, prepped$ecoregion_levels)

res$camera_level  <- this$camera_level
res$n_site_keep   <- cam$n_site_keep
res$site_replace  <- cam$site_replace
res$inat_level    <- this$inat_level
res$global_mult   <- inat$global_mult
res$info_by_region <- list(info_by_region)

OUT <- file.path(PROTO, "sim_results_effort",
                 paste0(this$camera_level, "_", this$inat_level))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
saveRDS(res, file.path(OUT, sprintf("rep_%03d.RDS", this$rep)))
cat(sprintf("status: %s | bias_all=%s rmse_all=%s frac_sign_correct=%s\n",
            res$status, round(res$bias_all, 4), round(res$rmse_all, 4), round(res$frac_sign_correct, 3)))
