#!/usr/bin/env Rscript
# =============================================================================
# 01i_run_estimator_sweep.R
# -----------------------------------------------------------------------------
# Estimator-comparison sweep: does moving the camera submodel to array level,
# or to Royle-Nichols at array level, change our ability to detect a TREND
# when integrated with iNaturalist data?
#
# DESIGN
#   estimator  x  abundance  x  scenario  x  replicate
#   3 arms        3 levels      2            n_rep
#
#   estimators: camera_occ (baseline, current production estimator)
#               array_occ  (cameras as replicates within array-year)
#               array_rn   (Royle-Nichols; abundance drives detection)
#
# WHY ABUNDANCE IS CROSSED IN, NOT FIXED
# The saturation problem RN is meant to solve IS abundance-dependent: camera
# information about lambda scales as 1 - psi, so occupancy degrades as a
# species becomes common while RN should not. A sweep at one abundance would
# answer the wrong question -- it could show RN offering nothing (at low
# abundance, where occupancy is fine) or everything (at high abundance) and
# both would be artifacts of the chosen level. The abundance ladder is the axis
# along which the estimators are expected to DIVERGE.
#
# WHY BOTH SCENARIOS
#   varying: true regional trend present -> gives POWER (and bias/coverage)
#   null:    no trend                    -> gives the FALSE-POSITIVE rate
# Reporting power without the null rate is how a method that "detects" trends
# everywhere looks good.
#
# RESEEDING -- READ THIS BEFORE EDITING
# build_reduced_constants() calls set.seed(seed) INTERNALLY (sim_helpers.R:97)
# and simulate_replicate_data() does not reseed. Without an explicit reseed
# AFTER build_reduced_constants() and BEFORE the fit, every array task
# simulates the IDENTICAL dataset and the sweep is silently n=1 per cell.
# This project has shipped that bug twice (job 436723, and again in
# 01e_run_abundance_sweep.R). The reseed on line marked RESEED below is not
# optional and must not be moved above build_reduced_constants().
#
# NO CHUNK-RESUME
# fit_replicate() uses a single Cmcmc$run() call. Do not introduce chunked
# resumption here: the project's checkpoint/restore mechanism is known broken
# (it restarts from inits at every boundary), so a chunked run would be
# uninterpretable.
#
# PREREQUISITES -- two additions are needed before this script can run
#
# 1. inputs$site_array  (in sim_inputs.RDS)
#    The real array label per camera site, from umflist.RDS siteCovs
#    subproject_name, spatially split per the agreed rule. 00_prep_sim_inputs.R
#    does not currently carry it. The array arms REFUSE to run without it
#    rather than silently substituting the spatial approximation in
#    sim_helpers_array.R -- that approximation exists for development only.
#
#    NOTE on array construction: subproject_name x year is necessary but NOT
#    sufficient. 27% of subprojects span multiple years (hence the x year), AND
#    the largest array-year has 163 cameras spanning 162 km, which is an
#    administrative aggregate rather than an array. Oversized array-years must
#    be spatially split before use. That construction requires PI sign-off on
#    the resulting unit sizes and extents (see estimator_sweep_spec.md).
#
# 2. constants$site_keep  (from build_reduced_constants)
#    The index of retained sites, needed to align site_array with the reduced
#    design. build_reduced_constants() computes site_keep internally but does
#    not currently return it; it must be added to its return list.
#
# USAGE
#   Rscript 01i_run_estimator_sweep.R <row_id>
# where row_id indexes the design_df below (1..nrow). On Slurm, row_id is
# SLURM_ARRAY_TASK_ID. Locally, loop or use parallel::mclapply with
# mc.cores <= 4 (each fit peaks ~8 GB; memory binds before cores).
# =============================================================================

suppressPackageStartupMessages({
  library(nimble)
  library(nimbleEcology)
})

args   <- commandArgs(trailingOnly = TRUE)
row_id <- if (length(args) >= 1) as.integer(args[1]) else 1L

BASE <- Sys.getenv("SIM_BASE", ".")
# calcIntensity_SVC() / calcIntensity_noSVC() are nimbleFunctions defined in
# integration_helper.R, NOT in any sim_helpers file. Every model_code here --
# all three arms share the byte-identical iNat intensity block -- calls
# calcIntensity_SVC, so without this source() nimbleModel() fails at "Defining
# model" with "R function 'calcIntensity_SVC' ... does not exist" for ALL
# THREE arms. 01e_run_abundance_sweep.R:50 sources the same file from the same
# per-species location; this driver simply omitted it.
PROJ_ROOT <- normalizePath(file.path(BASE, ".."), mustWork = FALSE)
source(file.path(PROJ_ROOT, "HPC", "bobcat", "integration_helper.R"))
source(file.path(BASE, "sim_helpers.R"))
source(file.path(BASE, "sim_helpers_abundance.R"))
source(file.path(BASE, "sim_helpers_array.R"))
source(file.path(BASE, "sim_helpers_estimator_metrics.R"))
source(file.path(BASE, "model_code_national_scalar.R"))
source(file.path(BASE, "model_code_array_occ.R"))
source(file.path(BASE, "model_code_array_rn.R"))

DESIGN_SEED <- 20260823L
N_REP       <- as.integer(Sys.getenv("N_REP", "20"))
N_BURNIN    <- as.integer(Sys.getenv("N_BURNIN", "2000"))
N_ITER      <- as.integer(Sys.getenv("N_ITER", "8000"))
ARRAY_RADIUS_KM   <- as.numeric(Sys.getenv("ARRAY_RADIUS_KM", "5.0"))
ARRAY_MAX_PER     <- as.integer(Sys.getenv("ARRAY_MAX_PER", "25"))
OUTDIR      <- file.path(BASE, "estimator_sweep_out")
dir.create(OUTDIR, showWarnings = FALSE)

# ---- design ----------------------------------------------------------------
#' @name build_design_df
#' @description One row per (estimator x abundance x scenario x replicate).
#' @param n_rep Integer, replicates per cell.
#' @return data.frame with row_id and the four design factors.
build_design_df <- function(n_rep = N_REP) {
  d <- expand.grid(
    rep_id    = seq_len(n_rep),
    scenario  = c("varying", "null"),
    abundance = c("bobcat_baseline", "moderate", "common_deerlike"),
    estimator = c("camera_occ", "array_occ", "array_rn"),
    stringsAsFactors = FALSE
  )
  d <- d[order(d$estimator, d$abundance, d$scenario, d$rep_id), ]
  d$row_id <- seq_len(nrow(d))
  rownames(d) <- NULL
  d
}

design_df <- build_design_df()
if (row_id > nrow(design_df)) {
  stop(sprintf("row_id %d exceeds design size %d", row_id, nrow(design_df)))
}
cfg <- design_df[row_id, ]
cat(sprintf("[row %d/%d] estimator=%s abundance=%s scenario=%s rep=%d\n",
            row_id, nrow(design_df), cfg$estimator, cfg$abundance,
            cfg$scenario, cfg$rep_id))

outfile <- file.path(OUTDIR, sprintf("row_%05d.rds", row_id))
if (file.exists(outfile)) {
  cat("already done, exiting\n"); quit(save = "no", status = 0)
}

# ---- one replicate ---------------------------------------------------------
#' @name run_estimator_replicate
#' @description Simulate one dataset at the given abundance and trend
#'   scenario, then fit it with the requested camera estimator. The SIMULATION
#'   is always camera-level (that is how the real data are generated); only the
#'   FIT changes, with array arms aggregating the simulated camera data first.
#'   This is deliberate: it mirrors the real decision, which is how to model
#'   data we already have, not what data to collect.
#' @param cfg One row of design_df.
#' @return A one-row data.frame of metrics, or an "Errored" row.
run_estimator_replicate <- function(cfg) {
  # FILE NAME -- the prep step writes "prepped_sim_inputs.RDS" (00_prep_sim_inputs.R:105)
  # and every working driver in this project reads that name (01e_run_abundance_sweep.R:56,
  # 01_run_sim_validation.R:71). "sim_inputs.RDS" does not exist on Hazel, so every
  # task died here before running. The _with_array variant is that same object plus
  # inputs$site_array (the 5 km array-year label per camera).
  inputs <- readRDS(file.path(BASE, "prepped_sim_inputs_with_array.RDS"))

  # build_reduced_constants() returns a WRAPPER list --
  # list(constants=, inat_effort=, y_ncol=, cell50_keep=, site_keep=) -- not
  # the NIMBLE constants themselves. Every other driver in this project
  # (01_run_sim_validation.R etc.) unwraps this via reduced$constants /
  # reduced$inat_effort; this script previously used the wrapper object
  # directly as `constants`, so constants$adj/$num/$J/$year_occ all resolved
  # to NULL (those live inside reduced$constants, not at the wrapper's top
  # level) -- silently breaking model construction for ALL THREE arms, not
  # just the array arms. Fixed to unwrap explicitly. Also: the reduced
  # inat_effort matrix (subset to cell50_keep, matching constants$ncell50)
  # must be used everywhere a full-size inputs$inat_effort was previously
  # passed -- the two have different row counts and passing the wrong one
  # would be a dimension mismatch against the reduced design.
  # Site reduction: sample WHOLE array-year units, not independent sites.
  # Independent-site sampling (build_reduced_constants()'s own default)
  # scatters cameras across array-year units and leaves almost none of them
  # intact -- confirmed on the real bobcat structure: n_site_keep=700 turns
  # 3,327 real array-years (median 4 cameras, mean 6.17) into 107 surviving
  # array-years at a degenerate median 2 / mean 2.36, using only 253 of the
  # 700 sampled sites. That would validate the array arms against arrays
  # that don't look like real arrays. sample_sites_by_array() (in
  # sim_helpers_array.R) instead selects complete array-year units up to
  # the same ~700-camera budget, so build_reduced_constants() is called
  # with a fixed site_keep_override rather than its own random sampling.
  # This changes camera_occ's retained-site SET (compared to the previous
  # buggy behavior) but not its estimator logic -- camera_occ fits whatever
  # sites are retained either way, so the comparison across arms stays
  # apples-to-apples for THIS driver's own design.
  array_sample <- sample_sites_by_array(
    array_field = inputs$site_array,
    site_year   = inputs$constants_list$year_occ,
    region_vec  = inputs$cell100_geo$region[match(inputs$constants_list$cell,
                                                  inputs$cell100_geo$cell100)],
    n_site_target = 700, drop_singletons = TRUE, seed = DESIGN_SEED
  )
  cat(sprintf("[site sampling] %d arrays kept, %d cameras kept\n",
              array_sample$n_arrays_kept, array_sample$n_sites_kept))

  reduced      <- build_reduced_constants(
    # FIELD NAMES -- prepped_sim_inputs.RDS carries constants_list /
    # inat_effort_real / real_y_template. inputs$cl, inputs$inat_effort and
    # inputs$y_template are all NULL, which build_reduced_constants() would
    # have received silently. Names verified against the RDS itself and
    # against 01e_run_abundance_sweep.R:93-94, which uses the same three.
    cl = inputs$constants_list, inat_effort_real = inputs$inat_effort_real,
    real_y_template = inputs$real_y_template, cell100_geo = inputs$cell100_geo,
    seed = DESIGN_SEED,
    site_keep_override = array_sample$site_keep
  )
  constants    <- reduced$constants
  inat_effort  <- reduced$inat_effort

  # RESEED -- must stay AFTER build_reduced_constants (see header)
  set.seed(DESIGN_SEED + row_id)

  year_effect_true <- make_true_year_effect(
    inputs$cell100_geo, constants$adj, constants$num,
    amplitude = if (cfg$scenario == "null") 0 else 0.3
  )
  truth <- true_param_list(inputs$real_post_means, year_effect_true)
  truth <- scale_truth_abundance(truth, label = cfg$abundance)

  # inputs$y_template does not exist (see the field-name fix above -- the
  # real field is inputs$real_y_template); build_reduced_constants() already
  # computes y_ncol correctly from it (reduced$y_ncol = ncol(real_y_template))
  # and returns it in the wrapper, so use that directly rather than
  # re-deriving it from a nonexistent field (which would silently produce
  # ncol(NULL) = NULL and fail deep inside simulate_replicate_data()'s
  # dimensions= argument).
  sim <- simulate_replicate_data(model_code_national_scalar, constants,
                                 inat_effort,
                                 y_ncol = reduced$y_ncol, truth)

  if (cfg$estimator == "camera_occ") {
    fit_constants <- constants
    fit_data      <- list(y = sim$y, y_inat = sim$y_inat)
    model_code    <- model_code_national_scalar
    array_diag    <- NULL
  } else {
    # Array structure comes from the REAL subproject_name x year field,
    # carried through 00_prep_sim_inputs.R as inputs$site_array (see
    # PREREQUISITE in the header). Real structure is preferable to a spatial
    # approximation here: the field exists (943 subprojects, 99.9% populated),
    # and using it means the simulated array sizes match what a real fit would
    # actually face. build_reduced_constants() samples sites, so subset the
    # array labels to the retained sites via the same index.
    if (is.null(inputs$site_array)) {
      stop("inputs$site_array missing -- re-run 00_prep_sim_inputs.R with the ",
           "array-structure addition (see estimator_sweep_spec.md). Refusing ",
           "to substitute a spatial approximation silently.")
    }
    array_id <- assign_arrays_from_field(
      # site_keep is a WRAPPER-level field on reduced (from
      # build_reduced_constants()'s return list), not inside
      # reduced$constants -- confirmed by inspection (constants$site_keep
      # is NULL; reduced$site_keep is the retained-site index).
      array_field = inputs$site_array[reduced$site_keep],
      site_year   = constants$year_occ,
      drop_singletons = TRUE
    )
    agg <- aggregate_to_array(sim$y, constants$J, array_id)
    fit_constants <- build_array_constants(constants, array_id, agg)
    fit_data      <- list(w = agg$w, y_inat = sim$y_inat)
    model_code    <- if (cfg$estimator == "array_occ") model_code_array_occ
                     else model_code_array_rn
    array_diag    <- summarize_array_structure(agg, constants, array_id)
    if (cfg$estimator == "array_rn") {
      # start N_a above the observed detections so the chain is feasible
      fit_data$N_a_init <- pmax(1L, rowSums(agg$w, na.rm = TRUE))
    }
  }

  t0 <- Sys.time()
  # None of the three sweep models (national_scalar, array_occ, array_rn)
  # has a CAR-based year_effect/tau_year trend field -- all three use the
  # scalar year_beta/year_var trend block directly (see the byte-identical
  # trend block, estimator_sweep_spec.md). fit_replicate()'s DEFAULT
  # trend_inits/extra_monitors are for the CAR ecoregion model and reference
  # nodes these three models do not have; passing them unmodified makes
  # configureMCMC()$addMonitors() error with "These variables are not in
  # model: year_effect,tau_year" for every arm. Override to match what
  # sim_helpers.R's own docstring recommends for the national-scalar model.
  samples <- fit_replicate(model_code, fit_constants, inat_effort,
                           fit_data, base_inits = inputs$base_inits,
                           n_burnin = N_BURNIN, n_iter = N_ITER,
                           trend_inits = list(year_beta = 0, year_var = 0),
                           extra_monitors = character(0))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  m <- compute_estimator_metrics(samples, truth, constants, inputs$cell100_geo)
  cbind(cfg, m, elapsed_sec = elapsed, status = "OK",
        if (is.null(array_diag)) data.frame(n_arrays = NA_integer_,
                                            mean_prop_detect = NA_real_,
                                            frac_saturated = NA_real_)
        else array_diag[, c("n_arrays", "mean_prop_detect", "frac_saturated")])
}

# ---- run, never letting one failure kill the batch --------------------------
res <- tryCatch(
  run_estimator_replicate(cfg),
  error = function(e) {
    cat("ERRORED:", conditionMessage(e), "\n")
    cbind(cfg, status = "Errored", error_msg = conditionMessage(e))
  }
)

saveRDS(res, outfile)
cat("wrote", outfile, "\n")
