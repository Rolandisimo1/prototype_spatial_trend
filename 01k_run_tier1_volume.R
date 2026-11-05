#!/usr/bin/env Rscript
# =============================================================================
# 01k_run_tier1_volume.R
# -----------------------------------------------------------------------------
# Tier 1 of SEED_data_volume_simulation.md: statistical power to detect a
# population trend, the range sizes at which it can be done, and what extra
# data would be needed. Production estimator only (camera_occ =
# model_code_national_scalar, unmodified); this study varies DATA, not the
# estimator. Design helpers are in sim_helpers_tier1.R.
#
# REDUCED DESIGN (the full A x B x C factorial at 3 effect sizes is too large)
#   Block RW  -- range size x camera window, fully crossed, at one mid-fleet
#                data volume (N_SITE_MID cameras, INAT_MID expected iNat
#                records; both totals held fixed across range sizes):
#                  range   25 / 50 / 100 / 300 / 900 cell100s
#                  window  100 / 60 / 40 / 20 % of the 18-yr iNat window,
#                          camera years at the END (18 / 11 / 7 / 4 yr)
#   Block VOL -- camera volume x iNat volume, fully crossed, at one range and
#                window (RANGE_VOL cells, WINDOW_VOL):
#                  cameras 1,000 / 3,000 / 9,000 / 25,000 sites
#                  iNat    1,000 / 5,000 / 20,000 / 100,000 expected records
#   Level anchors (fleet_volume.csv, all 34 v2u bundles): camera deployments
#   min 688 / Q1 2,933 / median 7,765 / max 24,463; iNat records min 211 /
#   Q1 7,343 / median 17,386 / max 188,941. Mid volume = 9,000 cameras and
#   20,000 records, which is also a VOL cell, so the blocks share a point.
#   Both blocks: effect null / -2 / -5 / -15 %/yr  x  N_REP replicates.
#
# THE NULL zeroes year_beta AND year_var (tier1_trend_truth(r = 0)); there is
# no spatial deviation field in this model to zero.
#
# SCORING: both year_beta (yb_*, the reported camera-anchored trend) and
# total_var_beta (tvb_*) -- see sim_helpers_estimator_metrics.R.
#
# RESEEDING: the design (which sites, cell50s, site years) is drawn under
# DESIGN_SEED inside build_tier1_constants() and is shared by all replicates
# of a design cell; the data are simulated after set.seed(RUN_SEED + row_id).
# The RESEED line must stay after build_tier1_constants() (see 01i header for
# the twice-shipped n=1 bug).
#
# USAGE  TIER1_BLOCK=RW|VOL  Rscript 01k_run_tier1_volume.R <row_id>
#        Rscript 01k_run_tier1_volume.R design   (writes the design CSVs, exits)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
BASE <- Sys.getenv("SIM_BASE", ".")
DESIGN_SEED <- 20261017L
RUN_SEED    <- 20261100L
N_REP       <- as.integer(Sys.getenv("N_REP", "30"))
N_BURNIN    <- as.integer(Sys.getenv("N_BURNIN", "2000"))
N_ITER      <- as.integer(Sys.getenv("N_ITER", "8000"))

RANGE_LEVELS  <- c(25L, 50L, 100L, 300L, 900L)
WINDOW_LEVELS <- c(1.0, 0.6, 0.4, 0.2)
CAM_LEVELS    <- c(1000L, 3000L, 9000L, 25000L)
INAT_LEVELS   <- c(1000, 5000, 20000, 100000)
EFFECT_LEVELS <- c(0, -0.02, -0.05, -0.15)
N_SITE_MID    <- 9000L
INAT_MID      <- 20000
RANGE_VOL     <- 300L
WINDOW_VOL    <- 0.4

#' @name build_tier1_design
#' @description One row per design cell x effect x replicate, for one block.
#'   Row order is fixed by an explicit sort, so a row_id always denotes the
#'   same cell for a given N_REP. A change of N_REP or of any level vector
#'   changes the mapping: use a fresh output directory.
#' @param block "RW" or "VOL".
#' @return data.frame with row_id and the design factors.
build_tier1_design <- function(block) {
  d <- if (block == "RW") {
    expand.grid(rep_id = seq_len(N_REP), effect = EFFECT_LEVELS,
                window_frac = WINDOW_LEVELS, range_cells = RANGE_LEVELS,
                n_site = N_SITE_MID, inat_records = INAT_MID)
  } else if (block == "VOL") {
    expand.grid(rep_id = seq_len(N_REP), effect = EFFECT_LEVELS,
                inat_records = INAT_LEVELS, n_site = CAM_LEVELS,
                range_cells = RANGE_VOL, window_frac = WINDOW_VOL)
  } else stop("TIER1_BLOCK must be RW or VOL")
  d <- d[order(d$range_cells, -d$window_frac, d$n_site, d$inat_records, -d$effect, d$rep_id), ]
  d$block <- block
  d$row_id <- seq_len(nrow(d)); rownames(d) <- NULL
  d[, c("block", "row_id", "range_cells", "window_frac", "n_site", "inat_records", "effect", "rep_id")]
}

if (length(args) && args[1] == "design") {
  for (b in c("RW", "VOL")) {
    d <- build_tier1_design(b); f <- sprintf("tier1_design_%s.csv", b)
    write.csv(d, f, row.names = FALSE); cat(b, nrow(d), "rows ->", f, "\n")
  }
  quit(save = "no")
}

# design mode above needs no NIMBLE; everything below does
suppressPackageStartupMessages({ library(nimble); library(nimbleEcology) })

PROJ_ROOT <- normalizePath(file.path(BASE, ".."), mustWork = FALSE)
.ih <- c(file.path(PROJ_ROOT, "HPC", "bobcat", "integration_helper.R"),
         file.path(BASE, "HPC", "bobcat", "integration_helper.R"))
.ih <- .ih[file.exists(.ih)]
if (!length(.ih)) stop("integration_helper.R not found (calcIntensity_SVC would be undefined)")
source(.ih[1])
source(file.path(BASE, "sim_helpers.R"))
source(file.path(BASE, "sim_helpers_effort.R"))
source(file.path(BASE, "sim_helpers_estimator_metrics.R"))
source(file.path(BASE, "sim_helpers_tier1.R"))
source(file.path(BASE, "model_code_national_scalar.R"))

BLOCK  <- Sys.getenv("TIER1_BLOCK", "RW")
row_id <- as.integer(args[1])
design <- build_tier1_design(BLOCK)
if (is.na(row_id) || row_id > nrow(design)) stop("row_id out of range 1..", nrow(design))
cfg <- design[row_id, ]
OUTDIR <- Sys.getenv("TIER1_OUTDIR", file.path(BASE, paste0("tier1_out_", BLOCK)))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(OUTDIR, sprintf("row_%05d.rds", row_id))
cat(sprintf("[%s row %d/%d] range=%d window=%.0f%% n_site=%d inat_records=%g effect=%+.0f%%/yr rep=%d\nOUTDIR %s\n",
            BLOCK, row_id, nrow(design), cfg$range_cells, 100 * cfg$window_frac, cfg$n_site,
            cfg$inat_records, 100 * cfg$effect, cfg$rep_id, OUTDIR))
if (file.exists(outfile)) { cat("already done, exiting\n"); quit(save = "no", status = 0) }

run_tier1_replicate <- function(cfg) {
  inputs <- readRDS(file.path(BASE, "prepped_sim_inputs_with_array.RDS"))
  cl <- inputs$constants_list
  seed_cell <- choose_range_seed(cl, min(RANGE_LEVELS), max(RANGE_LEVELS))

  des <- build_tier1_constants(inputs, cfg$range_cells, cfg$window_frac, cfg$n_site,
                               seed_cell, design_seed = DESIGN_SEED)
  constants <- des$constants

  truth0 <- true_param_list(inputs$real_post_means, rep(0, cl$ncell100))
  truth0$year_effect <- NULL; truth0$tau_year <- NULL   # not nodes of this model
  truth0 <- tier1_trend_truth(truth0, 0, des$s_occ, des$s_vals)
  cal <- calibrate_inat_effort(des, truth0, cfg$inat_records)
  des$diag$d_inat_mult <- cal$mult
  des$diag$d_inat_expected_null <- cal$expected_records

  # Axis-live checks: the design must actually carry the level it is labelled
  # with (the abundance ladder once recorded a label and applied nothing).
  stopifnot(constants$ncell100 == cl$ncell100,                   # full CAR grid, as in production
            all(constants$cell %in% des$range_keep), all(constants$inat_cell100 %in% des$range_keep),
            length(des$range_keep) == cfg$range_cells,
            constants$nsite == cfg$n_site,
            des$diag$d_cam_nyears == length(camera_window_years(cfg$window_frac)),
            des$diag$d_cam_year_max == 2025L,
            abs(cal$expected_records / cfg$inat_records - 1) < 1e-8)

  # RESEED -- must stay after build_tier1_constants()
  set.seed(RUN_SEED + row_id)

  truth <- tier1_trend_truth(truth0, cfg$effect, des$s_occ, des$s_vals)
  if (cfg$effect == 0) stopifnot(truth$year_beta == 0, truth$year_var == 0)
  cam_pct_true <- 100 * (exp(truth$year_beta * des$s_occ) - 1)
  stopifnot(abs(cam_pct_true - 100 * cfg$effect) < 1e-8)
  des$inat_effort <- cal$effort

  sim <- simulate_replicate_data(model_code_national_scalar, constants, des$inat_effort,
                                 y_ncol = des$y_ncol, truth)
  realized <- data.frame(
    r_cam_det_sites = sum(rowSums(sim$y, na.rm = TRUE) > 0),
    r_cam_naive_occ = mean(rowSums(sim$y, na.rm = TRUE) > 0),
    r_cam_det_total = sum(sim$y, na.rm = TRUE),
    r_inat_records  = sum(sim$y_inat, na.rm = TRUE))
  cat(sprintf("[design] CAR nodes=%d ncell50=%d nsite=%d (unique %d of pool %d) cam years %d-%d s_occ=%.4f\n",
              constants$ncell100, constants$ncell50, constants$nsite, des$diag$d_nsite_unique,
              des$diag$d_site_pool, des$diag$d_cam_year_min, des$diag$d_cam_year_max, des$s_occ))
  cat(sprintf("[range] %d cells: iNat in %d, cameras in %d; iNat effort x%.4g of real for %g expected records\n",
              cfg$range_cells, des$diag$d_inat_cells, des$diag$d_cam_cells, cal$mult, cal$expected_records))
  cat(sprintf("[truth] year_beta=%.5f year_var=%.5f camera %+.2f%%/yr\n", truth$year_beta, truth$year_var, cam_pct_true))
  cat(sprintf("[data] detecting sites=%d (naive occ %.3f) camera detections=%d iNat records=%d\n",
              realized$r_cam_det_sites, realized$r_cam_naive_occ, realized$r_cam_det_total, realized$r_inat_records))

  inits <- inputs$base_inits
  t0 <- Sys.time()
  samples <- fit_replicate(model_code_national_scalar, constants, des$inat_effort,
                           list(y = sim$y, y_inat = sim$y_inat), base_inits = inits,
                           n_burnin = N_BURNIN, n_iter = N_ITER,
                           trend_inits = list(year_beta = 0, year_var = 0,
                                              sigma_year_beta = 1, sigma_year_var = 1),
                           extra_monitors = character(0))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  m <- compute_estimator_metrics(samples, truth, constants, inputs$cell100_geo)
  cbind(cfg, des$diag, realized, cam_pct_true = cam_pct_true,
        yv_true = truth$year_var, m, elapsed_sec = elapsed, n_burnin = N_BURNIN, n_iter = N_ITER,
        status = "OK", stringsAsFactors = FALSE)
}

res <- tryCatch(run_tier1_replicate(cfg), error = function(e) {
  cat("ERRORED:", conditionMessage(e), "\n")
  cbind(cfg, status = "Errored", error_msg = conditionMessage(e))
})
saveRDS(res, outfile)
cat("wrote", outfile, "\n")
