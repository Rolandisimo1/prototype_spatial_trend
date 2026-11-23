#!/usr/bin/env Rscript
# =============================================================================
# 01l_run_camanch_regional.R
# -----------------------------------------------------------------------------
# Validation simulation for the camera-anchored regional trend model
# (model_code_ecoregion_camanchored.R; SPEC_camera_anchored_regional.md, Step 2).
# Data are generated FROM the new model, so a regional trend reaches both the
# camera and the iNaturalist stream.
#
# DESIGN (4 scenarios x 3 camera x 3 iNat x N_REP = 720 fits at N_REP = 20; Roland 2026-11-23)
#   scenario       year_region        year_region_var     tests
#   NULL           0                  0                   false-positive rate of the regional term
#   REAL           N(0, 0.2) draws    0                   power, bias, coverage of camera-anchored regional trends
#   SPURIOUS       0                  N(0, 0.2) draws     does iNat-only regional change leak into year_region?
#   REAL_SPURIOUS  N(0, 0.2)          N(0, 0.2), indep.   the same, with a real signal present
#   camera sites   1,000 / 3,000 / 9,000 DISTINCT real sites (no resampling with replacement)
#   iNat records   5,000 / 20,000 / 100,000 expected records
#   one abundance level (bobcat-like, the real posterior means), one range:
#   RANGE_CELLS = 900, the full grid, so all seven real EPA Level I ecoregions
#   carry data (the 300-cell Tier 1 range put 242 of its 300 cells in Eastern
#   Temperate Forests; Roland 2026-11-23: 900 cells at 20 replicates).
#   SAME-SPAN WINDOW (Roland 2026-11-23): camera AND iNat both cover the last
#   WINDOW_YEARS = 5 calendar years (2021-2025), as in the reported win5 fits.
#   The iNat side of the Tier 1 design is cut to those years and year_vals is
#   re-z-scored on them, exactly as the windowed prep does. Camera and iNat year
#   steps are then ~0.71 and ~0.63 per year (ratio ~0.89; fitted deer win5 0.91),
#   so REAL with year_region_var = 0 leaves only a ~10% difference between the
#   streams' implied %/yr. Both per-stream %/yr values are recorded per region.
# The national trend is null (year_beta = year_var = 0) so the regional terms
# are the only trend signal. Regional truths are drawn per replicate AFTER the
# reseed, only for ecoregions that hold data in the range; others are 0.
#
# Built on the Tier 1 design infrastructure (build_tier1_constants(),
# calibrate_inat_effort(), sim_helpers_tier1.R), unmodified. The no-replacement
# rule is enforced here: the run stops if any site would be drawn twice.
#
# RESEEDING: the design is drawn under DESIGN_SEED inside build_tier1_constants()
# and shared by all replicates of a design cell; truths and data are drawn after
# set.seed(RUN_SEED + row_id). The RESEED line must stay after
# build_tier1_constants() (twice-shipped n=1 bug, see 01i header).
#
# SMOKE MODE (CAMANCH_SMOKE=1): rows are read from camanch_smoke_rows.txt as
# "<scenario> <n_site> <inat_records> <arm>", arm = joint | camonly. camonly
# fits the SAME simulated replicate with inat_on = 0 (iNat likelihood removed).
#
# USAGE  Rscript 01l_run_camanch_regional.R <row_id>
#        Rscript 01l_run_camanch_regional.R design     (writes camanch_design.csv)
# Output: <CAMANCH_OUTDIR>/row_%05d.rds = list(rep = one-row data.frame,
#         region = one row per ecoregion present in the range)
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
BASE <- Sys.getenv("SIM_BASE", ".")
DESIGN_SEED <- 20261017L          # same as Tier 1: identical site / cell50 draws per design cell
RUN_SEED    <- 20261123L
N_REP       <- as.integer(Sys.getenv("N_REP", "20"))
N_BURNIN    <- as.integer(Sys.getenv("N_BURNIN", "2000"))
N_ITER      <- as.integer(Sys.getenv("N_ITER", "8000"))
SMOKE       <- Sys.getenv("CAMANCH_SMOKE", "0") == "1"

SCENARIOS   <- c("NULL", "REAL", "SPURIOUS", "REAL_SPURIOUS")
CAM_LEVELS  <- c(1000L, 3000L, 9000L)
INAT_LEVELS <- c(5000, 20000, 100000)
RANGE_CELLS  <- as.integer(Sys.getenv("RANGE_CELLS", "900"))
WINDOW_YEARS <- as.integer(Sys.getenv("WINDOW_YEARS", "5"))
INAT_YEARS   <- 2008:2025
WINDOW_FRAC  <- WINDOW_YEARS / length(INAT_YEARS)   # camera_window_years() then returns exactly the last WINDOW_YEARS
REGION_SD   <- 0.2

#' @name build_camanch_design
#' @description One row per scenario x camera level x iNat level x replicate.
#'   Sorted by camera level first so array tasks group by runtime class.
#' @return data.frame with row_id and design factors.
build_camanch_design <- function() {
  d <- expand.grid(rep_id = seq_len(N_REP), inat_records = INAT_LEVELS,
                   scenario = SCENARIOS, n_site = CAM_LEVELS, stringsAsFactors = FALSE)
  d <- d[order(d$n_site, match(d$scenario, SCENARIOS), d$inat_records, d$rep_id), ]
  d$arm <- "joint"; d$row_id <- seq_len(nrow(d)); rownames(d) <- NULL
  d[, c("row_id", "scenario", "n_site", "inat_records", "rep_id", "arm")]
}

if (length(args) && args[1] == "design") {
  d <- build_camanch_design(); write.csv(d, "camanch_design.csv", row.names = FALSE)
  cat(nrow(d), "rows -> camanch_design.csv\n"); quit(save = "no")
}

suppressPackageStartupMessages({ library(nimble); library(nimbleEcology) })
PROJ_ROOT <- normalizePath(file.path(BASE, ".."), mustWork = FALSE)
.ih <- c(file.path(PROJ_ROOT, "HPC", "bobcat", "integration_helper.R"),
         file.path(BASE, "HPC", "bobcat", "integration_helper.R"))
.ih <- .ih[file.exists(.ih)]
if (!length(.ih)) stop("integration_helper.R not found (calcIntensity_SVC would be undefined)")
source(.ih[1])
source(file.path(BASE, "sim_helpers.R"))
source(file.path(BASE, "sim_helpers_estimator_metrics.R"))   # .split_rhat
source(file.path(BASE, "sim_helpers_tier1.R"))
source(file.path(BASE, "model_code_ecoregion_camanchored.R"))

row_id <- as.integer(args[1])
if (SMOKE) {
  sm <- read.table(file.path(BASE, "camanch_smoke_rows.txt"), stringsAsFactors = FALSE,
                   col.names = c("scenario", "n_site", "inat_records", "arm"))
  if (is.na(row_id) || row_id > nrow(sm)) stop("smoke row out of range 1..", nrow(sm))
  cfg <- cbind(row_id = row_id, sm[row_id, ], rep_id = 1L)
} else {
  design <- build_camanch_design()
  if (is.na(row_id) || row_id > nrow(design)) stop("row_id out of range 1..", nrow(design))
  cfg <- design[row_id, ]
}
stopifnot(cfg$scenario %in% SCENARIOS, cfg$arm %in% c("joint", "camonly"))
OUTDIR <- Sys.getenv("CAMANCH_OUTDIR", file.path(BASE, if (SMOKE) "camanch_smoke" else "camanch_out"))
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
outfile <- file.path(OUTDIR, sprintf("row_%05d.rds", row_id))
cat(sprintf("[row %d] scenario=%s n_site=%d inat_records=%g rep=%d arm=%s range=%d iter=%d+%d\nOUTDIR %s\n",
            row_id, cfg$scenario, cfg$n_site, cfg$inat_records, cfg$rep_id, cfg$arm, RANGE_CELLS,
            N_BURNIN, N_ITER, OUTDIR))
if (file.exists(outfile)) { cat("already done, exiting\n"); quit(save = "no", status = 0) }

#' @name cut_inat_window
#' @description Cut the iNat side of a Tier 1 design to the camera window so both
#'   streams span the same calendar years, re-z-scoring year_vals on the kept
#'   years as the windowed prep (run_fleet_prep.R, year_vals = scale(years_in)) does.
#' @param des Output of build_tier1_constants().
#' @param keep_years Calendar years to keep (the camera window).
#' @return des with constants, inat_effort_base and s_vals updated.
cut_inat_window <- function(des, keep_years) {
  tk <- match(keep_years, INAT_YEARS); stopifnot(!anyNA(tk))
  k <- des$constants; n50 <- k$ncell50; stopifnot(k$nyear == length(INAT_YEARS))
  k$nyear <- length(tk)
  k$year_vals <- as.numeric(scale(keep_years))
  k$xdat_inat <- k$xdat_inat[, , tk, drop = FALSE]
  k$MWMT_inat <- k$MWMT_inat[, tk, drop = FALSE]
  k$MCMT_inat <- k$MCMT_inat[, tk, drop = FALSE]
  k$n_cells_year <- rep(n50, length(tk))
  k$inat_cells_by_year <- matrix(rep(seq_len(n50), length(tk)), nrow = n50, ncol = length(tk))
  des$constants <- k
  des$inat_effort_base <- des$inat_effort_base[, tk, drop = FALSE]
  des$s_vals <- 1 / stats::sd(keep_years)
  des$diag$d_s_vals <- des$s_vals; des$diag$d_inat_year_min <- min(keep_years); des$diag$d_inat_year_max <- max(keep_years)
  des
}

#' @name ess1
#' @description Effective sample size of one chain (coda), NA if constant.
ess1 <- function(x) { x <- x[is.finite(x)]; if (length(x) < 10 || stats::var(x) == 0) NA_real_ else unname(coda::effectiveSize(x)) }

#' @name score_param
#' @description Posterior summary of one draw vector against its truth.
score_param <- function(x, truth) {
  q <- stats::quantile(x, c(0.025, 0.975), names = FALSE)
  data.frame(true = truth, mean = mean(x), sd = stats::sd(x), q025 = q[1], q975 = q[2],
             excl0 = q[1] > 0 | q[2] < 0, sign_ok = sign(mean(x)) == sign(truth) & truth != 0,
             covered = truth >= q[1] & truth <= q[2], rhat = .split_rhat(x), ess = ess1(x))
}

run_camanch_replicate <- function(cfg) {
  inputs <- readRDS(file.path(BASE, "prepped_sim_inputs_with_array.RDS"))
  cl <- inputs$constants_list
  er <- as.integer(inputs$ecoregion_of_cell100); nregion <- as.integer(inputs$nregion)
  stopifnot(length(er) == cl$ncell100, !anyNA(er), all(er %in% seq_len(nregion)))
  seed_cell <- choose_range_seed(cl, 25L, 900L)            # the Tier 1 seed, so the range is Tier 1's

  des <- build_tier1_constants(inputs, RANGE_CELLS, WINDOW_FRAC, cfg$n_site, seed_cell,
                               design_seed = DESIGN_SEED)
  # Spec: distinct real camera sites, never resampled with replacement.
  if (des$diag$d_nsite_unique != cfg$n_site)
    stop(sprintf("site pool %d < %d: would resample with replacement", des$diag$d_site_pool, cfg$n_site))
  cam_years <- sort(unique(des$site_cal_year))
  stopifnot(identical(cam_years, tail(INAT_YEARS, WINDOW_YEARS)))
  des <- cut_inat_window(des, cam_years)                    # both streams: the same WINDOW_YEARS
  constants <- des$constants
  constants$ecoregion_of_cell100 <- er
  constants$nregion <- nregion
  constants$region_of_site <- er[constants$cell]
  constants$inat_on <- 1L

  # Scored regions: real ecoregions with data in the range. "WATER" is a mapping
  # artefact of the Level I layer (26 cells), not a region: truth 0, never scored.
  regions_in <- sort(unique(c(er[des$range_keep])))
  regions_in <- regions_in[toupper(as.character(inputs$ecoregion_levels[regions_in])) != "WATER"]
  if (RANGE_CELLS == 900L) stopifnot(length(regions_in) == 7L)
  truth0 <- true_param_list(inputs$real_post_means, rep(0, cl$ncell100))
  truth0$year_effect <- NULL; truth0$tau_year <- NULL
  truth0 <- tier1_trend_truth(truth0, 0, des$s_occ, des$s_vals)    # national null
  cal <- calibrate_inat_effort(des, truth0, cfg$inat_records)        # regional terms 0 at the null
  des$inat_effort <- cal$effort

  # Axis-live checks (the abundance ladder once recorded a label and applied nothing)
  stopifnot(constants$nsite == cfg$n_site, length(unique(des$site_keep)) == cfg$n_site,
            length(des$range_keep) == RANGE_CELLS,
            abs(cal$expected_records / cfg$inat_records - 1) < 1e-8,
            length(regions_in) >= 2)

  # RESEED -- must stay after build_tier1_constants()
  set.seed(RUN_SEED + row_id)

  yr <- rep(0, nregion); yrv <- rep(0, nregion)
  if (cfg$scenario %in% c("REAL", "REAL_SPURIOUS"))     yr[regions_in]  <- stats::rnorm(length(regions_in), 0, REGION_SD)
  if (cfg$scenario %in% c("SPURIOUS", "REAL_SPURIOUS")) yrv[regions_in] <- stats::rnorm(length(regions_in), 0, REGION_SD)
  truth <- truth0
  truth$year_region <- yr; truth$year_region_var <- yrv
  truth$sigma_region <- max(stats::sd(yr[regions_in]), 0.05)
  truth$sigma_region_var <- max(stats::sd(yrv[regions_in]), 0.05)
  if (cfg$scenario == "NULL") stopifnot(all(yr == 0), all(yrv == 0))
  if (cfg$scenario == "REAL") stopifnot(any(yr != 0), all(yrv == 0))
  if (cfg$scenario == "SPURIOUS") stopifnot(all(yr == 0), any(yrv != 0))

  sim <- simulate_replicate_data(model_code_ecoregion_camanchored, constants, des$inat_effort,
                                 y_ncol = des$y_ncol, truth)

  # realized data per region: camera detecting sites and iNat records
  det_site <- rowSums(sim$y, na.rm = TRUE) > 0
  inat_rec_g <- rowSums(sim$y_inat, na.rm = TRUE)
  reg_cells50 <- er[constants$inat_cell100]
  realized_region <- data.frame(
    region = seq_len(nregion),
    r_cam_sites = tabulate(constants$region_of_site, nregion),
    r_cam_det_sites = tabulate(constants$region_of_site[det_site], nregion),
    r_inat_cell50 = tabulate(reg_cells50, nregion),
    r_inat_records = vapply(seq_len(nregion), function(r) sum(inat_rec_g[reg_cells50 == r]), numeric(1)))

  fit_constants <- constants
  if (cfg$arm == "camonly") fit_constants$inat_on <- 0L
  t0 <- Sys.time()
  samples <- fit_replicate(model_code_ecoregion_camanchored, fit_constants, des$inat_effort,
                           list(y = sim$y, y_inat = sim$y_inat), base_inits = inputs$base_inits,
                           n_burnin = N_BURNIN, n_iter = N_ITER,
                           trend_inits = list(year_beta = 0, year_var = 0, sigma_year_beta = 1, sigma_year_var = 1,
                                              year_region = rep(0, nregion), year_region_var = rep(0, nregion),
                                              sigma_region = 1, sigma_region_var = 1),
                           extra_monitors = c("year_region", "year_region_var", "sigma_region", "sigma_region_var"))
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  yb <- samples[, "year_beta"]
  reg <- do.call(rbind, lapply(regions_in, function(r) {
    rep_draw <- yb + samples[, sprintf("year_region[%d]", r)]
    a <- score_param(rep_draw, truth$year_beta + yr[r]); names(a) <- paste0("rep_", names(a))
    b <- score_param(samples[, sprintf("year_region_var[%d]", r)], yrv[r]); names(b) <- paste0("yrv_", names(b))
    cbind(region = r, region_label = as.character(inputs$ecoregion_levels[r]),
          range_cells = sum(er[des$range_keep] == r),
          rep_pct_true = 100 * (exp((truth$year_beta + yr[r]) * des$s_occ) - 1),
          rep_pct_mean = mean(100 * (exp(rep_draw * des$s_occ) - 1)),
          # implied %/yr per stream (Roland 2026-11-23: make the residual stream gap visible)
          cam_pct_true  = 100 * (exp((truth$year_beta + yr[r]) * des$s_occ) - 1),
          inat_pct_true = 100 * (exp((truth$year_beta + truth$year_var + yr[r] + yrv[r]) * des$s_vals) - 1),
          cam_pct_mean  = mean(100 * (exp(rep_draw * des$s_occ) - 1)),
          inat_pct_mean = if (cfg$arm == "joint") mean(100 * (exp((rep_draw + samples[, "year_var"] +
                            samples[, sprintf("year_region_var[%d]", r)]) * des$s_vals) - 1)) else NA_real_,
          a, b, realized_region[r, -1], stringsAsFactors = FALSE)
  }))
  rep_row <- cbind(cfg, des$diag, range_cells = RANGE_CELLS, window_years = WINDOW_YEARS, n_regions_in = length(regions_in),
                   r_cam_det_sites = sum(det_site), r_inat_records = sum(sim$y_inat, na.rm = TRUE),
                   yb_mean = mean(yb), yb_rhat = .split_rhat(yb), yb_ess = ess1(yb),
                   sigma_region_mean = mean(samples[, "sigma_region"]),
                   sigma_region_var_mean = mean(samples[, "sigma_region_var"]),
                   s_occ = des$s_occ, s_vals = des$s_vals, elapsed_sec = elapsed,
                   n_burnin = N_BURNIN, n_iter = N_ITER, status = "OK", stringsAsFactors = FALSE)
  cat(sprintf("[design] nsite=%d (unique %d, pool %d) cam cells=%d ncell50=%d regions in range=%s\n",
              constants$nsite, des$diag$d_nsite_unique, des$diag$d_site_pool, des$diag$d_cam_cells,
              constants$ncell50, paste(regions_in, collapse = ",")))
  cat(sprintf("[data] detecting sites=%d iNat records=%d  fit %.1f min\n", sum(det_site),
              sum(sim$y_inat, na.rm = TRUE), elapsed / 60))
  print(reg[, c("region", "range_cells", "r_cam_sites", "r_cam_det_sites", "r_inat_records",
                "rep_true", "rep_mean", "rep_q025", "rep_q975", "rep_rhat", "yrv_true", "yrv_mean")], digits = 3)
  list(rep = rep_row, region = cbind(row_id = cfg$row_id, scenario = cfg$scenario, arm = cfg$arm,
                                     n_site = cfg$n_site, inat_records = cfg$inat_records,
                                     rep_id = cfg$rep_id, reg, stringsAsFactors = FALSE))
}

res <- tryCatch(run_camanch_replicate(cfg), error = function(e) {
  cat("ERRORED:", conditionMessage(e), "\n")
  list(rep = cbind(cfg, status = "Errored", error_msg = conditionMessage(e)), region = NULL)
})
saveRDS(res, outfile)
cat("wrote", outfile, "\n")
