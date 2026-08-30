#!/usr/bin/env Rscript
# =============================================================================
# 00z_make_synthetic_inputs.R
# -----------------------------------------------------------------------------
# Build a SYNTHETIC sim_inputs.RDS so the simulation code base runs end-to-end
# on any machine, with no access to Hazel and no real data.
#
# WHY THIS EXISTS
# The production inputs (`sim_inputs.RDS`) are built by 00_prep_sim_inputs.R
# from the real camera/iNaturalist bundle and a fitted bobcat posterior, both of
# which live on Hazel behind institutional access. Without them the pipeline
# cannot be executed by an independent reviewer at all. This script writes an
# object with the SAME STRUCTURE and plausible magnitudes, so every downstream
# script (01e/01f/01g/01i and their collectors) runs unmodified.
#
# WHAT THIS IS NOT
# The numbers produced from synthetic inputs are NOT our reported results and
# will not reproduce them. Real inputs carry the true spatial adjacency, the
# real camera effort distribution, and posterior-mean parameter magnitudes from
# an actual fit. This script exists to verify that the CODE does what the
# methods say -- that the estimators are implemented as described, that the
# metrics compute correctly, and that the design grid is what we claim. To
# reproduce the reported numbers you need the real bundle; see README.
#
# Sized to PRODUCTION SCALE (ncell100 = 600, nsite = 1200). This is deliberate
# and was learned the hard way: 01i_run_estimator_sweep.R hardcodes
# n_site_target = 700 in its sample_sites_by_array() call, and that function
# fails with "missing value where TRUE/FALSE needed" whenever a region's array
# pool cannot reach its share of the target. A small toy grid therefore cannot
# execute the array arms at all. Runtime is dominated by MCMC iterations, so
# use small N_BURNIN/N_ITER for a quick structural check rather than shrinking
# the grid.
#
# Usage:  Rscript 00z_make_synthetic_inputs.R [outfile]
# =============================================================================

set.seed(20260827)

OUTFILE <- if (length(commandArgs(TRUE)) > 0) commandArgs(TRUE)[1] else "sim_inputs.RDS"

# ---- grid geometry ----------------------------------------------------------
# A rectangular lattice standing in for the real 100 km CAR grid. Rook adjacency
# gives a valid dcar_normal structure (every cell has >= 2 neighbours).
NX <- 30L; NY <- 20L
ncell100 <- NX * NY
gx <- rep(seq_len(NX), times = NY)
gy <- rep(seq_len(NY), each  = NX)

adj <- integer(0); num <- integer(ncell100)
for (k in seq_len(ncell100)) {
  nb <- which((abs(gx - gx[k]) == 1 & gy == gy[k]) |
              (abs(gy - gy[k]) == 1 & gx == gx[k]))
  adj <- c(adj, nb); num[k] <- length(nb)
}
weights <- rep(1, length(adj))

# region labels: the NE/West/other split make_true_year_effect() expects
region <- ifelse(gx <= 3, "West", ifelse(gx >= 8 & gy >= 4, "NE", "other"))
cell100_geo <- data.frame(
  cell100 = seq_len(ncell100), lon = gx, lat = gy, region = region,
  stringsAsFactors = FALSE
)

# ---- 50 km cells, and the subcell ROWS nested inside each -------------------
# Two distinct nestings, easy to conflate and fatal if conflated:
#   1. each 100 km CAR cell contains SUBPER 50 km cells  -> inat_cell100
#   2. each 50 km cell owns a contiguous BLOCK OF ROWS in xdat_inat, indexed by
#      inat_cell50_start[g]:inat_cell50_end[g]
# The start/end pair indexes ROWS OF xdat_inat, NOT 50 km cell ids. Making them
# index cells (start = seq(1, ncell50, by = SUBPER)) produces indices that run
# past the end of the row space for high g and yields NA, which surfaces far
# downstream as "NA/NaN argument" inside build_reduced_constants(). Caught by
# running a real replicate against these inputs.
SUBPER   <- 4L                      # 50 km cells per 100 km cell
ROWS_PER <- 4L                      # xdat_inat rows per 50 km cell
ncell50 <- ncell100 * SUBPER
inat_cell100 <- rep(seq_len(ncell100), each = SUBPER)
subcell_n <- rep(SUBPER, ncell100)
n_inat_rows <- ncell50 * ROWS_PER
inat_cell50_start <- seq(1L, n_inat_rows, by = ROWS_PER)
inat_cell50_end   <- inat_cell50_start + ROWS_PER - 1L
stopifnot(length(inat_cell50_start) == ncell50,
          max(inat_cell50_end) == n_inat_rows)

# ---- years and covariates ---------------------------------------------------
nyear <- 18L
year_vals <- scale(seq_len(nyear))[, 1]
numOccCovars <- 10L

nsite <- 1200L
cell <- sample.int(ncell100, nsite, replace = TRUE)
J <- pmax(1L, rpois(nsite, 3.2))                    # matches measured mean J = 3.22
maxJ <- max(J)
yday <- matrix(rnorm(nsite * maxJ), nsite, maxJ)
year_occ <- year_vals[sample.int(nyear, nsite, replace = TRUE)]

cl <- list(
  ncell100 = ncell100, ncell50 = ncell50, nsite = nsite, nyear = nyear,
  numOccCovars = numOccCovars, year_vals = year_vals,
  adj = adj, num = num, weights = weights,
  nadj = length(adj), nnum = length(num),
  cell = cell, J = J, yday = yday, year_occ = year_occ,
  canopy_height = rnorm(nsite), log_roaddist = rnorm(nsite),
  MWMT = rnorm(nsite), MCMT = rnorm(nsite),
  occ_covars = matrix(rnorm(nsite * numOccCovars), nsite, numOccCovars),
  inat_cell100 = inat_cell100,
  inat_cell50_start = inat_cell50_start, inat_cell50_end = inat_cell50_end,
  # xdat_inat is indexed by ROW (n_inat_rows), not by 50 km cell -- the
  # start/end pair above slices into this first dimension.
  xdat_inat  = array(rnorm(n_inat_rows * numOccCovars * nyear),
                     dim = c(n_inat_rows, numOccCovars, nyear)),
  MWMT_inat  = matrix(rnorm(n_inat_rows * nyear), n_inat_rows, nyear),
  MCMT_inat  = matrix(rnorm(n_inat_rows * nyear), n_inat_rows, nyear),
  # Three build-time flags the model code branches on inside nimbleCode().
  # These are resolved by nimbleModel() at graph-construction time, not at run
  # time, so a missing or mis-cased value fails inside codeProcessIfThenElse()
  # with a bare "argument is of length zero" and no line number.
  #   has_SVC / hasSVC : BOTH are read (lines 27, 80, 121 of the baseline);
  #                      they are separate names, not aliases.
  #   prior_type       : compared to the literal "Normal" -- CASE-SENSITIVE.
  #                      "normal" matches neither branch and yields a
  #                      zero-length condition.
  has_SVC = TRUE, hasSVC = TRUE, prior_type = "Normal",
  interaction_group = NULL
)

# every 50 km cell observed every year, so inat_cells_by_year is a full matrix
cl$n_cells_year <- rep(ncell50, nyear)
cl$inat_cells_by_year <- matrix(rep(seq_len(ncell50), nyear), ncell50, nyear)

# ---- effort and templates ---------------------------------------------------
# lognormal effort, the shape the real iNat effort surface has
inat_effort_real <- matrix(rlnorm(ncell50 * nyear, meanlog = 2, sdlog = 1),
                           ncell50, nyear)
y_template     <- matrix(NA_integer_, nsite, maxJ)
real_y_template <- y_template

# ---- parameter magnitudes ---------------------------------------------------
# Plausible stand-ins for the bobcat posterior means the real prep supplies.
# link_occ_intercept is given genuine spatial structure so the CAR field and the
# SVC terms have something to recover.
real_post_means <- list(
  occ_beta           = rnorm(numOccCovars, 0, 0.3),
  link_occ_intercept = as.numeric(scale(gx + gy)) * 0.4 - 2.0,
  MWMT_effect        = as.numeric(scale(gx)) * 0.2,
  MCMT_effect        = as.numeric(scale(gy)) * 0.2,
  intercept_tau      = 1.5,
  MWMT_tau           = 2.0,
  MCMT_tau           = 2.0,
  det_intercept      = 0.25,
  p_beta             = rnorm(4, 0, 0.2),
  theta0             = -1.0,
  theta1             = 1.0,
  overdisp_inat      = 0.5,
  year_beta          = -0.12,   # near the real moose estimate
  year_var           = -0.06
)

base_inits <- list(
  occ_beta = real_post_means$occ_beta,
  link_occ_intercept = real_post_means$link_occ_intercept,
  MWMT_effect = real_post_means$MWMT_effect,
  MCMT_effect = real_post_means$MCMT_effect,
  intercept_tau = real_post_means$intercept_tau,
  MWMT_tau = real_post_means$MWMT_tau, MCMT_tau = real_post_means$MCMT_tau,
  det_intercept = real_post_means$det_intercept,
  p_beta = real_post_means$p_beta, q_beta = rep(0, 5),
  theta0 = real_post_means$theta0, theta1 = real_post_means$theta1,
  overdisp_inat = real_post_means$overdisp_inat,
  sigma_year_beta = 0.5, sigma_year_var = 0.5
)

# ---- array structure --------------------------------------------------------
# site_array is the array-level arms' prerequisite. Real runs take it from the
# 5 km spatial split of the camera_trap_array / subproject_name field
# (00c_build_site_array.R); here sites in a 100 km cell are split into groups of
# roughly 6, reproducing the measured median of 4 cameras per array-year.
site_array <- ave(seq_len(nsite), cell, FUN = function(ix)
  paste0("arr", ceiling(seq_along(ix) / 6)))
site_array <- paste0("c", cell, "_", site_array)

out <- list(
  cell100_geo = cell100_geo, constants_list = cl, cl = cl,
  subcell_n = subcell_n, real_post_means = real_post_means,
  base_inits = base_inits,
  numOccCovars = numOccCovars, nyear = nyear, year_vals = year_vals,
  inat_effort = inat_effort_real, inat_effort_real = inat_effort_real,
  y_template = y_template, real_y_template = real_y_template,
  site_array = site_array,
  SYNTHETIC = TRUE,
  SYNTHETIC_NOTE = paste(
    "Generated by 00z_make_synthetic_inputs.R. Structure matches the real",
    "prep output; VALUES ARE NOT REAL and will not reproduce reported results.",
    "Use to verify code behaviour only.")
)

saveRDS(out, OUTFILE)
cat(sprintf("wrote %s  (SYNTHETIC: ncell100=%d, ncell50=%d, inat_rows=%d, nsite=%d, nyear=%d, arrays=%d)\n",
            OUTFILE, ncell100, ncell50, n_inat_rows, nsite, nyear, length(unique(site_array))))
