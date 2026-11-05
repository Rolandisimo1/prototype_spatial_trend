# =============================================================================
# sim_helpers_tier1.R
# -----------------------------------------------------------------------------
# Design builder for Tier 1 of SEED_data_volume_simulation.md: how much camera
# and iNaturalist data a trend needs, over how large a range, with how long a
# camera record. Source AFTER sim_helpers.R.
#
# Varies three things the estimator sweep held fixed. The camera estimator does
# not change (production camera_occ, model_code_national_scalar, untouched):
#
#   RANGE SIZE  -- number of cell100s that hold data. A range of k cells is
#     the first k cells of a breadth-first walk over the real CAR adjacency
#     from one fixed seed cell, so ranges are compact, connected and NESTED
#     (25 inside 50 inside 100 ...). The CAR graph itself stays the full
#     908-node grid, exactly as in every production bundle (all 34 v2u
#     bundles have ncell100 = 908; the range mask only decides which cells
#     carry data). So a small range means a CAR field with data in few
#     nodes -- the mechanism behind the real range-size result. Cameras come
#     from real sites inside the range; iNat uses one cell50 per range cell100
#     (where one exists). Total camera count and total expected iNat records
#     are held fixed across range sizes, so range is not confounded with
#     data volume.
#   CAMERA WINDOW -- calendar years the cameras cover, always ENDING at the
#     last iNaturalist year (real structure: cameras 2019-2025 against iNat
#     2008-2025). Each camera site is given a calendar year from the window,
#     balanced across years, and year_occ is z-scored on those site-years
#     exactly as the production prep does it.
#   DATA VOLUME -- number of camera sites (resampled with replacement from
#     the real sites inside the range when the target exceeds the pool; y is
#     simulated fresh, so a repeated covariate row is two independent cameras
#     at the same habitat) and a target number of EXPECTED iNat records. The
#     real effort matrix keeps its shape (cells, growth over years) and is
#     rescaled so that sum(effort * mu) at the null truth equals the target
#     (calibrate_inat_effort()); the same effort is used at every effect size.
#
# EFFECT SIZE IS SET IN %/yr, NOT AS A COEFFICIENT. year_beta multiplies
# year_occ, whose step per calendar year depends on the camera window
# (1/sd(site years)); total_var_beta multiplies year_vals (fixed 18-year
# z-score, step ~0.1873). A population changing at r per year therefore has
#     year_beta      = log(1 + r) / s_occ
#     total_var_beta = log(1 + r) / s_vals
#     year_var       = total_var_beta - year_beta
# i.e. both streams see the same population change, and the true sink is only
# the units correction between the two year covariates. This is the same
# conversion trend_table_fleet.R uses on real fits (camera %/yr =
# 100*(exp(year_beta*s_occ)-1)), so simulated and real %/yr are comparable.
# =============================================================================

#' @name car_neighbors
#' @description Neighbor list from nimble dcar_normal adj/num vectors.
#' @param adj,num CAR adjacency vectors.
#' @return List, one integer vector of neighbors per node.
car_neighbors <- function(adj, num) {
  st <- cumsum(c(1, num))[seq_along(num)]; en <- cumsum(num)
  lapply(seq_along(num), function(i) if (num[i] > 0) adj[st[i]:en[i]] else integer(0))
}

#' @name bfs_order
#' @description Breadth-first order of nodes reachable from seed. Ties inside a
#'   ring are broken by node index, so the order is deterministic.
#' @param nbrs Neighbor list from car_neighbors().
#' @param seed_cell Starting node.
#' @return Integer vector of nodes in visit order.
bfs_order <- function(nbrs, seed_cell) {
  seen <- seed_cell; front <- seed_cell
  while (length(front)) {
    nx <- sort(setdiff(unique(unlist(nbrs[front])), seen))
    seen <- c(seen, nx); front <- nx
  }
  seen
}

#' @name choose_range_seed
#' @description The seed cell for the nested ranges: among cells whose
#'   connected component can reach the largest range requested, the one whose
#'   smallest range holds the most real camera sites. Chosen this way so the
#'   25-cell range is not starved of real covariate rows to resample.
#' @param cl Real constants_list.
#' @param k_small,k_large Smallest and largest range sizes in the design.
#' @return Integer seed cell.
choose_range_seed <- function(cl, k_small, k_large) {
  nbrs <- car_neighbors(cl$adj, cl$num)
  sites_per <- tabulate(cl$cell, cl$ncell100)
  sc <- vapply(seq_len(cl$ncell100), function(s) {
    b <- bfs_order(nbrs, s)
    if (length(b) < k_large) NA_real_ else sum(sites_per[b[seq_len(k_small)]])
  }, numeric(1))
  which.max(sc)
}

#' @name camera_window_years
#' @description Calendar years a camera window of the given fraction covers,
#'   ending at the last iNat year.
#' @param frac Fraction of the iNat window (1, 0.6, 0.4, 0.2).
#' @param inat_years Calendar years of the iNat record (default 2008:2025).
#' @return Integer vector of calendar years.
camera_window_years <- function(frac, inat_years = 2008:2025) {
  n <- max(2L, as.integer(round(frac * length(inat_years))))
  tail(inat_years, n)
}

#' @name build_tier1_constants
#' @description Build one Tier 1 design: range restriction, camera window,
#'   camera count and iNat cell choice, as a nimble constants list for
#'   model_code_national_scalar plus the UNSCALED effort matrix (scale it with
#'   calibrate_inat_effort()).
#' @param inputs prepped_sim_inputs_with_array.RDS contents.
#' @param range_cells Number of cell100s in the range.
#' @param window_frac Camera window as a fraction of the iNat window.
#' @param n_site Number of camera sites.
#' @param seed_cell Range seed cell (from choose_range_seed()).
#' @param max_subcell_per_cell Subcells kept per cell50 (15, as in
#'   build_reduced_constants()).
#' @param design_seed RNG seed for the design (site and cell50 draws, year
#'   assignment), shared by every replicate of a design cell; data are
#'   simulated afterwards under a per-row seed.
#' @return List(constants, inat_effort_base, y_ncol, range_keep, site_keep,
#'   cell50_keep, site_cal_year, s_occ, s_vals, diag).
build_tier1_constants <- function(inputs, range_cells, window_frac, n_site,
                                  seed_cell, max_subcell_per_cell = 15,
                                  design_seed = 20261017L) {
  cl <- inputs$constants_list
  set.seed(design_seed)
  ord <- bfs_order(car_neighbors(cl$adj, cl$num), seed_cell)
  if (length(ord) < range_cells) stop("seed component smaller than range_cells")
  range_keep <- ord[seq_len(range_cells)]

  # ---- iNat side: one cell50 per range cell100 ---------------------------
  cell50_keep <- sort(unlist(lapply(range_keep, function(c100) {
    pool <- which(cl$inat_cell100 == c100)
    if (length(pool) > 1) sample(pool, 1) else pool
  })))
  n50 <- length(cell50_keep)
  if (n50 < 3) stop("fewer than 3 cell50s in range")
  new_start <- integer(n50); new_end <- integer(n50); keep_rows <- integer(0); cursor <- 1L
  for (i in seq_len(n50)) {
    g <- cell50_keep[i]
    rr <- cl$inat_cell50_start[g]:cl$inat_cell50_end[g]
    take <- rr[seq_len(min(length(rr), max_subcell_per_cell))]
    keep_rows <- c(keep_rows, take)
    new_start[i] <- cursor; new_end[i] <- cursor + length(take) - 1L; cursor <- new_end[i] + 1L
  }

  # ---- camera side: real sites inside the range --------------------------
  pool_site <- which(cl$cell %in% range_keep)
  if (length(pool_site) < 10) stop("fewer than 10 real camera sites in range")
  site_keep <- if (n_site <= length(pool_site)) sort(sample(pool_site, n_site))
               else sort(c(pool_site, sample(pool_site, n_site - length(pool_site), replace = TRUE)))
  # NIMBLE sizes y[i, 1:J[i]] ~ dOcc_v from the LAST site; J[nsite] == 1
  # collapses it to a scalar and nimbleModel() fails (PIPELINE.md section 5).
  # Sites are exchangeable, so move the last J >= 2 site to the end.
  Jk <- cl$J[site_keep]
  if (Jk[length(Jk)] < 2) {
    k <- max(which(Jk >= 2)); site_keep <- c(site_keep[-k], site_keep[k])
  }

  # ---- camera window: balanced calendar years, z-scored as in prep --------
  win <- camera_window_years(window_frac)
  site_cal_year <- sample(rep(win, length.out = n_site))
  year_occ <- as.numeric(scale(site_cal_year))
  s_occ  <- 1 / stats::sd(site_cal_year)
  s_vals <- stats::median(diff(sort(unique(cl$year_vals))))

  constants <- list(
    has_SVC = cl$has_SVC, hasSVC = cl$hasSVC, prior_type = cl$prior_type,
    nyear = cl$nyear, year_vals = cl$year_vals,
    nsite = length(site_keep), J = cl$J[site_keep],
    numOccCovars = cl$numOccCovars,
    occ_covars = cl$occ_covars[site_keep, , drop = FALSE],
    year_occ = year_occ,
    yday = cl$yday[site_keep, , drop = FALSE],
    canopy_height = cl$canopy_height[site_keep],
    log_roaddist = cl$log_roaddist[site_keep],
    MWMT = cl$MWMT[site_keep], MCMT = cl$MCMT[site_keep],
    interaction_group = cl$interaction_group,
    adj = cl$adj, num = cl$num, nadj = cl$nadj, nnum = cl$nnum, ncell100 = cl$ncell100,
    cell = cl$cell[site_keep],
    ncell50 = n50,
    inat_cell100 = cl$inat_cell100[cell50_keep],
    xdat_inat = cl$xdat_inat[keep_rows, , , drop = FALSE],
    inat_cell50_start = new_start, inat_cell50_end = new_end,
    MCMT_inat = cl$MCMT_inat[keep_rows, , drop = FALSE],
    MWMT_inat = cl$MWMT_inat[keep_rows, , drop = FALSE],
    n_cells_year = rep(n50, cl$nyear),
    inat_cells_by_year = matrix(rep(seq_len(n50), cl$nyear), nrow = n50, ncol = cl$nyear)
  )
  stopifnot(constants$J[constants$nsite] >= 2, length(year_occ) == constants$nsite,
            all(constants$cell %in% range_keep), all(constants$inat_cell100 %in% range_keep))

  list(constants = constants, inat_effort_base = inputs$inat_effort_real[cell50_keep, , drop = FALSE],
       y_ncol = ncol(inputs$real_y_template),
       range_keep = range_keep, site_keep = site_keep, cell50_keep = cell50_keep,
       site_cal_year = site_cal_year, s_occ = s_occ, s_vals = s_vals,
       diag = data.frame(
         d_range_cells = range_cells,
         d_inat_cells = length(unique(constants$inat_cell100)), d_cam_cells = length(unique(constants$cell)),
         d_ncell50 = n50, d_nsubcell = length(keep_rows),
         d_nsite = length(site_keep), d_nsite_unique = length(unique(site_keep)),
         d_site_pool = length(pool_site),
         d_cam_year_min = min(site_cal_year), d_cam_year_max = max(site_cal_year),
         d_cam_nyears = length(unique(site_cal_year)),
         d_s_occ = s_occ, d_s_vals = s_vals, stringsAsFactors = FALSE))
}

#' @name inat_mu_r
#' @description iNat intensity mu[g, t] computed in plain R, identical in
#'   formula to calcIntensity_SVC() (integration_helper.R): log_lambda per
#'   subcell, summed on the natural scale, then exp(theta0 + theta1 * log(.)).
#'   Used only to calibrate effort; the data are still simulated by NIMBLE.
#' @param constants Tier 1 constants.
#' @param truth Truth list.
#' @return ncell50 x nyear matrix.
inat_mu_r <- function(constants, truth) {
  tvb <- truth$year_beta + truth$year_var
  mu <- matrix(NA_real_, constants$ncell50, constants$nyear)
  for (g in seq_len(constants$ncell50)) {
    rr <- constants$inat_cell50_start[g]:constants$inat_cell50_end[g]
    c100 <- constants$inat_cell100[g]
    for (t in seq_len(constants$nyear)) {
      ll <- truth$link_occ_intercept[c100] +
        as.vector(matrix(constants$xdat_inat[rr, , t], nrow = length(rr)) %*% truth$occ_beta) +
        constants$MWMT_inat[rr, t] * truth$MWMT_effect[c100] +
        constants$MCMT_inat[rr, t] * truth$MCMT_effect[c100] + tvb * constants$year_vals[t]
      mu[g, t] <- exp(truth$theta0 + truth$theta1 * log(sum(exp(ll))))
    }
  }
  mu
}

#' @name calibrate_inat_effort
#' @description Rescale the real effort matrix so the EXPECTED iNat record
#'   total, sum(effort * mu), equals target_records at the NULL truth (trend
#'   zeroed). The expected value of dnbinom(size = 1/od, prob = 1/(1 + od *
#'   effort * mu)) is effort * mu, so this sets the iNat volume in records,
#'   the unit the fleet is described in. Calibrated once per design cell at
#'   the null so every effect size gets the same effort.
#' @param des Output of build_tier1_constants().
#' @param truth0 Truth list with year_beta = year_var = 0.
#' @param target_records Expected record total.
#' @return List(effort, mult, expected_records).
calibrate_inat_effort <- function(des, truth0, target_records) {
  stopifnot(truth0$year_beta == 0, truth0$year_var == 0)
  mu0 <- inat_mu_r(des$constants, truth0)
  e0 <- sum(des$inat_effort_base * mu0)
  if (!is.finite(e0) || e0 <= 0) stop("expected iNat records at base effort is not positive")
  mult <- target_records / e0
  eff <- des$inat_effort_base * mult
  list(effort = eff, mult = mult, expected_records = sum(eff * mu0))
}

#' @name tier1_trend_truth
#' @description Set the trend truth from a %/yr effect (see file header).
#' @param truth Truth list.
#' @param r Annual proportional change (e.g. -0.05); 0 is the null.
#' @param s_occ,s_vals Per-calendar-year steps of year_occ and year_vals.
#' @return truth with year_beta and year_var set.
tier1_trend_truth <- function(truth, r, s_occ, s_vals) {
  if (r == 0) { truth$year_beta <- 0; truth$year_var <- 0; return(truth) }
  truth$year_beta <- log(1 + r) / s_occ
  truth$year_var  <- log(1 + r) / s_vals - truth$year_beta
  truth
}
