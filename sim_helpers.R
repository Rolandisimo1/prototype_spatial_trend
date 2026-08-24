#!/usr/bin/env Rscript
# =============================================================================
# sim_helpers.R
# Helper functions for the spatially-varying-trend simulation validation.
# Never run directly -- sourced by 01_run_sim_validation.R.
# =============================================================================

#' @name smooth_over_car_graph
#' @description Blend a piecewise value per CAR node with its neighbors'
#'   mean, a few passes, to turn a blocky region assignment into a smooth
#'   field over the real adjacency graph (used to build the TRUE year_effect
#'   field for simulation -- not part of the fitted model).
#' @param x Numeric vector, length ncell100, one value per CAR node.
#' @param adj,num Real CAR adjacency vectors (nimble dcar_normal convention:
#'   adj is the flattened neighbor list, num gives the neighbor count per node).
#' @param n_iter Number of smoothing passes.
#' @param w Blend weight toward the neighbor mean (0 = no smoothing, 1 = fully
#'   replace with neighbor mean each pass).
#' @return Smoothed numeric vector, same length as x.
smooth_over_car_graph <- function(x, adj, num, n_iter = 3, w = 0.5) {
  ncell <- length(x)
  nbr_start <- cumsum(c(1, num))[1:ncell]
  nbr_end   <- cumsum(num)
  out <- x
  for (iter in seq_len(n_iter)) {
    new_out <- out
    for (i in seq_len(ncell)) {
      if (num[i] == 0) next
      nbrs <- adj[nbr_start[i]:nbr_end[i]]
      new_out[i] <- (1 - w) * out[i] + w * mean(out[nbrs])
    }
    out <- new_out
  }
  out
}

#' @name make_true_year_effect
#' @description Build the KNOWN true spatially-varying trend deviation used
#'   to simulate data: positive in the Northeast, negative in the West,
#'   ~zero elsewhere, smoothed over the real CAR graph, then centered to
#'   exactly zero mean (matching the fitted model's zero_mean=1 constraint).
#' @param cell100_geo Data frame with columns cell100, region ("NE"/"West"/"other"),
#'   in cell100 order 1:ncell100 (from prepped_sim_inputs.RDS).
#' @param adj,num Real CAR adjacency vectors.
#' @param amplitude Magnitude of the raw regional deviation before smoothing/centering.
#' @return Numeric vector, length ncell100, zero-mean.
make_true_year_effect <- function(cell100_geo, adj, num, amplitude = 0.3) {
  raw <- ifelse(cell100_geo$region == "NE", amplitude,
         ifelse(cell100_geo$region == "West", -amplitude, 0))
  smoothed <- smooth_over_car_graph(raw, adj, num, n_iter = 3, w = 0.5)
  smoothed - mean(smoothed)
}

#' @name build_reduced_constants
#' @description Build a computationally tractable but structurally faithful
#'   constants_list for one simulation replicate design (same across all
#'   replicates -- the design is fixed, only the simulated data differ). The
#'   real 908-cell CAR adjacency graph (adj/num/nadj/nnum/ncell100) is kept
#'   EXACTLY as in the real bobcat fit -- it is cheap and it is literally the
#'   object year_effect is estimated on. What's reduced is DATA VOLUME:
#'   number of iNat grid cells (cell50), subcells per cell, and camera sites
#'   -- these are what drove the ~7.5 sec/iteration real-data cost.
#'   Cell50s and sites are sampled stratified by region so all three regions
#'   (NE/West/other) are represented; retained cells/sites keep their REAL
#'   covariate and effort values (only the response data are simulated).
#' @param cl Real constants_list (from input_data_bobcat.RDS).
#' @param inat_effort_real Real inat_effort matrix (ncell50 x nyear).
#' @param real_y_template Real real_data$y (nsite x J_max), used only for
#'   column count (occasions), not its values.
#' @param cell100_geo Region-labeled cell100 table (for stratified sampling).
#' @param n_cell50_keep,max_subcell_per_cell,n_site_keep Reduction targets.
#' @param seed RNG seed for the (fixed, shared-across-replicates) subsample.
#' @param stratify_by Column in cell100_geo to stratify on -- defaults to
#'   "region" (the original NE/West/other labeling used by the CAR-field
#'   prototype). Pass "ecoregion" to stratify on the EPA Level I ecoregion
#'   instead (see 00b_prep_ecoregion.R) so every region is represented in
#'   each replicate's subsample, generalized to however many levels
#'   stratify_by actually has (not hardcoded to 3).
#' @param ecoregion_of_cell100_full Optional, full length-ncell100 ecoregion
#'   lookup vector (from prepped_sim_inputs.RDS$ecoregion_of_cell100, added
#'   by 00b_prep_ecoregion.R) -- if supplied,
#'   attached to the returned constants as `ecoregion_of_cell100` for the
#'   ecoregion model (needs the FULL 908-length lookup, not a subsampled
#'   one, since it's indexed by inat_cell100[g] which ranges over all real
#'   CAR cells regardless of which cell50s are retained). NULL (default)
#'   omits it -- harmless for the CAR/national-scalar models, which don't
#'   reference it.
#' @param nregion Optional region count, attached alongside
#'   ecoregion_of_cell100_full when supplied.
#' @param site_replace FALSE (default, unchanged prior behavior): sites are
#'   subsampled WITHOUT replacement, capped at each stratum's real pool size
#'   (requesting more than the real count silently returns the whole pool).
#'   TRUE (effort sweep's "above real" camera levels, 01g_run_effort_sweep.R):
#'   sites are sampled WITH replacement per stratum, so n_site_keep CAN
#'   exceed the real total -- an explicit extrapolation (real covariate rows
#'   reused for the "extra" virtual sites; y is simulated fresh per
#'   replicate regardless, so a duplicated covariate row just means two
#'   sites share psi/p but draw independent y -- statistically fine, not a
#'   real-data duplication concern). Only affects site_keep, never
#'   cell50_keep (the iNat axis is controlled by effort scaling, not more
#'   cells -- see scale_inat_effort() in sim_helpers_effort.R).
#' @param site_jitter_sd Only used when site_replace=TRUE. Named list/vector
#'   of small SDs (as a FRACTION of each covariate's real population SD) to
#'   perturb canopy_height/log_roaddist/yday for DUPLICATE (2nd+ occurrence)
#'   resampled sites only -- the first (real) occurrence of every site is
#'   always left exactly as observed. 0 (default) = no jitter, pure
#'   duplication. Purely cosmetic/optional (see file header note in
#'   01g_run_effort_sweep.R) -- avoids literal duplicate covariate rows in
#'   the design without changing anything about model validity.
#' @param site_keep_override Optional pre-computed site index (integer,
#'   into 1:cl$nsite) to use INSTEAD of this function's own independent
#'   stratified random sampling -- e.g. from sample_sites_by_array() in
#'   sim_helpers_array.R, which selects whole array-year units so retained
#'   arrays keep a realistic camera-count distribution (see that function's
#'   docstring for why independent per-site sampling defeats the purpose of
#'   the array-estimator sweep). NULL (default) preserves exact prior
#'   behavior -- n_site_keep/site_replace/seed still control sampling as
#'   before. When supplied, n_site_keep and site_replace are ignored for
#'   site selection (seed still applies to cell50 sampling and any
#'   downstream jitter).
#' @return List with `constants` (nimble constants_list for the reduced
#'   design) and `inat_effort` (reduced-rows effort matrix, passed as data).
build_reduced_constants <- function(cl, inat_effort_real, real_y_template, cell100_geo,
                                    n_cell50_keep = 210, max_subcell_per_cell = 15,
                                    n_site_keep = 700, seed = 20260712,
                                    stratify_by = "region",
                                    ecoregion_of_cell100_full = NULL, nregion = NULL,
                                    site_replace = FALSE, site_jitter_sd = 0,
                                    site_keep_override = NULL) {
  set.seed(seed)

  region_vec <- cell100_geo[[stratify_by]]
  cell50_region <- region_vec[match(cl$inat_cell100, cell100_geo$cell100)]
  site_region   <- region_vec[match(cl$cell, cell100_geo$cell100)]
  # CRITICAL: fixed explicit order for the legacy "region" column, not
  # unique(region_vec)'s data-dependent first-occurrence order. set.seed()
  # + sequential sample() calls inside the region loop below means the
  # ITERATION ORDER over regions determines which random draws land in
  # which region, even for an identical seed -- silently reordering this
  # would change cell50_keep/site_keep for every EXISTING CAR-based script
  # (01_run_sim_validation.R, 01c/01d grain sweep) that relies on
  # stratify_by="region" + seed=20260712 reproducing the exact original
  # 210-cell/700-site subsample. New stratify_by values (e.g. "ecoregion")
  # have no such prior reproducibility commitment, so a data-driven order
  # is fine there.
  region_levels <- if (stratify_by == "region") c("NE", "West", "other") else unique(region_vec)

  stratified_sample <- function(idx, region, n_total, replace = FALSE) {
    n_per <- ceiling(n_total / length(region_levels))
    out <- unlist(lapply(region_levels, function(r) {
      pool <- idx[region == r]
      if (replace) sample(pool, size = n_per, replace = TRUE)
      else sample(pool, size = min(n_per, length(pool)))
    }))
    if (replace) sort(out) else sort(unique(out))
  }

  cell50_keep <- stratified_sample(seq_len(cl$ncell50), cell50_region, n_cell50_keep)
  # site_keep_override lets a caller supply a PRE-COMPUTED site index (e.g.
  # from sample_sites_by_array() in sim_helpers_array.R, which selects whole
  # array-year units rather than independent sites) instead of this
  # function's own independent-per-site random sampling. NULL (default)
  # preserves exact prior behavior for every existing caller.
  site_keep <- if (!is.null(site_keep_override)) {
    sort(unique(site_keep_override))
  } else {
    stratified_sample(seq_len(cl$nsite), site_region, n_site_keep,
                      replace = site_replace)
  }
  n_new_cell50 <- length(cell50_keep)
  n_new_site   <- length(site_keep)

  # ---- optional jitter for duplicate (resampled-with-replacement) sites --
  # only touches the 2nd+ occurrence of a repeated real site index; the
  # first occurrence (and everything when site_replace=FALSE, since
  # site_keep then has no duplicates) is untouched.
  canopy_height_v <- cl$canopy_height[site_keep]
  log_roaddist_v  <- cl$log_roaddist[site_keep]
  yday_m          <- cl$yday[site_keep, , drop = FALSE]
  if (site_replace && !identical(site_jitter_sd, 0)) {
    dup_mask <- duplicated(site_keep)
    if (any(dup_mask)) {
      jit <- if (is.list(site_jitter_sd)) site_jitter_sd else
        list(canopy_height = site_jitter_sd, log_roaddist = site_jitter_sd, yday = site_jitter_sd)
      n_dup <- sum(dup_mask)
      canopy_height_v[dup_mask] <- canopy_height_v[dup_mask] +
        rnorm(n_dup, 0, jit$canopy_height * sd(cl$canopy_height, na.rm = TRUE))
      log_roaddist_v[dup_mask] <- log_roaddist_v[dup_mask] +
        rnorm(n_dup, 0, jit$log_roaddist * sd(cl$log_roaddist, na.rm = TRUE))
      yday_sd <- sd(cl$yday, na.rm = TRUE)
      yday_m[dup_mask, ] <- yday_m[dup_mask, , drop = FALSE] +
        matrix(rnorm(n_dup * ncol(yday_m), 0, jit$yday * yday_sd), nrow = n_dup)
    }
  }

  # --- subcells: for each retained cell50, keep up to max_subcell_per_cell
  # rows from its real subcell range (first N rows -- deterministic, no extra
  # RNG draw needed since which subcells within a cell is arbitrary anyway).
  new_start <- integer(n_new_cell50); new_end <- integer(n_new_cell50)
  keep_rows <- integer(0)
  cursor <- 1L
  for (i in seq_len(n_new_cell50)) {
    g <- cell50_keep[i]
    real_rows <- cl$inat_cell50_start[g]:cl$inat_cell50_end[g]
    take <- real_rows[seq_len(min(length(real_rows), max_subcell_per_cell))]
    keep_rows <- c(keep_rows, take)
    new_start[i] <- cursor
    new_end[i]   <- cursor + length(take) - 1L
    cursor <- new_end[i] + 1L
  }

  xdat_inat_new <- cl$xdat_inat[keep_rows, , , drop = FALSE]
  MWMT_inat_new <- cl$MWMT_inat[keep_rows, , drop = FALSE]
  MCMT_inat_new <- cl$MCMT_inat[keep_rows, , drop = FALSE]

  # dense design: every retained cell50 observed every year (simplification
  # vs the real growing-coverage-over-time pattern; documented in README)
  n_cells_year_new <- rep(n_new_cell50, cl$nyear)
  inat_cells_by_year_new <- matrix(rep(seq_len(n_new_cell50), cl$nyear),
                                   nrow = n_new_cell50, ncol = cl$nyear)

  constants <- list(
    has_SVC = cl$has_SVC, hasSVC = cl$hasSVC,
    prior_type = cl$prior_type,
    nyear = cl$nyear, year_vals = cl$year_vals,
    nsite = n_new_site,
    J = cl$J[site_keep],
    numOccCovars = cl$numOccCovars,
    occ_covars = cl$occ_covars[site_keep, , drop = FALSE],
    year_occ = cl$year_occ[site_keep],
    yday = yday_m,
    canopy_height = canopy_height_v,
    log_roaddist = log_roaddist_v,
    MWMT = cl$MWMT[site_keep],
    MCMT = cl$MCMT[site_keep],
    interaction_group = cl$interaction_group,
    adj = cl$adj, num = cl$num, nadj = cl$nadj, nnum = cl$nnum, ncell100 = cl$ncell100,
    cell = cl$cell[site_keep],
    ncell50 = n_new_cell50,
    inat_cell100 = cl$inat_cell100[cell50_keep],
    xdat_inat = xdat_inat_new,
    inat_cell50_start = new_start,
    inat_cell50_end = new_end,
    MCMT_inat = MCMT_inat_new,
    MWMT_inat = MWMT_inat_new,
    n_cells_year = n_cells_year_new,
    inat_cells_by_year = inat_cells_by_year_new
  )

  if (!is.null(ecoregion_of_cell100_full)) {
    constants$ecoregion_of_cell100 <- ecoregion_of_cell100_full
    constants$nregion <- nregion
  }

  list(constants = constants,
       inat_effort = inat_effort_real[cell50_keep, , drop = FALSE],
       y_ncol = ncol(real_y_template),
       cell50_keep = cell50_keep, site_keep = site_keep)
}

#' @name true_param_list
#' @description Assemble the full set of TRUE parameter values used to
#'   simulate one replicate's data: every non-trend parameter fixed at its
#'   REAL posterior-mean value (from Arielle's actual converged fit -- keeps
#'   simulated magnitudes realistic), plus the controlled trend truth
#'   (year_beta/year_var at their real posterior means, "kept GLOBAL for
#'   identifiability" per the prototype's scope; year_effect the new known
#'   spatial field under test).
#' @param real_post_means List from prepped_sim_inputs.RDS.
#' @param year_effect_true Numeric vector, length ncell100.
#' @return Named list suitable as nimbleModel(inits = ...).
true_param_list <- function(real_post_means, year_effect_true) {
  rp <- real_post_means
  list(
    occ_beta = unname(rp$occ_beta),
    link_occ_intercept = unname(rp$link_occ_intercept),
    MWMT_effect = unname(rp$MWMT_effect),
    MCMT_effect = unname(rp$MCMT_effect),
    intercept_tau = rp$intercept_tau, MWMT_tau = rp$MWMT_tau, MCMT_tau = rp$MCMT_tau,
    det_intercept = rp$det_intercept,
    p_beta = unname(rp$p_beta),
    theta0 = rp$theta0, theta1 = rp$theta1,
    overdisp_inat = rp$overdisp_inat,
    year_beta = rp$year_beta, year_var = rp$year_var,
    sigma_year_beta = max(abs(rp$year_beta) * 2, 0.1),
    sigma_year_var  = max(abs(rp$year_var)  * 2, 0.1),
    year_effect = year_effect_true,
    tau_year = 1 / var(year_effect_true)
  )
}

#' @name make_true_year_region
#' @description Build the KNOWN true ecoregion-level trend deviations for the
#'   two-scenario ecoregion simulation study (see 01e_run_ecoregion_sim.R and
#'   01e_run_abundance_sweep.R).
#' @param cell100_geo Cell100 table carrying an `ecoregion` factor / integer
#'   `ecoregion_id` column (from 00b_prep_ecoregion.R, baked into
#'   prepped_sim_inputs.RDS) -- nregion is derived from it
#'   (nlevels(cell100_geo$ecoregion)), not passed directly, so callers only
#'   need to thread through the one table they already have.
#' @param scenario "varying": a deterministic, evenly-spaced spread from
#'   -target_sd*k to +target_sd*k (guarantees a clear mix of clearly-positive,
#'   clearly-negative, and near-zero regions regardless of K, rather than
#'   relying on a random draw landing well for a small number of regions) --
#'   tests recovery of real regional structure. "null": every region gets
#'   the SAME trend (all zero deviations) -- tests whether the model invents
#'   regional structure that isn't there.
#' @param target_sd SD target for "varying" (matches the CAR prototype's
#'   amplitude=0.3 field, SD~0.2, for comparable effect size).
#' @return Numeric vector, length nregion, zero-mean, ordered to match
#'   levels(cell100_geo$ecoregion) / cell100_geo$ecoregion_id (1:nregion).
make_true_year_region <- function(cell100_geo, scenario = c("varying", "null"), target_sd = 0.2) {
  scenario <- match.arg(scenario)
  nregion <- nlevels(cell100_geo$ecoregion)
  if (scenario == "null") return(rep(0, nregion))
  raw <- seq(-1, 1, length.out = nregion)
  centered <- raw - mean(raw)
  centered * (target_sd / sd(centered))
}

#' @name true_param_list_ecoregion
#' @description Like true_param_list(), but for the ecoregion model: swaps
#'   year_effect/tau_year (CAR field) for year_region/sigma_region. Every
#'   other parameter is identical -- fixed at its REAL posterior-mean value
#'   from Arielle's converged fit, including the climate SVCs
#'   (MWMT_effect/MCMT_effect), so the simulation tests whether the
#'   ecoregion trend is identifiable ALONGSIDE the temperature response.
#' @param real_post_means List from prepped_sim_inputs.RDS.
#' @param year_region_true Numeric vector, length nregion (from
#'   make_true_year_region()).
#' @return Named list suitable as nimbleModel(inits = ...) for
#'   model_code_ecoregion_trend.
true_param_list_ecoregion <- function(real_post_means, year_region_true) {
  rp <- real_post_means
  list(
    occ_beta = unname(rp$occ_beta),
    link_occ_intercept = unname(rp$link_occ_intercept),
    MWMT_effect = unname(rp$MWMT_effect),
    MCMT_effect = unname(rp$MCMT_effect),
    intercept_tau = rp$intercept_tau, MWMT_tau = rp$MWMT_tau, MCMT_tau = rp$MCMT_tau,
    det_intercept = rp$det_intercept,
    p_beta = unname(rp$p_beta),
    theta0 = rp$theta0, theta1 = rp$theta1,
    overdisp_inat = rp$overdisp_inat,
    year_beta = rp$year_beta, year_var = rp$year_var,
    sigma_year_beta = max(abs(rp$year_beta) * 2, 0.1),
    sigma_year_var  = max(abs(rp$year_var)  * 2, 0.1),
    year_region = year_region_true,
    sigma_region = max(sd(year_region_true), 0.05)
  )
}

#' @name simulate_replicate_data
#' @description Build the extended model with all parameters fixed at their
#'   TRUE values (data nodes left unassigned) and draw y / y_inat from the
#'   model's own likelihoods (dOcc_v, dnbinom) via nimble's simulate(). This
#'   is what makes the simulation faithful to the actual fitted model rather
#'   than a hand-rolled approximation of it.
#' @param model_code The (shared, forked) nimbleCode object.
#' @param constants Reduced constants_list from build_reduced_constants().
#' @param inat_effort Reduced inat_effort matrix (data, not simulated).
#' @param y_ncol Number of occasion columns for the camera y array.
#' @param truth Named list from true_param_list().
#' @return List with simulated y (nsite x y_ncol) and y_inat (ncell50 x nyear).
simulate_replicate_data <- function(model_code, constants, inat_effort, y_ncol, truth) {
  # y[i, 1:J[i]] ~ dOcc_v(...) is a RAGGED vectorized declaration -- J[i]
  # varies by site, and when the real data includes a J[i]==1 site, NIMBLE's
  # automatic node-dimension inference can disagree with the fixed
  # nsite x y_ncol shape the model actually needs (a documented NIMBLE
  # gotcha for vectorized declarations with a variable upper bound; see the
  # NIMBLE user manual section on specifying `dimensions`). y_ncol was
  # accepted as a parameter and documented ("Number of occasion columns for
  # the camera y array") but never actually passed to nimbleModel() -- this
  # is almost certainly the cause of the "Dimension of 'y[i, 1:J[i]]' does
  # not match required dimension for the distribution 'dOcc_v'. Necessary
  # dimension is 1" compile error hit on Hazel for ALL THREE arms (this
  # node is shared/byte-identical across estimators). Passing it explicitly
  # pins y's shape regardless of any individual J[i] value.
  # NOT YET VERIFIED: dyn.load()/compileNimble() is blocked in this sandbox,
  # so this fix is untested past nimbleModel() construction -- confirm on
  # Hazel before trusting it.
  truth_model <- nimbleModel(
    model_code, constants = constants,
    data = list(inat_effort = inat_effort),
    inits = truth, calculate = FALSE,
    dimensions = list(y = c(constants$nsite, y_ncol)))
  truth_model$calculate()
  simulate(truth_model, nodes = c("y_inat", "y"), includeData = TRUE)
  list(y = truth_model$y, y_inat = truth_model$y_inat)
}

#' @name fit_replicate
#' @description Fit the extended model to one replicate's simulated data:
#'   single chain (this is a recovery/coverage calibration across many
#'   replicates, not a per-replicate convergence certification -- see
#'   README), monitors include year_effect and tau_year.
#'   IMPORTANT: every stochastic node needs a real starting value here, not
#'   just the new trend nodes -- CAR fields (link_occ_intercept, MWMT_effect,
#'   MCMT_effect) depend on their tau, and if tau is left unset nimble's
#'   auto-init tries to simulate the CAR field from an undefined precision
#'   and produces NaN, which then cascades into every downstream
#'   deterministic node (psi, mu, ...). Caught via smoke-test job 435929.
#' @param model_code,constants,inat_effort As above.
#' @param sim_data List with y, y_inat (from simulate_replicate_data()).
#' @param base_inits Real inits_list from input_data_bobcat.RDS (occ_beta,
#'   p_beta, det_intercept, link_occ_intercept, intercept_tau, overdisp_inat,
#'   theta0, theta1, MWMT_effect, MCMT_effect, MWMT_tau, MCMT_tau all reused
#'   as-is -- their dimensions are untouched by the cell50/site reduction).
#' @param n_burnin,n_iter MCMC budget for this replicate.
#' @param trend_inits Named list of inits for whatever trend nodes this
#'   model_code actually has, overriding base_inits -- defaults to the CAR
#'   model's (year_beta, year_var, tau_year, year_effect), preserving exact
#'   prior behavior for existing callers. Pass e.g.
#'   list(year_beta=0, year_var=0, sigma_region=1, year_region=rep(0,nregion))
#'   for the ecoregion model, or list(year_beta=0, year_var=0) for the
#'   national-scalar model (no extra trend nodes at all).
#' @param extra_monitors Character vector of additional monitor names beyond
#'   the common set -- defaults to c("year_effect","tau_year") for the CAR
#'   model; pass c("year_region","sigma_region") for the ecoregion model, or
#'   character(0) for national-scalar.
#' @param compute_waic If TRUE, enables WAIC in configureMCMC and returns
#'   list(samples=..., waic=...) instead of the bare samples matrix (default
#'   FALSE preserves the original return type for existing callers).
#' @return Matrix of post-burnin posterior samples (all monitored nodes), or
#'   (if compute_waic) a list with that matrix plus the WAIC value.
fit_replicate <- function(model_code, constants, inat_effort, sim_data, base_inits,
                          n_burnin = 1000, n_iter = 3000,
                          trend_inits = list(year_beta = 0, year_var = 0, tau_year = 1,
                                             year_effect = rep(0, constants$ncell100)),
                          extra_monitors = c("year_effect", "tau_year"),
                          compute_waic = FALSE) {
  fit_inits <- base_inits
  fit_inits[names(trend_inits)] <- trend_inits

  # Data list is built from whatever camera-observation field sim_data
  # actually carries, not hardcoded to 'y'. The camera-level models (baseline,
  # national scalar) supply sim_data$y; the array-level arms (array_occ,
  # array_rn) supply sim_data$w instead (see aggregate_to_array() /
  # build_array_constants() in sim_helpers_array.R) and have NO 'y' node at
  # all. Hardcoding data=list(y=sim_data$y, ...) here silently dropped 'w' --
  # nimbleModel() does not error on an unused data-list entry ('y' is
  # provided in data but is not a variable in the model and is being
  # ignored"), so the array models would have sampled w from its OWN PRIOR
  # instead of conditioning on the simulated detections. Caught via local
  # compile testing before any array_occ/array_rn fit was attempted.
  cam_data <- list()
  if (!is.null(sim_data$y)) cam_data$y <- sim_data$y
  if (!is.null(sim_data$w)) cam_data$w <- sim_data$w
  if (length(cam_data) == 0) {
    stop("fit_replicate(): sim_data has neither $y nor $w -- no camera ",
         "observation data to condition on.")
  }

  # N_a is array_rn's latent discrete abundance node. It has no default
  # nimble auto-init that respects the observed detections, so an explicit
  # starting value (>= observed detecting-camera count per array, per the
  # N_a_init the 01i driver already computes) is required for the RN arm to
  # initialize in a non-zero-probability region.
  if (!is.null(sim_data$N_a_init)) fit_inits$N_a <- sim_data$N_a_init

  # y[i, 1:J[i]] ~ dOcc_v(...) (camera arms) and w[a, 1:n_a[a]] ~ dOcc_v(...)
  # (array_occ) are RAGGED vectorized declarations -- whenever some J[i] (or
  # n_a[a]) equals 1, "1:1" collapses to a scalar rather than a length-1
  # vector, and dOcc_v (a vector-valued distribution) rejects a scalar with
  # "Dimension of '...' does not match required dimension ... Necessary
  # dimension is 1" -- a compile-time-only failure (nimbleModel() build
  # succeeds; compileNimble() fails). This is the SAME mechanics as the
  # array-singleton (n_a[a]==1) bug already fixed by dropping singleton
  # arrays upstream -- but J[i]==1 camera sites are real, retained data, not
  # a droppable degenerate case, so this fixes it the general way: pin the
  # node's dimensions explicitly to the actual supplied matrix's shape
  # (which is already the correct fixed nsite/narray x max-occasions shape
  # regardless of any individual row's real length) so NIMBLE's inference
  # never collapses it.
  # NOT YET VERIFIED: dyn.load()/compileNimble() is blocked in this sandbox;
  # this fix is untested past nimbleModel() construction -- confirm on Hazel.
  fit_dims <- list()
  if (!is.null(cam_data$y)) fit_dims$y <- dim(cam_data$y)
  if (!is.null(cam_data$w)) fit_dims$w <- dim(cam_data$w)

  fit_model <- nimbleModel(
    model_code, constants = constants,
    data = c(cam_data, list(y_inat = sim_data$y_inat, inat_effort = inat_effort)),
    inits = fit_inits,
    dimensions = fit_dims,
    calculate = FALSE)

  Cmodel <- compileNimble(fit_model)
  conf <- configureMCMC(fit_model, enableWAIC = compute_waic)
  conf$addMonitors(c("occ_beta", "link_occ_intercept", "year_beta", "year_var",
                     "total_var_beta", "MWMT_effect", "MCMT_effect",
                     "trend_robust_indicator", extra_monitors))
  # enableWAIC only goes on configureMCMC() -- buildMCMC() rejects it when
  # given a config object (only allowed when building directly off a raw
  # model), confirmed via ecoregion smoke test job 444429.
  mcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(mcmc, project = fit_model)

  Cmcmc$run(n_burnin + n_iter)
  samples <- as.matrix(Cmcmc$mvSamples)
  post_samples <- samples[(n_burnin + 1):(n_burnin + n_iter), , drop = FALSE]

  if (!compute_waic) return(post_samples)
  waic_result <- calculateWAIC(Cmcmc, Cmodel)
  list(samples = post_samples, waic = waic_result$WAIC)
}

#' @name make_true_year_effect_grain
#' @description Build a KNOWN true spatially-varying trend deviation with
#'   CONTROLLABLE spatial grain, for the phase-1 grain sweep (how fine a
#'   patch size can the model actually detect, given current data density).
#'   Unlike make_true_year_effect() (a fixed NE+/West- regional contrast),
#'   this seeds white noise on the real CAR graph and diffuses it via the
#'   existing smooth_over_car_graph() for n_smooth_iter passes: fewer passes
#'   = patchier/finer grain, more passes = smoother/broader.
#'   n_smooth_iter = 3 reproduces the exact grain of make_true_year_effect()'s
#'   already-validated continental-scale field (same smoothing recipe), so
#'   that existing result can be reused as one point on the grain curve
#'   instead of re-simulated. Rescaled to match that same field's SD across
#'   every grain level, so grain is the only thing that varies -- not overall
#'   signal magnitude -- and centered to exact zero mean (the fitted model's
#'   zero_mean=1 constraint).
#' @param cell100_geo Cell100 table (only nrow() used here -- the truth is
#'   unstructured noise, not tied to region labels, unlike
#'   make_true_year_effect()).
#' @param adj,num Real CAR adjacency vectors.
#' @param n_smooth_iter Number of smoothing passes (grain control). 0 = pure
#'   white noise (finest possible grain).
#' @param target_sd SD to rescale to; defaults to make_true_year_effect()'s
#'   SD (amplitude=0.3, n_iter=3) for direct comparability across the sweep.
#' @param seed RNG seed for the white-noise draw -- fixed per grain level so
#'   the TRUTH field itself is reproducible/inspectable, independent of the
#'   per-replicate data-simulation seeding done by the caller afterward.
#' @return Numeric vector, length ncell100, zero-mean, SD == target_sd.
make_true_year_effect_grain <- function(cell100_geo, adj, num, n_smooth_iter,
                                        target_sd = NULL, seed = 20260713) {
  ncell <- nrow(cell100_geo)
  if (is.null(target_sd)) {
    reference <- make_true_year_effect(cell100_geo, adj, num, amplitude = 0.3)
    target_sd <- sd(reference)
  }
  set.seed(seed)
  raw <- rnorm(ncell)
  smoothed <- if (n_smooth_iter > 0) smooth_over_car_graph(raw, adj, num, n_iter = n_smooth_iter, w = 0.5) else raw
  centered <- smoothed - mean(smoothed)
  centered * (target_sd / sd(centered))
}

#' @name compute_grain_metrics
#' @description Recovery metrics for the grain sweep: unlike
#'   compute_replicate_metrics(), makes NO reference to NE/West region labels
#'   (not meaningful once the truth is unstructured, patchy noise rather than
#'   a two-region contrast) -- global bias/RMSE/coverage over informed cells,
#'   plus the Pearson correlation between the true and posterior-mean
#'   year_effect field, a shape-recovery metric that works for any true
#'   spatial pattern.
#' @param samples,year_effect_true,informed_cell100 As in compute_replicate_metrics().
#' @return One-row data frame of metrics.
compute_grain_metrics <- function(samples, year_effect_true, cell100_geo, informed_cell100) {
  ye_cols <- paste0("year_effect[", seq_along(year_effect_true), "]")
  est_mean <- colMeans(samples[, ye_cols, drop = FALSE])
  est_lb   <- apply(samples[, ye_cols, drop = FALSE], 2, quantile, probs = 0.025)
  est_ub   <- apply(samples[, ye_cols, drop = FALSE], 2, quantile, probs = 0.975)

  is_informed <- seq_along(year_effect_true) %in% informed_cell100
  err <- est_mean - year_effect_true
  covered <- (year_effect_true >= est_lb) & (year_effect_true <= est_ub)

  m <- is_informed
  spatial_cor <- if (sum(m) > 2) cor(est_mean[m], year_effect_true[m]) else NA

  data.frame(
    bias_all = mean(err[m]), rmse_all = sqrt(mean(err[m]^2)),
    coverage_all = mean(covered[m]), n_informed = sum(m),
    spatial_cor = spatial_cor,
    row.names = NULL
  )
}

#' @name compute_replicate_metrics
#' @description Compute per-replicate recovery metrics for year_effect:
#'   bias/RMSE against truth, 95% CI coverage, restricted to the "informed"
#'   CAR cells (those with at least one retained cell50 mapped into them --
#'   cells with no data are prior/neighbor-smoothed only and aren't a
#'   meaningful recoverability check). Also the NE-vs-West sign check.
#' @param samples Posterior samples matrix from fit_replicate().
#' @param year_effect_true Length-ncell100 truth vector.
#' @param cell100_geo Region table.
#' @param informed_cell100 Vector of cell100 ids with retained data this replicate.
#' @return One-row data frame of metrics.
compute_replicate_metrics <- function(samples, year_effect_true, cell100_geo, informed_cell100) {
  ye_cols <- paste0("year_effect[", seq_along(year_effect_true), "]")
  est_mean <- colMeans(samples[, ye_cols, drop = FALSE])
  est_lb   <- apply(samples[, ye_cols, drop = FALSE], 2, quantile, probs = 0.025)
  est_ub   <- apply(samples[, ye_cols, drop = FALSE], 2, quantile, probs = 0.975)

  is_informed <- seq_along(year_effect_true) %in% informed_cell100
  err <- est_mean - year_effect_true
  covered <- (year_effect_true >= est_lb) & (year_effect_true <= est_ub)

  region <- cell100_geo$region
  region_metrics <- function(mask) {
    m <- mask & is_informed
    if (sum(m) == 0) return(c(bias = NA, rmse = NA, coverage = NA, n = 0))
    c(bias = mean(err[m]), rmse = sqrt(mean(err[m]^2)),
      coverage = mean(covered[m]), n = sum(m))
  }
  ne_m    <- region_metrics(region == "NE")
  west_m  <- region_metrics(region == "West")
  other_m <- region_metrics(region == "other")
  all_m   <- region_metrics(rep(TRUE, length(region)))

  ne_minus_west_est  <- mean(est_mean[region == "NE" & is_informed]) -
                        mean(est_mean[region == "West" & is_informed])
  sign_recovered <- ne_minus_west_est > 0   # true NE-West is positive by construction

  data.frame(
    bias_all = all_m["bias"], rmse_all = all_m["rmse"], coverage_all = all_m["coverage"], n_informed = all_m["n"],
    bias_NE = ne_m["bias"], rmse_NE = ne_m["rmse"], coverage_NE = ne_m["coverage"], n_NE = ne_m["n"],
    bias_West = west_m["bias"], rmse_West = west_m["rmse"], coverage_West = west_m["coverage"], n_West = west_m["n"],
    bias_other = other_m["bias"], rmse_other = other_m["rmse"], coverage_other = other_m["coverage"], n_other = other_m["n"],
    ne_minus_west_est = ne_minus_west_est, sign_recovered = sign_recovered,
    row.names = NULL
  )
}

#' @name run_one_replicate
#' @description Simulate + fit + score one replicate, wrapped in tryCatch so
#'   one failed replicate returns an "Errored" row instead of killing the
#'   parallel batch (group's simulation-validation convention).
#' @param metrics_fn Metrics function to score the fit against truth. Called
#'   as metrics_fn(samples, ...) -- whatever extra NAMED arguments a
#'   particular metrics_fn needs (e.g. year_effect_true/cell100_geo/
#'   informed_cell100 for compute_replicate_metrics(), or
#'   year_region_true/ecoregion_levels for compute_ecoregion_metrics()) are
#'   passed straight through via this function's own `...` -- existing
#'   callers that already pass those as named args need no changes; R's `...`
#'   forwarding matches them by name.
#' @param trend_inits,extra_monitors Passed through to fit_replicate() for
#'   the PRIMARY model -- defaults match the CAR model's, unchanged from
#'   before this was generalized.
#' @param fit_null_scalar If TRUE, ALSO fits model_code_national_scalar (the
#'   national-scalar production structure with no regional trend term at
#'   all -- must already be sourced/in scope as a global, from
#'   model_code_national_scalar.R) to the SAME simulated data, and attaches
#'   both models' WAIC plus their difference -- "is the added complexity
#'   actually warranted" comparison (see 01e_run_abundance_sweep.R).
#' @param n_burnin,n_iter MCMC budget (both models, if fit_null_scalar).
#' @param ... Extra named args forwarded to metrics_fn().
#' @return One-row data frame: rep id, status, metrics (NA if errored),
#'   elapsed_sec, sim_data (list-column, the raw simulated y/y_inat -- e.g.
#'   for summarize_simulated_information() in sim_helpers_abundance.R), and
#'   (if fit_null_scalar) waic_primary/waic_national/waic_diff (positive =
#'   primary model preferred, i.e. lower/better WAIC).
run_one_replicate <- function(rep_id, model_code, constants, inat_effort, y_ncol,
                              truth, base_inits, n_burnin = 1000, n_iter = 3000,
                              metrics_fn = compute_replicate_metrics,
                              trend_inits = list(year_beta = 0, year_var = 0, tau_year = 1,
                                                 year_effect = rep(0, constants$ncell100)),
                              extra_monitors = c("year_effect", "tau_year"),
                              fit_null_scalar = FALSE, ...) {
  t0 <- Sys.time()
  result <- tryCatch({
    sim_data <- simulate_replicate_data(model_code, constants, inat_effort, y_ncol, truth)

    if (fit_null_scalar) {
      fit_primary <- fit_replicate(model_code, constants, inat_effort, sim_data, base_inits,
                                   n_burnin, n_iter, trend_inits = trend_inits,
                                   extra_monitors = extra_monitors, compute_waic = TRUE)
      fit_national <- fit_replicate(model_code_national_scalar, constants, inat_effort, sim_data,
                                    base_inits, n_burnin, n_iter,
                                    trend_inits = list(year_beta = 0, year_var = 0),
                                    extra_monitors = character(0), compute_waic = TRUE)
      metrics <- metrics_fn(fit_primary$samples, ...)
      out <- cbind(rep = rep_id, status = "ok", metrics,
                  waic_primary = fit_primary$waic, waic_national = fit_national$waic,
                  waic_diff = fit_national$waic - fit_primary$waic)
    } else {
      samples <- fit_replicate(model_code, constants, inat_effort, sim_data, base_inits,
                               n_burnin, n_iter, trend_inits = trend_inits,
                               extra_monitors = extra_monitors)
      metrics <- metrics_fn(samples, ...)
      out <- cbind(rep = rep_id, status = "ok", metrics)
    }
    out$sim_data <- list(sim_data)
    out
  }, error = function(e) {
    message("Rep ", rep_id, " ERRORED: ", conditionMessage(e))
    data.frame(rep = rep_id, status = "errored", detail = conditionMessage(e))
  })
  result$elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  result
}

#' @name compute_ecoregion_metrics
#' @description Per-region and overall recovery metrics for year_region,
#'   scored against whatever true value was used to simulate. Serves BOTH
#'   scenarios identically: under VARYING the truth is nonzero/region-
#'   specific, under NULL it's uniformly 0 -- "coverage of the true value" is
#'   exactly "does the CI cover 0" in that case, no separate branch needed.
#'   Also reports sigma_region's posterior summary (should concentrate near
#'   0 under NULL, away from 0 under VARYING) and frac_sign_correct: the
#'   fraction of regions where sign(posterior mean) matches sign(truth) --
#'   meaningless under NULL (truth is exactly 0 everywhere, sign undefined)
#'   but the region-recovery headline number under VARYING.
#' @param samples Posterior samples matrix from fit_replicate() (ecoregion model).
#' @param year_region_true Length-nregion truth vector.
#' @param ecoregion_levels Character vector of region names, length nregion
#'   (for readable per-region column names).
#' @return One-row data frame of metrics. Called via run_one_replicate()'s
#'   metrics_fn(samples, ...) -- year_region_true/ecoregion_levels are passed
#'   through run_one_replicate()'s own `...`.
compute_ecoregion_metrics <- function(samples, year_region_true, ecoregion_levels) {
  nregion <- length(year_region_true)
  yr_cols <- paste0("year_region[", seq_len(nregion), "]")
  est_mean <- colMeans(samples[, yr_cols, drop = FALSE])
  est_lb   <- apply(samples[, yr_cols, drop = FALSE], 2, quantile, probs = 0.025)
  est_ub   <- apply(samples[, yr_cols, drop = FALSE], 2, quantile, probs = 0.975)

  err <- est_mean - year_region_true
  covered <- (year_region_true >= est_lb) & (year_region_true <= est_ub)
  nonzero <- year_region_true != 0
  frac_sign_correct <- if (any(nonzero)) {
    mean(sign(est_mean[nonzero]) == sign(year_region_true[nonzero]))
  } else NA_real_

  region_names <- make.names(ecoregion_levels)
  per_region_est <- as.data.frame(t(setNames(est_mean, paste0("est_", region_names))))
  per_region_err <- as.data.frame(t(setNames(err, paste0("err_", region_names))))
  per_region_cov <- as.data.frame(t(setNames(as.numeric(covered), paste0("cov_", region_names))))

  cbind(
    data.frame(bias_all = mean(err), rmse_all = sqrt(mean(err^2)), coverage_all = mean(covered),
              frac_sign_correct = frac_sign_correct),
    sigma_region_mean = mean(samples[, "sigma_region"]),
    sigma_region_lb = unname(quantile(samples[, "sigma_region"], 0.025)),
    sigma_region_ub = unname(quantile(samples[, "sigma_region"], 0.975)),
    per_region_est, per_region_err, per_region_cov
  )
}

# =============================================================================
# RJMCMC indicator test (model_code_ecoregion_switch) -- ADDITIONS ONLY,
# nothing above this line is touched. See rjmcmc_indicator_test_plan.md and
# 01f_run_indicator_test.R. These functions are dedicated to the switch
# model: it needs >=2 MCMC chains (Gelman-Rubin needs multiple chains) and an
# explicit configureRJ() call, which fit_replicate()/run_one_replicate()
# don't do and shouldn't be complicated to do -- kept as a separate,
# parallel entry point instead of overloading the existing ones.
# =============================================================================

#' @name fit_replicate_switch
#' @description Fit model_code_ecoregion_switch to one replicate's simulated
#'   data with an explicitly-configured RJMCMC sampler on the gamma
#'   indicator (configureRJ(), NOT the default sampler NIMBLE would otherwise
#'   assign to a lone Bernoulli node). Runs >=2 chains from deliberately
#'   different gamma starting values (chain 1 starts gamma=1/"on", chain 2
#'   starts gamma=0/"off") specifically so convergence isn't an artifact of a
#'   lucky start.
#'   CONVERGENCE GOTCHA (rjmcmc_indicator_test_plan.md gotcha #1): do NOT
#'   judge convergence from Gelman-Rubin on gamma itself -- a Bernoulli
#'   indicator can sit at 0 or 1 for long stretches even while mixing
#'   perfectly well, so R-hat on gamma is not a meaningful diagnostic here.
#'   Instead this function computes Gelman-Rubin (coda::gelman.diag) and
#'   effective sample size (coda::effectiveSize) on total_data_logLik -- the
#'   model's own log-likelihood of the observed y/y_inat data at each
#'   retained draw (see the logLik_y/logLik_y_inat_year/total_data_logLik
#'   nodes added in model_code_ecoregion_switch.R) -- and reports THAT as the
#'   convergence readout. gamma is monitored separately, only to compute its
#'   posterior mean (the inclusion probability / support for the switch).
#' @param constants Must include pi_gamma (prior inclusion probability for
#'   gamma) and slab_sd (FIXED sd for year_region_raw's prior -- see model
#'   file header for why this can't be an estimated node under RJMCMC),
#'   beyond the usual ecoregion constants (ecoregion_of_cell100, nregion).
#' @param n_chains Number of chains (>=2 required for Gelman-Rubin; default 2).
#' @return List: gamma_mean (pooled posterior mean of gamma across chains,
#'   the P(gamma=1) readout), rhat_logLik, ess_logLik (convergence, judged on
#'   total_data_logLik, NOT gamma), gamma_by_chain (per-chain means, a cheap
#'   sanity check that both starting points agree), n_chains, n_iter_per_chain.
fit_replicate_switch <- function(constants, inat_effort, sim_data, base_inits,
                                 n_burnin = 1000, n_iter = 3000, n_chains = 2,
                                 pi_gamma = 0.5, slab_sd = 0.227) {
  stopifnot(!is.null(constants$pi_gamma), !is.null(constants$slab_sd))
  nregion <- constants$nregion

  fit_model <- nimbleModel(
    model_code_ecoregion_switch, constants = constants,
    data = list(y = sim_data$y, y_inat = sim_data$y_inat, inat_effort = inat_effort),
    inits = c(base_inits, list(gamma = 1, year_region_raw = rep(0, nregion))),
    calculate = FALSE)

  Cmodel <- compileNimble(fit_model)
  conf <- configureMCMC(fit_model)
  conf$addMonitors(c("occ_beta", "link_occ_intercept", "year_beta", "year_var",
                     "total_var_beta", "MWMT_effect", "MCMT_effect",
                     "trend_robust_indicator", "gamma", "year_region",
                     "year_region_raw", "total_data_logLik"))

  # Explicit RJMCMC configuration -- NOT the default sampler NIMBLE would
  # assign to a lone Bernoulli node (a plain binary sampler, which is the
  # Kuo & Mallick indicator-selection sampler, not Green's reversible jump).
  # ONE shared indicator (gamma) gates all nregion target nodes -- verified
  # feasible for a multi-target shared indicator via toy models
  # (rj_group_test.R/rj_group_test2.R); the real blocker turned out to be
  # that RJMCMC requires the target nodes' own prior to have a CONSTANT
  # hyperparameter (slab_sd, not an estimated sigma_region -- see model file
  # header), not the shared-indicator structure itself.
  configureRJ(conf,
             targetNodes    = paste0("year_region_raw[", seq_len(nregion), "]"),
             indicatorNodes = rep("gamma", nregion),
             control        = list(mean = 0, scale = slab_sd))

  mcmc <- buildMCMC(conf)
  Cmcmc <- compileNimble(mcmc, project = fit_model)

  # Deliberately different gamma starts per chain (base_inits/other nodes
  # shared) -- a convergence check that's insensitive to the RJ sampler's
  # starting regime is a stronger claim than one that only ever started "on".
  inits_list <- lapply(seq_len(n_chains), function(ch) {
    g0 <- if (ch == 1) 1 else 0
    c(base_inits, list(gamma = g0, year_region_raw = rep(0, nregion)))
  })

  samples <- runMCMC(Cmcmc, niter = n_burnin + n_iter, nburnin = n_burnin,
                     nchains = n_chains, inits = inits_list,
                     samplesAsCodaMCMC = TRUE, progressBar = FALSE)
  # samples is a coda::mcmc.list when n_chains > 1 (each element one chain's
  # post-burnin draws), or a single coda::mcmc object when n_chains == 1.
  chain_list <- if (n_chains > 1) samples else list(samples)

  gamma_by_chain <- sapply(chain_list, function(ch) mean(ch[, "gamma"]))
  gamma_pooled <- mean(unlist(lapply(chain_list, function(ch) ch[, "gamma"])))

  logLik_mcmc_list <- coda::mcmc.list(lapply(chain_list, function(ch)
    coda::mcmc(ch[, "total_data_logLik", drop = FALSE])))
  rhat_logLik <- if (n_chains > 1) {
    tryCatch(coda::gelman.diag(logLik_mcmc_list)$psrf[1, "Point est."],
            error = function(e) NA_real_)
  } else NA_real_
  ess_logLik <- sum(sapply(chain_list, function(ch) coda::effectiveSize(ch[, "total_data_logLik"])))

  list(gamma_mean = gamma_pooled, gamma_by_chain = gamma_by_chain,
      rhat_logLik = rhat_logLik, ess_logLik = ess_logLik,
      n_chains = n_chains, n_iter_per_chain = n_iter)
}

#' @name run_one_indicator_replicate
#' @description Simulate-or-reuse + fit_replicate_switch() + package one
#'   replicate's result, tryCatch-wrapped so one failure returns an
#'   "errored" row rather than killing the array task (same convention as
#'   run_one_replicate()). Unlike run_one_replicate(), does NOT re-simulate
#'   by default -- sim_data is REUSED from the corrected abundance sweep's
#'   stored rep_*.RDS files (the whole point of the indicator test is to
#'   read the SAME 180 datasets under a different model/method, not to draw
#'   new ones) -- pass sim_data directly. Only falls back to simulating fresh
#'   data if sim_data is NULL, using the identical per-replicate seed
#'   convention (set.seed(20260712 + row_id) by the caller, BEFORE calling
#'   this function) as 01e_run_abundance_sweep.R, so the fallback path
#'   produces byte-identical data to what's already stored.
#' @param sim_data List(y=, y_inat=) reloaded from a stored rep_*.RDS's
#'   res$sim_data[[1]], or NULL to simulate fresh (fallback path only).
#' @return One-row data frame: rep id, status, gamma_mean, rhat_logLik,
#'   ess_logLik, gamma_by_chain (semicolon-joined string, informational),
#'   elapsed_sec.
run_one_indicator_replicate <- function(rep_id, constants, inat_effort, y_ncol,
                                        truth, base_inits, sim_data = NULL,
                                        n_burnin = 1000, n_iter = 3000, n_chains = 2) {
  t0 <- Sys.time()
  result <- tryCatch({
    if (is.null(sim_data)) {
      sim_data <- simulate_replicate_data(model_code_ecoregion, constants, inat_effort, y_ncol, truth)
    }
    fit <- fit_replicate_switch(constants, inat_effort, sim_data, base_inits,
                                n_burnin = n_burnin, n_iter = n_iter, n_chains = n_chains,
                                pi_gamma = constants$pi_gamma, slab_sd = constants$slab_sd)
    data.frame(rep = rep_id, status = "ok",
              gamma_mean = fit$gamma_mean,
              rhat_logLik = fit$rhat_logLik, ess_logLik = fit$ess_logLik,
              gamma_by_chain = paste(round(fit$gamma_by_chain, 4), collapse = ";"),
              n_chains = fit$n_chains, n_iter_per_chain = fit$n_iter_per_chain)
  }, error = function(e) {
    message("Rep ", rep_id, " ERRORED: ", conditionMessage(e))
    data.frame(rep = rep_id, status = "errored", detail = conditionMessage(e))
  })
  result$elapsed_sec <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  result
}
