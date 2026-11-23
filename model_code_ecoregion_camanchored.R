#!/usr/bin/env Rscript
# =============================================================================
# model_code_ecoregion_camanchored.R
# -----------------------------------------------------------------------------
# Defines model_code_ecoregion_camanchored. FORKED 2026-11-23 from
# model_code_ecoregion_trend.R (left unedited) per SPEC_camera_anchored_regional.md.
#
# WHY THIS FORK: in model_code_ecoregion_trend.R the regional trend deviation
# year_region enters ONLY the iNaturalist intensity (line 71, via the
# year_effect pass-through at line 188); camera occupancy carries bare
# year_beta (lines 117-126). year_region was therefore identified by
# iNaturalist alone, and every sweep that simulated from that model (01e, 01f,
# 01g, 01h) generated camera data with no regional trend at all. Roland's
# decision: cameras are the gold standard for every reported trend, national
# AND regional.
#
# CHANGE: the national shared-parameter + sink structure is mirrored at the
# regional level.
#   camera occupancy   ... + (year_beta + year_region[region_of_site[i]]) * year_occ[i]
#   iNat intensity     total_var_beta + year_effect[inat_cell100[g]], where
#                      year_effect[c] <- year_region[e(c)] + year_region_var[e(c)]
#   so the iNat trend in cell c is
#                      year_beta + year_var + year_region[e(c)] + year_region_var[e(c)]
#   year_region[r]     ~ dnorm(0, sd = sigma_region)      SHARED, camera-anchored, REPORTABLE
#   year_region_var[r] ~ dnorm(0, sd = sigma_region_var)  iNat-only regional sink, NEVER reported
# year_region_var is the regional analogue of year_var: it absorbs iNat-specific
# regional change (observer growth, campaigns) so it cannot leak into the
# camera-anchored term. It is not optional.
#
# REPORTABLE regional trend: year_beta + year_region[r], converted to %/yr with
# the CAMERA year scaling, 100 * (exp((year_beta + year_region[r]) * s_occ) - 1).
#
# OTHER CHANGES (structural only, no change to any likelihood or prior):
#   * The regional nodes are declared outside the has_SVC block, and the
#     non-SVC branches carry them too (the parent excluded non-SVC species).
#   * inat_on (constant, 1/0): when 0 the iNat intensity and likelihood are not
#     built, giving the camera-only control fit from the SAME file. year_var,
#     year_region_var and theta are then prior-only and are not reported.
#
# UNCHANGED, byte-for-byte in substance: camera detection submodel, dOcc_v
# camera likelihood, iNat dnbinom likelihood, link_occ_intercept and climate
# SVC CAR fields, every prior inherited from model_code_ecoregion_trend.R,
# sigma_region ~ dexp(1), and the snr / trend_robust_indicator block.
#
# Requires (in constants): region_of_site (length nsite, integer 1:nregion,
# ecoregion of the site's cell100), ecoregion_of_cell100 (length ncell100),
# nregion, inat_on.
# =============================================================================

#' @name model_code_ecoregion_camanchored
#' @description Camera-anchored ecoregion trend model (see file header).
#'   Regional trend year_region is shared by the camera and iNaturalist streams;
#'   year_region_var is an iNaturalist-only regional sink.
#' @return nimbleCode object.
model_code_ecoregion_camanchored <- nimbleCode({

  if (inat_on) {
    for (g in 1:ncell50) {
      for (t in 1:nyear) {

        if (hasSVC) {
          mu[g,t] <- calcIntensity_SVC(
            intensity_intercept = link_occ_intercept[inat_cell100[g]],
            theta0 = theta0,
            theta1 = theta1,
            MWMT_effect = MWMT_effect[inat_cell100[g]],
            MCMT_effect = MCMT_effect[inat_cell100[g]],
            beta = occ_beta[1:numOccCovars],
            xdat = xdat_inat[inat_cell50_start[g]:inat_cell50_end[g], 1:numOccCovars, t],
            MWMT_dat = MWMT_inat[inat_cell50_start[g]:inat_cell50_end[g], t],
            MCMT_dat = MCMT_inat[inat_cell50_start[g]:inat_cell50_end[g], t],
            year_dat = year_vals[t],
            # call site as in the parent; year_effect now = shared + sink
            total_var_beta = total_var_beta + year_effect[inat_cell100[g]]
          )
        } else {
          mu[g,t] <- calcIntensity_noSVC(
            intensity_intercept = link_occ_intercept[inat_cell100[g]],
            theta0 = theta0,
            theta1 = theta1,
            beta = occ_beta[1:numOccCovars],
            xdat = xdat_inat[inat_cell50_start[g]:inat_cell50_end[g], 1:numOccCovars, t],
            year_dat = year_vals[t],
            total_var_beta = total_var_beta + year_effect[inat_cell100[g]]
          )
        }
      }
    }

    for (t in 1:nyear) {
      for(k in 1:n_cells_year[t]) {
        y_inat[inat_cells_by_year[k,t],t] ~ dnbinom(size = 1 / overdisp_inat,
                                                    prob = 1 / (1 + overdisp_inat * inat_effort[inat_cells_by_year[k,t],t] *
                                                                  mu[inat_cells_by_year[k,t],t]))
      }
    }

    # iNat regional trend deviation = shared camera-anchored term + iNat-only sink
    for (c in 1:ncell100) {
      year_effect[c] <- year_region[ecoregion_of_cell100[c]] + year_region_var[ecoregion_of_cell100[c]]
    }
  }

  # Camera detection submodel -- UNCHANGED
  for (i in 1:nsite) {
    for (j in 1:J[i]) {
      cloglog(p[i, j]) <- link_det_intercept +
        p_beta[1] * yday[i, j] +
        p_beta[2] * yday[i, j]^2 +
        p_beta[3] * canopy_height[i] +
        p_beta[4] * log_roaddist[i]
    }

    # CHANGED: the camera trend now carries the regional deviation of the
    # site's ecoregion, so year_region is identified by the cameras.
    if (has_SVC) {
      cloglog(psi[i]) <- link_occ_intercept[cell[i]] +
        MWMT[i] * MWMT_effect[cell[i]] +
        MCMT[i] * MCMT_effect[cell[i]] +
        inprod(occ_beta[1:numOccCovars], occ_covars[i, 1:numOccCovars]) +
        (year_beta + year_region[region_of_site[i]]) * year_occ[i]
    } else {
      cloglog(psi[i]) <- link_occ_intercept[cell[i]] +
        inprod(occ_beta[1:numOccCovars], occ_covars[i, 1:numOccCovars]) +
        (year_beta + year_region[region_of_site[i]]) * year_occ[i]
    }

    y[i, 1:J[i]] ~ dOcc_v(probOcc = psi[i], probDetect = p[i, 1:J[i]], len = J[i])
  }

  for (i in 1:numOccCovars) {
    if (prior_type == "Normal") {
      occ_beta[i] ~ dnorm(0, sd = 5)
    } else if (prior_type == "Laplace") {
      occ_beta[i] ~ ddexp(0, scale = lambda[interaction_group[i]])
    }
  }

  for (i in 1:4) {
    p_beta[i] ~ dnorm(0, sd = 5)
  }

  det_intercept ~ dunif(0, 1)
  link_det_intercept <- cloglog(det_intercept)

  link_occ_intercept[1:ncell100] ~ dcar_normal(
    adj = adj[1:nadj],
    num = num[1:nnum],
    tau = intercept_tau
  )
  intercept_tau ~ dgamma(0.1, 0.1)

  if (has_SVC) {
    MWMT_effect[1:ncell100] ~ dcar_normal(adj = adj[1:nadj], num = num[1:nnum], tau = MWMT_tau)
    MCMT_effect[1:ncell100] ~ dcar_normal(adj = adj[1:nadj], num = num[1:nnum], tau = MCMT_tau)
    MWMT_tau ~ dgamma(0.1, 0.1)
    MCMT_tau ~ dgamma(0.1, 0.1)
  }

  # ---- regional trend: shared (reportable) + iNat-only sink -------------------
  for (r in 1:nregion) {
    year_region[r] ~ dnorm(0, sd = sigma_region)          # camera-anchored, reported
    year_region_var[r] ~ dnorm(0, sd = sigma_region_var)  # iNat-only sink, never reported
  }
  sigma_region ~ dexp(1)
  sigma_region_var ~ dexp(1)

  theta0 ~ dnorm(0, sd = 10)
  theta1 ~ dnorm(1, sd = 1)
  overdisp_inat ~ dgamma(shape = 1, scale = 5)

  # National trend -- UNCHANGED
  year_beta ~ dnorm(0, sd = sigma_year_beta)
  year_var ~ dnorm(0, sd = sigma_year_var)
  total_var_beta <- year_beta + year_var
  sigma_year_beta ~ dunif(0, 2)
  sigma_year_var ~ dunif(0, 2)

  # Robustness indicator -- UNCHANGED (national terms only)
  snr <- year_beta / year_var
  trend_robust_indicator <- step(snr - 1)
})
