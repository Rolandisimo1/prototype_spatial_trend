#!/usr/bin/env Rscript
# =============================================================================
# model_code_national_scalar.R
# Exact fork of Arielle's original $PROJ/HPC/bobcat/HPC_run_model_chunks_chain1.R
# nimbleCode block (sha256 19f59aec394f940204cca1f37e7a7a3966c629892a49ba14592faf71cec3a86e,
# see README.md) -- byte-for-byte identical except the object name
# (model_code -> model_code_national_scalar, so it can be sourced alongside
# the other forks without collision) and this header. NO spatial/regional
# trend term of any kind: total_var_beta = year_beta + year_var, applied
# identically to every grid cell, exactly the CURRENT PRODUCTION structure.
#
# WHY THIS FILE EXISTS: the ecoregion-trend simulation study fits BOTH this
# model and model_code_ecoregion_trend.R to the SAME simulated data and
# compares via WAIC, for both a spatially-VARYING and a spatially-NULL truth
# scenario (see 01e_run_ecoregion_sim.R) -- this is the "is the added
# complexity of an ecoregion term actually warranted over the production
# model" comparison the task requires. Never edit this file to differ from
# the real production model; if Arielle's original changes, re-fork it.
# =============================================================================

model_code_national_scalar <- nimbleCode({

  for (g in 1:ncell50) {
    for (t in 1:nyear) {

      #Using original intensity formula from Ben's code, grid cell index only
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
          total_var_beta = total_var_beta
        )

      } else {

        mu[g,t] <- calcIntensity_noSVC(
          intensity_intercept = link_occ_intercept[inat_cell100[g]],
          theta0 = theta0,
          theta1 = theta1,
          beta = occ_beta[1:numOccCovars],
          xdat = xdat_inat[inat_cell50_start[g]:inat_cell50_end[g], 1:numOccCovars, t],
          year_dat = year_vals[t],
          total_var_beta = total_var_beta
        )

      }
    }
  }

  for (t in 1:nyear) {
    for(k in 1:n_cells_year[t]) {

      #iNat input data also indexed by year
      y_inat[inat_cells_by_year[k,t],t] ~ dnbinom(size = 1 / overdisp_inat,
                                                  prob = 1 / (1 + overdisp_inat * inat_effort[inat_cells_by_year[k,t],t] *
                                                                mu[inat_cells_by_year[k,t],t]))

    }
  }

  #Occ detection model does not change
  for (i in 1:nsite) {
    for (j in 1:J[i]) {
      cloglog(p[i, j]) <- link_det_intercept +
        p_beta[1] * yday[i, j] +
        p_beta[2] * yday[i, j]^2 +
        p_beta[3] * canopy_height[i] +
        p_beta[4] * log_roaddist[i]
    }

    #Year trend added to the occupancy submodel
    if (has_SVC) {
      cloglog(psi[i]) <- link_occ_intercept[cell[i]] +
        MWMT[i] * MWMT_effect[cell[i]] +
        MCMT[i] * MCMT_effect[cell[i]] +
        inprod(occ_beta[1:numOccCovars], occ_covars[i, 1:numOccCovars]) +
        year_beta * year_occ[i]

    } else {
      cloglog(psi[i]) <- link_occ_intercept[cell[i]] +
        inprod(occ_beta[1:numOccCovars], occ_covars[i, 1:numOccCovars]) +
        year_beta * year_occ[i]
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

    MWMT_effect[1:ncell100] ~ dcar_normal(
      adj = adj[1:nadj],
      num = num[1:nnum],
      tau = MWMT_tau
    )

    MCMT_effect[1:ncell100] ~ dcar_normal(
      adj = adj[1:nadj],
      num = num[1:nnum],
      tau = MCMT_tau
    )

    MWMT_tau ~ dgamma(0.1, 0.1)
    MCMT_tau ~ dgamma(0.1, 0.1)

  }

  theta0 ~ dnorm(0, sd = 10)
  theta1 ~ dnorm(1, sd = 1)

  overdisp_inat ~ dgamma(shape = 1, scale = 5)

  # Temporal trend priors
  year_beta ~ dnorm(0, sd = sigma_year_beta)
  year_var ~ dnorm(0, sd = sigma_year_var)
  total_var_beta <- year_beta + year_var

  sigma_year_beta ~ dunif(0, 2)
  sigma_year_var ~ dunif(0, 2)

  ###############
  # Robustness indicator
  ##############

  # Signal - to - noise ratio
  snr <- year_beta / year_var

  # Probabilistic “robustness” indicator based on snr
  # Posterior prob that annual occ trend dominates annual noise
  trend_robust_indicator <- step(snr - 1) # 1 if >1, 0 if <=1

})
