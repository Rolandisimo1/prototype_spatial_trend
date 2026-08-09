#!/usr/bin/env Rscript
# =============================================================================
# model_code_spatial_trend.R
# FORKED from Arielle's $PROJ/HPC/bobcat/HPC_run_model_chunks_chain1.R
# (nimbleCode block only; sha256 of the original recorded in DIFF_NOTES.md in
# this directory). This is the ONLY file that changes the model structure --
# both the sim-validation script and the (later, conditional) real bobcat
# refit source this shared definition, so they are guaranteed to fit the
# identical model.
#
# CHANGE: add a spatially-varying trend deviation, mirroring the existing
# MWMT_effect / MCMT_effect temperature CAR fields:
#   year_effect[1:ncell100] ~ dcar_normal(adj[], num[], tau_year, zero_mean=1)
#   tau_year ~ dgamma(0.1, 0.1)
# year_beta remains the GLOBAL mean trend (unchanged prior, unchanged use in
# the camera occupancy submodel). year_effect is a zero-mean spatial
# deviation added ONLY into the iNat intensity pathway (mu[g,t], via the
# total_var_beta argument passed to calcIntensity_SVC): the linear predictor
# term that was `total_var_beta * year_dat[t]` (where
# total_var_beta = year_beta + year_var) becomes
# `(total_var_beta + year_effect[cell100[g]]) * year_dat[t]`
# i.e. year_beta + year_var + year_effect[cell100[g]], all still added
# identically to Arielle's existing intensity linear predictor.
#
# SCOPE DECISIONS FOR THIS PROTOTYPE (both explicit in the task spec):
#   * year_var (the existing iNat-vs-camera global deviation) is left GLOBAL.
#     A second, separately-spatial iNat deviation field is a later extension,
#     not attempted here -- adding two competing spatial fields to the same
#     pathway in one step would be an identifiability nightmare to debug.
#   * The camera occupancy term `year_beta * year_occ[i]` in cloglog(psi[i])
#     is UNCHANGED -- the spatial trend deviation enters only the iNat
#     pathway in this prototype. calcIntensity_noSVC (the !hasSVC branch) is
#     likewise unchanged, since year_effect is only defined/meaningful when
#     hasSVC (bobcat, the validation species, has hasSVC = TRUE).
#   * iNat dnbinom likelihood, camera detection submodel, and the
#     link_occ_intercept CAR field are byte-for-byte unchanged.
# =============================================================================

model_code_spatial_trend <- nimbleCode({

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
          # SPATIAL TREND FIX: was `total_var_beta` (global year_beta+year_var
          # only); now adds the zero-mean spatial deviation for this grid
          # cell's CAR cell, i.e. year_beta + year_var + year_effect[cell100[g]]
          total_var_beta = total_var_beta + year_effect[inat_cell100[g]]
        )

      } else {

        # unchanged: no SVC infrastructure (incl. no year_effect) for
        # non-SVC species in this prototype -- out of scope, see header note
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

      #iNat input data also indexed by year -- UNCHANGED
      y_inat[inat_cells_by_year[k,t],t] ~ dnbinom(size = 1 / overdisp_inat,
                                                  prob = 1 / (1 + overdisp_inat * inat_effort[inat_cells_by_year[k,t],t] *
                                                                mu[inat_cells_by_year[k,t],t]))

    }
  }

  #Occ detection model does not change -- UNCHANGED
  for (i in 1:nsite) {
    for (j in 1:J[i]) {
      cloglog(p[i, j]) <- link_det_intercept +
        p_beta[1] * yday[i, j] +
        p_beta[2] * yday[i, j]^2 +
        p_beta[3] * canopy_height[i] +
        p_beta[4] * log_roaddist[i]
    }

    # Year trend added to the occupancy submodel -- UNCHANGED (global year_beta
    # only; the spatial trend deviation is scoped to the iNat pathway only in
    # this prototype, see header note)
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

    # ---- NEW: spatially-varying trend deviation, same CAR template as above.
    # zero_mean = 1 is required for identifiability: without it, year_effect
    # would be confounded with the global year_beta intercept-like term.
    year_effect[1:ncell100] ~ dcar_normal(
      adj = adj[1:nadj],
      num = num[1:nnum],
      tau = tau_year,
      zero_mean = 1
    )

    tau_year ~ dgamma(0.1, 0.1)

  }

  theta0 ~ dnorm(0, sd = 10)
  theta1 ~ dnorm(1, sd = 1)

  overdisp_inat ~ dgamma(shape = 1, scale = 5)

  # Temporal trend priors -- UNCHANGED. total_var_beta remains the global
  # (year_beta + year_var) iNat trend; year_effect is added on top of it only
  # at the calcIntensity_SVC call site above, not redefined here.
  year_beta ~ dnorm(0, sd = sigma_year_beta)
  year_var ~ dnorm(0, sd = sigma_year_var)
  total_var_beta <- year_beta + year_var

  sigma_year_beta ~ dunif(0, 2)
  sigma_year_var ~ dunif(0, 2)

  ###############
  # Robustness indicator -- UNCHANGED (based on the global year_beta/year_var
  # only; not respatialized in this prototype)
  ##############

  # Signal - to - noise ratio
  snr <- year_beta / year_var

  # Probabilistic "robustness" indicator based on snr
  # Posterior prob that annual occ trend dominates annual noise
  trend_robust_indicator <- step(snr - 1) # 1 if >1, 0 if <=1

})
