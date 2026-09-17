# model_code_national_scalar_camcount.R
#
# National-scalar integrated model with a COUNT camera submodel.
#
# Differs from model_code_national_scalar.R in the camera likelihood only. The
# iNaturalist likelihood, calcIntensity, the CAR fields, the climate SVCs,
# occ_beta, the trend block and all priors are unchanged, so the two arms are
# directly comparable and year_beta carries the same meaning in each.
#
# New data required (all per deployment i, length nsite):
#   y_cam[i]          integer count of detection events
#   log_effort_cam[i] log(survey_nights[i]); entered as an offset
#   yday_mean[i]      mean day-of-year over the deployment, standardised
# No longer required: y[i, 1:J[i]], J[i], yday[i, j]
#
# New parameter: overdisp_cam. Add to the monitor list alongside overdisp_inat.

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
# model" comparison the task requires. This file is a DELIBERATE fork: the
# camera likelihood differs from the production model by design. Everything
# outside it must stay byte-identical; if the original changes, re-fork and
# re-apply only the camera block below.
# =============================================================================

model_code_national_scalar_camcount <- nimbleCode({

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

  # Camera submodel: detection-event COUNTS with effort as an offset, replacing
  # the baseline's detection/non-detection occupancy submodel.
  #
  # Rationale. cloglog(psi) saturates: once a site is a near-certain detection,
  # presence/absence carries no further information about abundance, and for
  # white-tailed deer a majority of sites sit in that regime (median fitted psi
  # 0.952; among the 72% of deployments that detect, the detection RATE still
  # spans 29-fold between its 10th and 90th percentiles). A count model reads
  # that spread instead of discarding it.
  #
  # Why negative binomial and not Poisson: deer counts have a variance-to-mean
  # ratio of 99.4 across 26,689 deployments. overdisp_cam is the free parameter
  # that absorbs it; a Poisson has none and would understate every interval.
  #
  # NOTE the deliberate parameter sharing with the baseline. link_occ_intercept,
  # occ_beta, the climate SVCs and year_beta keep their names, their priors and
  # their roles, so this arm remains comparable with the occupancy arm and
  # year_beta means the same thing in both: the camera-anchored trend on the
  # log-abundance scale.
  for (i in 1:nsite) {

    # Detectability now multiplies a RATE, not a per-visit probability, so the
    # baseline's per-visit p[i, j] has no counterpart. yday enters as a
    # deployment-level mean rather than per-window.
    log_det_effect[i] <- p_beta[1] * yday_mean[i] +
      p_beta[2] * yday_mean[i]^2 +
      p_beta[3] * canopy_height[i] +
      p_beta[4] * log_roaddist[i]

    if (has_SVC) {
      log(lambda_cam[i]) <- link_occ_intercept[cell[i]] +
        MWMT[i] * MWMT_effect[cell[i]] +
        MCMT[i] * MCMT_effect[cell[i]] +
        inprod(occ_beta[1:numOccCovars], occ_covars[i, 1:numOccCovars]) +
        year_beta * year_occ[i] +
        log_det_effect[i] +
        log_effort_cam[i]

    } else {
      log(lambda_cam[i]) <- link_occ_intercept[cell[i]] +
        inprod(occ_beta[1:numOccCovars], occ_covars[i, 1:numOccCovars]) +
        year_beta * year_occ[i] +
        log_det_effect[i] +
        log_effort_cam[i]
    }

    # log_effort_cam[i] = log(survey_nights[i]) and is DATA, entered with a
    # fixed coefficient of 1 -- the offset. It converts the linear predictor
    # from a rate per camera-night into an expected count for this deployment's
    # actual effort, which ranges 1 to 2,688 nights (median 32).
    y_cam[i] ~ dnbinom(size = 1 / overdisp_cam,
                       prob = 1 / (1 + overdisp_cam * lambda_cam[i]))
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
  # det_intercept is DROPPED from the camera linear predictor, not merely
  # unused. With no repeat-visit structure the camera data see only the
  # product of effort and rate, so a constant c added to
  # link_occ_intercept can be exactly cancelled by subtracting c here, with
  # theta0 absorbing the difference on the iNat side -- the likelihood is
  # flat along that ridge and link_occ_intercept, det_intercept and theta0
  # would all fail to converge. Holding the camera-side detection intercept
  # at 0 restores identification; the baseline log detection rate is then
  # carried by link_occ_intercept, which is what the iNat side reads anyway.
  # det_intercept / link_det_intercept are therefore absent from this fork.

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

  # Camera-side overdispersion, same weakly-informative prior as the iNat side.
  overdisp_cam ~ dgamma(shape = 1, scale = 5)

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
