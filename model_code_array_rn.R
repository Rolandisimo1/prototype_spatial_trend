# =============================================================================
# model_code_array_rn.R
# -----------------------------------------------------------------------------
# Array-level Royle-Nichols variant of the integrated SDM, for the estimator
# comparison simulation.
#
# WHY RN, AND WHY IT IS NOT A DEPARTURE FROM THE EXISTING MODEL
# The whole integration rests on cloglog(psi) = log(lambda), which follows from
# psi = 1 - exp(-lambda) under a Poisson process of rate lambda -- i.e. "the
# site was used at least once". Camera information about lambda scales as
# dpsi/dlambda = 1 - psi, so as psi -> 1 the occupancy link goes FLAT and stops
# discriminating abundance. That is the saturation problem, and aggregating
# cameras into arrays makes it worse for abundant species.
#
# Royle-Nichols observes the SAME latent quantity less lossily:
#     N_a ~ Poisson(lambda_a)                  local abundance
#     per-camera detection = 1 - (1 - r)^{N_a}  more animals -> more detections
# The expected proportion of detecting cameras keeps changing as N_a grows, so
# it does not saturate the way psi does. RN is defined on
# detection/non-detection data -- which is exactly what this pipeline retains
# -- so no counts are required.
#
# Classic N-mixture was considered and is NOT implementable here: it needs
# counts per visit, and camera data is binarized to 0/1 per 10-day window
# before the model sees it. The usable content of that idea is this model: "k
# of n_a cameras detected" is a binomial count with abundance driving
# detection, which IS Royle-Nichols.
#
# WHAT DIFFERS FROM model_code_national_scalar.R
# ONLY the camera observation block (baseline lines 70-93). The iNat
# likelihood, calcIntensity_SVC, the CAR fields, the climate SVCs, occ_beta,
# and every trend parameter (year_beta, year_var, total_var_beta, snr,
# trend_robust_indicator) are byte-identical to the baseline.
#
# THE KEY MODELLING DECISION -- lambda_a on the same scale as the baseline
# The baseline's linear predictor IS log(lambda) (that is what
# cloglog(psi) = log(lambda) means). So here the SAME linear predictor, with
# the SAME parameters, is used directly as log(lambda_a) with no
# reparameterization:
#     log(lambda_a) = link_occ_intercept[cell] + MWMT*MWMT_effect + ...
#                     + year_beta * year_occ
# This is deliberate and is what makes the arms comparable: the trend
# parameters act on the identical quantity in all three models, differing only
# in how cameras observe it.
#
# HETEROGENEITY -- RN's DOCUMENTED WEAKNESS, carry into the writeup
# RN assumes cameras within an array share detection probability r given N_a.
# Real arrays have habitat variation, and unmodelled detection heterogeneity is
# RN's known failure mode: it inflates abundance estimates. Per-camera
# covariates on r (canopy, road distance, effort) mitigate but do not eliminate
# this. Any inference from this arm should be read with that caveat.
#
# IMPLEMENTATION NOTE
# N_a is a latent discrete parameter, so NIMBLE will assign it a slice or RW
# sampler and mixing may be slower than the occupancy arms. N_max caps the
# support for tractability; check the posterior of N_a is not piling up at
# N_max (which would mean the cap is binding and should be raised).
# =============================================================================

model_code_array_rn <- nimbleCode({

  # ---- iNat intensity + likelihood: UNCHANGED from baseline ----------------
  for (t in 1:nyear) {
    for (k in 1:n_cells_year[t]) {
      mu[inat_cells_by_year[k, t], t] <- calcIntensity_SVC(
        intensity_intercept = link_occ_intercept[inat_cell100[inat_cells_by_year[k, t]]],
        xdat         = xdat_inat[inat_cell50_start[inat_cells_by_year[k, t]]:inat_cell50_end[inat_cells_by_year[k, t]], 1:numOccCovars, t],
        beta         = occ_beta[1:numOccCovars],
        MWMT_dat     = MWMT_inat[inat_cell50_start[inat_cells_by_year[k, t]]:inat_cell50_end[inat_cells_by_year[k, t]], t],
        MWMT_effect  = MWMT_effect[inat_cell100[inat_cells_by_year[k, t]]],
        MCMT_dat     = MCMT_inat[inat_cell50_start[inat_cells_by_year[k, t]]:inat_cell50_end[inat_cells_by_year[k, t]], t],
        MCMT_effect  = MCMT_effect[inat_cell100[inat_cells_by_year[k, t]]],
        total_var_beta = total_var_beta,
        year_dat     = year_vals[t],
        theta0       = theta0,
        theta1       = theta1
      )
      y_inat[inat_cells_by_year[k, t], t] ~ dnbinom(
        size = 1 / overdisp_inat,
        prob = 1 / (1 + overdisp_inat * inat_effort[inat_cells_by_year[k, t], t] *
                      mu[inat_cells_by_year[k, t], t])
      )
    }
  }

  # ---- camera observation block: RN, THE ONLY CHANGE -----------------------
  for (a in 1:narray) {

    # log(lambda_a): the SAME linear predictor the baseline uses, unchanged.
    # In the baseline this quantity equals log(lambda) via
    # cloglog(psi) = log(lambda); here it drives abundance directly.
    log(lambda_a[a]) <- link_occ_intercept[array_cell[a]] +
      array_MWMT[a] * MWMT_effect[array_cell[a]] +
      array_MCMT[a] * MCMT_effect[array_cell[a]] +
      inprod(occ_beta[1:numOccCovars], array_occ_covars[a, 1:numOccCovars]) +
      year_beta * array_year_occ[a]

    # latent local abundance
    N_a[a] ~ dpois(lambda_a[a])

    # per-camera detection: rises with N_a -- this is the non-saturating part
    for (i in 1:n_a[a]) {
      cloglog(r[a, i]) <- link_det_intercept +
        q_beta[1] * array_yday[a, i] +
        q_beta[2] * array_yday[a, i]^2 +
        q_beta[3] * array_canopy[a, i] +
        q_beta[4] * array_roaddist[a, i] +
        q_beta[5] * log_effort[a, i]

      # P(camera i detects | N_a animals present)
      p_cam[a, i] <- 1 - pow(1 - r[a, i], N_a[a])
      w[a, i] ~ dbern(p_cam[a, i])
    }
  }

  # ---- priors --------------------------------------------------------------
  for (i in 1:5) {
    q_beta[i] ~ dnorm(0, sd = 5)
  }
  det_intercept ~ dunif(0, 1)
  link_det_intercept <- cloglog(det_intercept)

  # ---- occupancy / intensity priors: UNCHANGED from baseline ---------------
  for (i in 1:numOccCovars) {
    occ_beta[i] ~ dnorm(0, sd = 5)
  }

  link_occ_intercept[1:ncell100] ~ dcar_normal(adj[1:nadj], weights[1:nadj],
                                                num[1:nnum], intercept_tau,
                                                zero_mean = 0)
  intercept_tau ~ dgamma(1, 1)

  MWMT_effect[1:ncell100] ~ dcar_normal(adj[1:nadj], weights[1:nadj],
                                         num[1:nnum], MWMT_tau, zero_mean = 0)
  MCMT_effect[1:ncell100] ~ dcar_normal(adj[1:nadj], weights[1:nadj],
                                         num[1:nnum], MCMT_tau, zero_mean = 0)
  MWMT_tau ~ dgamma(1, 1)
  MCMT_tau ~ dgamma(1, 1)

  theta0 ~ dnorm(0, sd = 10)
  theta1 ~ dnorm(1, sd = 1)
  overdisp_inat ~ dgamma(1, 5)

  # ---- trend block: BYTE-IDENTICAL to baseline ----------------------------
  year_beta ~ dnorm(0, sd = sigma_year_beta)
  year_var ~ dnorm(0, sd = sigma_year_var)
  total_var_beta <- year_beta + year_var

  sigma_year_beta ~ dunif(0, 2)
  sigma_year_var ~ dunif(0, 2)

  snr <- year_beta / year_var

  # Probabilistic “robustness” indicator based on snr
  trend_robust_indicator <- step(snr - 1) # 1 if >1, 0 if <=1
})
