# =============================================================================
# model_code_array_occ.R
# -----------------------------------------------------------------------------
# Array-level occupancy variant of the integrated SDM, for the estimator
# comparison simulation.
#
# WHAT DIFFERS FROM model_code_national_scalar.R
# ONLY the camera observation block. Specifically, baseline lines 70-93 (the
# `for (i in 1:nsite)` loop over camera sites and 10-day windows) are replaced
# by a loop over ARRAYS with CAMERAS as the replicate dimension. Everything
# else -- the iNat likelihood, calcIntensity_SVC, link_occ_intercept's CAR
# field, MWMT_effect / MCMT_effect, occ_beta, and every trend parameter
# (year_beta, year_var, total_var_beta, snr, trend_robust_indicator) -- is
# byte-identical to the baseline. That invariant is what lets a difference in
# trend recovery be attributed to the estimator rather than to anything else.
#
# THE CHANGE, PRECISELY
#   baseline:  y[i, 1:J[i]] ~ dOcc_v(probOcc = psi[i],
#                                    probDetect = p[i, 1:J[i]], len = J[i])
#              replicates = 10-day windows within a camera
#
#   here:      w[a, 1:n_a[a]] ~ dOcc_v(probOcc = psi_a[a],
#                                      probDetect = q[a, 1:n_a[a]], len = n_a[a])
#              replicates = cameras within an array-year
#
# EFFORT -- THE CRITICAL ADDITION
# At camera level, effort enters AUTOMATICALLY: J[i] is the window count, so a
# long deployment contributes more Bernoulli trials and dOcc_v's `len`
# accounts for it exactly. Aggregating to arrays collapses the window
# dimension, so effort must re-enter as a per-camera detection covariate.
# q_beta[5] * log_effort[a, i] is that term. Measured J spread is 1-10
# occasions (median 3), a 10x range, so omitting it would confound deployment
# length with occupancy.
#
# INTERPRETATION SHIFT -- state this in any writeup
# psi_a is the probability the ARRAY was used, not that each camera's site was
# used. q absorbs within-array habitat patchiness in addition to detection.
# Occupancy covariates are array means (see build_array_constants); within-array
# covariate variation is therefore not represented on the occupancy side.
#
# yday is carried as the array's per-camera mean and its square (the PI-agreed
# choice); within-array seasonal spread is lost.
# =============================================================================

model_code_array_occ <- nimbleCode({

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

  # ---- camera observation block: THE ONLY CHANGE ---------------------------
  for (a in 1:narray) {

    # per-camera detection within the array; cameras are the replicates
    for (i in 1:n_a[a]) {
      cloglog(q[a, i]) <- link_det_intercept +
        q_beta[1] * array_yday[a, i] +
        q_beta[2] * array_yday[a, i]^2 +
        q_beta[3] * array_canopy[a, i] +
        q_beta[4] * array_roaddist[a, i] +
        q_beta[5] * log_effort[a, i]      # effort, no longer automatic
    }

    # array-level occupancy, on the SAME latent intensity scale as baseline
    cloglog(psi_a[a]) <- link_occ_intercept[array_cell[a]] +
      array_MWMT[a] * MWMT_effect[array_cell[a]] +
      array_MCMT[a] * MCMT_effect[array_cell[a]] +
      inprod(occ_beta[1:numOccCovars], array_occ_covars[a, 1:numOccCovars]) +
      year_beta * array_year_occ[a]

    w[a, 1:n_a[a]] ~ dOcc_v(probOcc = psi_a[a],
                            probDetect = q[a, 1:n_a[a]],
                            len = n_a[a])
  }

  # ---- priors --------------------------------------------------------------
  # detection: baseline's 4 p_beta plus the effort term, same prior
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
