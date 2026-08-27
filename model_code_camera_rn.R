# =============================================================================
# model_code_camera_rn.R
# -----------------------------------------------------------------------------
# Royle-Nichols at the CAMERA level -- the fourth arm, completing the 2x2:
#
#                    occupancy            Royle-Nichols
#   camera level     national_scalar      THIS FILE
#   array level      array_occ            array_rn
#
# With this cell filled, the estimator effect (occupancy vs RN) and the
# aggregation effect (camera vs array) become separately estimable instead of
# being confounded.
#
# STRUCTURE
#   N_i ~ Poisson(lambda_i)                    latent local abundance
#   per-window detection = 1 - (1 - r)^{N_i}   more animals -> more detections
# The replicate dimension is the 10-day window, as in the baseline. Effort is
# therefore handled automatically by the number of windows J[i]; no log_effort
# covariate is added (the array models need one only because aggregation
# collapses the window dimension).
#
# WHY THIS ARM IS EXPECTED TO HELP WHERE OCCUPANCY FAILS
# Fisher information about log(lambda), J=4 windows, r=0.1:
#   lambda:        1      3      5     10     20     50
#   occupancy:  2.33   1.89   0.68  0.018   ~0     0
#   RN:         0.38   1.02   1.52   2.27   2.39  0.63
# Occupancy's information dies by lambda~3-10 because psi saturates at 1. RN
# reads detection FREQUENCY, so its information keeps rising to lambda~20
# before its own (much later) saturation. Camera-level occupancy is therefore
# weakest exactly where camera-level RN is strongest -- the abundant-species
# case (white-tailed deer: 54.9% of detecting cameras fire in every window,
# where occupancy is uninformative but RN is not).
#
# TRADE-OFF vs ARRAY-LEVEL RN
# Camera level offers a median of 4 window-replicates per unit; array level
# offers a median of 4 cameras (mean 8.1) per unit. Camera level keeps far more
# units and better spatial resolution; array level has more replicates per unit
# and pools detection heterogeneity differently.
#
# HETEROGENEITY -- RN's DOCUMENTED WEAKNESS
# RN assumes windows within a camera share r given N_i. Unmodelled detection
# heterogeneity is RN's known failure mode: it inflates abundance. In
# simulation this is absent by construction, so this arm's performance here is
# OPTIMISTIC relative to real data. Carry this into any writeup.
#
# IMPLEMENTATION NOTE
# N_cam is a latent discrete parameter at nsite (~700) nodes -- more latent
# discrete nodes than array_rn (~97). Expect slower mixing; check that the
# N_cam posterior is not piling up against any cap.
# =============================================================================

model_code_camera_rn <- nimbleCode({


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

  # ---- camera observation block: RN at CAMERA level, THE ONLY CHANGE -------
  # Replicate dimension is the 10-day WINDOW within a camera (j in 1:J[i]),
  # exactly as in the baseline. The baseline's dOcc_v is replaced by a latent
  # local abundance N_i plus abundance-driven per-window detection.
  #
  # Effort enters AUTOMATICALLY here, as in the baseline: a longer deployment
  # has more windows, hence more Bernoulli trials. No log_effort term is added
  # (that term exists in the ARRAY models only because aggregating to arrays
  # collapses the window dimension).
  for (i in 1:nsite) {

    # log(lambda_i): the SAME linear predictor the baseline uses for
    # cloglog(psi_i). Identical parameters, so arms remain comparable.
    log(lambda_cam[i]) <- link_occ_intercept[cell[i]] +
      MWMT[i] * MWMT_effect[cell[i]] +
      MCMT[i] * MCMT_effect[cell[i]] +
      inprod(occ_beta[1:numOccCovars], occ_covars[i, 1:numOccCovars]) +
      year_beta * year_occ[i]

    # latent local abundance in the camera's viewshed
    N_cam[i] ~ dpois(lambda_cam[i])

    for (j in 1:J[i]) {
      # per-window, per-individual detection probability -- covariates are the
      # baseline's p_beta set, unchanged
      cloglog(r[i, j]) <- link_det_intercept +
        p_beta[1] * yday[i, j] +
        p_beta[2] * yday[i, j]^2 +
        p_beta[3] * canopy_height[i] +
        p_beta[4] * log_roaddist[i]

      # P(detect in window j | N_cam animals present)
      p_cam[i, j] <- 1 - pow(1 - r[i, j], N_cam[i])
      y[i, j] ~ dbern(p_cam[i, j])
    }
  }

  # ---- priors --------------------------------------------------------------
  for (i in 1:4) {
    p_beta[i] ~ dnorm(0, sd = 5)
  }
  det_intercept ~ dunif(0, 1)
  link_det_intercept <- cloglog(det_intercept)

  # ---- occupancy / intensity priors: UNCHANGED from baseline ---------------
  for (i in 1:numOccCovars) {
    occ_beta[i] ~ dnorm(0, sd = 5)
  }

  # CAR fields: NO `weights` argument, matching the BASELINE camera model.
  #
  # This is load-bearing and was a real bug when this file was first drafted
  # from model_code_array_rn.R. The array models pass weights[1:nadj]
  # explicitly, and that works there only because build_array_constants()
  # creates `weights` (a vector of 1s). Camera-level arms never call
  # build_array_constants(), and build_reduced_constants() does NOT supply
  # `weights` -- so an array-style call here leaves the node undefined and
  # buildMCMC() dies inside CAR_normal_processParams() with "missing value
  # where TRUE/FALSE needed", after the model has already built and samplers
  # have been configured. Caught by running a replicate end to end; a parse
  # check and a nimbleModel() build check both pass without detecting it.
  #
  # Omitting weights is equivalent: NIMBLE defaults dcar_normal to unit
  # weights, which is exactly what the array path constructs.
  # PRIORS AND zero_mean MATCH THE CAMERA BASELINE EXACTLY -- dgamma(0.1, 0.1)
  # and no zero_mean argument.
  #
  # This arm's whole purpose is the contrast camera_occ vs camera_rn, which is
  # interpretable as an ESTIMATOR effect only if nothing else differs. An
  # earlier draft of this file copied array_rn's dgamma(1, 1) + zero_mean = 0,
  # which would have confounded the estimator change with a CAR hyperprior
  # change. The two priors are far apart in practice, not a cosmetic
  # difference: median tau is 0.006 under dgamma(0.1, 0.1) versus 0.69 under
  # dgamma(1, 1), and P(tau < 0.1) is 0.66 versus 0.10.
  #
  # The direction of the fix is forced. model_code_national_scalar.R was
  # verified token-identical to Arielle's original (only the object name
  # differs), and dgamma(0.1, 0.1) is HER prior -- so the baseline cannot be
  # edited to match the array models, and the array models are the ones that
  # departed. Matching the baseline here yields:
  #     camera_occ vs camera_rn : estimator effect only          (clean)
  #     array_occ  vs array_rn  : estimator effect only          (clean)
  #     camera     vs array     : aggregation + prior            (confounded)
  # That confound is PRE-EXISTING -- it is already baked into the 540
  # completed rows -- and this choice confines it to the aggregation contrast
  # instead of also contaminating the camera-level estimator contrast.
  link_occ_intercept[1:ncell100] ~ dcar_normal(adj = adj[1:nadj],
                                               num = num[1:nnum],
                                               tau = intercept_tau)
  intercept_tau ~ dgamma(0.1, 0.1)

  MWMT_effect[1:ncell100] ~ dcar_normal(adj = adj[1:nadj], num = num[1:nnum],
                                        tau = MWMT_tau)
  MCMT_effect[1:ncell100] ~ dcar_normal(adj = adj[1:nadj], num = num[1:nnum],
                                        tau = MCMT_tau)
  MWMT_tau ~ dgamma(0.1, 0.1)
  MCMT_tau ~ dgamma(0.1, 0.1)

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
