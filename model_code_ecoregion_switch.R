#!/usr/bin/env Rscript
# =============================================================================
# model_code_ecoregion_switch.R
# Defines model_code_ecoregion_switch. FORKED from model_code_ecoregion_trend.R
# (which defines model_code_ecoregion). Replaces the WAIC comparison gate --
# shown by the corrected abundance sweep to be unusable, null vs varying WAIC
# distributions overlap almost completely, p = 0.80/0.91/0.61 at every
# abundance level, see abundance_sweep_evaluation.md -- with the model-
# comparison approach from Goldstein et al. (bioRxiv 2025.01.17.633640, same
# lab lineage, same camera+iNat iSDM, same NIMBLE stack): a Bernoulli "switch"
# indicator on the spatially-varying component, fit with reversible-jump MCMC
# (RJMCMC), read via posterior inclusion probability P(gamma=1). See
# rjmcmc_indicator_test_plan.md.
#
# SINGLE INDICATOR, MINIMAL SCOPE (explicit instruction): one gamma gates the
# ENTIRE ecoregion trend deviation as one group, not one gamma per region.
#   gamma = 1 -> year_region[r] ~ dnorm(0, slab_sd)   (spatially-varying trend)
#   gamma = 0 -> year_region[r] = 0 for all r          (national scalar trend only)
# Implemented as the standard NIMBLE indicator-selection "zeroing" pattern:
# year_region_raw[r] is the always-present continuous node RJMCMC adds/removes
# from the model; year_region[r] <- gamma * year_region_raw[r] is what the
# rest of the model actually sees, so the calcIntensity_SVC call site below
# is BYTE-IDENTICAL to model_code_ecoregion (same as that model was
# byte-identical to the CAR model at its call site) -- only how year_effect
# is defined upstream changes, never how it's consumed downstream.
#
# ONE DOCUMENTED, DELIBERATE DEVIATION from "everything else stays
# byte-identical" -- flagged rather than silently made, per empirical testing
# (rj_group_test.R / rj_group_test2.R on Hazel, nimble 1.4.2):
#   sigma_region is NO LONGER an estimated node in this model. NIMBLE's
#   configureRJ() explicitly refuses target nodes whose prior has a
#   non-constant (i.e. model-estimated) hyperparameter -- confirmed directly:
#   configureRJ() on year_region_raw[r] ~ dnorm(0, sd = sigma_region) with
#   sigma_region ~ dexp(1) errors with "Reversible jump target node ...
#   appears to have a non-constant hyper-parameter, which is not currently
#   supported for reversible jump MCMC." The identical model with a FIXED sd
#   instead of an estimated one lets configureRJ() run cleanly (verified,
#   including one shared indicator across all 8 region nodes -- that part
#   works fine and was NOT the blocker). This is not a NIMBLE-specific
#   workaround: RJMCMC birth/death moves need a fixed, known proposal
#   distribution to compute a tractable acceptance ratio, so slab variances
#   are conventionally fixed constants (not jointly estimated) in most
#   implementations of this method. slab_sd is supplied as a CONSTANT (not
#   hardcoded here -- see build_reduced_constants(..., slab_sd=) /
#   01f_run_indicator_test.R), default 0.227 (the mean posterior
#   sigma_region_mean from the corrected abundance sweep's varying scenario,
#   close to the 0.2 truth-generation target in make_true_year_region()).
#   pi_gamma (prior inclusion probability, default 0.5) is likewise a
#   constant, not hardcoded into dbern() here.
#
# UNCHANGED, byte-for-byte vs model_code_ecoregion: everything else --
# climate SVCs (MWMT_effect, MCMT_effect), the calcIntensity_SVC call site,
# both likelihoods (iNat dnbinom, camera dOcc_v), year_beta/year_var/
# total_var_beta, the camera occupancy trend term, the robustness indicator.
#
# Requires (in constants), beyond what model_code_ecoregion needs:
#   pi_gamma  (scalar, prior inclusion probability for gamma)
#   slab_sd   (scalar, FIXED sd for year_region_raw's prior -- replaces the
#             estimated sigma_region node)
# Still requires: ecoregion_of_cell100 (length ncell100), nregion.
# =============================================================================

model_code_ecoregion_switch <- nimbleCode({

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
          # UNCHANGED vs model_code_ecoregion's call site -- year_effect is
          # now gated by gamma upstream, but this line itself doesn't change.
          total_var_beta = total_var_beta + year_effect[inat_cell100[g]]
        )

      } else {

        # unchanged: no SVC infrastructure for non-SVC species -- out of scope
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

      # ---- NEW: deterministic re-expression of the SAME dnbinom density in
      # log form, for the data-log-likelihood convergence diagnostic (see
      # rjmcmc_indicator_test_plan.md gotcha #1 -- Gelman-Rubin is meaningless
      # on a Bernoulli indicator that sits at 0/1 for long stretches even
      # when mixing well; instead monitor total_data_logLik below and compute
      # R-hat/ESS on THAT). Pure bookkeeping: recomputes the density of the
      # already-declared y_inat node at its current value, has zero influence
      # on sampling of any other node.
      logLik_y_inat[k, t] <- dnbinom(y_inat[inat_cells_by_year[k,t],t],
                                     size = 1 / overdisp_inat,
                                     prob = 1 / (1 + overdisp_inat * inat_effort[inat_cells_by_year[k,t],t] *
                                                   mu[inat_cells_by_year[k,t],t]),
                                     log = 1)

    }
    logLik_y_inat_year[t] <- sum(logLik_y_inat[1:n_cells_year[t], t])
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
    # only; the ecoregion trend deviation is scoped to the iNat pathway only
    # in this prototype, same scope decision as model_code_ecoregion)
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

    # ---- NEW: same bookkeeping as logLik_y_inat above, for the camera
    # detection-history likelihood.
    logLik_y[i] <- dOcc_v(y[i, 1:J[i]], probOcc = psi[i], probDetect = p[i, 1:J[i]],
                          len = J[i], log = 1)
  }

  # ---- NEW: single scalar, the total log-likelihood of ALL observed data
  # (y and y_inat) at the current parameter draw -- this is the quantity
  # monitored for convergence (R-hat + ESS), NOT gamma. See fit_replicate_switch()
  # in sim_helpers.R.
  total_data_logLik <- sum(logLik_y[1:nsite]) + sum(logLik_y_inat_year[1:nyear])

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

    # ---- NEW: single Bernoulli switch on the whole ecoregion trend
    # deviation, replacing model_code_ecoregion's always-on partial-pooled
    # year_region. gamma is the ONE indicator (single-indicator scope).
    # year_region_raw[r] is the node RJMCMC adds/removes from the model (its
    # prior sd is the FIXED slab_sd constant, not an estimated hyperparameter
    # -- see file header). year_region[r] is a deterministic gate: identical
    # to year_region_raw when gamma=1, forced to exactly 0 when gamma=0.
    # year_effect[c] pass-through is otherwise UNCHANGED from model_code_ecoregion.
    gamma ~ dbern(pi_gamma)

    for (r in 1:nregion) {
      year_region_raw[r] ~ dnorm(0, sd = slab_sd)
      year_region[r] <- gamma * year_region_raw[r]
    }

    for (c in 1:ncell100) {
      year_effect[c] <- year_region[ecoregion_of_cell100[c]]
    }

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
