#!/usr/bin/env Rscript
# =============================================================================
# model_code_ecoregion_trend.R
# Defines model_code_ecoregion. FORKED from model_code_spatial_trend.R (the
# CAR-field spatial-trend prototype), which is itself forked from Arielle's
# $PROJ/HPC/bobcat/HPC_run_model_chunks_chain1.R. Diff vs both recorded in
# model_code_diff_ecoregion.txt in this directory.
#
# WHY THIS FORK: the CAR-field prototype was simulation-validated and found
# that at the ~900-cell CAR grid's grain, the fine-scale field shrinks
# under-informed cells toward zero (NE bias ~-0.3 at low data density,
# roughly halving when informed-cell density doubled) -- the fine grain is
# too data-hungry to resolve reliably. This pivots to a COARSER, structured
# estimand: an EPA Level I ecoregion partial-pooled random effect (K=8
# regions over the bobcat range, see 00b_prep_ecoregion.R), trading spatial
# resolution for identifiability.
#
# CHANGE vs model_code_spatial_trend.R: replace the year_effect CAR field
# with a partial-pooled ecoregion deviation. year_region[r] ~ dnorm(0,
# sigma_region) is the actual ESTIMATED/monitored node (one deviation per
# ecoregion, partial pooling toward the shared national mean -- year_beta --
# rather than toward each other via spatial adjacency like the CAR version).
# year_effect[c] is now a DETERMINISTIC pass-through,
# year_effect[c] <- year_region[ecoregion_of_cell100[c]], so the
# calcIntensity_SVC call site is BYTE-IDENTICAL to the CAR model
# (`total_var_beta + year_effect[inat_cell100[g]]`) -- only the definition
# of year_effect changed (CAR-distributed vs. deterministically inherited
# from a coarser regional random effect), not how it's used downstream. This
# also means true_param_list_ecoregion() (sim_helpers.R) only ever needs to
# set year_region as an init (year_effect can't take one directly -- it's
# deterministic, nimble computes it from its parents via calculate()).
#
# UNCHANGED, byte-for-byte (all inherited from model_code_spatial_trend.R):
#   * year_beta stays the global mean trend; year_var stays the global
#     iNat-only deviation; the camera occupancy trend term
#     (`year_beta * year_occ[i]` in cloglog(psi[i])) is untouched.
#   * iNat dnbinom likelihood, camera detection submodel,
#     link_occ_intercept CAR field.
#   * The climate SVCs (MWMT_effect, MCMT_effect) are UNCHANGED and remain
#     in the model -- this is deliberate: the simulation must test whether
#     the ecoregion trend is identifiable ALONGSIDE the temperature
#     response, not in isolation (they could in principle be confounded if
#     ecoregions correlate with temperature gradients, which they plausibly
#     do -- that's exactly the kind of confound this simulation checks for).
#
# Requires (in constants): ecoregion_of_cell100 (length ncell100, integer
# 1:nregion lookup) and nregion -- see build_reduced_constants(...,
# ecoregion_of_cell100_full=, nregion=) in sim_helpers.R.
# =============================================================================

model_code_ecoregion <- nimbleCode({

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
          # UNCHANGED vs the CAR model's call site -- year_effect is now
          # deterministically inherited from year_region (see below) instead
          # of CAR-distributed, but this line itself doesn't change.
          total_var_beta = total_var_beta + year_effect[inat_cell100[g]]
        )

      } else {

        # unchanged: no SVC infrastructure (incl. no year_region/year_effect)
        # for non-SVC species in this prototype -- out of scope, see header note
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
    # only; the ecoregion trend deviation is scoped to the iNat pathway only
    # in this prototype, same scope decision as the CAR-field version)
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

    # ---- NEW: ecoregion-level trend deviation, replacing the CAR
    # year_effect field. year_region[r] is the actual estimated/monitored
    # partial-pooled deviation, one per ecoregion, shrunk toward 0 (not
    # toward each other via spatial adjacency, unlike the CAR version) with
    # strength governed by sigma_region. year_effect[c] is a deterministic
    # pass-through so every grid cell in the same ecoregion shares exactly
    # the same trend deviation, and the calcIntensity_SVC call site above
    # doesn't need to know anything changed.
    for (r in 1:nregion) {
      year_region[r] ~ dnorm(0, sd = sigma_region)
    }
    sigma_region ~ dexp(1)

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
