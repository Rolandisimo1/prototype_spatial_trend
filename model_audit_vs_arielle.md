# Model audit — prototype model code vs Arielle's original

Verifies that the simulation prototype did NOT alter Arielle's core model.
Ground truth = Arielle's `model_code <- nimbleCode({...})` from
`Run_models_Hazel/HPC_scripts/HPC_run_model_chunks_chain1.R` (in
codebase.zip) and `calcIntensity_SVC` from her `integration_helper.R`
(Temporal_trend_final.zip). Comparison method: strip comments + variable
name + whitespace, then token-level char diff and md5.

## Result: core model is UNCHANGED

**`model_code_national_scalar.R` == Arielle's original, token-identical.**
- Normalized token stream md5: 95518deb4ade70bdf7b8d7c227b324c2 (both).
- The only diffs are (a) the nimbleCode object name
  (`model_code` -> `model_code_national_scalar`), (b) the closing `})`,
  (c) whitespace/blank lines. Zero logic difference.
- This is the baseline the simulations fit as the "national scalar" model,
  and it IS Arielle's model. Nothing was silently changed.

**`calcIntensity_SVC` is NOT redefined anywhere in the prototype** — it is
inherited from Arielle's `integration_helper.R` on Hazel. The intensity
assembly (log_lambda, theta0/theta1 log-sum) is hers, untouched.

## The two extended models differ from the original in EXACTLY the added term

**`model_code_ecoregion_trend.R`** — two changes, both additive:
1. Intensity call site (1 line): `total_var_beta = total_var_beta` becomes
   `total_var_beta = total_var_beta + year_effect[inat_cell100[g]]`
   (iNat pathway only; camera occupancy line unchanged).
2. New prior block appended:
   ```
   for (r in 1:nregion) { year_region[r] ~ dnorm(0, sd = sigma_region) }
   sigma_region ~ dexp(1)
   for (c in 1:ncell100) { year_effect[c] <- year_region[ecoregion_of_cell100[c]] }
   ```
Everything else (camera occupancy/detection submodel, iNat dnbinom, all
Arielle priors, CAR fields on intercept/MWMT/MCMT, year_beta/year_var/
total_var_beta, snr robustness indicator) is byte-identical to the original.

**`model_code_ecoregion_switch.R`** — the ecoregion change above, but with
the RJMCMC reparameterization (verified necessary, see
rjmcmc_indicator_test_plan.md):
- `gamma ~ dbern(pi_gamma)`; `year_region_raw[r] ~ dnorm(0, sd = slab_sd)`
  (FIXED slab, not estimated sigma_region — required by configureRJ);
  `year_region[r] <- gamma * year_region_raw[r]`.
- PLUS deterministic log-likelihood tracking nodes (logLik_y, logLik_y_inat,
  total_data_logLik) added ONLY so convergence can be judged on the data
  log-prob instead of on gamma (the paper's required diagnostic). These are
  monitoring-only; they do not enter any likelihood or prior and cannot
  change the fit.
All Arielle priors and both likelihoods otherwise byte-identical.

## Bottom line
No change was made to Arielle's model as received. The national-scalar
baseline reproduces it exactly; the two extensions add only the ecoregion
trend term (and, for the switch model, a fixed-slab RJMCMC reparameterization
+ monitoring-only logLik nodes). The `+ year_effect` deviation is confined to
the iNat intensity pathway in both, matching the documented scope decision.
