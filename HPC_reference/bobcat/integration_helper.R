# =============================================================================
# REFERENCE COPY -- not part of the runnable local pipeline. Added 2026-09-01.
# -----------------------------------------------------------------------------
# This is a Hazel-side forked copy of Arielle's `integration_helper.R` as it
# existed for the bobcat run, recovered from the untracked `isdm_sim_codebase/`
# snapshot bundle. It is committed for the historical record only.
#
# WHY IT MATTERS: the Fix 1 iNat effort/count sort-order fix (`names_sort=TRUE`)
# was originally applied BY HAND directly on Hazel, in per-species forked copies
# of this file, for the v1fix jobs (499944-499949). That left the actual code
# those fits ran with no local record -- a documented gap. This copy closes part
# of that gap: it is one of the real cluster-side files, not a reconstruction.
#
# Do NOT source this from local scripts. The local, versioned equivalent is
# `integration_helper_fix1.R` in the repo root. Diff the two to confirm whether
# the hand-applied cluster fix and the committed fix actually match.
# =============================================================================

# =============================================================================
# integration_helper.R  --  MINIMAL EXTRACT for the reproducible bundle
# -----------------------------------------------------------------------------
# Contains ONLY the two nimbleFunctions the model code calls:
# calcIntensity_SVC() and calcIntensity_noSVC(). Both are reproduced verbatim
# from Arielle's original
#   codebase/Run_models_Hazel/HPC_scripts/integration_helper.R
# and were verified byte-identical at the token level (whitespace/CRLF
# normalized) before extraction:
#   calcIntensity_SVC    md5 8cd666696ee435d8dddd910fe2fb0a70
#   calcIntensity_noSVC  md5 acfda1097358bf9080a4329118b2f415
#
# The full original additionally contains data-preparation helpers and
# hardcoded cluster paths (PROJ_DIR, DATA_DIR, GRID50_PATH pointing into
# /rsstu/...). Those are omitted here because they are not used by the
# simulation code and would fail on any machine without that filesystem. No
# model logic is modified or omitted.
# =============================================================================

calcIntensity_SVC <- nimbleFunction(run = function(
    intensity_intercept = double(0),
    theta0 = double(0),
    theta1 = double(0),
    MWMT_effect = double(0),
    MCMT_effect = double(0),
    beta   = double(1),
    xdat   = double(2),
    MWMT_dat = double(1),
    MCMT_dat = double(1),
    year_dat = double(0),
    total_var_beta = double(0)) {
  
  log_lambda <- intensity_intercept + 
    (xdat %*% matrix(beta, ncol = 1)) + 
    (MWMT_dat * MWMT_effect) + 
    (MCMT_dat * MCMT_effect) +
    (total_var_beta * year_dat)
  
  mu <- exp(theta0 + theta1 * log(sum(exp(log_lambda))))
  
  return(mu)
  returnType(double(0))
})

calcIntensity_noSVC <- nimbleFunction(run = function(
    intensity_intercept = double(0),
    theta0 = double(0),
    theta1 = double(0),
    beta   = double(1),
    xdat   = double(2),
    year_dat = double(0),
    total_var_beta = double(0)) {
  
  log_lambda <- intensity_intercept + 
    (xdat %*% matrix(beta, ncol = 1)) +
    (total_var_beta * year_dat)
  
  mu <- exp(theta0 + theta1 * log(sum(exp(log_lambda))))
  
  return(mu)
  returnType(double(0))
})
