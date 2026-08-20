# =============================================================================
# run_chain_chunk_fix_resume_burnin.R
# -----------------------------------------------------------------------------
# Fixes a fleet-wide sampler-restart bug in run_chain_chunk() (present,
# byte-identical in structure, in HPC_wtd_national_scalar_chain1.R,
# HPC_wtd_ecoregion_chain1.R, and the equivalent moose/bobcat chain scripts).
#
# BUG (found 2026-08-20, diagnosed via block-mean traces on real chain RDS):
# save_chain_state() persists only Cmcmc$mvSaved (current node VALUES) --
# it never saves the adaptive samplers' internal proposal-scale/covariance
# state. Every resumed chunk rebuilds nimbleModel -> configureMCMC ->
# buildMCMC -> compileNimble from scratch, so the NEW Cmcmc's samplers start
# at their default proposal scale regardless of what was learned in prior
# rounds. The ORIGINAL code's burn-in call (`Cmcmc$run(burnin_once)`) sits
# only in the fresh-start (`else`) branch -- a resumed chunk skips straight
# to `Cmcmc$run(chunk_iter)` with no re-adaptation burn-in at all, so the
# ~1000 draws where each sampler is quietly re-learning its proposal scale
# get recorded as if they were converged posterior draws. Confirmed
# reproducible across all three chains of multiple species/models at every
# round boundary (block 11 of each resumed chunk, ~7% of that chunk's
# draws).
#
# IMPACT: point estimates are essentially unaffected (confirmed on moose
# v2b, both models) -- the transient washes out under a mean/median. Credible
# intervals were inflated on one tail by ~30%, and R-hat was distorted
# (direction depends on whether the transient hit all chains transiently,
# which deflates R-hat and makes convergence look BETTER than it is, or
# permanently displaced one chain, which inflates it).
#
# FIX: give every resumed chunk its own re-adaptation burn-in, run under
# the SAME reset semantics used for a fresh chain's `Cmcmc$run(burnin_once)`
# call (full re-adaptation of the samplers), but WITHOUT wiping the
# restored mvSaved state or double-counting these draws in the recorded
# samples matrix. NIMBLE's exact reset/resetMV default behavior differs by
# version -- confirm the two calls below against `?runMCMC` and
# `?buildMCMC`'s "Reset" section for whatever NIMBLE version is loaded on
# Hazel BEFORE trusting this in the fleet-wide rollout. The intent is:
#
#   1. After restore_chain_state() sets Cmcmc$mvSaved from the checkpoint,
#      run burnin_once iterations that DO adapt the samplers but are
#      DISCARDED (never appended to the samples matrix returned to the
#      caller) -- do NOT reset the model to its initial values (must keep
#      the restored state), and do NOT clear the mvSamples history if that
#      history is being relied on elsewhere.
#   2. Only THEN run chunk_iter iterations and record those as the "real"
#      chunk, exactly as the fresh-start branch already does.
#
# VALIDATION PLAN (per 2026-08-20 decision -- do this before fleet rollout):
#   Apply this fix to ONE resumed chunk on an isolated/throwaway chain file
#   first. Re-run the same block-mean diagnostic used to find the bug
#   (per-1000-draw block means of theta0/theta1/overdisp_inat right after
#   the resume point) and confirm:
#     (a) the transient at the old block-11 boundary is gone -- block means
#         should be flat/stable from the first recorded block onward, and
#     (b) the restored mvSaved state is actually being used (i.e. the chain
#         resumes from where it left off, not from the model's initial
#         values) -- check the first post-resume draws are continuous with
#         the pre-resume chain's last draws, not a fresh start.
#   Only roll out fleet-wide once both (a) and (b) are confirmed empirically
#   on real output, not just via code review.
#
# Usage: source this file AFTER the model-specific chain script has defined
# `model_code_<...>`, `species`, `project_dir`, and BEFORE calling
# run_chain_chunk() -- it replaces the run_chain_chunk definition in-place
# via the same function name, so no call-site changes are needed.
# =============================================================================

run_chain_chunk <- function(chain_id,
                            input_data,
                            chunk_iter = 10000,
                            burnin_once = 5000,
                            resume_burnin = 1000) {

  chain_file <- file.path(project_dir,
                          paste0("chain_", species, "_", chain_id, ".RDS"))

  cat("Building model for chain", chain_id, "\n")

  model <- nimbleModel(model_code_national_scalar,  # NOTE: swap for the
                                                      # correct model_code_*
                                                      # object per script
                       constants = input_data$constants_list,
                       data = list(y = input_data$real_data$y,
                                   y_inat = input_data$inat_y,
                                   inat_effort = input_data$inat_effort),
                       inits = input_data$inits_list,
                       calculate = FALSE)

  Cmodel <- compileNimble(model)

  conf <- configureMCMC(model, enableWAIC = TRUE)

  conf$addMonitors("occ_beta",
                   "link_occ_intercept",
                   "year_beta",
                   "year_var",
                   "total_var_beta",
                   "MWMT_effect",
                   "MCMT_effect",
                   "trend_robust_indicator")

  mcmc  <- buildMCMC(conf)
  Cmcmc <- compileNimble(mcmc, project = model)

  if (file.exists(chain_file)) {

    cat("Restoring previous state...\n")

    obj <- restore_chain_state(chain_file, Cmcmc)

    samples_so_far <- obj$samples
    iter_total     <- obj$iter_total

    # FIX: re-adaptation burn-in for the resumed chunk, mirroring the
    # fresh-start branch's Cmcmc$run(burnin_once) below. VERIFY reset/resetMV
    # arguments against the loaded NIMBLE version's docs before trusting --
    # intent: adapt samplers from their DEFAULT (rebuilt) scale using the
    # RESTORED model/mvSaved state, discard these resume_burnin draws
    # entirely (never appended to samples_so_far or new_samples below).
    cat("Re-adaptation burn-in for resumed chunk (", resume_burnin,
        "iterations, discarded)...\n")
    Cmcmc$run(resume_burnin, reset = FALSE, resetMV = TRUE)
    # reset = FALSE  -> do NOT reinitialize model to initial values; keep
    #                   the restored mvSaved state as the starting point.
    # resetMV = TRUE -> DO clear mvSamples before this call, so these
    #                   resume_burnin draws are never mixed into
    #                   Cmcmc$mvSamples and never reach `new_samples` below.
    # CONFIRM: does reset=FALSE on THIS nimble version still allow the
    # adaptive samplers to re-adapt (i.e. adaptation is a property of the
    # sampler objects added fresh at buildMCMC() time this call, not of
    # `reset`), or does reset=FALSE also freeze adaptation? If the latter,
    # this call does nothing useful and needs a different argument.

  } else {

    cat("Starting new chain...\n")

    set.seed(1000 + chain_id)

    # Burn-in once only
    Cmcmc$run(burnin_once)

    samples_so_far <- NULL
    iter_total     <- 0
  }

  # Run next chunk -- resetMV = FALSE so this appends to (rather than wipes)
  # whatever mvSamples state exists after the branch above.
  cat("Running", chunk_iter, "iterations...\n")

  Cmcmc$run(chunk_iter, resetMV = (file.exists(chain_file)))
  # For a resumed chunk, mvSamples was just cleared by the burn-in call
  # above (resetMV = TRUE there), so this call should NOT also clear it --
  # CONFIRM this combination actually yields "only chunk_iter draws in
  # mvSamples after this call" empirically, per the validation plan above.

  new_samples <- as.matrix(Cmcmc$mvSamples)

  if (is.null(samples_so_far)) {
    combined_samples <- new_samples
  } else {
    combined_samples <- rbind(samples_so_far, new_samples)
  }

  iter_total <- iter_total + chunk_iter

  save_chain_state(chain_file,
                   Cmcmc,
                   combined_samples,
                   iter_total)

  cat("Chain", chain_id, "now at", iter_total, "iterations\n")
}
