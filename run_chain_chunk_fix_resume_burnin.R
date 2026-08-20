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
# NIMBLE SEMANTICS -- VERIFIED EMPIRICALLY, NOT ASSUMED.
# Probed on the Hazel build (NIMBLE 1.4.2, R 4.5.3) in jobs 603848/603879.
# run()'s formals are:
#   niter, reset, resetMV, time, progressBar, nburnin, thin, thin2,
#   resetWAIC, initializeModel, chain
# Findings that determine the correct calls:
#   * reset = TRUE WIPES every adaptive sampler's proposal scale back to its
#     default (measured: scale 1 -> 0.156961 after run(3000) -> back to 1
#     after a subsequent run(reset = TRUE)). reset = FALSE preserves it.
#   * reset = TRUE does NOT re-initialize model values; it resets sampler
#     adaptation and clears the sample buffer only.
#   * reset = TRUE clears mvSamples regardless of resetMV. resetMV is only
#     honoured when reset = FALSE.
#   * Assigning Cmcmc$mvSaved from a checkpoint DOES carry the chain over --
#     a fresh build inited at mu = 0 resumed at the donor's mu = 50. The
#     checkpoint/restore mechanism itself is sound.
#   * nburnin CANNOT be combined with reset = FALSE -- NIMBLE raises
#     "cannot specify nburnin when using reset = FALSE", so the tempting
#     single-call form run(b + n, reset = FALSE, nburnin = b) is illegal.
#
# CORRECTION TO THE FIRST DRAFT OF THIS FILE: it left `reset` at its TRUE
# default on the chunk call. Given the first finding above, that call would
# have thrown away exactly the adaptation the new resume burn-in buys, and
# reproduced the transient unchanged -- the fix would have silently done
# nothing. `reset = FALSE` on BOTH calls is the load-bearing part.
#
# FIX: give every resumed chunk its own re-adaptation burn-in, then run the
# recorded chunk WITHOUT resetting adaptation:
#   1. After restore_chain_state() sets Cmcmc$mvSaved from the checkpoint,
#      run resume_burnin iterations with reset = FALSE (keep the restored
#      state, let the samplers adapt) and resetMV = TRUE (discard the draws).
#   2. Then run chunk_iter with reset = FALSE (PRESERVE that adaptation) and
#      resetMV = TRUE (so exactly chunk_iter draws are recorded).
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
    # reset = FALSE  -> keep the restored mvSaved state as the starting point,
    #                   and let the samplers adapt (adaptation is NOT frozen by
    #                   reset = FALSE -- verified: scale kept its adapted value
    #                   0.156961 across a reset = FALSE call, and was reset to
    #                   the default 1 by a reset = TRUE one).
    # resetMV = TRUE -> clear mvSamples, so these resume_burnin draws never
    #                   reach `new_samples` below. Honoured because reset is
    #                   FALSE; with reset = TRUE it would be ignored.

  } else {

    cat("Starting new chain...\n")

    set.seed(1000 + chain_id)

    # Burn-in once only
    Cmcmc$run(burnin_once)

    samples_so_far <- NULL
    iter_total     <- 0
  }

  # Run the recorded chunk. Applies to BOTH branches: the fresh-start branch
  # has the same defect in the original code -- its Cmcmc$run(chunk_iter) with
  # reset = TRUE discards the adaptation that Cmcmc$run(burnin_once) just
  # bought, so round 1 also begins at the default proposal scale.
  cat("Running", chunk_iter, "iterations...\n")

  Cmcmc$run(chunk_iter, reset = FALSE, resetMV = TRUE)
  # reset = FALSE is the whole point: it PRESERVES the proposal scales adapted
  # by the branch above (resume_burnin on a resume, burnin_once on a fresh
  # start). Leaving reset at its TRUE default -- as the first draft did -- would
  # wipe them back to the build-time defaults and reproduce the transient.
  # resetMV = TRUE clears the buffer so exactly chunk_iter draws are recorded;
  # verified row counts: run(100); run(50, reset=FALSE) -> 150 (appends), then
  # run(30, reset=FALSE, resetMV=TRUE) -> 30 (clears, then records).

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
