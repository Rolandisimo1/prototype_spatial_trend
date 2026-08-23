# =============================================================================
# checkpoint_numeric_values.R
# -----------------------------------------------------------------------------
# Drop-in replacements for save_chain_state() / restore_chain_state() from
# HPC_run_model_chunks_chain{1,2,3}.R.
#
# WHY
# ---
# The original pair does:
#     save:    saveRDS(list(mcmc_state = Cmcmc$mvSaved, ...), file)
#     restore: Cmcmc$mvSaved <- obj$mcmc_state
#
# Cmcmc$mvSaved is a COMPILED CmodelValues object backed by an external
# pointer. saveRDS writes the R-level shell; the pointer is dead on read-back
# (touching it throws "Sextptr is not a valid external pointer"). The
# assignment therefore does not transfer any values -- and does not error --
# so the resumed chunk samples from the model's INITIAL VALUES.
#
# Audited 2026-08-21 across all 164 chain_*.RDS on Hazel (206 GB, 58 fit
# dirs): 1,234 / 1,234 resume boundaries restart from inits. Zero clean
# resumes. Confirmed by the post-resume draw landing exactly (<1e-9) on the
# hardcoded inits: det_intercept -> 0.5, overdisp_inat -> 0.1,
# intercept_tau -> 0.3.
#
# WHAT CHANGES
# ------------
# Persist PLAIN NUMERIC node values, not the compiled object. On resume,
# write them back into the freshly built model, then push model -> mvSaved
# through the compiled MCMC's own accessor so the sampler starts from the
# restored state.
#
# Two properties the original lacked:
#   1. The restore is VERIFIED, not assumed -- restore_chain_state_v2()
#      returns the restored values and the caller asserts continuity against
#      the last recorded draw. A silent no-op cannot pass.
#   2. It fails LOUDLY. A dead-pointer or missing-node condition raises
#      rather than quietly reverting to inits.
#
# WHAT THIS DOES **NOT** FIX
# --------------------------
# Adaptive sampler state (proposal scales / covariances) is still not carried
# across a rebuild -- NIMBLE has no supported way to serialize it, confirmed
# by the NIMBLE developers on the nimble-users list ("Pausing and Restarting
# MCMC iterations"), who also confirm the model build must be redone in each
# R session. So a resumed chunk still re-adapts its proposals from default
# scale. That is a real but much smaller effect than restarting from inits,
# and it is the ONLY thing the earlier resume_burnin draft addressed.
#
# Because of that, this file is NOT a licence to keep chunking. It exists so
# that IF chunk-resume is kept, it is at least correct. The cheaper structural
# option remains: run one chain per job and add chains rather than extending
# them.
#
# STATUS: UNVALIDATED. Do not roll out to production fits until the A/B test
# at the bottom of this file has been run and passed.
# =============================================================================

# --- which nodes define the chain's state -------------------------------------
# All non-data stochastic nodes: everything the sampler updates. Deliberately
# derived FROM THE MODEL rather than hardcoded, so a model change can't
# silently leave a node un-checkpointed.
state_node_names <- function(model) {
  nodes <- model$getNodeNames(stochOnly = TRUE, includeData = FALSE)
  if (!length(nodes)) stop("state_node_names(): no non-data stochastic nodes found")
  nodes
}

# --- SAVE ---------------------------------------------------------------------
save_chain_state_v2 <- function(file, model, Cmcmc, samples_matrix, iter_total) {

  nodes <- state_node_names(model)

  # values() pulls plain numerics out of the compiled model -- no pointers.
  vals <- values(model, nodes)

  if (!is.numeric(vals) || anyNA(vals) || any(!is.finite(vals))) {
    stop("save_chain_state_v2(): non-finite or NA node values; refusing to ",
         "write a corrupt checkpoint")
  }

  chain_list <- list(
    node_names  = nodes,
    node_values = vals,
    samples     = samples_matrix,
    iter_total  = iter_total,
    # provenance so a checkpoint can never be silently mixed across models
    n_nodes     = length(nodes),
    saved_at    = Sys.time(),
    format      = "numeric_values_v2"
  )

  saveRDS(chain_list, file = file)
}

# --- RESTORE ------------------------------------------------------------------
# Returns obj (as the original did) so the caller keeps samples/iter_total.
restore_chain_state_v2 <- function(file, model, Cmcmc) {

  obj <- readRDS(file)

  if (is.null(obj$format) || obj$format != "numeric_values_v2") {
    stop("restore_chain_state_v2(): checkpoint '", basename(file), "' is in the ",
         "OLD compiled-mvSaved format, which does not restore state. Refusing ",
         "to resume from it -- start this chain fresh, or the chain will ",
         "silently restart from inits.")
  }

  nodes <- state_node_names(model)

  if (!identical(nodes, obj$node_names)) {
    stop("restore_chain_state_v2(): model node structure does not match the ",
         "checkpoint (", length(nodes), " nodes now vs ", obj$n_nodes,
         " saved). Refusing to resume across a model change.")
  }

  # 1. plain numerics -> model
  values(model, nodes) <- obj$node_values

  # 2. recalculate deterministic nodes + logProb from the restored values
  lp <- model$calculate()
  if (!is.finite(lp)) {
    stop("restore_chain_state_v2(): model logProb is ", lp, " after restore; ",
         "the restored state is not usable")
  }

  # 3. model -> mvSaved, through the compiled MCMC's own accessor. This is the
  #    step the original tried to do by assigning the deserialized object.
  Cmcmc$mvSaved[["logProb"]]  # touch to confirm the compiled object is live
  nimCopy(from = model, to = Cmcmc$mvSaved, row = 1, logProb = TRUE)

  # 4. read back and confirm the copy actually landed
  check <- values(model, nodes)
  if (!isTRUE(all.equal(check, obj$node_values, tolerance = 1e-10))) {
    stop("restore_chain_state_v2(): node values did not survive the restore ",
         "(max abs diff ", max(abs(check - obj$node_values)), ")")
  }

  obj$restored_values <- setNames(obj$node_values, nodes)
  obj
}

# --- CALLER-SIDE CONTINUITY ASSERTION ----------------------------------------
# Call this immediately AFTER Cmcmc$run(chunk_iter) on a resumed chunk. It is
# the check whose absence let the original bug run undetected for months: it
# compares the first NEW draw against the last OLD draw for the monitored
# parameters, and fails if the chain jumped to a fresh-start configuration.
#
# `init_values` is a named numeric of the model's inits for a few sentinel
# monitored params (e.g. det_intercept = 0.5, overdisp_inat = 0.1). If the
# first post-resume draw sits exactly on those, we have the original bug back.
assert_resume_continuity <- function(samples_so_far, new_samples,
                                     init_values = NULL,
                                     sd_tol = 25) {

  if (is.null(samples_so_far) || !nrow(samples_so_far)) return(invisible(TRUE))

  last_old  <- samples_so_far[nrow(samples_so_far), , drop = TRUE]
  first_new <- new_samples[1, , drop = TRUE]
  common    <- intersect(names(last_old), names(first_new))

  # (a) exact-init detection -- the signature of the original defect
  if (!is.null(init_values)) {
    for (p in intersect(names(init_values), common)) {
      if (abs(first_new[[p]] - init_values[[p]]) < 1e-9 &&
          abs(last_old[[p]]  - init_values[[p]]) > 1e-6) {
        stop("assert_resume_continuity(): '", p, "' is exactly at its INIT (",
             init_values[[p]], ") on the first post-resume draw, but was ",
             last_old[[p]], " pre-resume. The chain restarted -- the ",
             "checkpoint restore did not work.")
      }
    }
  }

  # (b) implausible jump relative to the previous chunk's own scale
  ref <- tail(samples_so_far, min(2000L, nrow(samples_so_far)))
  sds <- apply(ref[, common, drop = FALSE], 2, sd)
  ok  <- sds > 0 & is.finite(sds)
  if (any(ok)) {
    z <- abs(first_new[common][ok] - last_old[common][ok]) / sds[ok]
    if (max(z) > sd_tol) {
      worst <- names(which.max(z))
      stop("assert_resume_continuity(): '", worst, "' jumped ",
           round(max(z), 1), " within-chunk SDs across the resume boundary ",
           "(", last_old[[worst]], " -> ", first_new[[worst]], "). ",
           "State did not carry over cleanly.")
    }
  }

  invisible(TRUE)
}

# =============================================================================
# VALIDATION PLAN -- run BEFORE any production use
# =============================================================================
# Two chunks, small (chunk_iter = 2000, burnin = 500), one real resume
# boundary, on a real species bundle. Two arms:
#
#   orig : unmodified save/restore
#   v2   : these functions + assert_resume_continuity()
#
# PASS requires ALL of:
#   1. orig reproduces the defect (first post-resume draw exactly on inits) --
#      confirms the test is sensitive enough to detect it.
#   2. v2's first post-resume draw is NOT on the inits, and its distance from
#      the last pre-resume draw is within a few within-chunk SDs.
#   3. v2's block means show no spike at the boundary block, on ALL monitored
#      params -- not just the one the bug was found on. The earlier
#      resume_burnin draft passed on overdisp_inat and made theta1 WORSE;
#      a per-parameter table is required, not a summary statistic.
#   4. assert_resume_continuity() fires on the orig arm and passes on v2 --
#      i.e. the guard itself is verified, not assumed.
#
# Note (3) is where the previous attempt failed. Do not accept a pass that is
# only demonstrated on one parameter.
#
# A residual, EXPECTED difference remains even on a pass: sampler proposal
# scales still reset (see WHAT THIS DOES NOT FIX). Judge the boundary against
# a single-run, never-resumed chain of the same total length -- that is the
# reference for "as good as not chunking," and the only fair benchmark.
# =============================================================================
