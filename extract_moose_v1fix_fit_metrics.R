#!/usr/bin/env Rscript
# extract_moose_v1fix_fit_metrics.R
#
# NOT YET RUN -- requires Hazel access (this session's ssh:hazel target was
# unreachable). Prepared for Claude Code / the user to execute on Hazel.
#
# Produces two overall model-fit metrics, requested as an "AUC-style" summary
# to complement the R-hat convergence diagnostics already extracted:
#
#   1. WAIC (both v1fix models) -- nimble's built-in predictive-fit score.
#      configureMCMC(..., enableWAIC = TRUE) was set at build time
#      (HPC_run_model_chunks_chain1.R L220 / moose chain scripts), but no
#      chunk-saving step in the checkpointed workflow ever called
#      Cmcmc$getWAIC() or saved a WAIC value -- obj$samples is all that's in
#      the chain_<species>_<n>.RDS files. WAIC must be computed post-hoc from
#      the compiled model + pooled samples using nimble::calculateWAIC(),
#      which needs the live nimble model object (cannot be done from CSVs
#      alone off Hazel).
#
#   2. Occupancy AUC (camera side only) -- a standard "how well does the
#      fitted psi/p separate detected from non-detected sites" metric.
#      Requires the real y[i,j] detection-history matrix from
#      input_data_moose_<model>.RDS, which is not available locally (only
#      posterior psi summaries were extracted, not the raw detection data
#      needed to score them against). Computed as posterior-mean AUC using
#      pROC::auc() on psi_mean vs. any(y[i,]==1) per site, averaged over a
#      subsample of posterior draws for a credible interval.
#
# No iNat-side AUC equivalent exists -- y_inat is a count, not binary
# detection/non-detection, so AUC does not apply there; a posterior
# predictive check (observed vs. simulated count distribution) would be the
# iNat-side analogue and is NOT attempted here (separate follow-up).

suppressPackageStartupMessages({
  library(nimble)
  library(coda)
})

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"

# ---------------------------------------------------------------------------
# Part 1: WAIC, computed post-hoc from the compiled model + full pooled samples
# ---------------------------------------------------------------------------
compute_waic <- function(species) {
  D <- file.path(PROJ, "HPC", species)
  input_data <- readRDS(file.path(D, paste0("input_data_", species, ".RDS")))

  model_code_file <- if (grepl("ecoregion", species)) "model_code_ecoregion.R" else "model_code_national_scalar.R"
  source(file.path(D, model_code_file))  # defines model_code_<...>

  model <- nimbleModel(get(ls(pattern = "^model_code_")[1]),
                        constants = input_data$constants_list,
                        data = input_data$data_list,
                        inits = input_data$inits_list,
                        calculate = FALSE)
  Cmodel <- compileNimble(model)

  chains <- list()
  for (chain_id in 1:3) {
    f <- file.path(D, paste0("chain_", species, "_", chain_id, ".RDS"))
    stopifnot(file.exists(f))
    obj <- readRDS(f)
    chains[[chain_id]] <- obj$samples
    rm(obj); gc(verbose = FALSE)
  }
  combined <- do.call(rbind, chains)

  # nimble::calculateWAIC requires monitored-node samples + the compiled model;
  # 'online = FALSE' style post-hoc WAIC (conditional, per-observation) --
  # see ?calculateWAIC for the exact signature in the installed nimble version.
  waic_result <- calculateWAIC(combined, Cmodel)

  cat(species, "WAIC:", waic_result$WAIC,
      " lppd:", waic_result$lppd, " pWAIC:", waic_result$pWAIC, "\n")
  data.frame(model = species, WAIC = waic_result$WAIC,
             lppd = waic_result$lppd, pWAIC = waic_result$pWAIC)
}

waic_ns  <- tryCatch(compute_waic("moose_v1fix_national_scalar"), error = function(e) {
  cat("WAIC FAILED for national_scalar:", conditionMessage(e), "\n"); NULL
})
waic_eco <- tryCatch(compute_waic("moose_v1fix_ecoregion"), error = function(e) {
  cat("WAIC FAILED for ecoregion:", conditionMessage(e), "\n"); NULL
})

waic_out <- do.call(rbind, list(waic_ns, waic_eco))
if (!is.null(waic_out)) {
  write.csv(waic_out, "moose_v1fix_waic.csv", row.names = FALSE)
  cat("wrote moose_v1fix_waic.csv\n")
}

# ---------------------------------------------------------------------------
# Part 2: Camera-side occupancy AUC (posterior mean + draw-level CI)
# ---------------------------------------------------------------------------
compute_occ_auc <- function(species, n_draws = 200) {
  D <- file.path(PROJ, "HPC", species)
  input_data <- readRDS(file.path(D, paste0("input_data_", species, ".RDS")))
  y <- input_data$data_list$y            # [nsite, maxJ] detection history
  J <- input_data$constants_list$J        # visits per site
  nsite <- input_data$constants_list$nsite

  # site-level "ever detected" outcome
  ever_detected <- vapply(seq_len(nsite), function(i) {
    any(y[i, seq_len(J[i])] == 1, na.rm = TRUE)
  }, logical(1))

  chains <- list()
  for (chain_id in 1:3) {
    f <- file.path(D, paste0("chain_", species, "_", chain_id, ".RDS"))
    obj <- readRDS(f)
    psi_cols <- grep("^psi\\[", colnames(obj$samples), value = TRUE)
    stopifnot(length(psi_cols) == nsite)  # psi must be monitored per site
    chains[[chain_id]] <- obj$samples[, psi_cols, drop = FALSE]
    rm(obj); gc(verbose = FALSE)
  }
  combined <- do.call(rbind, chains)

  if (!requireNamespace("pROC", quietly = TRUE)) {
    install.packages("pROC", repos = "https://cloud.r-project.org")
  }
  library(pROC)

  draw_idx <- sample(seq_len(nrow(combined)), min(n_draws, nrow(combined)))
  aucs <- vapply(draw_idx, function(d) {
    psi_draw <- combined[d, ]
    as.numeric(pROC::auc(pROC::roc(ever_detected, psi_draw, quiet = TRUE)))
  }, numeric(1))

  cat(species, "occupancy AUC: mean =", round(mean(aucs), 4),
      " 95% CI =", round(quantile(aucs, c(0.025, 0.975)), 4), "\n")
  data.frame(model = species, auc_mean = mean(aucs),
             auc_q025 = quantile(aucs, 0.025), auc_q975 = quantile(aucs, 0.975))
}

auc_ns  <- tryCatch(compute_occ_auc("moose_v1fix_national_scalar"), error = function(e) {
  cat("AUC FAILED for national_scalar:", conditionMessage(e), "\n"); NULL
})
auc_eco <- tryCatch(compute_occ_auc("moose_v1fix_ecoregion"), error = function(e) {
  cat("AUC FAILED for ecoregion:", conditionMessage(e), "\n"); NULL
})

auc_out <- do.call(rbind, list(auc_ns, auc_eco))
if (!is.null(auc_out)) {
  write.csv(auc_out, "moose_v1fix_occupancy_auc.csv", row.names = FALSE)
  cat("wrote moose_v1fix_occupancy_auc.csv\n")
}

cat("\nEXTRACTION DONE (check both CSVs; either can be NULL if psi/y not\n",
    "monitored/available in the checkpoint -- see FAILED messages above)\n")
