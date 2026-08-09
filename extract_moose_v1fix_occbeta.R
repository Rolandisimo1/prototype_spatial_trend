#!/usr/bin/env Rscript
# extract_moose_v1fix_occbeta.R
#
# NOT YET RUN -- requires Hazel access (this session's ssh:hazel target was
# unreachable). Prepared for Claude Code / the user to execute on Hazel.
#
# Pulls the full occ_beta[1:10] posterior (mean/median/CI/rhat) for both
# moose_v1fix models -- the 10 occupancy covariates retained per
# make_reduced_input.R's final decision:
#   Human_pop, NDVI_mean, Ag, Deciduous, Evergreen, Mixed,
#   terrain_ruggedness, soil_clay, soil_silt, soil_sand
# (order matches numOccCovars column order in occ_covars / xdat_inat -- must
# be verified against the actual input_data_moose_*.RDS colnames before
# trusting the label assignment below; this script does that check and
# fails loudly if the order can't be confirmed.)
#
# occ_beta is a MONITORED node (conf$addMonitors("occ_beta", ...) in the
# chain scripts) but was never extracted into a per-covariate CSV -- only
# its R-hat-by-family summary exists (moose_v1fix_*_rhat_by_family.csv,
# worst_param = "occ_beta[1]"), which gives no per-covariate breakdown.

suppressPackageStartupMessages(library(coda))

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"

COVAR_NAMES <- c("Human_pop", "NDVI_mean", "Ag", "Deciduous", "Evergreen",
                  "Mixed", "terrain_ruggedness", "soil_clay", "soil_silt", "soil_sand")

q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))

summarize_beta <- function(v, label, species) {
  data.frame(
    parameter  = label,
    mean       = mean(v, na.rm = TRUE),
    median     = median(v, na.rm = TRUE),
    sd         = sd(v, na.rm = TRUE),
    q025       = q(v, 0.025),
    q975       = q(v, 0.975),
    p_negative = mean(v < 0, na.rm = TRUE),
    model      = species,
    stringsAsFactors = FALSE
  )
}

run_occbeta <- function(species) {
  D <- file.path(PROJ, "HPC", species)

  # verify covariate order against the actual input bundle before labeling
  input_data <- readRDS(file.path(D, paste0("input_data_", species, ".RDS")))
  actual_cols <- colnames(input_data$constants_list$occ_covars)
  if (!is.null(actual_cols)) {
    stopifnot(
      "occ_covars column order does not match COVAR_NAMES -- fix the labels before trusting this output" =
        identical(actual_cols, COVAR_NAMES)
    )
    cat("  confirmed occ_covars column order matches COVAR_NAMES for", species, "\n")
  } else {
    cat("  WARNING: occ_covars has no column names in the input bundle for", species,
        "-- labels below are POSITIONAL ASSUMPTIONS from make_reduced_input.R,",
        "not independently verified. Cross-check before reporting.\n")
  }
  rm(input_data); gc(verbose = FALSE)

  chains <- list()
  for (chain_id in 1:3) {
    f <- file.path(D, paste0("chain_", species, "_", chain_id, ".RDS"))
    obj <- readRDS(f)
    s <- obj$samples
    beta_cols <- grep("^occ_beta\\[", colnames(s), value = TRUE)
    stopifnot(length(beta_cols) == length(COVAR_NAMES))
    # sort by numeric index to guarantee occ_beta[1] < occ_beta[2] < ... order
    idx <- as.integer(gsub("occ_beta\\[|\\]", "", beta_cols))
    beta_cols <- beta_cols[order(idx)]
    chains[[chain_id]] <- s[, beta_cols, drop = FALSE]
    rm(obj, s); gc(verbose = FALSE)
  }

  mcmc_list <- coda::mcmc.list(lapply(chains, coda::mcmc))
  gd <- coda::gelman.diag(mcmc_list, autoburnin = FALSE, multivariate = FALSE)
  rhat <- gd$psrf[, 1]

  combined <- do.call(rbind, chains)
  df <- do.call(rbind, lapply(seq_along(COVAR_NAMES), function(i)
    summarize_beta(combined[, i], COVAR_NAMES[i], species)))
  df$rhat <- unname(rhat)

  out <- file.path(D, paste0(species, "_occbeta_posterior.csv"))
  write.csv(df, out, row.names = FALSE)
  cat("wrote", out, "\n")
  print(df[, c("parameter", "mean", "q025", "q975", "rhat")], row.names = FALSE, digits = 4)
  invisible(df)
}

beta_ns  <- run_occbeta("moose_v1fix_national_scalar")
beta_eco <- run_occbeta("moose_v1fix_ecoregion")

combined_out <- rbind(beta_ns, beta_eco)
write.csv(combined_out, "moose_v1fix_occbeta_posterior_combined.csv", row.names = FALSE)
cat("\nwrote moose_v1fix_occbeta_posterior_combined.csv\n")
cat("\nEXTRACTION DONE\n")
