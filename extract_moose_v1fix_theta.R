#!/usr/bin/env Rscript
# extract_moose_v1fix_theta.R
#
# NOT YET RUN -- requires Hazel access (this session's ssh:hazel target was
# unreachable: "Permission denied (keyboard-interactive)"). Prepared for
# Claude Code / the user to execute on Hazel directly.
#
# Pulls theta0/theta1 full posterior summary (mean/median/sd/CI/rhat), both
# moose_v1fix models, per pending_hazel_extraction_notes.md item 7. Both
# nodes are directly monitored scalars in the same samples matrices already
# read by extract_moose_v1fix_posterior.R -- no new reconstruction needed.
#
# theta1 is Goldstein et al.'s (bioRxiv 2025.01.17.633640) primary camera/
# iNat congruence metric: log(mu) = theta0 + theta1*log(sum(lambda)).

suppressPackageStartupMessages(library(coda))

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"

THETA_NODES <- c("theta0", "theta1")

q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))

summarize_node <- function(v, label, species) {
  data.frame(
    parameter = label,
    mean      = mean(v, na.rm = TRUE),
    median    = median(v, na.rm = TRUE),
    sd        = sd(v, na.rm = TRUE),
    q025      = q(v, 0.025),
    q975      = q(v, 0.975),
    p_negative = mean(v < 0, na.rm = TRUE),
    n_na      = sum(!is.finite(v)),
    model     = species,
    stringsAsFactors = FALSE
  )
}

run_theta <- function(species) {
  D <- file.path(PROJ, "HPC", species)
  chains <- vector("list", 3)
  for (chain_id in 1:3) {
    f <- file.path(D, paste0("chain_", species, "_", chain_id, ".RDS"))
    stopifnot(file.exists(f))
    obj <- readRDS(f)
    s <- obj$samples
    have <- intersect(THETA_NODES, colnames(s))
    stopifnot(length(have) == length(THETA_NODES))  # both must be monitored
    chains[[chain_id]] <- s[, have, drop = FALSE]
    rm(obj, s); gc(verbose = FALSE)
  }

  mcmc_list <- coda::mcmc.list(lapply(chains, coda::mcmc))
  gd <- coda::gelman.diag(mcmc_list, autoburnin = FALSE, multivariate = FALSE)
  rhat <- gd$psrf[, 1]

  combined <- do.call(rbind, chains)
  df <- do.call(rbind, lapply(THETA_NODES, function(p)
    summarize_node(combined[, p], p, species)))
  df$rhat <- unname(rhat[df$parameter])

  out <- file.path(D, paste0(species, "_theta_posterior.csv"))
  write.csv(df, out, row.names = FALSE)
  cat("wrote", out, "\n")
  print(df, row.names = FALSE, digits = 4)
  invisible(df)
}

theta_ns  <- run_theta("moose_v1fix_national_scalar")
theta_eco <- run_theta("moose_v1fix_ecoregion")

combined_out <- rbind(theta_ns, theta_eco)
write.csv(combined_out, "moose_v1fix_theta_posterior_combined.csv", row.names = FALSE)
cat("\nwrote moose_v1fix_theta_posterior_combined.csv (combine both models' theta0/theta1)\n")
cat("\nEXTRACTION DONE\n")
