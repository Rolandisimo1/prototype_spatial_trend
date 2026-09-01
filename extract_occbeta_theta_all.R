#!/usr/bin/env Rscript
# =============================================================================
# extract_occbeta_theta_all.R
# -----------------------------------------------------------------------------
# occ_beta (per-covariate) + theta0/theta1 (Goldstein congruence) posteriors
# for every converged _single fit, for the report rebuild.
#
# Generalizes extract_moose_v1fix_occbeta.R and extract_moose_v1fix_theta.R,
# which are both hardcoded to the two moose_v1fix models AND read the CHUNKED
# chain_<sp>_<i>.RDS files. Those chunked files carry the resume-boundary
# defect, so they are not the fits any current diagnostic or posterior refers
# to. This reads chain_<sp>_<i>_single.RDS only and stops if one is missing.
#
# theta1 is Goldstein et al. (bioRxiv 2025.01.17.633640)'s primary camera/iNat
# congruence metric; our model gives it the identical structural role,
# log(mu) = theta0 + theta1 * log(sum(lambda)) inside calcIntensity_SVC.
#
# numOccCovars is DISCOVERED per model, not assumed to be moose_v1fix's 10 --
# the v2b bundles need not carry the same covariate set, and mislabeling a
# covariate is worse than not labeling it. Names come from the model's own
# occ_covars column names when present; when absent the CSV says so in a
# covariate_label column rather than inventing an order.
# =============================================================================
suppressPackageStartupMessages(library(coda))

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
MODELS <- c(
  "bobcat_v2b_national_scalar", "bobcat_v2b_ecoregion",
  "white-tailed_deer_v2b_national_scalar", "white-tailed_deer_v2b_ecoregion",
  "moose_v2b_national_scalar", "moose_v2b_ecoregion",
  "moose_v1fix9_national_scalar", "moose_v1fix9_ecoregion"
)

q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))
summ <- function(v, label) data.frame(
  parameter = label, mean = mean(v), median = median(v), sd = sd(v),
  q025 = q(v, 0.025), q975 = q(v, 0.975),
  p_negative = mean(v < 0), n_na = sum(!is.finite(v)),
  stringsAsFactors = FALSE)

all_occ <- list(); all_theta <- list()

for (sp in MODELS) {
  D <- file.path(PROJ, "HPC", sp)
  cat("\n=====", sp, "=====\n")
  files <- file.path(D, sprintf("chain_%s_%d_single.RDS", sp, 1:3))
  if (!all(file.exists(files))) { cat("  MISSING _single chains -- skipped\n"); next }

  # covariate labels from the model's own bundle, never assumed
  labs <- NULL
  ib <- file.path(D, paste0("input_data_", sp, ".RDS"))
  if (file.exists(ib)) {
    idat <- readRDS(ib); oc <- idat$constants_list$occ_covars
    if (!is.null(colnames(oc))) labs <- colnames(oc)
    cat("  numOccCovars =", idat$constants_list$numOccCovars,
        " labels:", if (is.null(labs)) "(none in bundle)" else paste(labs, collapse = ", "), "\n")
    rm(idat); gc(verbose = FALSE)
  }

  chains <- vector("list", 3)
  for (i in 1:3) {
    obj <- readRDS(files[i]); s <- obj$samples
    if (i == 1) {
      bc <- grep("^occ_beta\\[", colnames(s), value = TRUE)
      bc <- bc[order(as.integer(gsub("occ_beta\\[|\\]", "", bc)))]
      keep <<- c(bc, intersect(c("theta0", "theta1"), colnames(s)))
      cat("  occ_beta cols:", length(bc), "  theta present:",
          paste(intersect(c("theta0","theta1"), colnames(s)), collapse = ","), "\n")
    }
    chains[[i]] <- s[, keep, drop = FALSE]
    rm(obj, s); gc(verbose = FALSE)
  }
  rh <- gelman.diag(mcmc.list(lapply(chains, mcmc)), autoburnin = FALSE,
                    multivariate = FALSE)$psrf[, 1]
  comb <- do.call(rbind, chains)

  bcols <- grep("^occ_beta\\[", colnames(comb), value = TRUE)
  bcols <- bcols[order(as.integer(gsub("occ_beta\\[|\\]", "", bcols)))]
  od <- do.call(rbind, lapply(bcols, function(p) summ(comb[, p], p)))
  od$covariate_index <- as.integer(gsub("occ_beta\\[|\\]", "", bcols))
  od$covariate_label <- if (!is.null(labs) && length(labs) == nrow(od))
                          labs[od$covariate_index] else NA_character_
  od$rhat <- unname(rh[od$parameter]); od$model <- sp
  write.csv(od, file.path(D, paste0(sp, "_occbeta_posterior.csv")), row.names = FALSE)
  all_occ[[sp]] <- od
  print(od[, c("parameter","covariate_label","mean","q025","q975","rhat")],
        row.names = FALSE, digits = 4)

  tn <- intersect(c("theta0","theta1"), colnames(comb))
  if (length(tn)) {
    td <- do.call(rbind, lapply(tn, function(p) summ(comb[, p], p)))
    td$rhat <- unname(rh[td$parameter]); td$model <- sp
    write.csv(td, file.path(D, paste0(sp, "_theta_posterior.csv")), row.names = FALSE)
    all_theta[[sp]] <- td
    cat("  -- theta --\n"); print(td[, c("parameter","mean","q025","q975","rhat")],
                                  row.names = FALSE, digits = 4)
  }
  rm(comb, chains); gc(verbose = FALSE)
}

OUT <- file.path(PROJ, "prototype_spatial_trend")
write.csv(do.call(rbind, all_occ),   file.path(OUT, "occbeta_posterior_all_models.csv"), row.names = FALSE)
write.csv(do.call(rbind, all_theta), file.path(OUT, "theta_posterior_all_models.csv"),   row.names = FALSE)
cat("\n===== theta1 across models (Goldstein congruence) =====\n")
th <- do.call(rbind, all_theta); print(th[th$parameter == "theta1",
      c("model","mean","q025","q975","rhat")], row.names = FALSE, digits = 4)
cat("\nDONE\n")
