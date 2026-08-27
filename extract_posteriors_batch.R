#!/usr/bin/env Rscript
# =============================================================================
# extract_posteriors_batch.R
# -----------------------------------------------------------------------------
# Posterior extraction for the four converged _single fits that have R-hat CSVs
# but no posterior CSVs:
#     bobcat_v2b_national_scalar      (priority -- low-abundance validation
#                                      species, currently blocked from the report)
#     moose_v1fix9_national_scalar    } together these unblock the
#     moose_v1fix9_ecoregion          } mask-vs-covariate comparison vs v2b
#     bobcat_v2b_ecoregion            (trend params converged at R-hat <= 1.004;
#                                      13 CAR/MCMT_tau offenders -- see the
#                                      caveat written into its own CSV)
#
# Generalized from extract_moose_v1fix_posterior.R. Three differences:
#
#  1. READS THE _single CHAINS. The chunked chain_<species>_<i>.RDS files carry
#     the resume-boundary defect (chains restart from inits at every boundary,
#     REVISION_NOTE_resume_defect.md), so they are NOT the fits the 8/26
#     convergence diagnostics were run on. Every file this script reads is
#     chain_<species>_<i>_single.RDS. Asserted, not assumed: it stops if a
#     _single file is missing rather than silently falling back.
#  2. Output filenames are derived from the species, not hardcoded to
#     "moose_v1fix" -- the original would have overwritten the v1fix CSVs with
#     v1fix9 content under the old name.
#  3. Runs read-only against HPC/<species>/; writes only new CSVs.
#
# No refitting. Burn-in was run once and discarded before checkpointing
# (HPC_run_model_chunks_chain1.R::run_chain_chunk), so obj$samples already
# excludes it.
# =============================================================================

suppressPackageStartupMessages(library(coda))

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"

NATIONAL_NODES <- c("total_var_beta", "year_beta", "year_var", "snr",
                    "trend_robust_indicator")

q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))

title_case <- function(x) {
  vapply(strsplit(tolower(x), " ", fixed = TRUE), function(w) {
    paste(toupper(substring(w, 1, 1)), substring(w, 2), sep = "", collapse = " ")
  }, character(1))
}

summarize_node <- function(v, label) {
  data.frame(
    parameter = label,
    mean      = mean(v, na.rm = TRUE),
    median    = median(v, na.rm = TRUE),
    sd        = sd(v, na.rm = TRUE),
    q025      = q(v, 0.025),
    q975      = q(v, 0.975),
    p_negative = mean(v < 0, na.rm = TRUE),
    n_na      = sum(!is.finite(v)),
    stringsAsFactors = FALSE
  )
}

chain_file <- function(species, chain_id)
  file.path(PROJ, "HPC", species,
            sprintf("chain_%s_%d_single.RDS", species, chain_id))

load_chains <- function(species, keep_cols) {
  out   <- vector("list", 3)
  iters <- integer(3)
  for (chain_id in 1:3) {
    f <- chain_file(species, chain_id)
    if (!file.exists(f))
      stop("missing _single chain: ", f,
           "\nRefusing to fall back to the chunked chain_*.RDS -- those carry ",
           "the resume-boundary defect and are not what the 8/26 R-hat ",
           "diagnostics were run on.")
    obj <- readRDS(f)
    s   <- obj$samples
    have <- intersect(keep_cols, colnames(s))
    out[[chain_id]]  <- s[, have, drop = FALSE]
    iters[chain_id]  <- obj$iter_total
    cat("  chain", chain_id, ": iter_total =", obj$iter_total,
        " rows =", nrow(s), " kept cols =", length(have), "\n")
    rm(obj, s); gc(verbose = FALSE)
  }
  list(chains = out, iter_total = iters)
}

run_model <- function(species, is_ecoregion, caveat = NULL) {
  cat("\n\n############################################################\n")
  cat("##  ", species, if (is_ecoregion) "(ecoregion)" else "(national scalar)", "\n")
  cat("############################################################\n")
  D <- file.path(PROJ, "HPC", species)

  f1 <- chain_file(species, 1)
  if (!file.exists(f1)) stop("missing _single chain: ", f1)
  obj1 <- readRDS(f1)
  all_cols <- colnames(obj1$samples)
  cat("chain 1: ", ncol(obj1$samples), " monitored cols, ",
      nrow(obj1$samples), " draws\n", sep = "")
  rm(obj1); gc(verbose = FALSE)

  wanted <- NATIONAL_NODES
  nregion <- NA_integer_; region_name <- NULL
  n_camera_sites <- NULL; n_inat_cells <- NULL

  input_data <- readRDS(file.path(D, paste0("input_data_", species, ".RDS")))
  cl <- input_data$constants_list

  if (is_ecoregion) {
    yr_cols <- grep("^year_region\\[", all_cols, value = TRUE)
    nregion <- length(yr_cols)
    stopifnot(nregion > 0)
    wanted <- c(paste0("year_region[", seq_len(nregion), "]"),
                "sigma_region", NATIONAL_NODES)
    stopifnot(cl$nregion == nregion,
              length(cl$ecoregion_of_cell100) == cl$ncell100)
    prepped <- readRDS(file.path(PROJ, "prototype_spatial_trend",
                                 "prepped_sim_inputs.RDS"))
    lev <- prepped$ecoregion_levels
    stopifnot(length(lev) == nregion)
    region_name <- title_case(lev)
    rm(prepped); gc(verbose = FALSE)
    cat("nregion:", nregion, "\n")

    region_of_site <- cl$ecoregion_of_cell100[cl$cell]
    n_camera_sites <- as.integer(table(factor(region_of_site,
                                              levels = seq_len(nregion))))
    monitored_g <- unique(na.omit(as.vector(cl$inat_cells_by_year)))
    region_of_g <- cl$ecoregion_of_cell100[cl$inat_cell100[monitored_g]]
    n_inat_cells <- as.integer(table(factor(region_of_g,
                                            levels = seq_len(nregion))))
    cat("n_camera_sites:", paste(n_camera_sites, collapse = ", "), "\n")
    cat("n_inat_cells:  ", paste(n_inat_cells, collapse = ", "), "\n")
  }

  yv <- cl$year_vals
  cat("year_vals (n=", length(yv), "): ",
      paste(round(range(yv), 4), collapse = " .. "),
      "  step=", if (length(yv) > 1) round(diff(sort(yv))[1], 4) else NA,
      "\n", sep = "")
  rm(input_data); gc(verbose = FALSE)

  present <- intersect(wanted, all_cols)
  missing <- setdiff(wanted, all_cols)
  cat("\nrequested nodes PRESENT:", paste(present, collapse = ", "), "\n")
  cat("requested nodes MISSING:",
      if (length(missing)) paste(missing, collapse = ", ") else "(none)", "\n\n")

  ld <- load_chains(species, present)
  chain_full <- ld$chains

  mcmc_list <- mcmc.list(lapply(chain_full, mcmc))
  gd <- gelman.diag(mcmc_list, autoburnin = FALSE, multivariate = FALSE)
  rhat <- gd$psrf[, 1]
  cat("\n===== R-hat over reported params =====\n")
  print(round(sort(rhat, decreasing = TRUE), 4))
  rhat_gated <- rhat[!(names(rhat) %in% "trend_robust_indicator")]
  cat("MAX R-hat (reported params, excl. trend_robust_indicator):",
      round(suppressWarnings(max(rhat_gated, na.rm = TRUE)), 5), "\n")

  combined <- do.call(rbind, chain_full)
  cat("combined posterior draws:", nrow(combined), "\n")

  nat_present <- intersect(NATIONAL_NODES, colnames(combined))
  nat_df <- do.call(rbind, lapply(nat_present, function(p)
    summarize_node(combined[, p], p)))
  nat_df$rhat  <- unname(rhat[nat_df$parameter])
  nat_df$model <- species

  # snr: deterministic in the model code but NOT monitored, so derive it.
  # year_var straddles 0 -> ratio-of-normals, Cauchy-like, undefined mean.
  # Median and IQR only; trend_robust_indicator = P(snr > 1) is the usable form.
  if (all(c("year_beta", "year_var") %in% colnames(combined))) {
    snr_draws <- combined[, "year_beta"] / combined[, "year_var"]
    cat("\n----- derived snr = year_beta/year_var (NOT monitored) -----\n")
    cat("  median:", round(median(snr_draws), 4),
        "  IQR:", paste(round(q(snr_draws, c(0.25, 0.75)), 4), collapse = " .. "), "\n")
    cat("  P(snr > 1):", round(mean(snr_draws > 1), 4), "\n")
    cat("  share |snr| > 100 (tail blowup):",
        round(mean(abs(snr_draws) > 100), 4), "\n")
    nat_df <- rbind(nat_df, data.frame(
      parameter = "snr_derived", mean = NA_real_,
      median = median(snr_draws), sd = NA_real_,
      q025 = q(snr_draws, 0.25), q975 = q(snr_draws, 0.75),
      p_negative = mean(snr_draws < 0), n_na = sum(!is.finite(snr_draws)),
      rhat = NA_real_, model = species, stringsAsFactors = FALSE))
  }

  nat_df$caveat <- if (is.null(caveat)) "" else caveat

  if (is_ecoregion) {
    sig <- summarize_node(combined[, "sigma_region"], "sigma_region")
    sig$rhat <- unname(rhat["sigma_region"]); sig$model <- species
    sig$caveat <- nat_df$caveat[1]
    nat_df <- rbind(sig, nat_df)

    # year_region[r] is a mean-zero DEVIATION; the absolute regional iNat trend
    # is total_var_beta + year_region[r]. Report both.
    gr <- function(r) combined[, paste0("year_region[", r, "]")]
    abs_draws <- lapply(seq_len(nregion), function(r) combined[, "total_var_beta"] + gr(r))
    reg_df <- data.frame(
      region_index       = seq_len(nregion),
      region_name        = region_name,
      year_region_mean   = sapply(seq_len(nregion), function(r) mean(gr(r))),
      year_region_median = sapply(seq_len(nregion), function(r) median(gr(r))),
      year_region_sd     = sapply(seq_len(nregion), function(r) sd(gr(r))),
      year_region_q025   = sapply(seq_len(nregion), function(r) q(gr(r), 0.025)),
      year_region_q975   = sapply(seq_len(nregion), function(r) q(gr(r), 0.975)),
      p_negative         = sapply(seq_len(nregion), function(r) mean(gr(r) < 0)),
      rhat               = sapply(seq_len(nregion), function(r)
                                    unname(rhat[paste0("year_region[", r, "]")])),
      abs_trend_mean     = sapply(abs_draws, mean),
      abs_trend_median   = sapply(abs_draws, median),
      abs_trend_q025     = sapply(abs_draws, function(v) q(v, 0.025)),
      abs_trend_q975     = sapply(abs_draws, function(v) q(v, 0.975)),
      abs_p_negative     = sapply(abs_draws, function(v) mean(v < 0)),
      n_camera_sites     = n_camera_sites,
      n_inat_cells       = n_inat_cells,
      stringsAsFactors = FALSE
    )
    out_reg <- file.path(D, paste0(species, "_ecoregion_posterior.csv"))
    write.csv(reg_df, out_reg, row.names = FALSE); cat("\nwrote", out_reg, "\n")
    print(reg_df, row.names = FALSE, digits = 4)
    out_glob <- file.path(D, paste0(species, "_global.csv"))
    write.csv(nat_df, out_glob, row.names = FALSE); cat("\nwrote", out_glob, "\n")
  } else {
    out_nat <- file.path(D, paste0(species, "_posterior.csv"))
    write.csv(nat_df, out_nat, row.names = FALSE); cat("\nwrote", out_nat, "\n")
  }

  cat("\n===== national-level params (", species, ") =====\n")
  print(nat_df, row.names = FALSE, digits = 4)
  rm(combined, chain_full); gc(verbose = FALSE)
  invisible(NULL)
}

run_model("bobcat_v2b_national_scalar",   is_ecoregion = FALSE)
run_model("moose_v1fix9_national_scalar", is_ecoregion = FALSE)
run_model("moose_v1fix9_ecoregion",       is_ecoregion = TRUE)
run_model("bobcat_v2b_ecoregion",         is_ecoregion = TRUE,
          caveat = paste("Trend params converged (R-hat <= 1.004); the FIT as a whole did not:",
                         "13 gated params exceed 1.1 (MCMT_tau 1.182 plus link_occ_intercept",
                         "cells 854-855, 873-875, 883-885, 891-894). Trend numbers here are",
                         "usable; the CAR intercept field and MCMT SVC surface are NOT."))

cat("\nEXTRACTION DONE\n")
