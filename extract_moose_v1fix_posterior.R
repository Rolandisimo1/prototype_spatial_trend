#!/usr/bin/env Rscript
# extract_moose_v1fix_posterior.R
#
# Read-only extraction of the v1fix moose posteriors (moose_v1fix_national_scalar
# and moose_v1fix_ecoregion) for the model-vs-agency trend comparison.
#
# Adapted from HPC/moose_ecoregion/extract_moose_ecoregion_posterior.R, which
# produced the now-INVALIDATED pre-Fix1 CSVs. This version reads ONLY the
# _v1fix_ dirs (Fix 1 inat_effort column-sort + 10-covariate KEEP list) and
# never touches the old dirs.
#
# No refitting -- loads the checkpointed chain_<species>_{1,2,3}.RDS as-is.
# Burn-in (5000 iters) was run once and discarded before checkpointing began
# (HPC_run_model_chunks_chain1.R::run_chain_chunk), so obj$samples already
# excludes it; nothing further to drop.
#
# Node availability is DISCOVERED, not assumed: year_var / snr /
# trend_robust_indicator were not extracted by the old script, so we report
# which requested nodes are actually monitored rather than failing on absence.

suppressPackageStartupMessages(library(coda))

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"

# national-level scalars requested for BOTH models
NATIONAL_NODES <- c("total_var_beta", "year_beta", "year_var", "snr",
                    "trend_robust_indicator")

q <- function(x, p) unname(quantile(x, p, na.rm = TRUE))

# "EASTERN TEMPERATE FORESTS" -> "Eastern Temperate Forests"
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

# ---------------------------------------------------------------------------
# load the 3 chains of one model, keeping only the columns we need
# ---------------------------------------------------------------------------
load_chains <- function(species, keep_cols) {
  D <- file.path(PROJ, "HPC", species)
  out <- vector("list", 3)
  iters <- integer(3)
  for (chain_id in 1:3) {
    f <- file.path(D, paste0("chain_", species, "_", chain_id, ".RDS"))
    stopifnot(file.exists(f))
    obj <- readRDS(f)
    s <- obj$samples
    if (chain_id == 1) {
      cat("  all monitored nodes (", ncol(s), "cols ) -- scalar/global ones:\n   ")
      scal <- grep("\\[", colnames(s), invert = TRUE, value = TRUE)
      cat(paste(scal, collapse = ", "), "\n")
    }
    have <- intersect(keep_cols, colnames(s))
    out[[chain_id]] <- s[, have, drop = FALSE]
    iters[chain_id] <- obj$iter_total
    cat("  chain", chain_id, ": iter_total =", obj$iter_total,
        " rows =", nrow(s), " kept cols =", length(have), "\n")
    rm(obj, s); gc(verbose = FALSE)
  }
  list(chains = out, iter_total = iters)
}

# ---------------------------------------------------------------------------
# per-model driver
# ---------------------------------------------------------------------------
run_model <- function(species, is_ecoregion) {
  cat("\n\n############################################################\n")
  cat("##  ", species, "\n")
  cat("############################################################\n")
  D <- file.path(PROJ, "HPC", species)

  # ---- discover columns from chain 1 (cheap peek) ----
  f1 <- file.path(D, paste0("chain_", species, "_1.RDS"))
  obj1 <- readRDS(f1)
  all_cols <- colnames(obj1$samples)
  rm(obj1); gc(verbose = FALSE)

  wanted <- NATIONAL_NODES
  nregion <- NA_integer_
  region_name <- NULL

  if (is_ecoregion) {
    yr_cols <- grep("^year_region\\[", all_cols, value = TRUE)
    nregion <- length(yr_cols)
    stopifnot(nregion > 0)
    wanted <- c(paste0("year_region[", seq_len(nregion), "]"),
                "sigma_region", NATIONAL_NODES)

    # region names from the model's own bundle, cross-checked against the
    # sim-prep levels used for every earlier ecoregion figure.
    input_data <- readRDS(file.path(D, paste0("input_data_", species, ".RDS")))
    cl <- input_data$constants_list
    stopifnot(cl$nregion == nregion,
              length(cl$ecoregion_of_cell100) == cl$ncell100)

    prepped <- readRDS(file.path(PROJ, "prototype_spatial_trend",
                                 "prepped_sim_inputs.RDS"))
    lev <- prepped$ecoregion_levels
    stopifnot(length(lev) == nregion)
    region_name <- title_case(lev)
    cat("nregion:", nregion, "\n")
    cat("region levels:", paste(lev, collapse = " | "), "\n")

    # ---- data-support counts per region (same logic as the old script) ----
    region_of_site <- cl$ecoregion_of_cell100[cl$cell]
    n_camera_sites <- as.integer(table(factor(region_of_site,
                                              levels = seq_len(nregion))))
    monitored_g <- unique(na.omit(as.vector(cl$inat_cells_by_year)))
    region_of_g  <- cl$ecoregion_of_cell100[cl$inat_cell100[monitored_g]]
    n_inat_cells <- as.integer(table(factor(region_of_g,
                                            levels = seq_len(nregion))))
    cat("n_camera_sites:", paste(n_camera_sites, collapse = ", "), "\n")
    cat("n_inat_cells:  ", paste(n_inat_cells, collapse = ", "), "\n")

    # year_vals is the covariate total_var_beta/year_effect multiply, so the
    # trend's UNITS are "per one unit of year_vals" -- report its range so the
    # coefficient can be converted to a per-calendar-year rate.
    yv <- cl$year_vals
    cat("year_vals (n=", length(yv), "): ", paste(round(range(yv), 4), collapse = " .. "),
        "  step=", if (length(yv) > 1) round(diff(sort(yv))[1], 4) else NA, "\n", sep = "")
    rm(input_data); gc(verbose = FALSE)
  } else {
    input_data <- readRDS(file.path(D, paste0("input_data_", species, ".RDS")))
    yv <- input_data$constants_list$year_vals
    cat("year_vals (n=", length(yv), "): ", paste(round(range(yv), 4), collapse = " .. "),
        "  step=", if (length(yv) > 1) round(diff(sort(yv))[1], 4) else NA, "\n", sep = "")
    rm(input_data); gc(verbose = FALSE)
  }

  present <- intersect(wanted, all_cols)
  missing <- setdiff(wanted, all_cols)
  cat("\nrequested nodes PRESENT:", paste(present, collapse = ", "), "\n")
  cat("requested nodes MISSING:",
      if (length(missing)) paste(missing, collapse = ", ") else "(none)", "\n\n")

  ld <- load_chains(species, present)
  chain_full <- ld$chains

  # ---- R-hat restricted to the scientific params we report ----
  mcmc_list <- mcmc.list(lapply(chain_full, mcmc))
  gd <- gelman.diag(mcmc_list, autoburnin = FALSE, multivariate = FALSE)
  rhat <- gd$psrf[, 1]
  cat("\n===== R-hat over reported params =====\n")
  print(round(sort(rhat, decreasing = TRUE), 4))

  # trend_robust_indicator is a derived 0/1 step node -- excluded from the
  # convergence gate upstream, so exclude it from the headline max here too.
  rhat_gated <- rhat[!(names(rhat) %in% "trend_robust_indicator")]
  max_rhat <- suppressWarnings(max(rhat_gated, na.rm = TRUE))
  cat("MAX R-hat (reported params, excl. trend_robust_indicator):",
      round(max_rhat, 5), "\n")

  combined <- do.call(rbind, chain_full)
  cat("combined posterior draws:", nrow(combined), "\n")

  # ---- national-level table (all models) ----
  nat_present <- intersect(NATIONAL_NODES, colnames(combined))
  nat_df <- do.call(rbind, lapply(nat_present, function(p)
    summarize_node(combined[, p], p)))
  nat_df$rhat <- unname(rhat[nat_df$parameter])
  nat_df$model <- species

  # snr = year_beta/year_var is a deterministic node in the model code but is
  # NOT in the monitor list, so it cannot be read from samples -- derive it.
  # CAVEAT: year_var's posterior straddles 0, so this ratio is a ratio-of-
  # normals (Cauchy-like, undefined mean, exploding tails). Report median and
  # IQR only; the meaningful summary is trend_robust_indicator = P(snr > 1).
  if (all(c("year_beta", "year_var") %in% colnames(combined))) {
    snr_draws <- combined[, "year_beta"] / combined[, "year_var"]
    cat("\n----- derived snr = year_beta/year_var (NOT monitored) -----\n")
    cat("  median:", round(median(snr_draws), 4),
        "  IQR:", paste(round(q(snr_draws, c(0.25, 0.75)), 4), collapse = " .. "), "\n")
    cat("  P(snr > 1):", round(mean(snr_draws > 1), 4),
        "   [cross-check vs monitored trend_robust_indicator]\n")
    cat("  share of draws with |snr| > 100 (tail blowup):",
        round(mean(abs(snr_draws) > 100), 4), "\n")
    snr_row <- data.frame(
      parameter = "snr_derived", mean = NA_real_,
      median = median(snr_draws), sd = NA_real_,
      q025 = q(snr_draws, 0.25), q975 = q(snr_draws, 0.75),
      p_negative = mean(snr_draws < 0), n_na = sum(!is.finite(snr_draws)),
      rhat = NA_real_, model = species, stringsAsFactors = FALSE)
    nat_df <- rbind(nat_df, snr_row)
  }

  if (is_ecoregion) {
    sig <- summarize_node(combined[, "sigma_region"], "sigma_region")
    sig$rhat <- unname(rhat["sigma_region"]); sig$model <- species
    nat_df <- rbind(sig, nat_df)

    # year_region[r] is a mean-zero DEVIATION (year_region[r] ~ dnorm(0,
    # sigma_region)); the model applies it as
    # total_var_beta + year_effect[c], year_effect[c] = year_region[region].
    # So the ABSOLUTE regional iNat trend is total_var_beta + year_region[r].
    # The old (invalidated) CSV reported only the raw deviation -- for a
    # like-for-like comparison against agency absolute trend direction, the
    # absolute column is the one to map.
    abs_draws <- lapply(seq_len(nregion), function(r)
      combined[, "total_var_beta"] + combined[, paste0("year_region[", r, "]")])

    reg_df <- data.frame(
      region_index     = seq_len(nregion),
      region_name      = region_name,
      year_region_mean = sapply(seq_len(nregion), function(r) mean(combined[, paste0("year_region[", r, "]")])),
      year_region_median = sapply(seq_len(nregion), function(r) median(combined[, paste0("year_region[", r, "]")])),
      year_region_sd   = sapply(seq_len(nregion), function(r) sd(combined[, paste0("year_region[", r, "]")])),
      year_region_q025 = sapply(seq_len(nregion), function(r) q(combined[, paste0("year_region[", r, "]")], 0.025)),
      year_region_q975 = sapply(seq_len(nregion), function(r) q(combined[, paste0("year_region[", r, "]")], 0.975)),
      p_negative       = sapply(seq_len(nregion), function(r) mean(combined[, paste0("year_region[", r, "]")] < 0)),
      rhat             = sapply(seq_len(nregion), function(r) unname(rhat[paste0("year_region[", r, "]")])),
      abs_trend_mean   = sapply(abs_draws, mean),
      abs_trend_median = sapply(abs_draws, median),
      abs_trend_q025   = sapply(abs_draws, function(v) q(v, 0.025)),
      abs_trend_q975   = sapply(abs_draws, function(v) q(v, 0.975)),
      abs_p_negative   = sapply(abs_draws, function(v) mean(v < 0)),
      n_camera_sites   = n_camera_sites,
      n_inat_cells     = n_inat_cells,
      stringsAsFactors = FALSE
    )
    out_reg <- file.path(D, "moose_v1fix_ecoregion_posterior.csv")
    write.csv(reg_df, out_reg, row.names = FALSE)
    cat("\nwrote", out_reg, "\n")
    cat("\n===== ecoregion year_region =====\n")
    print(reg_df, row.names = FALSE, digits = 4)

    out_glob <- file.path(D, "moose_v1fix_ecoregion_global.csv")
    write.csv(nat_df, out_glob, row.names = FALSE)
    cat("\nwrote", out_glob, "\n")
  } else {
    out_nat <- file.path(D, "moose_v1fix_national_scalar_posterior.csv")
    write.csv(nat_df, out_nat, row.names = FALSE)
    cat("\nwrote", out_nat, "\n")
  }

  cat("\n===== national-level params (", species, ") =====\n")
  print(nat_df, row.names = FALSE, digits = 4)

  rm(combined, chain_full); gc(verbose = FALSE)
  invisible(NULL)
}

run_model("moose_v1fix_national_scalar", is_ecoregion = FALSE)
run_model("moose_v1fix_ecoregion",       is_ecoregion = TRUE)

cat("\nEXTRACTION DONE\n")
