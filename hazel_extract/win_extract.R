#!/usr/bin/env Rscript
# Posterior extraction (2026-11-01): summaries, convergence and year scaling for ONE run dir (env RUN).
# Chains (~1 GiB each; 1.2 TB windowed) cannot leave Hazel, so each run is reduced to small CSVs for the laptop.
# Chains are read ONE AT A TIME and cut to the scalar parameters before the next is read.
# Windowed runs:   <stem>_national_scalar_win{5,10}[_camonly], 12 single-shot chains chain_<RUN>_<k>_single.RDS
# Full-record runs: <stem>_national_scalar[_camonly], chain files given by env CHAIN_GLOB (default chain_<RUN>_*_single.RDS)
# Outputs (OUT, default ~/posteriors_20261101); posterior names put the window suffix LAST to match reportlib/inputs.py:
#   <stem>[_camonly]_posterior[_win5|_win10].csv   pooled summary over ALL chains (no chain is excluded)
#   convergence_<RUN>.csv                          per parameter: rank-normalised split R-hat, bulk/tail ESS, flags
#   perchain_<RUN>.csv                             per chain x parameter: mean, sd, that chain's bulk/tail ESS
#   year_scaling_<RUN>.csv                         year standardisation per submodel, from the input bundle
# camonly runs have no iNat likelihood: theta0, theta1, year_var, overdisp_inat (and total_var_beta, sigma_year_var,
# which are functions of / hyperpriors on year_var, plus trend_robust_indicator) are prior draws -> prior_only = TRUE.
# AFFECTED (env, TRUE/FALSE/NA) is copied into every table as affected_by_resume_defect; it is decided by the log check
# outside this script, never here.
suppressPackageStartupMessages(library(posterior))
H   <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/HPC"
OUT <- Sys.getenv("OUT", "/home/rwkays/posteriors_20261101"); dir.create(OUT, showWarnings = FALSE)
run <- Sys.getenv("RUN"); stopifnot(nzchar(run))
AFFECTED <- as.logical(Sys.getenv("AFFECTED", "FALSE"))
REQUIRED <- c("year_beta", "year_var", "total_var_beta", "theta0", "theta1")
PRIOR_ONLY_CAMONLY <- c("theta0", "theta1", "year_var", "overdisp_inat", "total_var_beta", "sigma_year_var",
                        "trend_robust_indicator")

m <- regmatches(run, regexec("^(.+_national_scalar)(_win[0-9]+)?(_camonly)?$", run))[[1]]
if (!length(m)) stop("run name does not parse as <stem>_national_scalar[_win<N>][_camonly]: ", run)
stem <- paste0(m[2], m[4]); wtag <- if (nzchar(m[3])) sub("^_", "", m[3]) else "full"; camonly <- nzchar(m[4])
species <- sub("_v2[a-z]*_national_scalar$", "", m[2]); mask <- sub("^.*_(v2[a-z]*)_national_scalar$", "\\1", m[2])
post_file <- file.path(OUT, sprintf("%s_posterior%s.csv", stem, if (wtag == "full") "" else paste0("_", wtag)))
cat("RUN:", run, " window:", wtag, " camonly:", camonly, " affected:", AFFECTED, "\nposterior ->", post_file, "\n")

# ---- chains, one at a time ---------------------------------------------------
glob <- Sys.getenv("CHAIN_GLOB", sprintf("chain_%s_*_single.RDS", run))
f <- Sys.glob(file.path(H, run, glob))
f <- f[order(as.integer(sub(".*_([0-9]+)_[^_]*\\.RDS$", "\\1", f)))]
if (wtag != "full" && length(f) != 12) cat("WARNING: windowed run has", length(f), "chain files, expected 12\n")
if (length(f) < 2) stop("fewer than 2 chain files matching ", glob)
keep <- NULL; draws <- list()
for (x in f) {
  s <- readRDS(x)$samples
  if (is.null(keep)) {
    scal <- colnames(s)[!grepl("\\[", colnames(s))]          # every non-indexed (scalar) monitored node
    miss <- setdiff(REQUIRED, scal)
    if (length(miss)) stop("required parameters not in chains: ", paste(miss, collapse = ", "),
                           "\nscalar columns present: ", paste(scal, collapse = ", "))
    keep <- c(REQUIRED, setdiff(scal, REQUIRED))
  }
  stopifnot(all(keep %in% colnames(s)))
  draws[[length(draws) + 1]] <- s[, keep, drop = FALSE]
  cat(sprintf("%s: %d x %d read, kept %d\n", basename(x), nrow(s), ncol(s), length(keep)))
  rm(s); invisible(gc())
}
n_it <- vapply(draws, nrow, 1L)
if (length(unique(n_it)) != 1) stop("chains differ in length: ", paste(n_it, collapse = "/"))
arr <- array(NA_real_, dim = c(n_it[1], length(draws), length(keep)), dimnames = list(NULL, NULL, keep))
for (k in seq_along(draws)) arr[, k, ] <- draws[[k]]
rm(draws); invisible(gc())
n_chain <- dim(arr)[2]; n_draw <- prod(dim(arr)[1:2])
chain_ids <- as.integer(sub(".*_([0-9]+)_[^_]*\\.RDS$", "\\1", f))
prior_only <- camonly & keep %in% PRIOR_ONLY_CAMONLY
lab <- data.frame(run = run, species = species, mask = mask, window = wtag, camonly = camonly,
                  affected_by_resume_defect = AFFECTED)

# ---- A. pooled posterior summary (every chain in) ------------------------------
probs <- c(.025, .25, .5, .75, .975)
post <- do.call(rbind, lapply(keep, function(p) {
  v <- as.vector(arr[, , p]); q <- quantile(v, probs, names = FALSE)
  data.frame(lab, param = p, prior_only = prior_only[keep == p], mean = mean(v), median = median(v), sd = sd(v),
             q2.5 = q[1], q25 = q[2], q50 = q[3], q75 = q[4], q97.5 = q[5], n_chains = n_chain, n_draws = n_draw)
}))
write.csv(post, post_file, row.names = FALSE)   # write.csv keeps 15 significant digits; nothing is rounded here

# ---- B. convergence + per-chain ------------------------------------------------
safe <- function(fun, x) if (var(as.vector(x)) == 0) NA_real_ else fun(x)
conv <- do.call(rbind, lapply(keep, function(p) {
  x <- arr[, , p, drop = TRUE]; rh <- safe(rhat, x); eb <- safe(ess_bulk, x); et <- safe(ess_tail, x)
  data.frame(lab, param = p, prior_only = prior_only[keep == p], rhat = rh, ess_bulk = eb, ess_tail = et,
             n_chains = n_chain, n_draws = n_draw, constant = var(as.vector(x)) == 0,
             flag_rhat = isTRUE(rh > 1.01), flag_ess = isTRUE(min(eb, et) < 400))
}))
conv$flag <- conv$flag_rhat | conv$flag_ess
write.csv(conv, file.path(OUT, sprintf("convergence_%s.csv", run)), row.names = FALSE)
pc <- do.call(rbind, lapply(keep, function(p) do.call(rbind, lapply(seq_len(n_chain), function(k) {
  v <- arr[, k, p]; cst <- var(v) == 0
  data.frame(lab, param = p, prior_only = prior_only[keep == p], chain = chain_ids[k], n_iter = length(v),
             mean = mean(v), sd = sd(v), ess_bulk = if (cst) NA_real_ else ess_bulk(v),
             ess_tail = if (cst) NA_real_ else ess_tail(v))
}))))
write.csv(pc, file.path(OUT, sprintf("perchain_%s.csv", run)), row.names = FALSE)

# ---- C. year scaling, from the input bundle --------------------------------------
b  <- readRDS(file.path(H, run, sprintf("input_data_%s.RDS", run)))
cl <- b$constants_list
if (wtag != "full") {
  if (is.null(b$window$years)) stop("windowed run but bundle has no $window$years")
  yrs <- b$window$years; w0 <- b$window$year_start; w1 <- b$window$year_end
} else {
  yrs <- if (!is.null(b$window$years)) b$window$years else seq(2008, 2025)   # run_fleet_prep.R full-record years_in
  w0 <- min(yrs); w1 <- max(yrs)
}
ycam <- as.numeric(as.character(b$real_data$occ.covs$year))    # camera deployment years (year_occ = scale(ycam))
# occ.covs columns that were factors became integer codes in prep (deployment_id/project_id trap); a code vector would
# still pass the scale() check below, so require real calendar years inside the window.
if (!all(ycam >= w0 & ycam <= w1)) stop("occ.covs$year is not calendar years inside ", w0, "-", w1,
                                        ": range ", paste(range(ycam), collapse = "-"))
chk <- function(raw, std) if (length(raw) != length(std)) NA_real_ else max(abs((raw - mean(raw)) / sd(raw) - std))
ys <- data.frame(
  lab, submodel = c("camera_occupancy", "inat_intensity"),
  covariate = c("year_occ", "year_vals"), coefficient = c("year_beta", "total_var_beta"),
  standardised_over = c("camera deployment rows (one value per site)", "calendar years in window (one value per year)"),
  n_values = c(length(ycam), length(yrs)),
  year_mean = c(mean(ycam), mean(yrs)), year_sd = c(sd(ycam), sd(yrs)),
  data_year_min = c(min(ycam), min(yrs)), data_year_max = c(max(ycam), max(yrs)),
  window_start = w0, window_end = w1,
  max_abs_check = c(chk(ycam, cl$year_occ), chk(yrs, cl$year_vals)))
if (any(is.na(ys$max_abs_check) | ys$max_abs_check > 1e-8))
  stop("bundle year covariates are not scale() of the raw years:\n", paste(capture.output(print(ys)), collapse = "\n"))
write.csv(ys, file.path(OUT, sprintf("year_scaling_%s.csv", run)), row.names = FALSE)

options(width = 220)
print(conv[, c("param", "prior_only", "rhat", "ess_bulk", "ess_tail")], row.names = FALSE, digits = 5)
cat("DONE\n")
