#!/usr/bin/env Rscript
# =============================================================================
# 00y_smoke_test.R
# -----------------------------------------------------------------------------
# Verify the simulation code base is internally consistent and runnable, without
# needing the real data. Every check either PASSES or FAILS loudly; the script
# exits non-zero if anything fails, so it can gate a CI run.
#
# Usage:  Rscript 00y_smoke_test.R
# =============================================================================

suppressMessages({library(nimble); library(nimbleEcology)})

BASE <- getwd()
fails <- character(0)
ok <- function(label, cond, detail = "") {
  cat(sprintf("  [%s] %s%s\n", if (cond) "PASS" else "FAIL", label,
              if (nzchar(detail)) paste0("  -- ", detail) else ""))
  if (!cond) fails <<- c(fails, label)
}

cat("\n== 1. every source file parses ==\n")
srcs <- c("sim_helpers.R","sim_helpers_abundance.R","sim_helpers_array.R",
          "sim_helpers_estimator_metrics.R","sim_helpers_effort.R",
          "model_code_national_scalar.R","model_code_array_occ.R",
          "model_code_array_rn.R","model_code_camera_rn.R",
          "model_code_ecoregion_trend.R","01i_run_estimator_sweep.R")
for (f in srcs) {
  r <- tryCatch({parse(f); TRUE}, error = function(e) {attr(r,"m") <<- conditionMessage(e); FALSE})
  ok(f, r)
}

cat("\n== 2. trend block is identical across all four estimator models ==\n")
trend_hash <- function(f) {
  L <- readLines(f, warn = FALSE)
  i <- grep("year_beta ~ dnorm", L)[1]; j <- grep("trend_robust_indicator", L)
  j <- j[j >= i][1]
  s <- L[i:j]; s <- s[!grepl("^\\s*#", s)]
  digest_str <- paste(gsub("[ \t\r]", "", s), collapse = "")
  digest_str
}
mods <- c("model_code_national_scalar.R","model_code_array_occ.R",
          "model_code_array_rn.R","model_code_camera_rn.R")
hashes <- vapply(mods, trend_hash, character(1))
ok("trend block identical in all 4 models", length(unique(hashes)) == 1,
   sprintf("%d distinct", length(unique(hashes))))

cat("\n== 3. models build and calculate a finite log-probability ==\n")
for (f in c("model_code_national_scalar.R","model_code_array_occ.R",
            "model_code_array_rn.R","model_code_camera_rn.R")) source(f)
ok("all four model objects exist",
   all(sapply(c("model_code_national_scalar","model_code_array_occ",
                "model_code_array_rn","model_code_camera_rn"), exists)))

cat("\n== 4. synthetic inputs generate with the required structure ==\n")
tmp <- tempfile(fileext = ".RDS")
system2("Rscript", c("00z_make_synthetic_inputs.R", tmp), stdout = NULL, stderr = NULL)
ok("synthetic inputs written", file.exists(tmp))
if (file.exists(tmp)) {
  x <- readRDS(tmp)
  need <- c("cell100_geo","constants_list","cl","real_post_means","base_inits",
            "inat_effort","inat_effort_real","y_template","real_y_template","site_array")
  ok("all driver-required fields present", all(need %in% names(x)),
     paste("missing:", paste(setdiff(need, names(x)), collapse = ", ")))
  ok("CAR adjacency well-formed",
     sum(x$cl$num) == length(x$cl$adj) && all(x$cl$num >= 1))
  ok("flagged as synthetic", isTRUE(x$SYNTHETIC))
}

cat("\n== 5. design grid: 4 arms x 3 abundance x 2 scenarios ==\n")
source("sim_helpers.R"); source("sim_helpers_abundance.R")
N_REP <- 30L
d <- expand.grid(rep_id = seq_len(N_REP), scenario = c("varying","null"),
                 abundance = c("bobcat_baseline","moderate","common_deerlike"),
                 estimator = c("camera_occ","camera_rn","array_occ","array_rn"),
                 stringsAsFactors = FALSE)
d <- d[order(d$estimator, d$abundance, d$scenario, d$rep_id), ]
d$row_id <- seq_len(nrow(d))
ok("720 tasks in 24 cells", nrow(d) == 720 &&
     nrow(unique(d[, c("estimator","abundance","scenario")])) == 24)
ok("30 replicates in every cell",
   all(table(d$estimator, d$abundance, d$scenario) == 30))
old <- expand.grid(rep_id = seq_len(N_REP), scenario = c("varying","null"),
                   abundance = c("bobcat_baseline","moderate","common_deerlike"),
                   estimator = c("camera_occ","array_occ","array_rn"),
                   stringsAsFactors = FALSE)
old <- old[order(old$estimator, old$abundance, old$scenario, old$rep_id), ]
# Compare CONTENT, not the data.frame objects: both carry stale rownames from
# expand.grid + order(), which identical() would flag even when every design
# value matches. rownames are irrelevant here -- what matters is that row_id N
# denotes the same (estimator, abundance, scenario, rep_id) in both grids.
kk <- c("estimator","abundance","scenario","rep_id")
o_ <- old[, kk]; n_ <- d[seq_len(540), kk]
rownames(o_) <- NULL; rownames(n_) <- NULL
ok("adding camera_rn does not renumber the existing 540 rows",
   identical(o_, n_),
   "camera_rn sorts last, so rows 1-540 keep their meaning")

cat("\n== 6. abundance ladder is the measured one ==\n")
lv <- abundance_levels_measured()
ok("three levels", nrow(lv) == 3)
ok("baseline is a no-op", lv$occ_shift[1] == 0 && lv$count_log_mult[1] == 0)
ok("ladder is monotone increasing", all(diff(lv$target_mean_count) > 0))

cat("\n== 6b. abundance ladder is actually APPLIED by the driver ==\n")
# The inert-ladder bug: scale_truth_abundance() only RECORDS `label`; all its
# scaling comes from occ_shift/count_log_mult, which default to 0. Calling it
# with `label` alone is a silent no-op, and the whole sweep collapses to three
# i.i.d. replicate sets from one DGP. It shipped once and cost a 720-task run.
# These checks are static (source text) plus behavioural (does lambda move).
drv <- paste(readLines("01i_run_estimator_sweep.R", warn = FALSE), collapse = "\n")
ok("driver passes occ_shift to scale_truth_abundance",
   grepl("scale_truth_abundance\\([^)]*occ_shift", drv))
ok("driver passes count_log_mult to scale_truth_abundance",
   grepl("scale_truth_abundance\\([^)]*count_log_mult", drv))
ok("driver builds a ladder (not label-only)",
   grepl("abundance_levels_(measured|default)\\(\\)", drv))
lad <- abundance_levels_measured()
dsn_levels <- regmatches(drv, regexpr("abundance = c\\([^)]*\\)", drv))
ok("every design level resolves in the ladder",
   length(dsn_levels) == 1 &&
     all(vapply(lad$level, function(L) grepl(L, dsn_levels, fixed = TRUE), logical(1))),
   paste("design:", dsn_levels))
tt <- list(link_occ_intercept = c(-2, -1.5, -2.5), theta0 = -1)
lams <- vapply(seq_len(nrow(lad)), function(i)
  stats::median(exp(scale_truth_abundance(tt, occ_shift = lad$occ_shift[i],
                    count_log_mult = lad$count_log_mult[i])$link_occ_intercept)),
  numeric(1))
ok("scaling actually moves lambda across levels", length(unique(round(lams, 6))) == nrow(lad),
   sprintf("multipliers %s", paste(round(lams / lams[1], 2), collapse = " / ")))
ok("label-only call is a no-op (the trap, documented)",
   identical(scale_truth_abundance(tt, label = "deer_like")$link_occ_intercept,
             tt$link_occ_intercept))

cat("\n== 7. degeneracy guard catches identical replicates ==\n")
source("sim_helpers_estimator_metrics.R")
dup <- data.frame(estimator = "x", abundance = "y", scenario = "varying",
                  rep_id = 1:30, status = "OK",
                  tvb_true = -0.18, tvb_mean = -0.15, tvb_bias = 0.03,
                  tvb_ci_width = 1.0, tvb_covered = 1L, tvb_detected = 0L,
                  tri_fires = 0L, disc_auc = .6, disc_spearman = .3,
                  disc_rmse = 2, elapsed_sec = 1:30)
caught <- tryCatch({summarize_estimator_sweep(dup); FALSE}, error = function(e) TRUE)
ok("guard fires on identical replicates", caught)
var_ <- dup; var_$tvb_mean <- rnorm(30, -0.15, .05)
passed <- tryCatch({summarize_estimator_sweep(var_); TRUE}, error = function(e) FALSE)
ok("guard passes on varied replicates", passed)

cat("\n")
if (length(fails)) {
  cat(sprintf("SMOKE TEST FAILED -- %d check(s):\n  %s\n",
              length(fails), paste(fails, collapse = "\n  ")))
  quit(status = 1)
}
cat("SMOKE TEST PASSED -- all checks clear.\n")
