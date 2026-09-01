# =============================================================================
# RECOVERED FILE -- provenance note added 2026-09-01
# -----------------------------------------------------------------------------
# This file existed ONLY inside the untracked `isdm_sim_codebase/` snapshot
# bundle and was never in version control. It is recovered here so the sweep
# collector is tracked alongside `01i_run_estimator_sweep.R`, which it consumes.
#
# TWO CAUTIONS before trusting this copy:
#
# 1. STALENESS. The bundle is an August snapshot. For at least one sibling file
#    (`01i_camera_rn_sbatch.sh`) the bundle copy is an OLDER, broken draft and
#    the repo-root copy is correct. Verify this collector against whatever is
#    live on Hazel before treating it as authoritative.
#
# 2. KNOWN DEFECT -- null-scenario mislabelling. Around the `scenario == "null"`
#    block below, this script asserts that `detect_rate` under the null arm is
#    a FALSE-POSITIVE rate. That is incorrect: the null scenario zeroed only the
#    spatial deviation amplitude, never `year_beta`/`year_var`, so the true
#    national trend was present in BOTH scenarios. The null column is therefore
#    power-under-a-spatially-flat-truth, NOT Type I error. The earlier
#    "false-positive > power" finding derived from it was retracted. A corrected
#    sweep (with the trend genuinely zeroed) was requested; if that fix has
#    landed on Hazel, THAT version supersedes this file.
# =============================================================================

#!/usr/bin/env Rscript
# =============================================================================
# 01j_collect_estimator_sweep.R
# -----------------------------------------------------------------------------
# Collector for the estimator sweep. Stacks estimator_sweep_out/row_*.rds and
# runs summarize_estimator_sweep() -- which is where the degeneracy guard
# lives, so this is also the reseed check.
#
# The guard is deliberately called on the RAW stack before anything else is
# printed: a summary table that got as far as being formatted is a summary
# somebody will read, and the whole point of the guard is that a degenerate
# sweep must not produce a readable table.
# =============================================================================
BASE <- Sys.getenv("SIM_BASE", ".")
setwd(BASE)
source("sim_helpers_estimator_metrics.R")

# Must match the driver's SWEEP_OUTDIR for the run being collected -- the 3-rep
# pilot and the 30-rep sweep write to different directories on purpose (see the
# OUTDIR note in 01i_run_estimator_sweep.R).
OUTDIR <- Sys.getenv("SWEEP_OUTDIR", file.path(BASE, "estimator_sweep_out"))
TAG    <- Sys.getenv("SWEEP_TAG", "pilot")
cat("OUTDIR:", OUTDIR, "  tag:", TAG, "\n")
files  <- sort(list.files(OUTDIR, pattern = "^row_\\d+\\.rds$", full.names = TRUE))
cat("row files found:", length(files), "\n")

rows <- do.call(rbind, lapply(files, function(f) {
  r <- readRDS(f)
  if (is.list(r) && !is.data.frame(r)) r <- as.data.frame(r, stringsAsFactors = FALSE)
  r
}))
cat("stacked rows:", nrow(rows), " cols:", paste(colnames(rows), collapse = ", "), "\n\n")

cat("===== STATUS =====\n")
print(table(rows$status, useNA = "ifany"))
cat("\n")
if (any(rows$status != "OK")) {
  cat("non-OK rows:\n")
  print(rows[rows$status != "OK",
             intersect(c("row_id","estimator","abundance","scenario","rep_id","status","message"),
                       colnames(rows))], row.names = FALSE)
  cat("\n")
}

cat("===== DESIGN COVERAGE =====\n")
print(with(rows[rows$status == "OK", , drop = FALSE],
           table(estimator, abundance, scenario)))
cat("\n")

# ---- reseed / degeneracy check, explicitly ---------------------------------
# summarize_estimator_sweep() stops on degeneracy, but it only inspects
# tvb_mean. Check every numeric metric column: a reseed failure makes the
# ENTIRE replicate byte-identical, so identical tvb_mean with differing
# elapsed_sec would mean something else is wrong.
cat("===== RESEED / DEGENERACY CHECK =====\n")
ok  <- rows[rows$status == "OK", , drop = FALSE]
key <- paste(ok$estimator, ok$abundance, ok$scenario, sep = "|")
num_cols <- names(which(vapply(ok, is.numeric, logical(1))))
num_cols <- setdiff(num_cols, c("row_id", "rep_id", "tvb_true"))
chk <- do.call(rbind, lapply(sort(unique(key)), function(k) {
  sub <- ok[key == k, num_cols, drop = FALSE]
  n_ident <- sum(vapply(sub, function(v) length(unique(round(v, 12))) == 1L, logical(1)))
  data.frame(cell = k, n_rep = nrow(sub),
             n_metrics = length(num_cols),
             n_identical_across_reps = n_ident,
             tvb_mean_distinct = length(unique(round(sub$tvb_mean, 12))),
             stringsAsFactors = FALSE)
}))
print(chk, row.names = FALSE)
cat("\ntvb_mean_distinct should equal n_rep in every cell. n_identical_across_reps\n",
    "counts metric columns that are constant within a cell -- a handful can be\n",
    "legitimate (tvb_true is excluded; binary flags like tvb_detected can\n",
    "genuinely agree), but ALL of them constant is the reseed bug.\n\n", sep = "")

# ---- the summary itself -----------------------------------------------------
cat("===== summarize_estimator_sweep() =====\n")
smry <- summarize_estimator_sweep(rows)
smry <- smry[order(smry$scenario, smry$abundance, smry$estimator), ]
print(smry, row.names = FALSE, digits = 4)
cat("\n")

write.csv(rows, file.path(BASE, sprintf("estimator_sweep_%s_rows.csv", TAG)), row.names = FALSE)
write.csv(smry, file.path(BASE, sprintf("estimator_sweep_%s_summary.csv", TAG)), row.names = FALSE)
cat(sprintf("wrote estimator_sweep_%s_rows.csv and estimator_sweep_%s_summary.csv\n\n", TAG, TAG))

# ---- against the stated baseline -------------------------------------------
# Baseline to beat: 66% sign recovery, 70% false-negative rate. The false
# negative rate is 1 - detect_rate under scenario == "varying"; the null arm
# gives the false-POSITIVE rate, which is a separate number and must not be
# read as the same thing.
cat("===== VS BASELINE (66% sign recovery / 70% false negative) =====\n")
v <- smry[smry$scenario == "varying", ]
if (nrow(v)) {
  v$false_negative <- 1 - v$detect_rate
  print(v[, c("estimator","abundance","n_rep","detect_rate","false_negative",
              "bias","coverage","ci_width","tri_rate","auc","spearman")],
        row.names = FALSE, digits = 4)
}
n <- smry[smry$scenario == "null", ]
if (nrow(n)) {
  cat("\nnull scenario (detect_rate here is the FALSE-POSITIVE rate):\n")
  print(n[, c("estimator","abundance","n_rep","detect_rate","bias","coverage",
              "tri_rate","auc")], row.names = FALSE, digits = 4)
}
nr <- if (nrow(smry)) max(smry$n_rep) else 0
cat(sprintf("\nNOTE: n_rep = %d per cell. ", nr))
if (nr < 10) cat("A rate from so few draws has a 95% interval of roughly\n",
  "+/- 0.5 -- these numbers say whether the machinery runs and roughly where it\n",
  "lands, NOT whether any arm beats the baseline.\n", sep = "") else
  cat(sprintf("A rate from %d draws has a 95%% interval of about +/- %.2f;\n", nr, 1.96*0.5/sqrt(nr)),
      "differences smaller than that are not resolvable by this design.\n", sep="")
cat("\nDONE\n")
