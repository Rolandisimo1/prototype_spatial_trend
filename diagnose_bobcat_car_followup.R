#!/usr/bin/env Rscript
# Follow-up to diagnose_bobcat_v2b_ecoregion_car.R.
#
# The first pass REFUTED both stated hypotheses: the 12 offending cells are
# data-RICH (median 25 camera sites vs 0 for the rest) and INTERIOR (degree 4 =
# field median; fewer empty neighbours than average), not sparse and not on a
# mask edge. Two things were left open:
#   1. n_cam_det came back NA -- the detection matrix is real_data$y, not
#      data$y. Site COUNT is not detection count; closing that properly.
#   2. The signal that did separate them was low within-cell MCMT spread.
#      That is the signature of an additive-identifiability problem, so test
#      it directly rather than by proxy.
#
# MECHANISM UNDER TEST. The occupancy linear predictor contains
#     link_occ_intercept[cell[i]] + MCMT[i] * MCMT_effect[cell[i]]
# Both terms are indexed by the SAME cell. If MCMT is near-constant across the
# sites within a cell, then MCMT[i]*MCMT_effect[c] is itself near-constant
# within that cell and is not separable from the intercept -- the two trade off
# along a ridge. That predicts exactly what is observed: a data-rich block
# (enough information to support both surfaces) with low within-cell MCMT
# spread (not enough to tell them apart), and MCMT_tau -- the precision
# governing the competing field -- as the single worst R-hat in the fit.
PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
D <- file.path(PROJ, "HPC", "bobcat_v2b_ecoregion")
OFF <- c(854,855,873,874,875,883,884,885,891,892,893,894)

d  <- readRDS(file.path(D, "input_data_bobcat_v2b_ecoregion.RDS"))
cl <- d$constants_list
nc <- cl$ncell100
y  <- d$real_data$y
stopifnot(!is.null(y), nrow(y) == length(cl$cell))

f <- factor(cl$cell, levels = seq_len(nc))
site_det <- rowSums(y, na.rm = TRUE)
site_occ <- rowSums(!is.na(y))

n_sites <- as.integer(table(f))
n_det   <- as.integer(tapply(site_det, f, sum, default = 0)); n_det[is.na(n_det)] <- 0L
n_occ   <- as.integer(tapply(site_occ, f, sum, default = 0)); n_occ[is.na(n_occ)] <- 0L
det_rate <- ifelse(n_occ > 0, n_det / n_occ, NA_real_)

wsd <- function(v) tapply(v, f, function(x) if (length(x) >= 2) sd(x) else NA_real_)
mcmt_sd <- wsd(cl$MCMT); mwmt_sd <- wsd(cl$MWMT)

rest <- setdiff(seq_len(nc), OFF)
has  <- function(idx) idx[n_sites[idx] >= 2]        # comparable subset
pctl <- function(v, x) round(100 * mean(v <= x, na.rm = TRUE), 1)

cat("=== 1. DATA SPARSITY, using real detections ===\n")
cat(sprintf("%-26s %12s %12s\n", "", "offenders", "rest"))
row <- function(lab, v, r = rest) cat(sprintf("%-26s %12.3f %12.3f\n", lab,
        median(v[OFF], na.rm = TRUE), median(v[r], na.rm = TRUE)))
row("median n_cam_sites",   as.numeric(n_sites))
row("median n_detections",  as.numeric(n_det))
row("median n_occasions",   as.numeric(n_occ))
row("median naive det rate", det_rate, has(rest))
cat(sprintf("\ncells with ZERO detections: offenders %d/12, rest %d/%d (%.1f%%)\n",
            sum(n_det[OFF] == 0), sum(n_det[rest] == 0), length(rest),
            100 * mean(n_det[rest] == 0)))
cat("offender detection counts:", paste(n_det[OFF], collapse = ", "), "\n")
cat("percentile of each offender in the all-cell detection distribution:\n  ",
    paste(sprintf("%d:%.0f%%", OFF, vapply(OFF, function(i) pctl(n_det, n_det[i]), numeric(1))),
          collapse = "  "), "\n")

cat("\n=== 2. WITHIN-CELL MCMT SPREAD (the identifiability test) ===\n")
cmp <- has(rest)
cat(sprintf("offenders median within-cell MCMT sd : %.5f (n=%d cells)\n",
            median(mcmt_sd[OFF], na.rm = TRUE), sum(!is.na(mcmt_sd[OFF]))))
cat(sprintf("rest      median within-cell MCMT sd : %.5f (n=%d cells)\n",
            median(mcmt_sd[cmp], na.rm = TRUE), sum(!is.na(mcmt_sd[cmp]))))
cat(sprintf("offenders median within-cell MWMT sd : %.5f\n", median(mwmt_sd[OFF], na.rm = TRUE)))
cat(sprintf("rest      median within-cell MWMT sd : %.5f\n", median(mwmt_sd[cmp], na.rm = TRUE)))
cat("\nper-offender MCMT sd and its percentile among cells with >=2 sites:\n")
for (i in OFF) cat(sprintf("  cell %d: n=%4d  MCMT sd=%.5f  pctl=%.0f%%  MWMT sd=%.5f\n",
                           i, n_sites[i], mcmt_sd[i], pctl(mcmt_sd[cmp], mcmt_sd[i]), mwmt_sd[i]))
w <- suppressWarnings(wilcox.test(mcmt_sd[OFF], mcmt_sd[cmp]))
cat(sprintf("\nWilcoxon offenders vs rest (within-cell MCMT sd): W=%.0f, p=%.4g\n",
            w$statistic, w$p.value))

cat("\n=== 3. THE COMBINATION ===\n")
cat("Cells that are BOTH data-rich (>= median offender sites) AND MCMT-flat\n")
cat("(<= median offender MCMT sd) -- i.e. share the proposed signature:\n")
thr_n <- median(n_sites[OFF]); thr_s <- median(mcmt_sd[OFF], na.rm = TRUE)
sig <- which(n_sites >= thr_n & mcmt_sd <= thr_s)
cat(sprintf("  %d cells field-wide; %d of them (%.0f%%) are among the 12 offenders\n",
            length(sig), length(intersect(sig, OFF)),
            100 * length(intersect(sig, OFF)) / max(1, length(sig))))
cat("  the signature cells:", paste(sort(sig), collapse = ", "), "\n")
cat("\n(If this set is small and mostly offenders, the signature is specific.\n",
    " If it is large, data-rich + MCMT-flat is common and something else\n",
    " distinguishes the 12.)\n", sep = "")
cat("\nFOLLOWUP DONE\n")
