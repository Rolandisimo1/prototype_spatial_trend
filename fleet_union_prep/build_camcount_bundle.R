#!/usr/bin/env Rscript
# =============================================================================
# build_camcount_bundle.R -- add the count-camera fields to a v2b bundle.
#
# Input : HPC/<SRC_TOKEN>/input_data_<SRC_TOKEN>.RDS  (converged occupancy
#         arm's bundle, 9-cov, presence mask -- read only, never modified)
# Output: HPC/<OUT_TOKEN>/input_data_<OUT_TOKEN>.RDS  = the same bundle plus
#         constants_list$log_effort_cam, constants_list$yday_mean,
#         y_cam (top level, data), inits_list$overdisp_cam, and without
#         inits_list$det_intercept (removed from model_code_national_scalar_camcount).
#
# y_cam[i]: INDEPENDENT detection events of the species at deployment i.
#   1. sequences file rows -> one row per (uid, sequence_id): the table carries
#      one row per species x age x sex class, so sequence_id is not a key.
#   2. keep sequences whose start date lies in [start_date, end_date], the
#      window survey_nights describes.
#   3. 30-MINUTE INDEPENDENCE FILTER, same rule as 00d_prep_rn_counts.R: walk
#      each deployment's sequences in start_time order and keep one only if it
#      starts >= 30 min after the last KEPT event. Sequences in this file are
#      split by a 60-second gap, not an independence interval (measured
#      2026-09-16: 30.1% of consecutive deer sequences are < 30 min apart).
# log_effort_cam[i] = log(survey_nights[i]) from combined_deployments_all.csv.
# yday_mean[i] = mean of the row's yday_scaled windows: the same linear
#   standardisation as the occupancy arm's yday, averaged over the deployment.
#
# uid = paste(project_id, deployment_id, sep = "___") -- the key 01_Prep_UMFs.R
# used to build the y rownames. Row order of the bundle is never changed.
# =============================================================================
suppressPackageStartupMessages(library(data.table))

PROJ  <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
RAW   <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data"
SRC   <- Sys.getenv("SRC_TOKEN", "white-tailed_deer_v2b_national_scalar")
OUT   <- Sys.getenv("OUT_TOKEN", "white-tailed_deer_v2b_camcount_national_scalar")
SCI   <- Sys.getenv("SCINAME",   "Odocoileus virginianus")
GAP_S <- as.numeric(Sys.getenv("INDEP_MIN", "30")) * 60

in_f  <- file.path(PROJ, "HPC", SRC, paste0("input_data_", SRC, ".RDS"))
out_d <- file.path(PROJ, "HPC", OUT); out_f <- file.path(out_d, paste0("input_data_", OUT, ".RDS"))
if (file.exists(out_f)) stop("refusing to overwrite ", out_f)
b  <- readRDS(in_f); cl <- b$constants_list
uid_in <- rownames(b$real_data$y)   # never call this `uid`: inside dt[...] that name is the column
stopifnot(!is.null(uid_in), !anyDuplicated(uid_in), length(uid_in) == cl$nsite)
cat("source bundle:", in_f, "\n  nsite", cl$nsite, " ncell50", cl$ncell50,
    " numOccCovars", cl$numOccCovars, "\n")
stopifnot(cl$numOccCovars == 9)

# ---- deployments / effort ----------------------------------------------------
dep <- fread(file.path(RAW, "combined_deployments_all.csv"), colClasses = list(character = c("project_id", "deployment_id")))
dep <- unique(dep)
dep[, uid := paste(project_id, deployment_id, sep = "___")]
dup <- dep[, .N, by = uid][N > 1]
cat("deployments rows:", nrow(dep), " duplicated uid (after distinct):", nrow(dup), "\n")
dep_b <- dep[uid %in% uid_in]
cat("bundle uids found in deployments file:", uniqueN(dep_b$uid), "of", length(uid_in), "\n")
if (nrow(dup[uid %in% uid_in])) {
  cat("  bundle uids with >1 deployment row:", nrow(dup[uid %in% uid_in]), " -- effort ambiguous\n")
  print(head(dep[uid %in% dup[uid %in% uid_in, uid]], 10))
}
stopifnot("every bundle deployment must have exactly one effort row" =
            uniqueN(dep_b$uid) == length(uid_in) && nrow(dep_b) == length(uid_in))
dep_b <- dep_b[match(uid_in, dep_b$uid)]
stopifnot(identical(dep_b$uid, uid_in))
sn <- dep_b$survey_nights
cat("survey_nights: min", min(sn), " median", median(sn), " max", max(sn),
    " n<=0:", sum(sn <= 0), " NA:", sum(is.na(sn)), "\n")
stopifnot("survey_nights must be positive and finite for log offset" = all(is.finite(sn) & sn > 0))
dep_b[, `:=`(sd = as.Date(substr(start_date, 1, 10)), ed = as.Date(substr(end_date, 1, 10)))]
span <- as.integer(dep_b$ed - dep_b$sd)
cat("survey_nights - (end-start): quantiles", paste(quantile(sn - span, c(0, .01, .5, .99, 1)), collapse = " "), "\n")
win_days <- 10 * cl$J
cat("windowed occupancy days (10*J) / survey_nights: median", round(median(win_days / sn), 3),
    " share of deployments with survey_nights > 10*J + 12:", round(mean(sn > win_days + 12), 3), "\n")

# ---- sequences -> independent events ---------------------------------------
sq <- fread(file.path(RAW, "combined_sequences_all.csv"),
            select = c("project_id", "deployment_id", "sequence_id", "genus", "species", "start_time"),
            colClasses = list(character = c("project_id", "deployment_id", "sequence_id")))
sq <- sq[paste(genus, species) == SCI]
cat("\n", SCI, "sequence-class rows:", nrow(sq), "\n")
sq[, uid := paste(project_id, deployment_id, sep = "___")]
sq <- unique(sq, by = c("uid", "sequence_id"))
cat("distinct (uid, sequence):", nrow(sq), "\n")
sq <- sq[uid %in% uid_in]
cat("  in bundle deployments:", nrow(sq), "\n")
sq[, ts := as.POSIXct(start_time, format = "%Y-%m-%dT%H:%M:%OS", tz = "UTC")]
# 31 deer rows in the whole file have an empty start_time; they cannot be
# placed in the deployment window or thinned, so they are dropped and counted.
# Anything unparseable that is NOT blank is a format problem and stops.
stopifnot("non-blank start_time failed to parse" = !any(is.na(sq$ts) & nzchar(sq$start_time) & !is.na(sq$start_time)))
cat("  blank start_time (dropped):", sum(is.na(sq$ts)), "\n")
sq <- sq[!is.na(ts)]
sq[dep_b, on = "uid", `:=`(sd = i.sd, ed = i.ed)]
inside <- as.Date(sq$ts) >= sq$sd & as.Date(sq$ts) <= sq$ed
cat("  outside [start_date, end_date] (dropped):", sum(!inside), "\n")
sq <- sq[inside]
setorder(sq, uid, ts)
thin <- function(t) {                     # keep if >= GAP_S after last KEPT event
  t <- as.numeric(t); keep <- logical(length(t)); last <- -Inf
  for (k in seq_along(t)) if (t[k] - last >= GAP_S) { keep[k] <- TRUE; last <- t[k] }
  keep
}
sq[, keep := thin(ts), by = uid]
cat("  30-min filter: kept", sum(sq$keep), "of", nrow(sq), sprintf("(removed %.1f%%)\n", 100 * mean(!sq$keep)))
ev <- sq[keep == TRUE, .N, by = uid]
y_cam <- ev$N[match(uid_in, ev$uid)]; y_cam[is.na(y_cam)] <- 0L
y_cam <- as.integer(y_cam)

# ---- consistency with the occupancy arm's y --------------------------------
occ_any <- rowSums(b$real_data$y, na.rm = TRUE) > 0
tab <- table(occ_detect = occ_any, count_pos = y_cam > 0)
cat("\noccupancy detection vs count>0:\n"); print(tab)
cat("  occ-detected but zero count (must be 0):", tab["TRUE", "FALSE"], "\n")
stopifnot("an occupancy detection has no counted event -- join or window logic wrong" = tab["TRUE", "FALSE"] == 0)
cat("  count>0 but no occupancy detection (events outside full 10-day windows):", tab["FALSE", "TRUE"], "\n")

rate <- y_cam / sn
det  <- y_cam > 0
cat("\ny_cam: detecting", sum(det), "of", length(y_cam), sprintf("(%.1f%%)", 100 * mean(det)),
    " total events", sum(y_cam), " max", max(y_cam), "\n")
cat("  variance/mean:", round(var(y_cam) / mean(y_cam), 1), "\n")
q <- quantile(rate[det], c(.1, .9))
cat("  rate among detecting, p10 / p90:", signif(q[1], 3), "/", signif(q[2], 3),
    " -> fold", round(q[2] / q[1], 1), "\n")

# ---- yday_mean ---------------------------------------------------------------
yday_mean <- rowMeans(cl$yday, na.rm = TRUE)
stopifnot(!anyNA(yday_mean))

# ---- write -------------------------------------------------------------------
b$constants_list$log_effort_cam <- log(sn)
b$constants_list$yday_mean      <- as.numeric(yday_mean)
b$y_cam <- y_cam
b$inits_list$det_intercept <- NULL
b$inits_list$overdisp_cam  <- 0.1
b$camcount_notes <- sprintf(paste(
  "Built %s from %s. y_cam = %s sequences collapsed to one per (project___deployment, sequence_id),",
  "restricted to [start_date,end_date], thinned to >= %g min between kept events (same rule as",
  "00d_prep_rn_counts.R); log_effort_cam = log(survey_nights); yday_mean = rowMeans(yday).",
  "Row order identical to source. For model_code_national_scalar_camcount."),
  format(Sys.time()), in_f, SCI, GAP_S / 60)
dir.create(out_d, showWarnings = FALSE)
saveRDS(b, out_f)
fwrite(data.table(uid = uid_in, y_cam, survey_nights = sn, log_effort_cam = log(sn), yday_mean, occ_any, J = cl$J),
       file.path(out_d, "camcount_site_table.csv"))
cat("\nwrote", out_f, "\nDONE\n")
