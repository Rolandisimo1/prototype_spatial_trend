#!/usr/bin/env Rscript
# =============================================================================
# 00d_prep_rn_counts.R
# -----------------------------------------------------------------------------
# Build the per-window DETECTION COUNT matrix that the Royle-Nichols arm needs,
# from the raw camera data. Writes rn_counts_<species>.RDS for use by
# 00_prep_sim_inputs.R.
#
# WHAT THIS PRODUCES AND WHY
# RN exploits detection FREQUENCY -- "more animals produce more detections" --
# so its response is the number of independent detections per camera-window,
# not a binary indicator and NOT group size.
#
# GROUP SIZE IS DELIBERATELY IGNORED. It was evaluated and rejected:
# moose mean group 1.08, bobcat 1.04 (both solitary; 93% / 97% of detections
# are a single animal), and group_size is 100% missing for 2 of 27 project
# groups -- structured whole-project absence, not random, so dropping NAs would
# drop two projects and their geography. Detection frequency carries the
# abundance signal for solitary species; group size does not.
#
# THREE QC STEPS, EACH LOAD-BEARING (see camera-trap-qc skill)
#
# 1. sequence_id IS NOT A UNIQUE KEY.
#    The sequences table carries one row per (species x age x sex) class within
#    an event: a doe with two fawns is legitimately two rows on one sequence_id.
#    Measured here: 2,517,928 rows -> 2,248,930 distinct sequence_ids;
#    122,369 multi-row sequences (72,858 genuinely multi-species,
#    49,511 same-species split by class).
#    Counting raw rows inflates per-deployment detections (moose max 146 -> 73,
#    a 2x error). Collapse to one row per (deployment, sequence, species)
#    FIRST. For an events/frequency response that is the correct collapse;
#    summing group_size would be correct only for an animals-counted response,
#    which this is not.
#
# 2. 30-MINUTE INDEPENDENCE FILTER.
#    Sequences are algorithmically grouped triggers, so one animal lingering
#    produces several. RN reads temporal clustering as more animals, which
#    biases abundance upward exactly where it matters. Thinning to >= 30 min
#    between consecutive detections of the same species at the same deployment
#    removes 22.6% of moose events, 7.7% of bobcat, and 30.0% of deer -- the
#    deer figure shows how large the artifact would have been.
#
# 3. EFFORT SANITY BOUNDS.
#    survey_nights has 87 zero-night deployments, 19 NA-year rows, and a max of
#    2,688 days (~7.4 years -- a deployment never closed). Bounded to
#    (0, 365] days here, dropping 115 of 26,798 deployments (0.43%). Effort now
#    enters the array arms as a covariate, so an unbounded outlier would
#    propagate into the detection model.
#
# WINDOW ALIGNMENT CHECK
# J = ceil(survey_nights / 10) gives median 4, mean 4.08 windows per
# deployment, against the existing prep pipeline's measured mean J of 3.22.
# The same order, slightly higher -- consistent with the ~7% of raw
# deployments that prep filters out. VERIFY per-deployment against the existing
# binary y before use: every deployment-window with count >= 1 must have y = 1.
# Any mismatch means the filtering differs and must be understood, not patched.
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
})

RAW  <- Sys.getenv("RAW_CAM", "../raw_cam_data")
OUT  <- Sys.getenv("SIM_BASE", ".")
INDEP_MIN <- as.numeric(Sys.getenv("INDEP_MIN", "30"))
WIN_DAYS  <- as.numeric(Sys.getenv("WIN_DAYS", "10"))
MAX_NIGHTS <- as.numeric(Sys.getenv("MAX_NIGHTS", "365"))

#' @name prep_rn_counts
#' @description Build a per-(deployment, window) independent-detection count
#'   table for one species, applying the three QC steps in the header.
#' @param species_regex Regex matched case-insensitively against common_name.
#' @param label Short species label used in the output filename.
#' @return Invisibly, the count data frame; also written to
#'   rn_counts_<label>.RDS with an attached QC summary.
prep_rn_counts <- function(species_regex, label) {

  dep <- read.csv(file.path(RAW, "combined_deployments_all.csv"),
                  stringsAsFactors = FALSE)
  seq <- read.csv(file.path(RAW, "combined_sequences_all.csv"),
                  stringsAsFactors = FALSE)

  n_dep_raw <- nrow(dep); n_seq_raw <- nrow(seq)

  # ---- QC 3: effort bounds ------------------------------------------------
  dep <- dep %>%
    filter(!is.na(year), !is.na(survey_nights),
           survey_nights > 0, survey_nights <= MAX_NIGHTS) %>%
    mutate(start_ts = as.POSIXct(start_date, tz = "UTC"),
           J = ceiling(survey_nights / WIN_DAYS))

  # ---- QC 1: collapse age/sex class rows to one row per event -------------
  # For a frequency response the correct collapse is distinct
  # (deployment, sequence, species) -- NOT summing group_size.
  ev <- seq %>%
    filter(grepl(species_regex, common_name, ignore.case = TRUE)) %>%
    mutate(ts = as.POSIXct(start_time, tz = "UTC")) %>%
    filter(!is.na(ts)) %>%
    distinct(deployment_id, sequence_id, common_name, .keep_all = TRUE)
  n_events <- nrow(ev)

  # ---- QC 2: 30-minute independence filter --------------------------------
  ev <- ev %>%
    arrange(deployment_id, common_name, ts) %>%
    group_by(deployment_id, common_name) %>%
    mutate(gap_min = as.numeric(difftime(ts, lag(ts), units = "mins"))) %>%
    filter(is.na(gap_min) | gap_min >= INDEP_MIN) %>%
    ungroup()
  n_indep <- nrow(ev)

  # ---- assign to 10-day windows within the deployment --------------------
  counts <- ev %>%
    inner_join(dep[, c("deployment_id", "start_ts", "J", "survey_nights")],
               by = "deployment_id") %>%
    mutate(win = as.integer(as.numeric(difftime(ts, start_ts, units = "days")) / WIN_DAYS)) %>%
    filter(win >= 0, win < J) %>%
    count(deployment_id, win, name = "count")

  qc <- list(
    label = label,
    deployments_raw = n_dep_raw, deployments_kept = nrow(dep),
    sequence_rows_raw = n_seq_raw,
    events_after_class_collapse = n_events,
    events_after_independence = n_indep,
    pct_removed_by_independence = round(100 * (1 - n_indep / max(n_events, 1)), 1),
    occupied_windows = nrow(counts),
    mean_count = round(mean(counts$count), 3),
    median_count = stats::median(counts$count),
    max_count = max(counts$count),
    frac_gt1 = round(mean(counts$count > 1), 4),
    indep_min = INDEP_MIN, win_days = WIN_DAYS
  )
  attr(counts, "qc") <- qc
  str(qc)

  outfile <- file.path(OUT, sprintf("rn_counts_%s.RDS", label))
  saveRDS(counts, outfile)
  cat("wrote", outfile, "\n")
  invisible(counts)
}

#' @name verify_counts_against_binary
#' @description Consistency check: every deployment-window with an independent
#'   detection must be a 1 in the existing binary y. Run BEFORE using these
#'   counts; a mismatch means the raw and prep filtering differ.
#' @param counts Output of prep_rn_counts().
#' @param y Existing binary detection matrix (nsite x max_J).
#' @param site_ids Character vector of deployment ids aligned to y's rows.
#' @return A list of mismatch counts; stops if any count>0 maps to y==0.
verify_counts_against_binary <- function(counts, y, site_ids) {
  idx <- match(counts$deployment_id, site_ids)
  ok <- !is.na(idx)
  yy <- vapply(seq_len(sum(ok)), function(k) {
    i <- idx[ok][k]; j <- counts$win[ok][k] + 1L
    if (j <= ncol(y)) y[i, j] else NA_real_
  }, numeric(1))
  n_bad <- sum(yy == 0, na.rm = TRUE)
  res <- list(matched = sum(ok), unmatched = sum(!ok),
              count_pos_but_y_zero = n_bad)
  if (n_bad > 0) {
    stop(sprintf(paste("INCONSISTENT: %d deployment-windows have an independent",
                       "detection but y == 0. Raw and prep filtering disagree;",
                       "resolve before using these counts."), n_bad))
  }
  res
}

if (!interactive()) {
  prep_rn_counts("Bobcat", "bobcat")
  prep_rn_counts("Moose", "moose")
  prep_rn_counts("White-tailed Deer", "white-tailed_deer")
}
