#!/usr/bin/env Rscript
# agency_check.R -- reproduce the v1fix agency sanity check, then apply to v2b.
#
# The check that originally caught the inat_effort sort bug: Pearson r between
# the model's per-region year_region DEVIATION and the agency trend score,
# across ecoregions. v1fix scored +0.971 (the pre-Fix1 bugged fit was -0.707).
#
# Two variants are reported because the memory records two numbers:
#   ALL     = every region with an agency score (drops Mediterranean California)
#   MODELED = additionally drops prior-only regions (0 camera AND 0 iNat)
#   CAMERA  = only regions with camera sites  (memory: v1fix r = +0.957)
# Reproducing v1fix's +0.971 / +0.957 validates the row set before v2b is read.

DIR <- "/Users/rwkays/claude_code/data_integration_arielle/prototype_spatial_trend"

ag <- read.csv(file.path(DIR, "moose_agency_by_ecoregion.csv"),
               stringsAsFactors = FALSE)
ag$key <- toupper(trimws(ag$ecoregion))

load_fit <- function(tag) {
  f <- file.path(DIR, sprintf("moose_%s_ecoregion_posterior.csv", tag))
  d <- read.csv(f, stringsAsFactors = FALSE)
  d$key <- toupper(trimws(d$region_name))
  d
}

report <- function(tag) {
  d <- load_fit(tag)
  m <- merge(d, ag, by = "key")
  stopifnot(nrow(m) == nrow(d))

  has_ag  <- !is.na(m$agency_score)
  modeled <- m$n_camera_sites > 0 | m$n_inat_cells > 0
  camera  <- m$n_camera_sites > 0

  sets <- list(
    ALL     = has_ag,
    MODELED = has_ag & modeled,
    CAMERA  = has_ag & camera
  )

  cat("\n=====================================================\n")
  cat("  ", tag, "\n")
  cat("=====================================================\n")
  for (nmset in names(sets)) {
    s <- sets[[nmset]]
    r_dev <- cor(m$year_region_mean[s], m$agency_score[s])
    r_abs <- cor(m$abs_trend_mean[s],   m$agency_score[s])
    rho   <- suppressWarnings(cor(m$year_region_mean[s], m$agency_score[s],
                                 method = "spearman"))
    cat(sprintf("\n  %-8s n=%d regions\n", nmset, sum(s)))
    cat(sprintf("     Pearson r  (deviation year_region vs agency) = %+.4f\n", r_dev))
    cat(sprintf("     Spearman   (deviation year_region vs agency) = %+.4f\n", rho))
    cat(sprintf("     Pearson r  (absolute abs_trend  vs agency)   = %+.4f\n", r_abs))
  }

  cat("\n  per-region detail (sorted by agency score):\n")
  mm <- m[order(m$agency_score), ]
  cat(sprintf("    %-34s %9s %9s %9s %6s %6s\n",
              "region", "agency", "year_reg", "abs_trend", "cam", "inat"))
  for (i in seq_len(nrow(mm))) {
    cat(sprintf("    %-34s %9s %9.4f %9.4f %6d %6d\n",
                mm$region_name[i],
                ifelse(is.na(mm$agency_score[i]), "NA",
                       sprintf("%.4f", mm$agency_score[i])),
                mm$year_region_mean[i], mm$abs_trend_mean[i],
                mm$n_camera_sites[i], mm$n_inat_cells[i]))
  }

  # sign agreement on the deviation, over regions with an agency score
  s <- sets$ALL
  agree <- sign(m$year_region_mean[s]) == sign(m$agency_score[s])
  cat(sprintf("\n  sign agreement (deviation vs agency, ALL set): %d / %d\n",
              sum(agree), sum(s)))

  invisible(m)
}

v1 <- report("v1fix")
v2 <- report("v2b")

cat("\n\n=====================================================\n")
cat("   SIDE BY SIDE: deviation year_region_mean\n")
cat("=====================================================\n")
cmp <- merge(v1[, c("key", "region_name", "year_region_mean", "abs_trend_mean",
                    "agency_score")],
             v2[, c("key", "year_region_mean", "abs_trend_mean")],
             by = "key", suffixes = c("_v1fix", "_v2b"))
cmp <- cmp[order(cmp$agency_score), ]
cat(sprintf("%-34s %8s %9s %9s | %9s %9s\n", "region", "agency",
            "dev_v1fix", "dev_v2b", "abs_v1fix", "abs_v2b"))
for (i in seq_len(nrow(cmp))) {
  cat(sprintf("%-34s %8s %9.4f %9.4f | %9.4f %9.4f\n",
              cmp$region_name[i],
              ifelse(is.na(cmp$agency_score[i]), "NA",
                     sprintf("%.4f", cmp$agency_score[i])),
              cmp$year_region_mean_v1fix[i], cmp$year_region_mean_v2b[i],
              cmp$abs_trend_mean_v1fix[i],   cmp$abs_trend_mean_v2b[i]))
}
cat(sprintf("\ncor(dev_v1fix, dev_v2b) across all 8 regions = %+.4f\n",
            cor(cmp$year_region_mean_v1fix, cmp$year_region_mean_v2b)))

# --- write the comparison table for the report ---------------------------
out <- merge(v2[, c("key","region_name","year_region_mean","year_region_q025",
                    "year_region_q975","abs_trend_mean","abs_trend_q025",
                    "abs_trend_q975","abs_p_negative","n_camera_sites",
                    "n_inat_cells","agency_score","dominant_dir","n_states",
                    "covered_area_frac")],
             v1[, c("key","year_region_mean","abs_trend_mean")],
             by = "key", suffixes = c("", "_v1fix"))
names(out)[names(out)=="year_region_mean_v1fix"] <- "v1fix_year_region_mean"
names(out)[names(out)=="abs_trend_mean_v1fix"]   <- "v1fix_abs_trend_mean"
out$key <- NULL
out <- out[order(-out$n_camera_sites), ]
write.csv(out, file.path(DIR, "moose_v2b_agency_comparison.csv"), row.names = FALSE)
cat("\nwrote moose_v2b_agency_comparison.csv\n")
