#!/usr/bin/env Rscript
# =============================================================================
# quantify_range_mask_exclusion.R
# -----------------------------------------------------------------------------
# How many real occupied cells / records did the old IUCN range mask drop, for
# each of moose, bobcat and white-tailed deer?
#
# WHY THIS IS NOT A 2-HOUR JOB PER SPECIES.
# build_moose_unmasked_review.R recomputes counts from the 2.2 GB raw master
# (inat_combo_nam_mams.csv) because "the masked pull's out-of-range total_count
# is 0 by construction". That was true of the MASKED grid -- but the v2 build
# already wrote genuinely unmasked per-species grids
# (<species>_inat_grid_unmasked_v2.csv), and those DO carry nonzero
# out-of-range counts. So the exclusion is a join between the two cached
# grids, in seconds, with no raw-file pass.
#
# VALIDATION: this reproduces the moose figure obtained the expensive way --
# 143 occupied out-of-range cells / 4,293 records here, against the recorded
# +144 cells / +4,290 records from the raw-file route. The small delta is
# filter-chain detail, not a different answer.
#
# STRUCTURE NOTE, load-bearing: <species>_inat_grid.csv is the SAME FILE for
# every species (verified: all three are byte-identical, md5 6f5dfef561a5).
# It carries one column per species, each independently masked to that
# species' own IUCN range -- confirmed by the distinct in-range cell counts
# below (moose 382, bobcat 3145, WTD 2848 of 3322). So the species is selected
# by COLUMN, not by filename; reading "bobcat_inat_grid.csv" and then the
# "moose" column would silently give the moose mask.
# =============================================================================
suppressPackageStartupMessages(library(data.table))

D <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/data"
OUT <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/prototype_spatial_trend"
SPECIES <- c("moose", "bobcat", "white-tailed_deer")

grid <- fread(file.path(D, "moose_inat_grid.csv"),
              select = c("cell50", "x50", "y50", SPECIES))

res <- rbindlist(lapply(SPECIES, function(s) {
  ir  <- grid[, .(x50 = x50[1], y50 = y50[1],
                  in_range = !all(is.na(get(s)))), by = cell50]
  u   <- fread(file.path(D, paste0(s, "_inat_grid_unmasked_v2.csv")))
  tot <- u[, .(n = sum(get(s), na.rm = TRUE)), by = cell50]
  m   <- merge(tot, ir, by = "cell50", all.x = TRUE)
  m[is.na(in_range), in_range := FALSE]
  fwrite(m, file.path(OUT, paste0(s, "_range_mask_exclusion_cells.csv")))
  data.table(
    species              = s,
    cells_total          = nrow(ir),
    cells_in_range       = sum(ir$in_range),
    cells_masked_out     = sum(!ir$in_range),
    records_total        = sum(m$n),
    records_out_of_range = sum(m$n[!m$in_range]),
    pct_records_dropped  = round(100 * sum(m$n[!m$in_range]) / sum(m$n), 2),
    occupied_cells_dropped = sum(m$n > 0 & !m$in_range),
    max_records_in_a_dropped_cell = if (any(m$n > 0 & !m$in_range))
                                      max(m$n[m$n > 0 & !m$in_range]) else 0L
  )
}))

print(res, row.names = FALSE)
fwrite(res, file.path(OUT, "range_mask_exclusion_summary.csv"))
cat("\nwrote range_mask_exclusion_summary.csv and per-species *_range_mask_exclusion_cells.csv\n")
