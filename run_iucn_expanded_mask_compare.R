#!/usr/bin/env Rscript
# =============================================================================
# run_iucn_expanded_mask_compare.R
# -----------------------------------------------------------------------------
# Wires build_iucn_expanded_mask.R's flood_fill_expand_mask() to real data and
# compares THREE mask definitions for moose / bobcat / WTD. Read-only,
# evaluation-only: writes comparison CSVs, wires nothing into any run script.
#
# THE TWO THINGS THE DRIVER'S stop() ASKED FOR, BOTH CONFIRMED ON HAZEL:
#
# (a) RANGE_DIR naming is "<Genus>_<species>_range.shp" -- UNDERSCORES.
#     The placeholder used list.files(pattern = sci_name) with sci_name as
#     "Alces alces"; the space form matches ZERO of the 512 shapefiles.
#     Second trap in the same line: range_files[1] off an unanchored pattern
#     returns the .dbf (dbf < prj < shp < shx < zip alphabetically), not the
#     .shp -- vect() on the .dbf gives attributes with no geometry. Both fixed
#     here by constructing the exact filename and asserting it exists.
#
# (b) Fix 2b's isolation/captive-filtered counts are ALREADY CACHED as
#     data/<species>_inat_grid_unmasked_v2.csv. integration_helper_v2.R:412
#     writes that file only after applying Fix 2(b) captive_cultivated
#     exclusion (:451-452) and Fix 2(c) spatial-isolation 100km/3yr filter
#     (:489-508). So the >2-detection threshold below reads exactly what Fix
#     2b calls a valid detection. Nothing is recomputed from the raw master.
#
# GRID GEOMETRY: flood_fill_expand_mask() wants INTEGER cell indices, not
# projected metres. x100/y100 are metres on a confirmed exact 100,000 m
# lattice, so index = coord / 1e5. Feeding metres directly would make every
# +/-1 neighbour lookup miss and the fill would expand nothing -- it would
# fail silently as "no expansion", which looks like a result rather than a bug.
#
# ENV: needs terra -> plotting_env, NOT nimble_env (nimble_env has no sf/terra).
# =============================================================================
# PROJ DATA PATH -- load-bearing. plotting_env ships share/proj/proj.db but
# neither PROJ_DATA nor PROJ_LIB is set, so PROJ cannot open its database and
# project() dies with "[project] Cannot do this transformation" after dumping
# the whole CRS as JSON (which buries the one-line cause). Set it explicitly
# rather than relying on the env's activation hooks, which do not run under a
# bare Rscript invocation.
.penv <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/HPC/conda_envs/plotting_env/share/proj"
if (dir.exists(.penv)) Sys.setenv(PROJ_DATA = .penv, PROJ_LIB = .penv)
suppressPackageStartupMessages({library(data.table); library(terra)})

PROJ      <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
DATA_DIR  <- file.path(PROJ, "data")
RAW_DIR   <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data"
RANGE_DIR <- "/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/iucn_ranges_extantonly"
OUT_DIR   <- file.path(PROJ, "prototype_spatial_trend")
CELL100_M <- 100000
MIN_DET   <- 2      # strict >, i.e. a cell needs >=3 filtered detections

source(file.path(OUT_DIR, "build_iucn_expanded_mask.R"))   # flood_fill_expand_mask()

SPECIES <- list(
  list(sci = "Alces alces",            lab = "moose"),
  list(sci = "Lynx rufus",             lab = "bobcat"),
  list(sci = "Odocoileus virginianus", lab = "white-tailed_deer")
)

grid100 <- rast(file.path(RAW_DIR, "grid100.tif"))
# The ORIGINAL IUCN mask as actually applied: one shared grid file, one column
# per species, each masked to that species' own range. Species selects the
# COLUMN, not the filename (all <species>_inat_grid.csv are byte-identical).
masked <- fread(file.path(DATA_DIR, "moose_inat_grid.csv"),
                select = c("cell50", "cell100", vapply(SPECIES, `[[`, "", "lab")))

run_one <- function(sci, lab) {
  cat("\n=====", lab, "(", sci, ") =====\n")

  # ---- (a) IUCN polygon -> cell100 base mask ------------------------------
  shp <- file.path(RANGE_DIR, paste0(gsub(" ", "_", sci), "_range.shp"))
  stopifnot("IUCN shapefile not found" = file.exists(shp))
  poly <- project(vect(shp), crs(grid100))
  base_ids <- unique(terra::extract(grid100, poly, ID = FALSE)[, 1])
  base_ids <- base_ids[!is.na(base_ids)]
  cat("  IUCN polygon covers", length(base_ids), "cell100\n")

  # ---- (b) Fix 2b filtered detection counts, aggregated to cell100 --------
  u <- fread(file.path(DATA_DIR, paste0(lab, "_inat_grid_unmasked_v2.csv")))
  c100 <- u[, .(count = sum(get(lab), na.rm = TRUE),
                x100 = x100[1], y100 = y100[1]), by = cell100]
  c100[, `:=`(xi = as.integer(round(x100 / CELL100_M)),
              yi = as.integer(round(y100 / CELL100_M)))]
  stopifnot("grid index collision" = !anyDuplicated(paste(c100$xi, c100$yi)))
  c100[, in_iucn := cell100 %in% base_ids]

  # ---- flood fill ---------------------------------------------------------
  c100[, in_expanded := flood_fill_expand_mask(
    x = xi, y = yi, count = count, in_range_base = in_iucn,
    min_detections = MIN_DET)]
  c100[, in_presence := count >= 1]

  # ---- join down to cell50 ------------------------------------------------
  m50 <- unique(masked[, .(cell50, cell100)])
  m50 <- merge(m50, c100[, .(cell100, in_iucn, in_expanded, in_presence)],
               by = "cell100", all.x = TRUE)
  for (v in c("in_iucn", "in_expanded", "in_presence"))
    m50[is.na(get(v)), (v) := FALSE]

  # cell50-native masks, for a like-for-like against how each was really used
  orig50 <- masked[, .(iucn50 = !all(is.na(get(lab)))), by = cell50]
  u50    <- u[, .(n50 = sum(get(lab), na.rm = TRUE)), by = cell50]
  n50    <- merge(orig50, u50, by = "cell50", all.x = TRUE)
  n50[is.na(n50), n50 := 0]

  fwrite(c100, file.path(OUT_DIR, paste0(lab, "_mask_compare_cell100.csv")))

  data.table(
    species = lab,
    cell100_total = nrow(c100),
    cell100_iucn = sum(c100$in_iucn),
    cell100_expanded = sum(c100$in_expanded),
    cell100_presence = sum(c100$in_presence),
    cell100_added_by_fill = sum(c100$in_expanded & !c100$in_iucn),
    cell50_total = nrow(n50),
    cell50_iucn_orig = sum(n50$iucn50),          # the mask as actually applied
    cell50_presence  = sum(n50$n50 >= 1),        # Fix 2b, its native cell50 rule
    cell50_expanded  = sum(m50$in_expanded)      # this candidate, via parent cell100
  )
}

res <- rbindlist(Map(function(s) run_one(s$sci, s$lab), SPECIES))
cat("\n\n================ MASK COMPARISON ================\n")
print(res, row.names = FALSE)
fwrite(res, file.path(OUT_DIR, "mask_comparison_summary.csv"))
cat("\nwrote mask_comparison_summary.csv + per-species *_mask_compare_cell100.csv\n")
