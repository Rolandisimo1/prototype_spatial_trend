#!/usr/bin/env Rscript
# =============================================================================
# 00b_prep_ecoregion.R
# Read-only w.r.t. Arielle's data. Assigns an EPA/CEC Level I ecoregion to
# every cell100 (the same 908-cell CAR grid used throughout this prototype),
# for the pivot from a fine CAR field (found to shrink under-informed cells
# toward zero at the ~900-cell grain) to a coarser, structured ecoregion
# random effect.
#
# Shapefile: EPA "Ecoregions of North America" / CEC NA_CEC_Eco_Level1
# (downloaded here on the LOGIN node -- compute nodes have no internet --
# from EPA's current S3-hosted data commons; the gaftp.epa.gov path referenced
# in older docs 404s now). Verified: 2140 features, fields NA_L1CODE/
# NA_L1NAME/NA_L1KEY, CRS = Lambert Azimuthal Equal Area (na_cec2010.prj).
#
# RUN (login node, needs internet for the download + plotting_env for sf/terra):
#   $PROJ/HPC/conda_envs/plotting_env/bin/Rscript 00b_prep_ecoregion.R
# Writes (UPDATES prepped_sim_inputs.RDS in place, does not create a separate
# file -- the abundance-sweep driver expects everything bundled in one self-
# contained RDS): adds cell100_geo$ecoregion / $ecoregion_id,
# ecoregion_of_cell100, nregion, ecoregion_levels, region_table, plus
# inat_effort_real / real_y_template / base_inits (pulled from
# input_data_bobcat.RDS, not previously bundled -- every prior CAR-based
# script loaded these separately; the abundance sweep expects them already
# in prepped_sim_inputs.RDS so it's a single self-contained load).
# =============================================================================

suppressPackageStartupMessages({ library(sf); library(dplyr) })

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
ECO_DIR <- paste0(PROTO_DIR, "/ecoregion_data")
setwd(PROTO_DIR)

SHP_URL <- "https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/cec_na/na_cec_eco_l1.zip"
SHP_ZIP <- paste0(ECO_DIR, "/na_cec_eco_l1.zip")
SHP_FILE <- paste0(ECO_DIR, "/NA_CEC_Eco_Level1.shp")

dir.create(ECO_DIR, showWarnings = FALSE, recursive = TRUE)
if (!file.exists(SHP_FILE)) {
  cat("downloading EPA Level I ecoregions shapefile (login node only -- compute nodes have no internet)...\n")
  ok <- download.file(SHP_URL, SHP_ZIP, mode = "wb", quiet = TRUE, timeout = 120)
  stopifnot(ok == 0)
  unzip(SHP_ZIP, exdir = ECO_DIR)
}
stopifnot(file.exists(SHP_FILE))

eco <- st_read(SHP_FILE, quiet = TRUE)
cat("ecoregion shapefile: n features=", nrow(eco), " CRS=", st_crs(eco)$input, "\n")

# ------------------------------ spatial join -----------------------------------
prepped <- readRDS("prepped_sim_inputs.RDS")
cell100_geo <- prepped$cell100_geo   # cell100, lon, lat (EPSG:4326), region (old NE/West -- untouched)

pts <- st_as_sf(cell100_geo, coords = c("lon", "lat"), crs = 4326)
pts_eco_crs <- st_transform(pts, st_crs(eco))

joined <- st_join(pts_eco_crs, eco[, "NA_L1NAME"], join = st_intersects)
n_missing <- sum(is.na(joined$NA_L1NAME))
cat("cells with no direct polygon match (ocean/edge):", n_missing, "/", nrow(joined), "\n")

if (n_missing > 0) {
  missing_idx <- which(is.na(joined$NA_L1NAME))
  nearest <- st_nearest_feature(pts_eco_crs[missing_idx, ], eco)
  joined$NA_L1NAME[missing_idx] <- eco$NA_L1NAME[nearest]
}
stopifnot(!any(is.na(joined$NA_L1NAME)))

cell100_geo$ecoregion_raw <- as.character(joined$NA_L1NAME)

cat("\nraw region table (cell100 counts, bobcat range only):\n")
print(sort(table(cell100_geo$ecoregion_raw), decreasing = TRUE))

# ------------------------------ fold small regions -------------------------------
MIN_REGION_SIZE <- 15
counts <- table(cell100_geo$ecoregion_raw)
small_regions <- names(counts)[counts < MIN_REGION_SIZE]

if (length(small_regions) > 0) {
  cat("\nfolding", length(small_regions), "region(s) with <", MIN_REGION_SIZE,
      "cells into their nearest large-region neighbor:", paste(small_regions, collapse=", "), "\n")

  large_idx <- which(!cell100_geo$ecoregion_raw %in% small_regions)
  # need projected (planar) coordinates for a sane nearest-neighbor distance --
  # reuse the eco CRS points already computed above
  coords_all <- st_coordinates(pts_eco_crs)

  for (r in small_regions) {
    small_idx <- which(cell100_geo$ecoregion_raw == r)
    for (i in small_idx) {
      d <- sqrt((coords_all[large_idx, 1] - coords_all[i, 1])^2 +
                (coords_all[large_idx, 2] - coords_all[i, 2])^2)
      nearest_large <- large_idx[which.min(d)]
      cell100_geo$ecoregion_raw[i] <- cell100_geo$ecoregion_raw[nearest_large]
    }
  }
}

cell100_geo$ecoregion <- factor(cell100_geo$ecoregion_raw)
cell100_geo$ecoregion_id <- as.integer(cell100_geo$ecoregion)
nregion <- nlevels(cell100_geo$ecoregion)
ecoregion_of_cell100 <- cell100_geo$ecoregion_id   # length ncell100, indexed 1:ncell100 by construction

cat("\n=== FINAL region table (K =", nregion, ") ===\n")
final_table <- as.data.frame(table(cell100_geo$ecoregion))
names(final_table) <- c("ecoregion", "n_cells")
final_table <- final_table[order(-final_table$n_cells), ]
print(final_table, row.names = FALSE)
stopifnot(all(final_table$n_cells >= MIN_REGION_SIZE))

# ------------------------------ pull the fields every prior script loaded ------
# separately from input_data_bobcat.RDS, so prepped_sim_inputs.RDS becomes a
# single self-contained bundle for the abundance sweep.
inp <- readRDS(paste0(PROJ, "/HPC/bobcat/input_data_bobcat.RDS"))

# ------------------------------ save (update prepped_sim_inputs.RDS in place) --
prepped$cell100_geo <- cell100_geo
prepped$ecoregion_of_cell100 <- ecoregion_of_cell100
prepped$nregion <- nregion
prepped$region_table <- final_table
prepped$ecoregion_levels <- levels(cell100_geo$ecoregion)
prepped$inat_effort_real <- inp$inat_effort
prepped$real_y_template <- inp$real_data$y
prepped$base_inits <- inp$inits_list

saveRDS(prepped, "prepped_sim_inputs.RDS")
cat("\nupdated prepped_sim_inputs.RDS with ecoregion fields + inat_effort_real/real_y_template/base_inits\n")
