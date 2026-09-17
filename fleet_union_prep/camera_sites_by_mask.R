#!/usr/bin/env Rscript
# =============================================================================
# camera_sites_by_mask.R -- camera deployments per species under each range
# filter, by the REAL site-level test (no 100 km cell approximation).
#
# Universe: umflist[[sci_name]] -- every CONUS deployment (01_Prep_UMFs.R).
# Filters (site in / out):
#   iucn_old_poly   point in old extantonly polygon: site_pts2[range_poly, ],
#                   exactly prep_spoccupancy_data() lines 101-105 (CURRENT)
#   iucn_new_poly   point in fleet34 present_day polygon
#   presence_cell   site's cell50 is a presence-mask cell (cell50 are the mask's
#                   own polygons, so this IS the point-in-polygon test for it)
#   union_cell      site's cell50 is in the grid's union mask
#                   (fleet34 IUCN rasterised touches=TRUE, OR presence) --
#                   the definition that makes camera and iNat masks identical
#   union_point     point in fleet34 polygon OR site's cell50 is a presence cell
# For each: sites in mask, and sites after the standard non-range filters
# prep_spoccupancy_data() applies next (complete siteCovs, complete static +
# temp covariates at the point) -- the rule the cached 4SPO files actually used.
# SELF-CHECK: iucn_old_poly after filters must equal the fitted bundle nsite.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(terra); library(sf) })
W    <- "/home/rwkays/isdm/fleet_v2u"
PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
source(file.path(W, "iucn_source_lib.R"))
NEW  <- file.path(W, "fleet34_isdm", "fleet34_iucn_ranges_present_day.geojson")
SPECIES <- strsplit(Sys.getenv("SPECIES", paste(fread(file.path(W, "table_fleet_launch_plan.csv"))$species, collapse = " ")), "[ ,]+")[[1]]

tk <- fread("/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/target_species.csv")
tk[, common_name_clean := gsub("[ -']", "_", tolower(common_name))]
rangefile <- function(sci) {                      # as prep_data_for_spoccupancy.R
  s <- gsub(" ", "_", sci); if (s == "Pekania_pennanti") s <- "Martes_pennanti"
  list.files("/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/iucn_ranges_extantonly/",
             pattern = paste0(s, ".*.shp"), full.names = TRUE)
}
static_covars <- c("terrain_ruggedness", "soil_clay", "soil_silt", "soil_sand", "elevation")
temp_covars   <- c("MWMT", "MCMT")
brick <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/covar_raster_brick.tif")
g50 <- rast(GRID50_PATH); g100 <- rast(GRID100_PATH)
pres_all <- fread(file.path(PROJ, "output", "inat_grids_v2", "bobcat_inat_grid_v2.csv"))
U <- unique(pres_all$cell50)
lyr <- vect(NEW)
um <- readRDS("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/umflist.RDS")
bundle_nsite <- c(moose = 1654L, bobcat = 20531L, "white-tailed_deer" = 21559L, coyote = 24113L)

rows <- list()
for (sp in SPECIES) tryCatch({
  sci <- tk[common_name_clean == sp, sci_name]
  if (!sp %in% names(pres_all)) stop("no presence column")
  if (all(is.na(pres_all[[sp]]))) cat("  WARNING", sp, ": presence column is ENTIRELY NA -- presence half is empty\n")
  u <- um[[sci]]; sc <- as.data.frame(attr(u, "siteCovs")); y <- attr(u, "y")
  n <- nrow(sc)
  pts2 <- vect(sc, geom = c("longitude", "latitude"), crs = "+proj=longlat")

  # CURRENT camera filter, verbatim logic
  range_old <- vect(st_read(rangefile(sci), quiet = TRUE))
  in_old <- as.character(pts2$deployment_id) %in% as.character(pts2[range_old, ]$deployment_id)
  # (deployment_id alone is how the prep maps back; report if it is not unique)
  dup_dep <- anyDuplicated(as.character(sc$deployment_id)) > 0

  poly_new <- lyr[lyr$species == sp, ]
  in_new <- is.related(pts2, project(poly_new, crs(pts2)), "intersects")

  cell50 <- extract(g50, project(pts2, crs(g100)))[, 2]
  pres_cells <- unique(pres_all[!is.na(get(sp)), cell50])
  iu <- iucn_cells50(sp, NEW, key = "species", universe = U)
  union_cells <- union(iu[in_iucn == TRUE, cell50], pres_cells)
  in_pres  <- !is.na(cell50) & cell50 %in% pres_cells
  in_union <- !is.na(cell50) & cell50 %in% union_cells
  in_upt   <- in_new | in_pres

  # standard non-range filters, evaluated for every site once
  cov <- extract(brick[[c(static_covars, temp_covars)]], pts2)[, -1]
  ok_cov <- rowSums(is.na(sc)) == 0 & rowSums(is.na(cov)) == 0
  ok_occ <- rowSums(is.na(y)) < (ncol(y) - 1)
  det <- rowSums(y, na.rm = TRUE) > 0
  # EXACT rule of the cached 4SPO files behind every fitted bundle (verified by
  # uid-set equality, job 858979): NA siteCovs + NA raster covariates only. The
  # single-occasion filter in today's prep code was never applied to them.
  keep <- ok_cov

  f <- list(iucn_old_poly = in_old, iucn_new_poly = in_new, presence_cell = in_pres,
            union_cell = in_union, union_point = in_upt)
  for (nm in names(f)) {
    m <- f[[nm]]
    rows[[length(rows) + 1]] <- data.table(species = sp, filter = nm,
      sites_in_mask = sum(m), sites_after_filters = sum(m & keep),
      detecting_after_filters = sum(m & keep & det))
  }
  cur <- sum(in_old & keep)
  cat(sprintf("%-18s universe %d | current-filter nsite %d vs bundle %s -> %s | deployment_id dup: %s | off-grid sites %d\n",
              sp, n, cur, bundle_nsite[sp],
              if (is.na(bundle_nsite[sp])) "n/a" else if (cur == bundle_nsite[sp]) "REPRODUCED" else "MISMATCH",
              dup_dep, sum(is.na(cell50))))
  cat(sprintf("   old-poly-in but union_cell-out: %d   union_cell-in but old-poly-out: %d   (after filters: %d / %d)\n",
              sum(in_old & !in_union), sum(in_union & !in_old),
              sum(in_old & !in_union & keep), sum(in_union & !in_old & keep)))
  cat("   iucn cells dropped outside universe:", attr(iu, "n_outside_universe"), "\n")
}, error = function(e) cat("  SKIPPED", sp, ":", conditionMessage(e), "\n"))
res <- rbindlist(rows)
base <- res[filter == "iucn_old_poly", .(species, base = sites_after_filters)]
res <- merge(res, base, by = "species", sort = FALSE)
res[, change_vs_current := sites_after_filters - base]
res[, pct_change_vs_current := round(100 * change_vs_current / base, 1)]
res[, base := NULL]
print(res, nrows = 500)
fwrite(res, file.path(W, "camera_sites_by_mask.csv"))
cat("DONE\n")
