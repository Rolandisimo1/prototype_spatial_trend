#!/usr/bin/env Rscript
# validate_iucn_layer.R -- gate before any grid is built from the fleet34 layer.
#  (1) RASTERISER VALIDATION: iucn_cells50() on the OLD extantonly shapefiles
#      must reproduce the IUCN cell sets embedded in data/moose_inat_grid.csv
#      (written by Arielle's prep_inat_data_grid) exactly. Hard stop otherwise.
#  (2) FLEET34 LAYER: every fleet species resolves to one feature touching > 0
#      grid cells; report old vs new IUCN cell50 counts side by side.
#  (3) ARMADILLO / PECCARY KEYS: which binomials the pipeline's own inputs use.
#  (4) CROSS-CHECK vs rwkays's 100 km grids (EPSG:5070 lattice, different cell
#      IDs): counts and geography only -- never cell IDs.
suppressPackageStartupMessages({ library(data.table); library(terra) })
W    <- "/home/rwkays/isdm/fleet_v2u"
PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
source(file.path(W, "iucn_source_lib.R"))
tk <- fread("/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/target_species.csv")
tk[, common_name_clean := gsub("[ -']", "_", tolower(common_name))]
OLD_CSV <- file.path(PROJ, "data", "moose_inat_grid.csv")
OLD_SHP <- "/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/iucn_ranges_extantonly/{sci_}_range.shp"
NEW     <- file.path(W, "fleet34_isdm", "fleet34_iucn_ranges_present_day.geojson")
g50 <- rast(GRID50_PATH); g100 <- rast(GRID100_PATH)
U <- unique(fread(OLD_CSV, select = "cell50")$cell50)
cat("pipeline grid universe:", length(U), "cell50\n")
cat("grid50 cells:", sum(!is.na(values(g50))), "  grid100 CRS:\n"); cat(crs(g100, describe = TRUE)$name, "\n")

cat("\n==== (1) rasteriser validation on old shapefiles ====\n")
ok <- TRUE
for (sp in c("moose", "bobcat", "white-tailed_deer", "coyote")) {
  a <- iucn_cells50(sp, OLD_CSV, universe = U)
  b <- iucn_cells50(sp, OLD_SHP, taxon_key = tk, universe = U)
  A <- sort(as.integer(a[in_iucn == TRUE, cell50])); B <- sort(as.integer(b[in_iucn == TRUE, cell50]))
  same <- identical(A, B)
  cat(sprintf("  %-18s csv %4d  shp-rasterised %4d  only_csv %d  only_shp %d  -> %s\n",
              sp, length(A), length(B), length(setdiff(A, B)), length(setdiff(B, A)),
              if (same) "IDENTICAL" else "DIFFERS"))
  ok <- ok && same
}
if (!ok) stop("RASTERISER VALIDATION FAILED -- not building from the new layer")
cat("  RASTERISER VALIDATED\n")

cat("\n==== (2) fleet34 present_day layer, all 34 plan species ====\n")
plan <- fread(file.path(W, "table_fleet_launch_plan.csv"))
lyr  <- vect(NEW)
cat("  features:", nrow(lyr), " crs:", crs(lyr, describe = TRUE)$name, "\n")
old_cols <- names(fread(OLD_CSV, nrows = 1))
res <- rbindlist(lapply(plan$species, function(sp) {
  n_new <- tryCatch(sum(iucn_cells50(sp, NEW, key = "species", universe = U)$in_iucn), error = function(e) { cat("  ERROR", sp, conditionMessage(e), "\n"); NA_integer_ })
  n_old <- if (sp %in% old_cols) sum(iucn_cells50(sp, OLD_CSV, universe = U)$in_iucn) else NA_integer_
  data.table(species = sp, tier = substr(plan[species == sp, tier], 1, 1),
             binomial_layer = lyr$binomial[lyr$species == sp],
             sci_name_pipeline = tk[common_name_clean == sp, sci_name][1],
             iucn50_old = n_old, iucn50_new = n_new)
}))
res[, ratio_new_old := round(iucn50_new / iucn50_old, 3)]
print(res, nrows = 100)
fwrite(res, file.path(W, "fleet34_iucn_old_vs_new.csv"))
if (anyNA(res$iucn50_new)) stop("some species did not resolve in the fleet34 layer")

cat("\n==== (3) binomial keys the pipeline itself uses ====\n")
print(res[binomial_layer != sci_name_pipeline | is.na(sci_name_pipeline)])
inat <- fread("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/inat_combo_nam_mams.csv",
              select = c("taxon_species_name", "public_positional_accuracy", "observed_on"))
inat <- inat[public_positional_accuracy < 1000]
inat[, Year := as.integer(substr(observed_on, 1, 4))]
inat <- inat[Year >= 2008 & Year <= 2025]
cat("  iNat pull (accuracy<1000, 2008-2025) records by name, Dasypus / Pecari / Dicotyles:\n")
print(inat[grepl("^(Dasypus|Pecari|Dicotyles)", taxon_species_name), .N, by = taxon_species_name])
tryCatch({
  um <- readRDS("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/umflist.RDS")
  cat("  umflist (camera) names matching Dasypus/Pecari:", grep("Dasypus|Pecari|Dicotyles", names(um), value = TRUE), "\n")
  for (nm in grep("Dasypus|Pecari", names(um), value = TRUE)) {
    y <- methods::slot(um[[nm]], "y")
    cat("   ", nm, " sites", nrow(y), " detecting", sum(rowSums(y, na.rm = TRUE) > 0), "\n")
  }
}, error = function(e) cat("  umflist check skipped:", conditionMessage(e), "\n"))
old_arm <- vect(sub("{sci_}", "Dasypus_novemcinctus", OLD_SHP, fixed = TRUE))
cat("  old extantonly Dasypus_novemcinctus polygon extent (camera-side filter):", as.vector(ext(old_arm)), "\n")

cat("\n==== (4) cross-check vs rwkays 100 km grids (counts + geography) ====\n")
pres_all <- fread(file.path(PROJ, "output", "inat_grids_v2", "bobcat_inat_grid_v2.csv"))
xc <- list()
for (sp in c("moose", "bobcat", "white-tailed_deer", "coyote")) {
  cells <- pres_all[, .(cell100 = cell100[1], x100 = x100[1], y100 = y100[1], x50 = x50[1], y50 = y50[1],
                        pres = any(!is.na(get(sp)))), by = cell50]
  new_iucn <- iucn_cells50(sp, NEW, key = "species", universe = U)
  cells <- merge(cells, new_iucn, by = "cell50")
  cells[, uni := pres | in_iucn]
  pip_pres100 <- uniqueN(cells[pres == TRUE, cell100]); pip_uni100 <- uniqueN(cells[uni == TRUE, cell100])
  rk_p <- vect(file.path(W, "fleet34_isdm", "grids_100km", paste0(sp, "_presence_100km.geojson")))
  rk_u <- vect(file.path(W, "fleet34_isdm", "grids_100km", paste0(sp, "_union_100km.geojson")))
  # geography: pipeline 50 km presence-cell centroids falling inside rwkays cells, and
  # rwkays presence-cell centroids falling on a pipeline presence cell50
  pts <- project(vect(as.matrix(cells[pres == TRUE, .(x50, y50)]), crs = crs(g100)), crs(rk_p))
  in_rk <- mean(!is.na(extract(rk_p, pts)[, 2]))
  rk_cent <- project(centroids(rk_p), crs(g50))
  on_pip <- extract(g50, rk_cent)[, 2]
  rk_on_pip <- mean(on_pip %in% cells[pres == TRUE, cell50])
  rk_off_grid <- mean(is.na(on_pip))
  xc[[sp]] <- data.table(species = sp, pipeline_c100_pres = pip_pres100, rk_pres100 = nrow(rk_p),
                         ratio_pres = round(nrow(rk_p) / pip_pres100, 2),
                         pipeline_c100_union_newIUCN = pip_uni100, rk_union100 = nrow(rk_u),
                         ratio_union = round(nrow(rk_u) / pip_uni100, 2),
                         pip_pres50_inside_rk_pres = round(in_rk, 3),
                         rk_pres_centroid_on_pip_pres50 = round(rk_on_pip, 3),
                         rk_pres_centroid_off_pipeline_grid = round(rk_off_grid, 3))
}
xc <- rbindlist(xc); print(xc)
fwrite(xc, file.path(W, "crosscheck_rk100km.csv"))
cat("\nVALIDATION DONE\n")
