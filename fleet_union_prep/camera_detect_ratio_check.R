#!/usr/bin/env Rscript
# =============================================================================
# camera_detect_ratio_check.R -- STANDING CHECK on the camera record (rwkays 2026-10-08).
# Per fleet species: deployments with >= 1 detection in the camera record (umflist / umflist_v2) divided by raw
# detecting deployments (distinct project___deployment with >= 1 sequence under ANY of the species' camera names
# from camera_name_map.csv -- MDD-resolved names plus recorded decisions; never keyed on common name).
# Denominator is restricted to deployments that fall on the modelling grid (grid100.tif, CONUS): her umflist
# drops non-CONUS deployments, so an unrestricted denominator flagged brown bear (21 of 103 raw detecting
# deployments are in Alaska/Canada) for a reason that is not a loss.
# Fleet median ~0.9 (the 10-100-inner-day and full-window filters). Flags ratio < RATIO_MIN: any loss mechanism
# shows here, not only renames (peccary was 0.511 under one name).
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(terra) })
D <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data"
W <- "/home/rwkays/isdm/fleet_v2u"
RATIO_MIN <- as.numeric(Sys.getenv("RATIO_MIN", "0.75"))
options(width = 220)
cw <- fread(file.path(W, "species_name_crosswalk.csv")); map <- fread(file.path(W, "camera_name_map.csv"))
recs <- list(umflist = readRDS(file.path(D, "umflist.RDS")))
if (file.exists(file.path(D, "umflist_v2.RDS"))) recs$umflist_v2 <- readRDS(file.path(D, "umflist_v2.RDS"))

sq <- fread(file.path(D, "combined_sequences_all.csv"), select = c("project_id", "deployment_id", "genus", "species"),
            colClasses = list(character = c("project_id", "deployment_id")))
sq[, sciname := paste(genus, species)][, uid := paste(project_id, deployment_id, sep = "___")]
dep <- unique(fread(file.path(D, "combined_deployments_all.csv"), select = c("project_id", "deployment_id", "latitude", "longitude"),
                    colClasses = list(character = c("project_id", "deployment_id")))[!is.na(latitude) & !is.na(longitude)])
dep[, uid := paste(project_id, deployment_id, sep = "___")]
g <- rast(file.path(D, "grid100.tif"))
dep[, on_grid := !is.na(extract(g, project(vect(dep, geom = c("longitude", "latitude"), crs = "+proj=longlat"), crs(g)))[, 2])]
grid_uids <- unique(dep[on_grid == TRUE, uid])
cat("raw deployments with coordinates:", uniqueN(dep$uid), " on the modelling grid:", length(grid_uids), "\n")

rat <- rbindlist(lapply(seq_len(nrow(cw)), function(i) {
  sp <- cw$species[i]; k <- cw$wi_name[i]; nms <- map[species == sp, sciname]
  if (!k %in% names(recs$umflist)) return(NULL)
  r_all <- uniqueN(sq[sciname %in% nms, uid]); r_grid <- uniqueN(sq[sciname %in% nms & uid %in% grid_uids, uid])
  o <- data.table(species = sp, camera_names = paste(nms, collapse = " | "), raw_detecting_all = r_all, raw_detecting_on_grid = r_grid)
  for (nmr in names(recs)) o[, (paste0(nmr, "_detecting")) := sum(rowSums(attr(recs[[nmr]][[k]], "y"), na.rm = TRUE) > 0)]
  o
}))
for (nmr in names(recs)) {
  rc <- paste0("ratio_", nmr)
  rat[, (rc) := round(get(paste0(nmr, "_detecting")) / raw_detecting_on_grid, 3)]
  rat[, (paste0("flag_", nmr)) := ifelse(get(rc) < RATIO_MIN, "LOW", "")]
  cat(sprintf("%s: median ratio %.3f, IQR %.3f-%.3f, flagged < %.2f: %s\n", nmr, median(rat[[rc]]),
              quantile(rat[[rc]], .25), quantile(rat[[rc]], .75), RATIO_MIN,
              paste(rat[get(rc) < RATIO_MIN, species], collapse = ", ")))
}
print(rat[order(get(paste0("ratio_", names(recs)[length(recs)])))][, !"camera_names"])
fwrite(rat, file.path(W, "camera_detect_ratio_check.csv"))
cat("DONE\n")
