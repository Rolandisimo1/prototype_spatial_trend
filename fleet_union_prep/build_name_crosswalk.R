#!/usr/bin/env Rscript
# =============================================================================
# build_name_crosswalk.R -- per-source species-name crosswalk, VERIFIED against
# each source, plus the deployment_id-vs-site range-test delta for all species.
#
# Columns: species (pipeline key), inat_name, iucn_name, wi_name, and the
# evidence that each resolves: n_inat_records (accuracy<1000, 2008-2025),
# iucn_layer_feature (fleet34 feature found and its iucn_mapped_name),
# n_wi_detecting (umflist detection histories). A proposal row comes from
# target_species.csv sci_name; the only deviations are the ones the SOURCES
# force, and each is printed with the evidence:
#   inat_name: if sci_name has 0 iNat records, take the single same-genus name
#              that has records (reported; stops if not unique).
#   iucn_name: fleet34 layer's iucn_mapped_name (the layer is authoritative).
#   wi_name:   sci_name if umflist has it; otherwise stop.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(terra); library(sf) })
W <- "/home/rwkays/isdm/fleet_v2u"
tk <- fread("/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/target_species.csv")
tk[, species := gsub("[ -']", "_", tolower(common_name))]
plan <- fread(file.path(W, "table_fleet_launch_plan.csv"))
lyr <- vect(file.path(W, "fleet34_isdm", "fleet34_iucn_ranges_present_day.geojson"))
lyr_dt <- as.data.table(as.data.frame(lyr))

inat <- fread("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/inat_combo_nam_mams.csv",
              select = c("taxon_species_name", "public_positional_accuracy", "observed_on"))
inat <- inat[public_positional_accuracy < 1000]
inat[, Year := as.integer(substr(observed_on, 1, 4))]
inat_n <- inat[Year >= 2008 & Year <= 2025, .N, by = taxon_species_name]

um <- readRDS("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/umflist.RDS")

cw <- rbindlist(lapply(plan$species, function(sp) {
  sci <- tk[species == sp, sci_name]
  if (length(sci) != 1) stop("no unique target_species row for ", sp)
  genus <- strsplit(sci, " ")[[1]][1]
  # iNaturalist
  n_sci <- inat_n[taxon_species_name == sci, N]; n_sci <- if (length(n_sci)) n_sci else 0L
  inat_name <- sci; note <- character(0)
  if (n_sci == 0) {
    alt <- inat_n[startsWith(taxon_species_name, paste0(genus, " "))][order(-N)]
    cat(sprintf("  %s: 0 iNat records as '%s'; same-genus names with records:\n", sp, sci)); print(alt)
    if (nrow(alt) != 1) stop("cannot resolve iNat name for ", sp, " unambiguously -- set it by hand")
    inat_name <- alt$taxon_species_name; note <- c(note, sprintf("iNat uses %s (0 records as %s)", inat_name, sci))
  }
  n_inat <- inat_n[taxon_species_name == inat_name, N]
  # IUCN (fleet34 layer)
  f <- lyr_dt[species == sp]
  if (nrow(f) != 1) stop("fleet34 layer has ", nrow(f), " features for ", sp)
  iucn_name <- f$iucn_mapped_name
  if (iucn_name != sci) note <- c(note, sprintf("IUCN maps as %s (layer binomial %s)", iucn_name, f$binomial))
  # Wildlife Insights / umflist
  if (!sci %in% names(um)) stop("umflist has no detection histories under '", sci, "' for ", sp)
  y <- attr(um[[sci]], "y")
  data.table(species = sp, inat_name = inat_name, iucn_name = iucn_name, wi_name = sci,
             mdd_binomial = f$binomial, n_inat_records = n_inat,
             n_wi_detecting = sum(rowSums(y, na.rm = TRUE) > 0),
             note = paste(note, collapse = "; "))
}))
cat("\n=== CROSSWALK ===\n"); print(cw, nrows = 100)
cat("\nnon-identity rows:\n"); print(cw[inat_name != wi_name | iucn_name != wi_name])
fwrite(cw, file.path(W, "species_name_crosswalk.csv"))

# ---- deployment_id match vs per-site test, current IUCN polygon, cached filters
cat("\n=== deployment_id-match vs per-site polygon test (current IUCN polygon, cached 4SPO filters) ===\n")
brick <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/covar_raster_brick.tif")
cv <- c("terrain_ruggedness", "soil_clay", "soil_silt", "soil_sand", "elevation", "MWMT", "MCMT")
rangefile <- function(sci) { s <- gsub(" ", "_", sci); if (s == "Pekania_pennanti") s <- "Martes_pennanti"
  list.files("/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/iucn_ranges_extantonly/",
             pattern = paste0(s, ".*.shp"), full.names = TRUE) }
sc0 <- as.data.frame(attr(um[[cw$wi_name[1]]], "siteCovs"))
pts <- vect(sc0, geom = c("longitude", "latitude"), crs = "+proj=longlat")
cov_ok <- rowSums(is.na(sc0)) == 0 & rowSums(is.na(extract(brick[[cv]], pts)[, -1])) == 0
idd <- rbindlist(lapply(seq_len(nrow(cw)), function(k) {
  sci <- cw$wi_name[k]
  sc <- as.data.frame(attr(um[[sci]], "siteCovs"))
  stopifnot(identical(sc$deployment_id, sc0$deployment_id), identical(sc$longitude, sc0$longitude))
  rp <- vect(st_read(rangefile(sci), quiet = TRUE))
  by_id <- as.character(pts$deployment_id) %in% as.character(pts[rp, ]$deployment_id)
  by_pt <- is.related(pts, rp, "intersects")
  data.table(species = cw$species[k], sites_depid_match = sum(by_id & cov_ok),
             sites_per_site = sum(by_pt & cov_ok), added = sum(by_pt & !by_id & cov_ok),
             removed = sum(by_id & !by_pt & cov_ok))
}))
idd[, delta := sites_per_site - sites_depid_match]
print(idd, nrows = 100)
fwrite(idd, file.path(W, "deployment_id_fix_delta.csv"))
cat("DONE\n")
