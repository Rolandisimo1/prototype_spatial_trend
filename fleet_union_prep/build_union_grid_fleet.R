#!/usr/bin/env Rscript
# =============================================================================
# build_union_grid_fleet.R -- UNION-mask iNat grid for any fleet species.
#
# Generalises ~/isdm/windowed_v2b/build_union_grid.R (2026-09-04), which was
# hardwired to moose/bobcat/WTD, to their cache files, and to the shared
# IUCN columns of data/moose_inat_grid.csv. The union RULE is unchanged:
# a cell50 is in range if it is inside the IUCN range OR clears the presence
# threshold (>= 1 post-isolation-filter record, PRESENCE_MIN_RECORDS = 1).
# Counts come from the presence grid wherever it has a value; the mask is
# carried as NA in the species column, never by dropping rows.
#
# ---- parameters (env) -------------------------------------------------------
# SPECIES       space- or comma-separated species tokens (common_name_clean).
# IUCN_SOURCE   where in-IUCN-range cell50s come from. Either
#                 - a CSV with a column per species (NA outside range), read
#                   once and selected by COLUMN -- today's source,
#                   data/moose_inat_grid.csv; or
#                 - a template containing {species}, resolved per species --
#                   for the consolidated per-species layer being built in
#                   geospatialdata/north_american_mammal_ranges/.
#               A species absent from the source is REFUSED (stop), not
#               special-cased. Vector layers are rasterised exactly as
#               Arielle's prep_inat_data_grid() did (iucn_source_lib.R);
#               validated 2026-09-17 by reproducing the shared-column cell sets
#               from the old extantonly shapefiles.
# PRES_GRID     presence-mask grid, path or {species} template. The v2 grid
#               writer emits every species' column in one file, so all
#               output/inat_grids_v2/*_inat_grid_v2.csv are byte-identical
#               (md5 63c163e8293b...) and any one of them serves a new species.
# UNMASKED_CACHE optional {species} template for the 2026-07-27 unmasked
#               cache. Used ONLY for IUCN-only cells (outside the presence
#               mask), reproducing the 09-04 build for moose/bobcat/WTD. When
#               absent, IUCN-only cells get 0: by construction they hold no
#               post-isolation-filter record in the presence grid's vintage.
# PINS_CSV      table of expected union/presence numbers. A pinned species
#               that deviates STOPS. An unpinned species is built, but a
#               REVIEW_REQUIRED flag is written next to its grid and the
#               downstream bundle build and fit runner refuse to proceed
#               while it exists.
# UNION_OUT_DIR output dir (default $PROJ/output/inat_grids_union).
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(terra) })

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
env  <- function(k, d = "") { v <- Sys.getenv(k, d); if (!nzchar(v)) stop("env ", k, " is required"); v }

SPECIES  <- strsplit(env("SPECIES"), "[ ,]+")[[1]]
IUCN_SRC <- env("IUCN_SOURCE", file.path(PROJ, "data", "moose_inat_grid.csv"))
PRES_T   <- env("PRES_GRID", file.path(PROJ, "output", "inat_grids_v2", "{species}_inat_grid_v2.csv"))
CACHE_T  <- Sys.getenv("UNMASKED_CACHE", "")
PINS_CSV <- env("PINS_CSV", "/home/rwkays/isdm/fleet_v2u/union_pins.csv")
OUT      <- env("UNION_OUT_DIR", file.path(PROJ, "output", "inat_grids_union"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

resolve <- function(tmpl, sp) gsub("{species}", sp, tmpl, fixed = TRUE)
pins <- if (file.exists(PINS_CSV)) fread(PINS_CSV) else data.table(species = character())
cat("IUCN_SOURCE :", IUCN_SRC, "\nPRES_GRID   :", PRES_T, "\nUNMASKED_CACHE:",
    if (nzchar(CACHE_T)) CACHE_T else "(none)", "\nPINS_CSV    :", PINS_CSV,
    "(", nrow(pins), "pinned )\nOUT         :", OUT, "\n")

# ---- IUCN source: see iucn_source_lib.R (csv columns, per-species layer, or file template)
source(Sys.getenv("IUCN_LIB", "/home/rwkays/isdm/fleet_v2u/iucn_source_lib.R"))
IUCN_KEY <- Sys.getenv("IUCN_KEY", "species")
taxon_key <- fread("/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/target_species.csv")
taxon_key[, common_name_clean := gsub("[ -']", "_", tolower(common_name))]
CW_PATH <- Sys.getenv("NAME_CROSSWALK", "")
if (nzchar(CW_PATH)) {                       # layer's IUCN name must agree with the crosswalk
  cw <- fread(CW_PATH); lyr_names <- as.data.table(as.data.frame(vect(IUCN_SRC)))
  for (sp in SPECIES) {
    a <- cw[species == sp, iucn_name]; b <- lyr_names[species == sp, iucn_mapped_name]
    if (length(a) != 1 || length(b) != 1 || a != b) stop("IUCN name mismatch for ", sp, ": crosswalk '", a, "' vs layer '", b, "'")
  }
  cat("crosswalk IUCN names agree with layer for:", SPECIES, "\n")
}
iucn_cells <- function(sp, universe) iucn_cells50(sp, IUCN_SRC, key = IUCN_KEY, taxon_key = taxon_key, universe = universe)

# ---- grid geometry for the edge report --------------------------------------
grid100 <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/grid100.tif")
to_lonlat <- function(x, y) {
  v <- project(vect(cbind(x, y), crs = crs(grid100)), "EPSG:4326")
  crds(v)
}

summ <- list()
for (sp in SPECIES) {
  cat("\n===== ", sp, " =====\n", sep = "")
  presence <- fread(resolve(PRES_T, sp))
  if (!sp %in% names(presence)) stop("presence grid has no column '", sp, "'")
  key <- c("Year", "cell50", "cell100")
  setkeyv(presence, key)

  iucn50 <- iucn_cells(sp, unique(presence$cell50))
  pres_cells <- unique(presence[!is.na(get(sp)), cell50])
  m <- merge(iucn50, data.table(cell50 = pres_cells, in_pres = TRUE), by = "cell50", all = TRUE)
  m[is.na(in_iucn), in_iucn := FALSE]; m[is.na(in_pres), in_pres := FALSE]
  m[, in_union := in_iucn | in_pres]
  union_cells <- m[in_union == TRUE, cell50]

  # ---- counts: presence grid first; IUCN-only cells from cache or 0 --------
  vp <- presence[[sp]]
  cache_f <- if (nzchar(CACHE_T)) resolve(CACHE_T, sp) else ""
  if (nzchar(cache_f) && file.exists(cache_f)) {
    unmasked <- fread(cache_f); setkeyv(unmasked, key)
    stopifnot("row keys differ" = identical(unmasked[, ..key], presence[, ..key]))
    vu <- unmasked[[sp]]; vu[is.na(vu)] <- 0L
    fill_src <- "07-27 unmasked cache"
  } else {
    vu <- rep(0L, length(vp)); fill_src <- "zero (no cache; IUCN-only cells hold no post-filter record)"
  }
  vals <- ifelse(is.na(vp), vu, vp)
  n_from_fill <- sum(is.na(vp) & presence$cell50 %in% union_cells & vu > 0)
  vals[!(presence$cell50 %in% union_cells)] <- NA
  out <- copy(presence); set(out, j = sp, value = vals)

  # ---- metrics ---------------------------------------------------------------
  grid_cells   <- presence[, .(cell100 = cell100[1], x50 = x50[1], y50 = y50[1]), by = cell50]
  n_grid50     <- nrow(grid_cells); n_grid100 <- uniqueN(grid_cells$cell100)
  c100 <- function(cells) uniqueN(grid_cells[cell50 %in% cells, cell100])
  got_cells    <- uniqueN(out[!is.na(get(sp)), cell50])
  union_recs   <- sum(out[[sp]], na.rm = TRUE)
  pres_recs    <- sum(presence[[sp]], na.rm = TRUE)
  stopifnot("union is not a superset of the presence data" = union_recs >= pres_recs,
            "emitted cell count != computed union"          = got_cells == length(union_cells))

  row <- data.table(
    species = sp,
    grid_cell50 = n_grid50, grid_cell100 = n_grid100,
    cell50_iucn = sum(m$in_iucn), cell50_presence = sum(m$in_pres), cell50_union = length(union_cells),
    cell100_iucn = c100(m[in_iucn == TRUE, cell50]), cell100_presence = c100(pres_cells),
    cell100_union = c100(union_cells),
    union_added_over_presence = length(union_cells) - sum(m$in_pres),
    presence_outside_iucn = sum(m$in_pres & !m$in_iucn),
    records_presence = pres_recs, records_union = union_recs,
    records_from_fill = n_from_fill, fill_source = fill_src,
    union_share_of_grid = round(length(union_cells) / n_grid50, 4),
    iucn_cells_dropped_outside_universe = if (is.null(attr(iucn50, "n_outside_universe"))) NA_integer_ else attr(iucn50, "n_outside_universe"))

  # ---- edge report: where does the union NOT reach? ------------------------
  gc <- merge(grid_cells, m, by = "cell50", all.x = TRUE)
  for (c in c("in_iucn", "in_pres", "in_union")) set(gc, which(is.na(gc[[c]])), c, FALSE)
  ll <- to_lonlat(gc$x50, gc$y50); gc[, `:=`(lon = ll[, 1], lat = ll[, 2])]
  gc[, cls := fifelse(in_iucn & in_pres, "both", fifelse(in_iucn, "IUCN only",
             fifelse(in_pres, "presence only", "outside union")))]
  fwrite(gc[, .(cell50, cell100, lon, lat, in_iucn, in_pres, in_union, cls)],
         file.path(OUT, paste0(sp, "_union_cells.csv")))
  outside <- gc[in_union == FALSE]
  row[, `:=`(outside_n = nrow(outside),
             lon_range_union = paste(round(range(gc[in_union == TRUE, lon]), 1), collapse = ".."),
             lat_range_union = paste(round(range(gc[in_union == TRUE, lat]), 1), collapse = ".."),
             lon_range_grid  = paste(round(range(gc$lon), 1), collapse = ".."),
             lat_range_grid  = paste(round(range(gc$lat), 1), collapse = ".."))]
  if (nrow(outside)) {
    outside[, `:=`(lat_band = cut(lat, seq(20, 55, 5)), lon_band = cut(lon, seq(-130, -60, 10)))]
    cat("  cells OUTSIDE the union, by 5-deg lat x 10-deg lon band:\n")
    print(dcast(outside, lat_band ~ lon_band, fun.aggregate = length, value.var = "cell50"))
  }
  png(file.path(OUT, paste0(sp, "_union_map.png")), width = 1400, height = 900, res = 130)
  cols <- c("both" = "#1b7837", "IUCN only" = "#9970ab", "presence only" = "#e08214", "outside union" = "#d9d9d9")
  plot(gc$lon, gc$lat, pch = 15, cex = 0.55, col = cols[gc$cls], asp = 1.25,
       xlab = "longitude", ylab = "latitude",
       main = sprintf("%s union mask: %d of %d cell50 (IUCN %d, presence %d)", sp,
                      length(union_cells), n_grid50, sum(m$in_iucn), sum(m$in_pres)))
  legend("bottomleft", legend = names(cols), col = cols, pch = 15, bty = "n")
  dev.off()

  print(t(row))

  # ---- pins: hard stop if pinned and different; flag if unpinned ------------
  flag <- file.path(OUT, paste0(sp, "_REVIEW_REQUIRED"))
  pin <- pins[species == sp]
  if (nrow(pin)) {
    bad <- c(cell50_union = pin$cell50_union != row$cell50_union,
             records_union = pin$records_union != row$records_union,
             cell50_presence = pin$cell50_presence != row$cell50_presence,
             records_presence = pin$records_presence != row$records_presence)
    if (any(bad)) stop("UNION PIN FAILED for ", sp, ": ", paste(names(bad)[bad], collapse = ", "),
                       " differ from ", PINS_CSV, " -- refusing to write the grid")
    cat("  PIN PASSED (", paste(names(bad), collapse = ", "), ")\n")
    if (file.exists(flag)) file.remove(flag)
  } else {
    writeLines(c(paste("species:", sp), paste("built:", format(Sys.time())),
                 "No pinned expectation exists. Numbers above must be reviewed and",
                 "added to the pins table before a bundle is built from this grid."), flag)
    cat("  UNPINNED -- wrote", flag, "; downstream bundle build will refuse until reviewed\n")
  }

  f <- file.path(OUT, paste0(sp, "_inat_grid_union.csv"))
  fwrite(out, f)
  cat("  wrote", f, "\n")
  summ[[sp]] <- row
}
res <- rbindlist(summ, fill = TRUE)
sf <- file.path(OUT, paste0("union_grid_report_", format(Sys.time(), "%Y%m%d_%H%M"), ".csv"))
fwrite(res, sf); cat("\nreport:", sf, "\nDONE\n")
