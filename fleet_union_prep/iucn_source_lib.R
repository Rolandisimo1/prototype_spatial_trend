# =============================================================================
# iucn_source_lib.R -- in-IUCN-range cell50 sets, from any supported source.
#
# iucn_cells50(sp, src, key = "species", taxon_key = NULL) returns a
# data.table(cell50, in_iucn) over every cell50 of the national grid.
#
# Supported sources:
#   *.csv            shared grid CSV with one column per species, NA outside
#                    range (data/moose_inat_grid.csv). Selected by COLUMN.
#   *.geojson/*.gpkg one layer, one feature per species, matched on attribute
#                    `key` (fleet34_iucn_ranges_present_day.geojson: key
#                    "species" = pipeline name, so taxonomy overrides such as
#                    nine-banded_armadillo -> Dasypus mexicanus are already
#                    resolved in the layer and never re-keyed here).
#   *.shp template   per-species files; "{species}" -> pipeline name,
#                    "{sci_}" -> taxon_key sci_name with "_" (old extantonly set).
# A species with no match is an ERROR, never an empty mask.
# universe: cell50 IDs of the grid the result will be joined to.
#
# Rasterisation replicates Arielle's prep_inat_data_grid() exactly:
#   vect() -> simplifyGeom(0.01) -> project(crs(grid100)) ->
#   rasterize(range, grid50, touches = TRUE); a cell50 is in range if the
#   rasterised layer is non-NA there.
# =============================================================================
suppressPackageStartupMessages({ library(data.table); library(terra) })

GRID50_PATH  <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/grid50.tif"
GRID100_PATH <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/grid100.tif"

.iucn_cache <- new.env()

iucn_cells50 <- function(sp, src, key = "species", taxon_key = NULL, universe = NULL) {
  grid50 <- rast(GRID50_PATH); grid100 <- rast(GRID100_PATH)
  all50 <- values(grid50)[, 1]; all50 <- sort(all50[!is.na(all50)])
  # grid50.tif has 3,364 cells; the pipeline grids only carry the 3,322 with any
  # iNat effort (same IDs and coordinates, verified 2026-09-17). Restrict to the
  # caller's universe so a polygon cannot add cells that have no grid rows.
  if (!is.null(universe)) all50 <- all50[all50 %in% universe]

  if (grepl("\\.csv$", src, ignore.case = TRUE)) {
    if (is.null(.iucn_cache[[src]])) .iucn_cache[[src]] <- fread(src)
    dt <- .iucn_cache[[src]]
    if (!sp %in% names(dt)) stop("species '", sp, "' has no column in IUCN source ", src)
    out <- dt[, .(in_iucn = !all(is.na(get(sp)))), by = cell50]
    return(if (is.null(universe)) out else out[cell50 %in% universe])
  }

  if (grepl("\\{(species|sci_)\\}", src)) {
    f <- gsub("{species}", sp, src, fixed = TRUE)
    if (grepl("{sci_}", f, fixed = TRUE)) {
      stopifnot(!is.null(taxon_key))
      sci <- taxon_key$sci_name[taxon_key$common_name_clean == sp]
      if (length(sci) != 1) stop("no unique sci_name for ", sp)
      f <- gsub("{sci_}", gsub(" ", "_", sci), f, fixed = TRUE)
    }
    if (!file.exists(f)) stop("IUCN file not found for ", sp, ": ", f)
    range <- vect(f)
    label <- f
  } else {
    if (!file.exists(src)) stop("IUCN source not found: ", src)
    if (is.null(.iucn_cache[[src]])) .iucn_cache[[src]] <- vect(src)
    lyr <- .iucn_cache[[src]]
    if (!key %in% names(lyr)) stop("IUCN layer has no attribute '", key, "'")
    range <- lyr[lyr[[key]][, 1] == sp, ]
    if (nrow(range) != 1) stop("expected exactly 1 feature with ", key, " == '", sp,
                               "' in ", src, ", found ", nrow(range))
    label <- paste0(src, " [", key, "=", sp, "]")
  }
  range <- project(simplifyGeom(range, 0.01), crs(grid100))
  r <- rasterize(range, grid50, touches = TRUE)
  in_ids <- values(grid50)[!is.na(values(r)[, 1]), 1]
  in_ids <- in_ids[!is.na(in_ids)]
  n_outside_universe <- if (is.null(universe)) 0L else sum(!in_ids %in% universe)
  if (!length(in_ids)) stop("IUCN polygon for ", sp, " touches 0 grid cells (", label, ") -- refusing an empty mask")
  out <- data.table(cell50 = all50, in_iucn = all50 %in% in_ids)
  attr(out, "n_outside_universe") <- n_outside_universe
  message(sprintf("  [iucn_cells50] %s: %d in-range raster cells, %d dropped as outside the %s-cell grid universe",
                  sp, length(in_ids), n_outside_universe, if (is.null(universe)) "full" else length(universe)))
  out
}
