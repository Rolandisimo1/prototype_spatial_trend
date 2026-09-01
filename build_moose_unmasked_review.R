#!/usr/bin/env Rscript
# Read-only: TRUE (unmasked) per-cell50 moose observation counts, regardless
# of range-mask status -- the masked pull's out-of-range total_count is 0 by
# construction (moose_inat_grid.csv already has NA baked in for out-of-range
# cells before it's written), so it can't show how much real moose activity
# falls outside the mask. This recomputes counts directly from the master
# raw file (same method as the fleet35 script's fresh-extraction path),
# skipping the NA-masking step entirely, then attaches in_range as a label
# only (not a filter). No pipeline change, no bundle touched.

t_start <- Sys.time()
log_step <- function(msg) { cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg)); flush(stdout()) }
log_step("script start")

suppressPackageStartupMessages({
  library(data.table)
  library(terra)
})
log_step("packages loaded")

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
DATA_DIR <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data"
CACHE_DIR <- file.path(PROJ, "data")
RANGE_DIR <- "/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/iucn_ranges_extantonly"
RAW_FILE <- file.path(DATA_DIR, "inat_combo_nam_mams.csv")
OUT_DIR <- file.path(PROJ, "prototype_spatial_trend")

# ------------------------------ base grid skeleton (all 3322 active cells) ------
log_step("reading moose_inat_grid.csv for the base cell50/x50/y50 skeleton + in_range flag")
base <- fread(file.path(CACHE_DIR, "moose_inat_grid.csv"), select = c("cell50","x50","y50","moose"))
skeleton <- base[, .(x50 = x50[1], y50 = y50[1],
                       in_range = !all(is.na(moose))), by = cell50]
log_step(paste("skeleton:", nrow(skeleton), "distinct cell50"))

# ------------------------------ TRUE unmasked moose counts, fresh from raw file -
log_step("loading grid50.tif / grid100.tif")
grid50 <- rast(file.path(DATA_DIR, "grid50.tif"))
grid100 <- rast(file.path(DATA_DIR, "grid100.tif"))
log_step("grids loaded")

log_step("reading + filtering master file for Alces alces")
dt <- fread(RAW_FILE, select = c("observed_on","public_positional_accuracy",
                                    "taxon_species_name","latitude","longitude"))
dt[, observed_on := as.Date(observed_on)]
dt <- dt[public_positional_accuracy < 1000]
dt[, Year := year(observed_on)]
dt <- dt[Year >= 2008 & Year <= 2025]
dt <- dt[taxon_species_name == "Alces alces"]
log_step(paste(nrow(dt), "raw Alces alces observations (accuracy+year filtered, NO range mask applied)"))

t0 <- Sys.time()
pts <- vect(dt, geom = c("longitude","latitude"), crs = "+proj=longlat")
pts <- project(pts, crs(grid100))
log_step(paste("point reprojection took", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec"))

t0 <- Sys.time()
cell50_vals <- terra::extract(grid50, pts)[, 2]
log_step(paste("grid50 extract took", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "sec"))

dt[, cell50 := cell50_vals]
dt_matched <- dt[!is.na(cell50)]
log_step(paste(nrow(dt_matched), "of", nrow(dt), "observations matched to a grid50 cell"))

true_counts <- dt_matched[, .(total_count = .N), by = cell50]

# ------------------------------ merge onto full skeleton (0-fill unobserved) ----
out <- merge(skeleton, true_counts, by = "cell50", all.x = TRUE)
out$total_count[is.na(out$total_count)] <- 0

# ------------------------------ lon/lat ------------------------------------------
log_step("reprojecting cell50 centroids to EPSG:4326")
centroid_pts <- vect(out, geom = c("x50", "y50"), crs = crs(grid50))
centroid_ll <- project(centroid_pts, "EPSG:4326")
ll <- geom(centroid_ll)[, c("x", "y")]
out$lon <- ll[, "x"]
out$lat <- ll[, "y"]

out_final <- out[, .(x50, y50, lon, lat, total_count, in_range)]
OUT_PATH <- file.path(OUT_DIR, "moose_cell_counts_unmasked.csv")
fwrite(out_final, OUT_PATH)
log_step(paste("wrote", OUT_PATH, "(", nrow(out_final), "rows )"))

# ------------------------------ summary, in-range vs out-of-range, TRUE counts --
in_r <- out_final[in_range == TRUE]
out_r <- out_final[in_range == FALSE]
total_obs <- sum(out_final$total_count)

cat("\n=== TRUE (unmasked) summary ===\n")
cat("total raw Alces alces observations matched to grid:", total_obs, "\n")
cat("in-range: ", sum(in_r$total_count), sprintf(" (%.2f%%)\n", 100*sum(in_r$total_count)/total_obs))
cat("out-of-range: ", sum(out_r$total_count), sprintf(" (%.2f%%)\n", 100*sum(out_r$total_count)/total_obs))
cat("\ndistinct cells in-range:", nrow(in_r), " with >=1 obs:", sum(in_r$total_count > 0), "\n")
cat("distinct cells out-of-range:", nrow(out_r), " with >=1 obs:", sum(out_r$total_count > 0), "\n")

summarize_counts <- function(x) c(min=min(x), median=median(x), mean=round(mean(x),3), max=max(x), n=length(x))
cat("\nper-cell total_count, IN-RANGE occupied cells:\n"); print(summarize_counts(in_r$total_count[in_r$total_count>0]))
cat("per-cell total_count, OUT-OF-RANGE occupied cells:\n"); print(summarize_counts(out_r$total_count[out_r$total_count>0]))

log_step(paste("TOTAL script time:", round(as.numeric(Sys.time() - t_start, units = "secs"), 1), "sec"))
cat("\ndone.\n")
