#!/usr/bin/env Rscript
# =============================================================================
# 00_prep_sim_inputs.R
# Read-only prep step for the spatially-varying-trend simulation validation.
# Pulls from Arielle's REAL, converged bobcat fit (never modifies it) the
# pieces needed to build a reduced-scale-but-structurally-faithful simulation:
#   (1) cell100 centroid lon/lat (from grid100.tif), to assign a NE/West/other
#       region label per CAR cell for the known true year_effect field.
#   (2) the real 908-cell CAR adjacency graph (adj/num/nadj/nnum) -- kept
#       EXACTLY as-is in the simulation (this is cheap; it is what year_effect
#       is actually estimated on).
#   (3) subcell-count-per-cell50 distribution, to size a tractable subsample.
#   (4) posterior-mean values for every NON-trend parameter from the real
#       chain_bobcat_{1,2,3}.RDS (occ_beta, MWMT/MCMT effects, intercepts,
#       detection params, theta0/theta1, overdisp_inat) so the simulation is
#       realistic in magnitude everywhere except the trend field under test.
#
# Run with plotting_env's Rscript (terra + coda/MCMCvis, no nimble needed):
#   $PROJ/HPC/conda_envs/plotting_env/bin/Rscript 00_prep_sim_inputs.R
# Writes: prepped_sim_inputs.RDS (in this directory)
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(coda); library(MCMCvis); library(dplyr)
})

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
DATA_DIR <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data"
BOBCAT_DIR <- paste0(PROJ, "/HPC/bobcat")

# --- (1) cell100 centroids + region label ------------------------------------
grid100 <- rast(paste0(DATA_DIR, "/grid100.tif"))
pts <- as.points(grid100)               # one point per non-NA cell, value = cell100 id
pts_ll <- project(pts, "EPSG:4326")
cell100_geo <- data.frame(cell100 = values(pts)[,1],
                          lon = crds(pts_ll)[,1], lat = crds(pts_ll)[,2])
cell100_geo <- cell100_geo[order(cell100_geo$cell100), ]
stopifnot(all(diff(cell100_geo$cell100) == 1), cell100_geo$cell100[1] == 1)

# Region assignment for the TRUE simulated trend field:
#   NE: lon > -85 & lat > 37   (positive trend)
#   West: lon < -100            (negative trend)
#   other: everything else      (~0, smoothed by CAR neighbors in between)
cell100_geo$region <- with(cell100_geo, ifelse(lon > -85 & lat > 37, "NE",
                                        ifelse(lon < -100, "West", "other")))
cat("region counts (of", nrow(cell100_geo), "cell100 cells):\n")
print(table(cell100_geo$region))

# --- (2) real CAR graph + cell50->cell100 map + subcell distribution --------
# NOTE: input_data_bobcat.RDS (NOT _reduced) is the one that matches the real
# converged chain_bobcat_{1,2,3}.RDS -- confirmed by occ_beta count (17, not
# 10). input_data_bobcat_reduced.RDS appears to be an orphaned artifact from
# unrelated earlier work and does NOT correspond to these chains.
inp <- readRDS(paste0(BOBCAT_DIR, "/input_data_bobcat.RDS"))
cl  <- inp$constants_list

subcell_n <- cl$inat_cell50_end - cl$inat_cell50_start + 1
cat("\nsubcells per cell50: n=", length(subcell_n),
    " median=", median(subcell_n), " mean=", round(mean(subcell_n),1),
    " range=", paste(range(subcell_n), collapse="-"), "\n")

cell_map <- readRDS(paste0(PROJ, "/cell_maps/cell_map_50_100_bobcat.RDS"))
stopifnot(nrow(cell_map) == cl$ncell50)
# sanity: cell_map's cell100 per cell50 should match constants_list$inat_cell100
stopifnot(all(cell_map$cell100 == cl$inat_cell100))

cat("\nsites per cell100 (camera): \n")
print(summary(as.numeric(table(cl$cell))))

# --- (3) real posterior means for non-trend params ---------------------------
chain_paths <- paste0(BOBCAT_DIR, "/chain_bobcat_", 1:3, ".RDS")
stopifnot(all(file.exists(chain_paths)))
chains <- lapply(chain_paths, readRDS)
samples <- do.call(rbind, lapply(chains, function(x) x$samples))
cat("\ncombined real posterior draws:", nrow(samples), "x", ncol(samples), "\n")

pmean <- function(pattern) {
  cols <- grep(pattern, colnames(samples), value = TRUE)
  colMeans(samples[, cols, drop = FALSE])
}
real_post_means <- list(
  occ_beta            = pmean("^occ_beta\\["),
  link_occ_intercept  = pmean("^link_occ_intercept\\["),
  MWMT_effect         = pmean("^MWMT_effect\\["),
  MCMT_effect         = pmean("^MCMT_effect\\["),
  intercept_tau       = mean(samples[, "intercept_tau"]),
  MWMT_tau            = mean(samples[, "MWMT_tau"]),
  MCMT_tau            = mean(samples[, "MCMT_tau"]),
  det_intercept       = if ("det_intercept" %in% colnames(samples)) mean(samples[, "det_intercept"]) else NA,
  p_beta              = pmean("^p_beta\\["),
  theta0              = mean(samples[, "theta0"]),
  theta1              = mean(samples[, "theta1"]),
  overdisp_inat       = mean(samples[, "overdisp_inat"]),
  year_beta           = mean(samples[, "year_beta"]),
  year_var            = mean(samples[, "year_var"])
)
cat("\nreal posterior means (non-trend params, for realistic sim magnitudes):\n")
str(real_post_means)

stopifnot(cl$has_SVC || cl$hasSVC)  # sim/model fork assumes SVC branch, per task scope

out <- list(cell100_geo = cell100_geo, constants_list = cl,
           subcell_n = subcell_n, real_post_means = real_post_means,
           numOccCovars = cl$numOccCovars, nyear = cl$nyear, year_vals = cl$year_vals)
saveRDS(out, "prepped_sim_inputs.RDS")
cat("\nwrote prepped_sim_inputs.RDS\n")
