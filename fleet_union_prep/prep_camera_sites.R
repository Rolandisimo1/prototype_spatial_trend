# =============================================================================
# prep_camera_sites.R -- camera-side prep with a PARAMETERISED range test.
#
# Fork of prep_spoccupancy_data() in ~/isdm/v2_prep/prep_data_for_spoccupancy.R.
# That function's cached outputs (data/forSPO/<species>_4SPO.RDS, April-June
# 2026) are what every fitted bundle used, and they were NOT produced by its
# current body: an exact uid-set comparison (job 858979, moose/bobcat/WTD/
# coyote) shows the caches kept   in-polygon (matched by deployment_id)
#                                 & complete siteCovs & complete raster covs
# with NO single-occasion filter. Rerunning today's function would drop ~17%
# more sites. So this fork:
#   - keeps every step of the body verbatim EXCEPT
#   - the single-occasion filter is OFF (matches the caches), and
#   - the range test is `range_mode`:
#       "legacy_depid"  old extantonly polygon, hits mapped back by
#                       deployment_id (non-unique) -- reproduces the caches;
#                       used only for validation
#       "iucn_site"     old polygon, tested per site
#       "union_cell"    site's cell50 is in the union mask read from the SAME
#                       union grid CSV the iNat side is built from, so camera
#                       and iNat range definitions are identical (approved 9/17)
#   - species names come from species_name_crosswalk.csv (wi_name keys umflist)
#   - output goes to data/forSPO_<tag>/, never over data/forSPO/
# =============================================================================
suppressPackageStartupMessages({ library(tidyverse); library(terra); library(sf); library(data.table) })

CAM_PROJ  <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
CAM_DATA  <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data"
CROSSWALK <- Sys.getenv("NAME_CROSSWALK", "/home/rwkays/isdm/fleet_v2u/species_name_crosswalk.csv")

cam_static_covars <- c("terrain_ruggedness", "soil_clay", "soil_silt", "soil_sand", "elevation")
cam_temp_covars   <- c("MWMT", "MCMT")
cam_annual_covars <- c("CWD","Human_pop","NDVI_mean","NDVI_sd","Ag","Deciduous",
                       "Evergreen","Impervious","Mixed","PDSI","PPT","Temp")

cam_add_latitude_layer <- function(raster_brick) {        # verbatim add_latitude_layer
  temp_ras <- raster_brick[[1]]
  terra::values(temp_ras) <- 0
  pts <- vect(crds(temp_ras, df = TRUE), geom = c("x", "y"), crs = crs(temp_ras))
  pts_ll <- project(pts, "EPSG:4326")
  lat <- as.data.frame(pts_ll, geom = "XY")$y
  terra::values(temp_ras) <- lat
  names(temp_ras) <- "latitude"
  temp_ras
}

cam_names <- function(species) {
  key <- species                     # never compare against bare `species` inside cw[...]: that is the column
  cw <- fread(CROSSWALK)
  r <- cw[cw[["species"]] == key]
  if (nrow(r) != 1) stop("species '", species, "' not in name crosswalk ", CROSSWALK)
  r
}

cam_old_rangefile <- function(sci) {                     # as taxon_key$rangefile
  s <- gsub(" ", "_", sci); if (s == "Pekania_pennanti") s <- "Martes_pennanti"
  f <- list.files(file.path("/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/iucn_ranges_extantonly"),
                  pattern = paste0(s, ".*.shp"), full.names = TRUE)
  if (length(f) != 1) stop("expected 1 old range shapefile for ", sci, ", found ", length(f))
  f
}

prep_camera_sites <- function(species, range_mode, union_grid = NULL, tag = range_mode,
                              out_dir = file.path(CAM_PROJ, "data", paste0("forSPO_", tag))) {
  stopifnot(range_mode %in% c("legacy_depid", "iucn_site", "union_cell"))
  nm  <- cam_names(species)
  data_outfile <- file.path(out_dir, paste0(species, "_4SPO.RDS"))
  if (file.exists(data_outfile)) stop("refusing to overwrite ", data_outfile)

  umflist <- readRDS(file.path(CAM_DATA, "umflist.RDS"))
  dat <- umflist[[nm$wi_name]]
  if (is.null(dat)) stop("umflist has no entry '", nm$wi_name, "'")

  reshape_to_3d_array <- function(x, y) {
    nrow_y <- nrow(y); ncol_y <- ncol(y)
    result_array <- array(NA, dim = c(nrow_y, ncol_y, ncol(x)))
    for (i in seq_len(ncol(x))) result_array[, , i] <- matrix(x[, i], nrow = nrow_y, ncol = ncol_y, byrow = TRUE)
    dimnames(result_array)[[3]] <- colnames(x)
    result_array
  }
  # Slots read as attributes: plotting_env may not have `unmarked`, and S4 slot
  # assignment needs the class definition. Values are identical.
  siteCovs <- attr(dat, "siteCovs"); y_all <- attr(dat, "y"); obsCovs <- attr(dat, "obsCovs")
  obs_dat2  <- reshape_to_3d_array(obsCovs, y_all)
  site_pts2 <- vect(siteCovs, geom = c("longitude", "latitude"), crs = "+proj=longlat")

  # ---- THE ONLY CHANGED STEP: range test --------------------------------------
  grid100 <- rast(file.path(CAM_DATA, "grid100.tif"))
  if (range_mode == "legacy_depid") {
    range_poly <- terra::vect(st_read(cam_old_rangefile(nm$wi_name), quiet = TRUE))
    site_pts <- site_pts2[range_poly, ]
    deps_in_range <- as.character(site_pts$deployment_id)
    ind <- which(site_pts2$deployment_id %in% deps_in_range)
  } else if (range_mode == "iucn_site") {
    range_poly <- terra::vect(st_read(cam_old_rangefile(nm$wi_name), quiet = TRUE))
    ind <- which(is.related(site_pts2, range_poly, "intersects"))
  } else {
    stopifnot(!is.null(union_grid), file.exists(union_grid))
    ug <- fread(union_grid, select = c("cell50", species))
    union_cells <- unique(ug[!is.na(ug[[species]]), cell50])
    grid50 <- rast(file.path(CAM_DATA, "grid50.tif"))
    site_cell50 <- extract(grid50, project(site_pts2, crs(grid100)))[, 2]
    ind <- which(!is.na(site_cell50) & site_cell50 %in% union_cells)
    cat("  union_cell range test:", length(ind), "of", length(site_cell50), "sites in",
        length(union_cells), "union cell50;", sum(is.na(site_cell50)), "sites off the grid\n")
  }
  range_n <- length(ind)
  # -----------------------------------------------------------------------------

  obs_dat <- obs_dat2[ind, , ]
  siteCovs <- siteCovs[ind, ]
  y_all <- y_all[ind, ]
  site_pts <- site_pts2[ind, ]

  raster_brick <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/covar_raster_brick.tif")
  raster_brick[["latitude_cov"]] <- cam_add_latitude_layer(raster_brick)
  covar_mtx <- extract(raster_brick, site_pts) %>% select(-ID)
  covar_mtx <- covar_mtx %>% select(all_of(c(cam_static_covars, cam_temp_covars)))

  inds1 <- which(rowSums(is.na(siteCovs)) == 0 & rowSums(is.na(covar_mtx)) == 0)
  # Single-occasion filter DELIBERATELY OFF: the cached 4SPO files behind every
  # fitted bundle never applied it (job 858979).
  inds_to_keep <- inds1

  coords_as_unique <- siteCovs[inds_to_keep, ] %>% ungroup() %>%
    distinct(longitude, latitude) %>% mutate(loc_ID = row_number())
  occ.covs.static <- siteCovs %>% bind_cols(covar_mtx) %>% .[inds_to_keep, ] %>%
    ungroup() %>% left_join(coords_as_unique) %>% mutate(order = seq(1, n()))
  occ.covs.split <- split(occ.covs.static, occ.covs.static$year)
  occ_covs_df2 <- list()
  for (i in 1:length(names(occ.covs.split))) {
    occ_covs_df <- occ.covs.split[[i]]
    year_i <- names(occ.covs.split)[i]
    site_pts_year <- vect(occ_covs_df, geom = c("longitude", "latitude"), crs = "+proj=longlat")
    rast_i <- rast(paste0(CAM_DATA, "/Raster_bricks/raster_brick_", year_i, ".tif"))
    covar_mtx2 <- extract(rast_i, site_pts_year) %>% select(-ID)
    occ_covs_df2[[i]] <- occ_covs_df %>% bind_cols(covar_mtx2)
  }
  occ.covs <- do.call("rbind", occ_covs_df2)
  occ.covs <- occ.covs %>%
    mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .))) %>%
    arrange(order)

  covars_toscale <- c(cam_static_covars, cam_temp_covars, cam_annual_covars)
  covars_toscale <- covars_toscale[!grepl("_sq$", covars_toscale)]
  scaling_factors <- data.frame(covar = covars_toscale,
                                mean = colMeans(occ.covs[, covars_toscale]),
                                sd = apply(occ.covs[, covars_toscale], 2, sd), type = "site")
  for (i in 1:length(covars_toscale)) {
    occ.covs[, scaling_factors$covar[i]] <-
      (occ.covs[, scaling_factors$covar[i]] - scaling_factors$mean[i]) / scaling_factors$sd[i]
  }
  occ.covs$MCMT_sq <- occ.covs$MCMT^2
  occ.covs$MWMT_sq <- occ.covs$MWMT^2

  det.covs <- list(
    yday_scaled = obs_dat[inds_to_keep, , "yday_scaled"],
    canopy_height_scaled = obs_dat[inds_to_keep, 1, "canopy_height_scaled"],
    log_roaddist_scaled = obs_dat[inds_to_keep, 1, "log_roaddist_scaled"])
  det.covs$yday_scaled_sq <- det.covs$yday_scaled^2

  spocc_dat <- list(y = y_all[inds_to_keep, ], occ.covs = occ.covs, det.covs = det.covs,
                    coords = as.matrix(coords_as_unique[, 1:2]), grid.index = occ.covs$loc_ID,
                    species = species, scaling_factors = scaling_factors)
  spocc_dat$y[spocc_dat$y > 0] <- 1

  adj <- numeric()
  num <- numeric(max(terra::values(grid100), na.rm = T))
  for (i in 1:max(terra::values(grid100), na.rm = T)) {
    cell_id <- which(terra::values(grid100) == i)
    neighbors <- as.numeric(adjacent(grid100, cell_id))
    neighbor_cell_IDs <- terra::values(grid100)[neighbors]
    neighbor_cell_IDs <- neighbor_cell_IDs[!is.na(neighbor_cell_IDs)]
    num[i] <- length(neighbor_cell_IDs)
    adj <- c(adj, neighbor_cell_IDs)
  }
  spocc_dat$continent_adj_info <- list(adj = adj, num = num)
  spocc_dat$occ.covs$spatcell <- extract(grid100, vect(spocc_dat$occ.covs[, c("longitude", "latitude")],
                                                       geom = c("longitude", "latitude"), crs = "+proj=longlat"))[, 2]
  spocc_dat$grid100 <- grid100
  spocc_dat$camera_prep <- list(range_mode = range_mode, union_grid = union_grid, names = as.list(nm),
                                sites_in_range = range_n, sites_kept = length(inds_to_keep),
                                built = format(Sys.time()))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  saveRDS(spocc_dat, data_outfile)
  cat("  camera prep", species, range_mode, ": in range", range_n, " kept", length(inds_to_keep),
      " detecting", sum(rowSums(spocc_dat$y, na.rm = TRUE) > 0), " ->", data_outfile, "\n")
  spocc_dat
}
