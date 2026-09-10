# integration_helper.R -- PATCHED for Hazel (/rsstu paths); iNat grid writes -> our space
if (!exists("PROJ_DIR"))    PROJ_DIR    <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
if (!exists("DATA_DIR"))    DATA_DIR    <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data"
if (!exists("GRID50_PATH")) GRID50_PATH <- paste0(DATA_DIR, "/grid50.tif")

# Functions needed for running the Nimble model and prepping iNat dataset

#Set the spatial covariate names
annual_covars <- c("CWD","Human_pop","NDVI_mean","NDVI_sd","Ag","Deciduous",
                   "Evergreen","Impervious","Mixed","PDSI","PPT","Temp")

static_covars <- c("terrain_ruggedness",
                   "soil_clay",
                   "soil_silt",
                   "soil_sand",
                   "elevation")

#Set the time-varying spatial covariates
temp_covars <- c("MWMT", "MCMT")

# Nimble functions for intensity using SVCs
calcIntensity_SVC <- nimbleFunction(run = function(
    intensity_intercept = double(0),
    theta0 = double(0),
    theta1 = double(0),
    MWMT_effect = double(0),
    MCMT_effect = double(0),
    beta   = double(1),
    xdat   = double(2),
    MWMT_dat = double(1),
    MCMT_dat = double(1),
    year_dat = double(0),
    total_var_beta = double(0)) {
  
  log_lambda <- intensity_intercept + 
    (xdat %*% matrix(beta, ncol = 1)) + 
    (MWMT_dat * MWMT_effect) + 
    (MCMT_dat * MCMT_effect) +
    (total_var_beta * year_dat)
  
  mu <- exp(theta0 + theta1 * log(sum(exp(log_lambda))))
  
  return(mu)
  returnType(double(0))
})

calcIntensity_noSVC <- nimbleFunction(run = function(
    intensity_intercept = double(0),
    theta0 = double(0),
    theta1 = double(0),
    beta   = double(1),
    xdat   = double(2),
    year_dat = double(0),
    total_var_beta = double(0)) {
  
  log_lambda <- intensity_intercept + 
    (xdat %*% matrix(beta, ncol = 1)) +
    (total_var_beta * year_dat)
  
  mu <- exp(theta0 + theta1 * log(sum(exp(log_lambda))))
  
  return(mu)
  returnType(double(0))
})

################## Functions for iNat data prep ################################

# Make the cell x year matrix for model input
make_inat_cell_year_matrix <- function(df, species) {
  
  mat_df <- df %>%
    select(cell50, Year, all_of(species)) %>%
    pivot_wider(
      names_from = Year,
      values_from = all_of(species),
      values_fill = NA
    ) %>%
    arrange(cell50)
  
  mat <- as.matrix(mat_df[,-1])
  rownames(mat) <- mat_df$cell50
  
  return(mat)
}

# Make the effort matrix for model input
make_inat_effort_matrix <- function(df) {
  # df must have columns: cell50, Year, effort
  
  effort_df <- df %>%
    group_by(cell50, Year) %>%
    summarise(effort_sum = sum(effort, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = Year,
      values_from = effort_sum,
      values_fill = 0  # fill zeros where no effort
    ) %>%
    arrange(cell50)
  
  # convert to matrix
  mat <- as.matrix(effort_df[,-1])
  rownames(mat) <- effort_df$cell50
  
  return(mat)
}

# Make the annual rangewide covariate matrix for model input - this is also used for prediction
make_inat_xdat <- function(inat_df, real_data, species) {
  
  # Using Ben's orig raster brick for grid scaling and static covs
  raster_brick <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/covar_raster_brick.tif")
  raster_brick[["latitude_cov"]] <- add_latitude_layer(raster_brick)
  
  # Just get the layers needed
  selected_layers <- c(static_covars, temp_covars)
  raster_brick <- raster_brick[[selected_layers]]
  
  grid100 <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/grid100.tif")
  grid50  <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/grid50.tif")
  
  # Get array dims
  rast1<-rast("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/Raster_bricks/raster_brick_2008.tif")
  rast1 <- resample(rast1, raster_brick, method = "bilinear")
  
  #Make a raster mask 
  mask <- mask(raster_brick["soil_clay"], rast1["NDVI_mean"])
  
  #Apply the mask across the bricks
  masked_brick <- mask(raster_brick, mask)
  masked_rast1 <- mask(rast1, mask)
  
  # Combine the rasters
  combo_brick <- c(masked_brick, masked_rast1)
  
  # Figure out which cell100 each raster_brick cell is in. 
  grid50_downscaled <- resample(grid50, combo_brick, method = "mode")
  terra::values(grid50_downscaled)[is.na(terra::values(grid50_downscaled))] <- 0
  terra::values(grid50_downscaled) <- ifelse(
    terra::values(grid50_downscaled) %in% inat_df$cell50,
    terra::values(grid50_downscaled),
    NA)
  
  combo_brick[["cell50"]] <- grid50_downscaled
  x <- as.data.frame(combo_brick)
  dim(x)
  
  # Build out the xdat matrix for inat data. Need one row per
  # 5 km brick cell contained inside a 50 km inat cell
  # for each year
  years <- 2008:2025

  xdat_inat_df_yr <- array(NA,
                           dim = c(dim(x)[1], dim(x)[2]+2, length(years)),
                           dimnames = list(1:dim(x)[1],
                                           c(colnames(x),"MCMT_sq","MWMT_sq"),
                                           years))
  
  for(i in 1:length(years)){
    
    year_i <- years[i]
    
    rast_i<-rast(paste0("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/Raster_bricks/raster_brick_", year_i, ".tif"))
    
    # Align with raster_brick
    rast_i2 <- resample(rast_i, raster_brick, method = "bilinear")
    
    #Apply the mask across the brick
    masked_rasti2 <- mask(rast_i2, mask)
    
    # Combine the two raster bricks
    combined_brick <- c(masked_brick, masked_rasti2)
    
    # Figure out which cell100 each raster_brick cell is in. 
    grid50_downscaled <- resample(grid50, combined_brick, method = "mode")
    terra::values(grid50_downscaled)[is.na(terra::values(grid50_downscaled))] <- 0
    terra::values(grid50_downscaled) <- ifelse(
      terra::values(grid50_downscaled) %in% inat_df$cell50,
      terra::values(grid50_downscaled),
      NA)
    
    #stopifnot(all(table(as.numeric(terra::values(grid50_downscaled))) == 100))
    
    combined_brick[["cell50"]] <- grid50_downscaled
    
    xdat_inat_df <- as.data.frame(combined_brick)
    
    xdat_mat <- as.matrix(xdat_inat_df)
    
    xdat_inat_df_yr[,1:ncol(xdat_mat),i] <- xdat_mat

  }  

  # Filter across the array to remove cells with NAs in cell50 or NDVI_mean in any year
  cell50_col <- which(colnames(xdat_inat_df_yr) == "cell50")
  ndvi_col   <- which(colnames(xdat_inat_df_yr) == "NDVI_mean")
  
  # Create a logical vector for rows to keep:
  keep_rows <- apply(xdat_inat_df_yr[, c(cell50_col, ndvi_col), ], 1, function(x) all(!is.na(x)))
  
  # Subset the array
  xdat_inat_df_yr_filtered <- xdat_inat_df_yr[keep_rows, , ]
  dim(xdat_inat_df_yr_filtered)
  
  # Scale the covariates
  covars_toscale <- c(static_covars, temp_covars, annual_covars)
  covars_toscale <- covars_toscale[!covars_toscale %in% c('soil_clay','soil_silt','soil_sand')]
  xdat_inat_df_yr_scaled <- xdat_inat_df_yr_filtered
  for (ii in seq_along(covars_toscale)) {

    covar_name <- covars_toscale[ii]
    col_idx <- which(colnames(xdat_inat_df_yr_filtered)==covar_name)
    
    # Loop over years
    for (yr in 1:dim(xdat_inat_df_yr_filtered)[3]) {
      
      # Extract the column for this year
      covar_values <- xdat_inat_df_yr_filtered[, col_idx, yr]
      
      # Interpolate NAs within each cell50 group
      cells_wNA <- unique(xdat_inat_df_yr_filtered[, which(colnames(xdat_inat_df_yr_filtered) == "cell50"), yr][is.na(covar_values)])
      
      for (g in cells_wNA) {
        inds <- which(xdat_inat_df_yr_filtered[, which(colnames(xdat_inat_df_yr_filtered) == "cell50"), yr] == g)
        covar_values[inds] <- ifelse(
          is.na(covar_values[inds]),
          mean(covar_values[inds], na.rm = TRUE),
          covar_values[inds])
      }
      
      # Set any final NAs to the mean
      covar_values[which(is.na(covar_values))] <- mean(covar_values, na.rm = TRUE)
      
      # Apply scaling
      scaling <- real_data$scaling_factors[which(real_data$scaling_factors$covar == covar_name),]
      covar_values2 <- (covar_values - scaling$mean) / scaling$sd
      
      # Put back into array
      xdat_inat_df_yr_scaled[, col_idx, yr] <- covar_values2
      
    }
  }
  
  # ===========================================================================
  # From here Claude has fixed the broken indexing
  # ===========================================================================
  
  # (2) NA GUARD -------------------------------------------------------------
  # Fill any residual NAs in covariate columns (e.g. the unscaled soil layers,
  # which are not handled by the scaling loop above). This prevents the
  # downstream model.matrix() calls in prep from dropping rows and thereby
  # breaking the correspondence with inat_cell50_start/end.
  all_cols  <- dimnames(xdat_inat_df_yr_scaled)[[2]]
  fill_cols <- setdiff(all_cols, "cell50")  # never overwrite the cell ID column
  for (yr in 1:dim(xdat_inat_df_yr_scaled)[3]) {
    for (cc in fill_cols) {
      v <- xdat_inat_df_yr_scaled[, cc, yr]
      if (any(is.na(v))) {
        v[is.na(v)] <- mean(v, na.rm = TRUE)
        xdat_inat_df_yr_scaled[, cc, yr] <- v
      }
    }
  }
  
  # (1) ROW ORDER FIX --------------------------------------------------------
  # Reorder the array rows by cell50 so that inat_cell50_start[g]:inat_cell50_end[g]
  # is a contiguous block containing exactly the fine cells of model cell g.
  # cell50 is identical across year slices, so slice 1 defines the ordering.
  ord <- order(xdat_inat_df_yr_scaled[, "cell50", 1])
  xdat_inat_df_yr_scaled <- xdat_inat_df_yr_scaled[ord, , , drop = FALSE]
  
  # Compute the index ranges ON THE REORDERED ARRAY (now contiguous by cell50)
  cell50_sorted   <- xdat_inat_df_yr_scaled[, "cell50", 1]
  cell50_forModel <- as.integer(as.factor(cell50_sorted))  # 1..G, contiguous
  n_groups        <- max(cell50_forModel)
  
  inat_cell50_start <- inat_cell50_end <- integer(n_groups)
  for (i in seq_len(n_groups)) {
    w <- which(cell50_forModel == i)
    inat_cell50_start[i] <- min(w)
    inat_cell50_end[i]   <- max(w)
  }
  
  # Sanity checks: ranges are contiguous and each maps to a single cell50
  stopifnot(all(inat_cell50_end >= inat_cell50_start))
  stopifnot(all(vapply(seq_len(n_groups),
                       function(i) length(unique(cell50_sorted[inat_cell50_start[i]:inat_cell50_end[i]])) == 1L,
                       logical(1))))
  
  cell50_model <- sort(unique(cell50_sorted))  # cell50 for group i is cell50_model[i]
  
  return(list(
    xdat_inat_df      = xdat_inat_df_yr_scaled,
    inat_cell50_start = inat_cell50_start,
    inat_cell50_end   = inat_cell50_end,
    cell50_model      = cell50_model))
}

# Prep the iNat data into a grid
prep_inat_data_grid <- function(taxon_key, species, redo = FALSE) {
  
  # PATCHED: write iNat grid to OUR space; reuse Arielle's precomputed file
  # (read-only) if present so we don't regenerate the slow all-species summary.
  if (!exists("INAT_GRID_DIR")) INAT_GRID_DIR <- paste0(PROJ_DIR, "/output/inat_grids")
  dir.create(INAT_GRID_DIR, showWarnings = FALSE, recursive = TRUE)
  outfile       <- paste0(INAT_GRID_DIR, "/", species, "_inat_grid.csv")
  arielle_file  <- paste0(PROJ_DIR, "/data/", species, "_inat_grid.csv")
  read_target   <- if (file.exists(outfile)) outfile else if (file.exists(arielle_file)) arielle_file else outfile

  sciname <- taxon_key$sci_name[taxon_key$common_name_clean == species]

  if (file.exists(read_target) && file.exists(GRID50_PATH) && !redo) {
    cat("  reusing precomputed iNat grid:", read_target, "\n")
    res <- read_csv(read_target)
    return(res)
  }
  
  # Using Ben's orig raster brick here to do the grid downscaling
  raster_brick <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/context_dependence_everything/data/covar_raster_brick.tif")
  
  grid100 <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/grid100.tif")
  names(grid100) <- "cell100"
  
  grid50 <- terra::aggregate(raster_brick[["MAP"]] * 
                               raster_brick[["agriculture_pct"]] * 
                               raster_brick[["soil_clay"]], 
                             fact = 10, na.rm = T)
  
  names(grid50) <- "cell50"
  terra::values(grid50)[!is.na(terra::values(grid50))] <- 1:sum(!is.na(terra::values(grid50)))
  
  if (!file.exists(GRID50_PATH)) {
    terra::writeRaster(grid50, paste0(INAT_GRID_DIR, "/grid50.tif"))  # never overwrite Arielle's
  }
  
  coords50 <- as.data.frame(grid50, xy = TRUE) %>% rename(x50 = x, y50 = y)
  coords100 <- as.data.frame(grid100, xy = TRUE) %>% rename(x100 = x, y100 = y)
  
  # iNat download
  inat_dat_all <- read_csv("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/inat_combo_nam_mams.csv") %>% 
    filter(public_positional_accuracy < 1000) 
  
  #Filter to just the years 2008-2025 to match the occupancy dataset
  inat_dat_all <- inat_dat_all %>% mutate(Year=year(observed_on))%>%
    filter(Year>=2008&Year<=2025)
  
  inat_pts <- inat_dat_all %>% 
    select(longitude, latitude) %>% 
    vect(geom = c("longitude", "latitude"), crs = "+proj=longlat") %>% 
    project(crs(grid100))
  
  inat_cells <- bind_cols(
    extract(grid50, inat_pts)[, 2],
    extract(grid100, inat_pts)[, 2])
  colnames(inat_cells) <- c("cell50", "cell100")
  
  inat_cells$sciname <- inat_dat_all$taxon_species_name
  
  # Add the year
  inat_cells$Year <- inat_dat_all$Year
  
  # Aggregate all obs, obs of each target species
  inat_cell_summary <- inat_cells %>% 
    group_by(Year) %>%
    count(cell50, cell100) %>% 
    rename(effort = n) %>% 
    filter(!is.na(cell50)) %>% 
    left_join(coords50) %>% 
    left_join(coords100)
  
  for (i in 1:nrow(taxon_key)) {
    this_spec_dat <- inat_cells %>% 
      filter(sciname == taxon_key$sci_name[i]) %>% 
      group_by(Year) %>%
      count(cell50, cell100)
    colnames(this_spec_dat)[4] <- taxon_key$common_name_clean[i]
    
    # NA out effort for cells that are outside the range
    range <- vect(taxon_key$rangefile[i]) %>% simplifyGeom(0.01) %>% project(crs(grid100))
    range_to50 <- rasterize(range, grid50, touches = TRUE) %>% as.data.frame(xy = TRUE)
    
    inat_cell_summary <- left_join(inat_cell_summary, this_spec_dat, by = c("Year","cell50", "cell100"))
    inat_cell_summary[[taxon_key$common_name_clean[i]]][is.na(inat_cell_summary[[taxon_key$common_name_clean[i]]])] <- 0
    
    inat_cell_summary[[taxon_key$common_name_clean[i]]][!paste0(inat_cell_summary$x50, inat_cell_summary$y50) %in%
                                                          paste0(range_to50$x, range_to50$y)] <- NA
  }
  
  write_csv(inat_cell_summary, outfile)
  
  inat_cell_summary
}

