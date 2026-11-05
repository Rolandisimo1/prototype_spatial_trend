#!/usr/bin/env Rscript
# =============================================================================
# run_fleet_prep.R -- GENERALISED fork of ~/isdm/v2_prep/run_v2b_prep.R
# (2026-09-16) so any fleet species can be built on any mask grid.
# With every env var at its default this script reproduces run_v2b_prep.R's
# bundle byte-for-byte (validated by rebuilding moose_v2b). Parameters:
#   SPECIES          taxon token (common_name_clean); must exist in target_species.csv
#   MASK_TAG         "v2b" (presence mask, default) or e.g. "v2u" (union mask);
#                    out dir HPC/<species>_<MASK_TAG>, cell map suffix _<MASK_TAG>
#   INAT_GRID        iNat grid CSV path; {species} template allowed.
#                    default output/inat_grids_v2/<species>_inat_grid_v2.csv
#   PINS_CSV         optional pins table (species, cell50_<kind>, records_<kind>);
#   PIN_KIND         "presence" or "union" -- which pinned columns apply
#   OUT_DIR          override output dir (validation rebuilds)
#   CELLMAP_DIR      override cell map dir (validation rebuilds)
# PIN LOGIC: pinned species that deviates -> stop. Unpinned species -> the
# bundle is built so it can be reviewed, but BUNDLE_REVIEW_REQUIRED is written
# into its directory and every fit runner refuses to start while it exists.
# If the grid itself carries <species>_REVIEW_REQUIRED next to it, the bundle
# inherits the flag regardless of pins.
# =============================================================================
#
# ---- original header (run_v2b_prep.R) ----
# v2b data prep driver -- FIX 1 + FIX 2b (data-driven PRESENCE MASK).
#
# Forked from run_v1fix_prep.R (2026-08-11). EXACTLY ONE thing differs from
# that driver: where the iNat grid comes from. v1fix used the original
# IUCN range-polygon mask via prep_inat_data_grid(); this uses the
# presence-mask grid (a cell is "in range" if it has >= PRESENCE_MIN_RECORDS
# total records, currently 1) produced by run_data_prep_v2.R. For moose that
# moves ncell50 from 382 to 422 and recovers +4,290 iNat records, including
# all of Colorado, which the IUCN polygon excluded entirely.
#
# Everything downstream of the iNat pull -- the 10-name KEEP covariate list,
# make_inat_xdat(), the Fix 1 matrix builders and their alignment guardrail,
# constants_list, inits_list -- is byte-identical to run_v1fix_prep.R, so the
# only intended difference between a v1fix fit and a v2b fit is the mask.
#
# Writes to NEW directories, HPC/<species>_v2b/, and never touches the
# v1fix bundles or the (unmasked, ncell50=3322, unfittable) _v2 bundles from
# 2026-07-27. NOTE the token is "_v2b", deliberately distinct from "_v2".
#
# Species is passed via the SPECIES env var (one of: bobcat, moose,
# white-tailed_deer).

if (!requireNamespace("nimbleEcology", quietly=TRUE)) message("NOTE: nimbleEcology skipped (unused in prep)") else suppressMessages(library(nimbleEcology))
library(tidyverse)
library(sf)
if (!requireNamespace("MCMCvis", quietly=TRUE)) message("NOTE: MCMCvis not installed - skipping (unused in prep)") else suppressMessages(library(MCMCvis))
library(parallel)
if (!requireNamespace("coda", quietly=TRUE)) message("NOTE: coda skipped (unused in prep)") else suppressMessages(library(coda))

setwd("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final")
project_dir <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"

source('/home/rwkays/isdm/v2_prep/prep_data_for_spoccupancy.R')
source('/home/rwkays/isdm/v2_prep/pipeline_helper.R')
source('/home/rwkays/isdm/v2_prep/integration_helper_v2.R')
# CAMERA_CACHE_DIR (optional): read the camera-side prep from this directory
# instead of data/forSPO/. Used for the union-mask camera filter built by
# prep_camera_sites.R. Unset = original behaviour (regression-validated).
CAMERA_CACHE_DIR <- Sys.getenv("CAMERA_CACHE_DIR", "")
if (nzchar(CAMERA_CACHE_DIR)) {
  prep_spoccupancy_data <- function(species, redo = FALSE) {
    f <- file.path(CAMERA_CACHE_DIR, paste0(species, "_4SPO.RDS"))
    if (!file.exists(f)) stop("camera cache not found: ", f)
    spocc_dat <- readRDS(f)
    spocc_dat$grid100 <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/grid100.tif")
    cat("camera cache:", f, " range_mode:", spocc_dat$camera_prep$range_mode, " nsite:", nrow(spocc_dat$y), "\n")
    spocc_dat
  }
}  # FIX 1 lives here; prep_inat_data_grid() (masked) is untouched in this same file

species <- Sys.getenv("SPECIES", unset = "")
if (species == "") stop("SPECIES env var must be set")
if (!species %in% taxon_key$common_name_clean) stop("SPECIES '", species, "' not in target_species.csv")
MASK_TAG <- Sys.getenv("MASK_TAG", "v2b")
INAT_GRID_T <- Sys.getenv("INAT_GRID", file.path(project_dir, "output", "inat_grids_v2", "{species}_inat_grid_v2.csv"))
PINS_CSV <- Sys.getenv("PINS_CSV", ""); PIN_KIND <- Sys.getenv("PIN_KIND", "presence")
review_reasons <- character(0)
cat("=== fleet prep for species:", species, " MASK_TAG:", MASK_TAG, " grid:", INAT_GRID_T, "===\n")
cat("start:", format(Sys.time()), "\n")

data_prep_v2b <- function(species) {

  annual_covars <- c("CWD","Human_pop","NDVI_mean","NDVI_sd","Ag","Deciduous",
                     "Evergreen","Impervious","Mixed","PDSI","PPT","Temp")
  static_covars <- c("terrain_ruggedness","soil_clay","soil_silt","soil_sand","elevation")
  temp_covars <- c("MWMT", "MCMT")

  cases <- case <- "SVC"
  years_in <- seq(2008, 2025, by = 1)

  real_data <- prep_spoccupancy_data(species)

  row <- as.data.frame(t(rep(1, length(c(annual_covars, static_covars, temp_covars)))))
  colnames(row) <- c(annual_covars, static_covars, temp_covars)
  row$temp_type <- 6

  # EXPLICIT covariate list, per make_reduced_input.R's documented final
  # decision (2026-07-28 review) -- NOT the is.nan() auto-selection filter.
  # That filter silently changes the design matrix (10 -> 17 covariates seen
  # today, 7 extra ones that pass only because the underlying covariate
  # raster source has drifted since the original bundles were built) and
  # was the proximate trigger for moose's NA initial logProb. MWMT/MCMT are
  # deliberately absent here -- they enter via their own SVC/CAR field, not
  # occ_covars.
  KEEP_COVARS <- c("Human_pop", "NDVI_mean", "Ag", "Deciduous", "Evergreen", "Mixed",
                     "terrain_ruggedness", "soil_clay", "soil_silt", "soil_sand")
  DROP_COVARS <- c("CWD", "NDVI_sd", "Impervious", "PDSI", "PPT", "elevation", "Temp")  # documented collinear set, for reference only

  candidate_covars <- c(annual_covars, static_covars, temp_covars)
  na_covars <- colSums(is.nan(as.matrix(real_data$occ.covs[, candidate_covars])))
  names(na_covars) <- candidate_covars

  missing_keep <- setdiff(KEEP_COVARS, candidate_covars)
  if (length(missing_keep) > 0) {
    stop("KEEP_COVARS names not found among candidate covariates for species '", species,
         "': ", paste(missing_keep, collapse = ", "))
  }
  nan_in_keep <- KEEP_COVARS[na_covars[KEEP_COVARS] > 0]
  if (length(nan_in_keep) > 0) {
    stop("Explicit KEEP_COVARS have NaN values in this fresh pull for species '", species,
         "' -- stopping rather than silently substituting: ", paste(nan_in_keep, collapse = ", "))
  }

  occ_covars <- KEEP_COVARS
  cat("\n--- occ_covars: EXPLICIT 10-name KEEP list (make_reduced_input.R decision) ---\n")
  cat(" ", paste(occ_covars, collapse = ", "), "\n")
  auto_would_have_picked <- candidate_covars[na_covars == 0]
  extra_available <- setdiff(auto_would_have_picked, occ_covars)
  if (length(extra_available) > 0) {
    cat("(", length(extra_available), "additional covariate(s) are NaN-free in this fresh pull",
        "but excluded per the documented collinearity decision:", paste(extra_available, collapse = ", "), ")\n")
  }
  stopifnot(
    "occ_covars does not exactly match the documented 10-covariate KEEP list" =
      identical(sort(occ_covars), sort(KEEP_COVARS))
  )
  occ_formula <- make_formula(row, occ_covars, occ_only = TRUE)

  occ_cov_mtx <- model.matrix(occ_formula, real_data$occ.covs)[, -1]
  occ_cov_mtx <- occ_cov_mtx[, !colnames(occ_cov_mtx) %in% c("MCMT", "MWMT", "MCMT_sq", "MWMT_sq")]

  occ_cov_terms <- colnames(occ_cov_mtx)
  occ_cov_names_df <- data.frame(name = occ_cov_terms) %>%
    mutate(param = paste0("occ_beta[", row_number(), "]"))
  occ_cov_names_df$cat <- "Main"
  occ_cov_names_df$cat[grepl(":", occ_cov_names_df$name)] <- "interaction"
  occ_cov_names_df$temp_cat <- "None"
  occ_cov_names_df$temp_cat[grepl("MWMT", occ_cov_names_df$name)] <- "MWMT"
  occ_cov_names_df$temp_cat[grepl("MCMT", occ_cov_names_df$name)] <- "MCMT"
  occ_cov_names_df$plot_cat <- ifelse(occ_cov_names_df$cat == "Main", "Main",
                                      paste0(occ_cov_names_df$temp_cat, " interaction"))

  # ---- iNat: FIX 2b path -- data-driven PRESENCE MASK, not the IUCN polygon ----
  # Reads the grid written by run_data_prep_v2.R (job 526880, 2026-08-10
  # 23:25) rather than re-invoking prep_inat_data_grid_v2(). That artifact is
  # the one that was traced end-to-end and signed off on 2026-08-11
  # (PRESENCE_MIN_RECORDS = 1); reading it directly pins this bundle to
  # exactly those numbers instead of re-deriving them, and avoids repeating
  # the ~2h15m prep for a result that should be identical.
  inat_grid_path <- gsub("{species}", species, INAT_GRID_T, fixed = TRUE)
  grid_flag <- file.path(dirname(inat_grid_path), paste0(species, "_REVIEW_REQUIRED"))
  if (file.exists(grid_flag)) review_reasons <<- c(review_reasons, paste("grid carries", grid_flag))
  cat("\n--- FIX2b (presence-mask) iNat grid for", species, "---\n")
  cat("reading:", inat_grid_path, "\n")
  stopifnot("presence-mask iNat grid not found" = file.exists(inat_grid_path))
  inat_df <- readr::read_csv(inat_grid_path, show_col_types = FALSE)

  # Same NA-filter the v1fix driver applies. The presence mask marks
  # out-of-range cells by NA-ing the species column WITHOUT dropping the rows,
  # so skipping this would silently build matrices over the full 3,322-cell
  # national grid and make the mask a no-op.
  inat_df <- inat_df[!is.na(inat_df[[species]]), ]
  cat(species, "in-range rows after presence mask:", nrow(inat_df),
      " distinct cell50:", length(unique(inat_df$cell50)),
      " total records:", sum(inat_df[[species]], na.rm = TRUE), "\n")

  # Pin to the traced numbers. If the grid is ever regenerated at a different
  # threshold, this stops the build rather than silently fitting different
  # data -- which is exactly the failure mode behind the >=3 / >=1 artifact
  # mixup found on 2026-08-11.
  EXPECTED <- if (MASK_TAG == "v2b") list(moose = list(ncell50 = 422L, records = 11582)) else list()
  if (nzchar(PINS_CSV) && file.exists(PINS_CSV)) {
    pt <- read.csv(PINS_CSV, stringsAsFactors = FALSE)
    pr <- pt[pt$species == species, ]
    if (nrow(pr) == 1) EXPECTED[[species]] <- list(ncell50 = as.integer(pr[[paste0("cell50_", PIN_KIND)]]),
                                                   records = as.numeric(pr[[paste0("records_", PIN_KIND)]]))
  }
  if (!is.null(EXPECTED[[species]])) {
    e <- EXPECTED[[species]]
    got_cells <- length(unique(inat_df$cell50))
    got_recs  <- sum(inat_df[[species]], na.rm = TRUE)
    if (got_cells != e$ncell50 || got_recs != e$records) {
      stop("PRESENCE-MASK GUARDRAIL FAILED for ", species,
           ": expected ncell50 = ", e$ncell50, " and records = ", e$records,
           ", got ncell50 = ", got_cells, " and records = ", got_recs,
           ". The iNat grid at ", inat_grid_path,
           " is not the traced 2026-08-10 presence-mask artifact -- refusing",
           " to build a bundle from an unverified grid.")
    }
    cat("PRESENCE-MASK GUARDRAIL PASSED: ncell50 =", e$ncell50,
        "  records =", e$records, "\n")
  } else {
    cat("NOTE: no pinned expectation for", species, "-- bundle will carry BUNDLE_REVIEW_REQUIRED\n")
    review_reasons <<- c(review_reasons, paste0("no pinned ", PIN_KIND, " expectation (ncell50 = ",
      length(unique(inat_df$cell50)), ", records = ", sum(inat_df[[species]], na.rm = TRUE), ")"))
  }

  inat_xdat_list <- make_inat_xdat(inat_df, real_data, species)

  for (i in 1:dim(inat_xdat_list$xdat_inat_df)[3]) {
    inat_xdat_list$xdat_inat_df[,'MCMT_sq',i] <- inat_xdat_list$xdat_inat_df[,'MCMT',i]^2
    inat_xdat_list$xdat_inat_df[,'MWMT_sq',i] <- inat_xdat_list$xdat_inat_df[,'MWMT',i]^2
  }

  mm <- model.matrix(occ_formula, as.data.frame(inat_xdat_list$xdat_inat_df[,,1]))[, -1]
  xdat_array <- array(NA, dim = c(dim(mm)[1], dim(mm)[2], length(years_in)),
                      dimnames = list(rownames(mm), colnames(mm), years_in))
  for (i in 1:dim(inat_xdat_list$xdat_inat_df)[3]) {
    xdat_array[,,i] <- model.matrix(occ_formula, as.data.frame(inat_xdat_list$xdat_inat_df[,,i]))[, -1]
  }

  inat_y <- make_inat_cell_year_matrix(inat_df, species)
  inat_effort <- make_inat_effort_matrix(inat_df)

  # FIX 1 GUARDRAIL -- must hold or we stop here rather than bake a silent
  # misalignment into a converged bundle again.
  stopifnot(
    "inat_y and inat_effort column (year) order diverged" =
      identical(colnames(inat_y), colnames(inat_effort)),
    "inat_y and inat_effort row (cell50) order diverged" =
      identical(rownames(inat_y), rownames(inat_effort))
  )
  cat("FIX1 GUARDRAIL PASSED: inat_y/inat_effort row+column alignment confirmed for", species, "\n")

  inat_year_index <- vector("list", dim(inat_y)[2])
  for (t in 1:dim(inat_y)[2]) inat_year_index[[t]] <- which(inat_effort[, t] > 0)

  max_cells <- max(sapply(inat_year_index, length))
  inat_cells_by_year <- matrix(NA, max_cells, dim(inat_y)[2])
  n_cells_year <- integer(dim(inat_y)[2])
  for (t in 1:dim(inat_y)[2]) {
    idx <- inat_year_index[[t]]
    n_cells_year[t] <- length(idx)
    inat_cells_by_year[1:n_cells_year[t], t] <- idx
  }
  cat("n_cells_year range:", paste(range(n_cells_year), collapse=" to "), "\n")

  cells <- inat_df %>% group_by(cell50) %>% select(cell50, cell100) %>% slice(1) %>% arrange(cell50)
  cellmap_dir <- Sys.getenv("CELLMAP_DIR", file.path(project_dir, "cell_maps"))
  dir.create(cellmap_dir, showWarnings = FALSE)
  saveRDS(cells, file = paste0(cellmap_dir, "/cell_map_50_100_", species, "_", MASK_TAG, ".RDS"))

  prior_type <- "Normal"

  constants_list <- list(
    has_SVC = case %in% c("SVC", "Both"),
    prior_type = prior_type,
    nyear = length(years_in),
    year_vals = scale(years_in)[,1],
    nsite = nrow(real_data$y),
    J = as.numeric(rowSums(!is.na(real_data$y))),
    numOccCovars = ncol(occ_cov_mtx),
    occ_covars = occ_cov_mtx,
    year_occ = scale(real_data$occ.covs$year)[,1],
    yday = real_data$det.covs$yday_scaled,
    canopy_height = real_data$det.covs$canopy_height_scaled,
    log_roaddist = real_data$det.covs$log_roaddist_scaled,
    MWMT = real_data$occ.covs$MWMT,
    MCMT = real_data$occ.covs$MCMT,
    interaction_group = as.numeric(as.factor(occ_cov_names_df$plot_cat)),
    adj = real_data$continent_adj_info$adj,
    num = real_data$continent_adj_info$num,
    nadj = length(real_data$continent_adj_info$adj),
    nnum = length(real_data$continent_adj_info$num),
    ncell100 = length(real_data$continent_adj_info$num),
    cell = real_data$occ.covs$spatcell,
    ncell50 = nrow(cells),
    inat_cell100 = cells$cell100,
    xdat_inat = xdat_array[, colnames(occ_cov_mtx), ],
    inat_cell50_start = inat_xdat_list$inat_cell50_start,
    inat_cell50_end = inat_xdat_list$inat_cell50_end,
    MCMT_inat = xdat_array[, "MCMT", ],
    MWMT_inat = xdat_array[, "MWMT", ],
    n_cells_year = n_cells_year,
    inat_cells_by_year = inat_cells_by_year,
    hasSVC = (case %in% c("SVC", "Both"))
  )

  inits_list <- list(
    occ_beta = rep(0, ncol(occ_cov_mtx)),
    p_beta = rep(0, 4),
    det_intercept = 0.5,
    link_occ_intercept = rep(-1, length(real_data$continent_adj_info$num)),
    intercept_tau = 0.3,
    overdisp_inat = 0.1,
    theta0 = -5,
    theta1 = 1,
    year_beta = 0,
    year_var = 0
  )
  if (case %in% c("SVC", "Both")) {
    inits_list$MWMT_effect <- rep(0, length(real_data$continent_adj_info$num))
    inits_list$MCMT_effect <- rep(0, length(real_data$continent_adj_info$num))
    inits_list$MWMT_tau <- 0.3
    inits_list$MCMT_tau <- 0.3
  }

  input_data <- list(real_data = real_data, inat_y = inat_y, inat_effort = inat_effort,
                     constants_list = constants_list, inits_list = inits_list,
                     occ_cov_mtx = occ_cov_mtx, annual_covars = annual_covars,
                     static_covars = static_covars, temp_covars = temp_covars,
                     species = species, prior_type = prior_type,
                     v2_notes = if (MASK_TAG != "v2b") paste0("Fleet build via run_fleet_prep.R, MASK_TAG=", MASK_TAG, ", iNat grid ", inat_grid_path, ", built ", format(Sys.time()), ". Downstream of the iNat grid identical to run_v2b_prep.R (10-name KEEP list; reduce to 9 cov with make_reduced_input_v2.R). Camera sites come from data/forSPO/<species>_4SPO.RDS, which prep_spoccupancy_data() restricts to the IUCN extant-range polygon.") else "Fix1 (names_sort pivot fix + alignment guardrail) + Fix2b (data-driven PRESENCE MASK, PRESENCE_MIN_RECORDS = 1) -- the IUCN range-polygon mask is REPLACED, not applied. iNat grid read from output/inat_grids_v2/<species>_inat_grid_v2.csv as written by run_data_prep_v2.R job 526880 on 2026-08-10 23:25, traced and signed off 2026-08-11 (moose: ncell50 422, 11,582 records, +4,290 vs IUCN, all of Colorado recovered). NOT the unmasked Fix2 v1 bundle (ncell50 3322, 2026-07-27, never fittable). Everything downstream of the iNat pull is byte-identical to run_v1fix_prep.R. Occupancy covariates use the EXPLICIT 10-name KEEP list from make_reduced_input.R's documented collinearity decision (2026-07-28), not the is.nan() auto-selection filter that produced a drifted 17-covariate design matrix and triggered moose's NA initial logProb.")

  out_dir <- Sys.getenv("OUT_DIR", paste0(project_dir, "/HPC/", species, "_", MASK_TAG))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  out_file <- paste0(out_dir, "/input_data_", species, "_", MASK_TAG, ".RDS")
  if (file.exists(out_file)) stop("refusing to overwrite ", out_file)
  saveRDS(input_data, file = out_file)
  cat("wrote", out_file, "\n")
  flag <- file.path(out_dir, "BUNDLE_REVIEW_REQUIRED")
  if (length(review_reasons)) { writeLines(c(format(Sys.time()), review_reasons), flag); cat("WROTE", flag, ":", review_reasons, sep = "\n  ") }
  cat("BUNDLE SUMMARY:", species, MASK_TAG, " nsite", constants_list$nsite,
      " detecting_sites", sum(rowSums(real_data$y, na.rm = TRUE) > 0),
      " ncell50", constants_list$ncell50, " inat_ncell100", length(unique(cells$cell100)),
      " cam_ncell100", length(unique(constants_list$cell)), " CAR_ncell100", constants_list$ncell100,
      " iNat_records", sum(inat_y, na.rm = TRUE), " nyear", constants_list$nyear, "\n")
  invisible(input_data)
}

data_prep_v2b(species)
cat("=== DONE v2b prep for", species, "at", format(Sys.time()), "===\n")
