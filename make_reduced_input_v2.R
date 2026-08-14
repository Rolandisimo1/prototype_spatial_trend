# =============================================================================
# make_reduced_input_v2.R
# -----------------------------------------------------------------------------
# Forked from make_reduced_input.R (10-covariate KEEP list, 2026-08 decision).
# Produces a further-reduced-covariate copy of input_data_<species>.RDS /
# input_data_<species>_v2.RDS, dropping ONE MORE covariate: soil_silt.
#
# WHY: soil_clay, soil_silt, and soil_sand are compositional data -- by
# construction they sum to ~100% of a soil sample, so any two of the three
# fully determine the third. Including all three as free linear coefficients
# in occ_beta builds an exactly-unidentified redundancy into the design
# matrix: the model can trade weight between them along the "sum" direction
# with zero change in fit. This is a latent design flaw in the ORIGINAL
# 10-covariate KEEP list, not something introduced by the presence-mask fix.
#
# HOW IT WAS FOUND (2026-08-14): moose_v2b_national_scalar and
# moose_v2b_ecoregion both plateaued with occ_beta[8] (soil_clay) and
# occ_beta[10] (soil_sand) R-hat matching to four decimal places -- the
# classic signature of two parameters drifting along an unidentified ridge.
# Direct correlation check on the bundles:
#   moose   nsite=1654   cor(soil_clay, soil_sand) = -0.9966
#   bobcat  nsite=20531  cor(soil_clay, soil_sand) = -0.5662
#   WTD     nsite=21559  cor(soil_clay, soil_sand) = -0.5664
# The ridge is present in ALL THREE species' v1fix bundles (identical
# correlation, since it depends only on camera-site soil data, and the
# presence mask only changes iNat cell membership, never camera sites). It
# was not BINDING for bobcat/WTD because their camera sites span enough
# geography for clay/sand to vary somewhat independently; moose's narrow
# northern camera-site band is what pushed the correlation to -0.997 and
# made the ridge bite.
#
# FIX: drop soil_silt (keep soil_clay + soil_sand as the two explicit
# texture-axis effects; silt becomes the implicit reference level). This is
# a real correction to the design matrix, not a per-species workaround --
# apply it to any v2/v2b build for any species, not just moose.
#
# Usage:  Rscript make_reduced_input_v2.R <species> <input_rds_path>
#   e.g.  Rscript make_reduced_input_v2.R moose \
#           $PROJ/output/bundles/input_data_moose_v2b_national_scalar.RDS
# Writes: <input_rds_path> with _v2cov appended before .RDS
#
# Base R only -- no nimble needed; runs under `module load R`.
# =============================================================================

args    <- commandArgs(trailingOnly = TRUE)
species  <- if (length(args) >= 1) args[1] else stop("species required")
fin      <- if (length(args) >= 2) args[2] else stop("input RDS path required")

# Additional drop, layered on top of the original 10-covariate KEEP list
# (Human_pop, NDVI_mean, Ag, Deciduous, Evergreen, Mixed, terrain_ruggedness,
#  soil_clay, soil_silt, soil_sand -- see make_reduced_input.R for the first
#  reduction, DROP = CWD/NDVI_sd/Impervious/PDSI/PPT/elevation/Temp).
DROP2 <- c("soil_silt")

fout <- sub("\\.RDS$", "_v2cov.RDS", fin)

cat("Reading :", fin, "\n")
stopifnot(file.exists(fin))
d  <- readRDS(fin)
cl <- d$constants_list

old_names <- colnames(cl$occ_covars)
stopifnot(!is.null(old_names))
keep_idx   <- which(!(old_names %in% DROP2))
keep_names <- old_names[keep_idx]

cat("\nIncoming occ covariates (", length(old_names), "):\n", sep=""); print(old_names)
cat("\nDropping (compositional redundancy fix):\n"); print(intersect(DROP2, old_names))
missing_drop <- setdiff(DROP2, old_names)
if (length(missing_drop)) {
  stop("soil_silt not found in incoming covariate list -- was this bundle already reduced? Names: ",
       paste(old_names, collapse=", "))
}
cat("\nRetained (", length(keep_names), "):\n", sep=""); print(keep_names)

old_ncov <- cl$numOccCovars
stopifnot(old_ncov == length(old_names))

cl$occ_covars <- cl$occ_covars[, keep_idx, drop = FALSE]
cl$xdat_inat  <- cl$xdat_inat[, keep_idx, , drop = FALSE]
cl$numOccCovars      <- length(keep_idx)
cl$interaction_group <- cl$interaction_group[keep_idx]

# SVC (MWMT/MCMT) fields are untouched -- same invariant as make_reduced_input.R.

if (!is.null(d$inits_list)) {
  for (nm in names(d$inits_list)) {
    v <- d$inits_list[[nm]]
    if (is.numeric(v) && is.null(dim(v)) && length(v) == old_ncov) {
      d$inits_list[[nm]] <- v[keep_idx]
      cat("Subset init '", nm, "' from ", old_ncov, " -> ", length(keep_idx), "\n", sep="")
    }
  }
}

d$constants_list <- cl
if (!is.null(d$occ_cov_mtx) && !is.null(colnames(d$occ_cov_mtx))) {
  cm <- colnames(d$occ_cov_mtx)
  d$occ_cov_mtx <- d$occ_cov_mtx[, which(!(cm %in% DROP2)), drop = FALSE]
}
if (!is.null(d$annual_covars)) d$annual_covars <- setdiff(d$annual_covars, DROP2)
if (!is.null(d$static_covars)) d$static_covars <- setdiff(d$static_covars, DROP2)

cat("\n=== POST-REDUCTION CHECK ===\n")
cat("numOccCovars       :", cl$numOccCovars, "(expected 9)\n")
cat("occ_covars dim     :", paste(dim(cl$occ_covars), collapse=" x "), "\n")
cat("xdat_inat dim      :", paste(dim(cl$xdat_inat), collapse=" x "), "\n")
cat("interaction_group  :", length(cl$interaction_group), "\n")
cat("has_SVC / hasSVC   :", cl$has_SVC, "/", cl$hasSVC, "\n")
stopifnot(ncol(cl$occ_covars) == cl$numOccCovars,
          dim(cl$xdat_inat)[2] == cl$numOccCovars,
          length(cl$interaction_group) == cl$numOccCovars,
          cl$numOccCovars == 9,
          isTRUE(cl$has_SVC), isTRUE(cl$hasSVC))

saveRDS(d, fout)
cat("\nWROTE:", fout, "\n")
cat("Size:", round(file.size(fout)/1e6, 1), "MB\n")
