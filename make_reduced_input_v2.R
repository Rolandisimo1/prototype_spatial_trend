# =============================================================================
# make_reduced_input_v2.R
# -----------------------------------------------------------------------------
# Forked from make_reduced_input.R (10-covariate KEEP list, 2026-08 decision).
# Produces a further-reduced-covariate copy of input_data_<species>.RDS /
# input_data_<species>_v2.RDS, dropping ONE MORE covariate: soil_sand.
#
# REVISION (2026-08-14, corrects the first version of this file, which
# dropped soil_silt): the first cut was wrong on TWO counts, found by a
# direct VIF/correlation check across all three species' actual bundles
# rather than assumed from the "compositional data" framing alone.
#
# 1. NOT a universal compositional identity. R^2 of sand ~ clay + silt is
#    0.9995 for moose (nsite=1654, narrow northern camera-site footprint --
#    a genuine near-collinearity specific to that footprint) but only 0.606
#    (bobcat, nsite=20531) and 0.617 (WTD, nsite=21559). Max VIF with all
#    10 original covariates: moose 2267 (severe), bobcat 3.0, WTD 3.1
#    (both completely benign). The "sums to ~100% therefore exactly
#    unidentified for every species" framing in the memo's Issue 3 is NOT
#    supported by bobcat/WTD data -- this is a moose-specific collinearity,
#    not a fleet-wide structural identity, even though a fleet-wide 9-cov
#    design is still the right call for report comparability.
#
# 2. Dropping soil_silt was the WRONG fraction. It leaves soil_clay +
#    soil_sand retained together -- exactly the -0.9966-correlated pair
#    (moose) that was occ_beta[8]/occ_beta[10], the two parameters that
#    actually plateaued both moose v2b fits. Max VIF over the RETAINED
#    pair, moose, by drop choice:
#      drop none        : max VIF 2267  (both moose v2b fits plateaued here)
#      drop soil_silt    : max VIF  152  (first version of this file --
#                                         removes the exact singularity but
#                                         leaves a severe ridge; would very
#                                         likely plateau again)
#      drop soil_sand    : max VIF  4.5  cor(clay,silt) = -0.2974
#      drop soil_clay    : max VIF  4.5  cor(silt,sand) = +0.2202
#
# FIX (corrected): drop soil_sand, keep soil_clay + soil_silt. Best
# conditioning across all three species (moose max VIF 4.5; bobcat/WTD both
# 1.9), and it resolves the ridge that blocked moose rather than just
# shrinking it. Applying one uniform 9-covariate design across all five
# builds is still the right call for report comparability -- it costs
# bobcat/WTD nothing, since their soil covariates were never near-collinear
# to begin with.
#
# HOW THE ORIGINAL PROBLEM WAS FOUND (2026-08-14): moose_v2b_national_scalar
# and moose_v2b_ecoregion both plateaued with occ_beta[8] (soil_clay) and
# occ_beta[10] (soil_sand) R-hat matching to four decimal places -- the
# classic signature of two parameters drifting along a shared ridge.
#
# FLEET-WIDE ADDENDUM (2026-08-14, later the same day -- supersedes the
# "moose-specific" wording in point 1 above): the same check run across the
# cached <species>_4SPO.RDS files for five additional narrow-range species
# shows the ridge is NOT confined to moose. Max VIF over the KEEP design
# (all 10 / drop silt / drop sand / drop clay):
#     moose      2267.4 / 152.2 / 4.5 / 4.5
#     grey_wolf  1732.7 /  19.9 / 2.8 / 2.8
#     fisher      242.9 /  30.9 / 3.2 / 3.2
#     elk          23.3 /  13.7 / 2.7 / 2.7
#     kit_fox       4.5 /   4.3 / 1.9 / 1.7
#     bobcat        3.0 /   1.9 / 1.9 / 2.1
#     WTD           3.1 /   1.9 / 1.9 / 2.1
#     pronghorn     2.7 /   2.7 / 2.7 / 2.7
# grey_wolf and fisher are in the same failure class as moose and elk is over
# the conventional VIF>10 line, so all three would have plateaued the same way
# on the uncorrected design. Dropping soil_silt is never sufficient wherever a
# real problem exists (moose 152, fisher 31, grey_wolf 20, elk 14). Dropping
# soil_sand is at or tied-with best for all eight species, fleet-wide worst
# case 4.5 -- so ONE uniform 9-covariate design is safe and the fix does not
# need to be species-specific. Any future fleet v2b build must use it.
#
# Underlying cause, measured rather than assumed: the raw fractions are
# MINERAL soil and their sum maxes at exactly 1.0 but is bimodal, tailing to
# ~0.002 where ground is rock, wetland/organic or urban fill. Share of sites
# within 1 percentage point of 1.0: moose 98.3%, grey_wolf 98.6%, fisher
# 92.5%, elk 79.6%, pronghorn 63.3%, kit_fox 55.0%, WTD 50.6%, bobcat 45.8%.
# So the composition closes in narrow northern/forested ranges and not in
# geographically diverse ones -- hence "implicit reference level" is strictly
# accurate only for the former; elsewhere the two retained coefficients are
# just two texture axes, not deviations from a fixed whole.
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
DROP2 <- c("soil_sand")

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
cat("\nDropping (moose-specific soil_clay/soil_sand near-collinearity fix):\n"); print(intersect(DROP2, old_names))
missing_drop <- setdiff(DROP2, old_names)
if (length(missing_drop)) {
  stop("soil_sand not found in incoming covariate list -- was this bundle already reduced? Names: ",
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
