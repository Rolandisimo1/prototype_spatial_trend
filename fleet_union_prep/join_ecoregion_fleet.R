#!/usr/bin/env Rscript
# Adapted from 00b2_join_ecoregion_real.R for the v2 (Fix1+Fix2) bundles.
# Read-only w.r.t. everything except the new HPC/<token>_v2b_{ecoregion,
# national_scalar}/ directories it creates. Never touches Arielle's
# originals or the existing converged bobcat/moose/WTD bundles.
#
# Usage: Rscript join_ecoregion_v2b.R <species_token>
# Example: Rscript join_ecoregion_v2b.R white-tailed_deer
# Reads:  HPC/<species_token>_v2b/input_data_<species_token>_v2b.RDS
# Writes: HPC/<species_token>_v2b_ecoregion/input_data_<species_token>_v2b_ecoregion.RDS
#         HPC/<species_token>_v2b_national_scalar/input_data_<species_token>_v2b_national_scalar.RDS

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: Rscript join_ecoregion_fleet.R <species_token> <mask_tag>")
species_token <- args[1]; TAG <- args[2]

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")

real_path <- paste0(PROJ, "/HPC/", species_token, "_", TAG, "/input_data_", species_token, "_", TAG, ".RDS")
stopifnot(file.exists(real_path))
real <- readRDS(real_path)
cl <- real$constants_list

prepped <- readRDS(paste0(PROTO_DIR, "/prepped_sim_inputs.RDS"))
stopifnot(!is.null(prepped$ecoregion_of_cell100), !is.null(prepped$nregion))

cat("=== identity check: v2 bundle's cell100 CAR graph vs prepped_sim_inputs' ===\n")
cat("v2 ncell100:", cl$ncell100, " prepped ncell100:", prepped$constants_list$ncell100, "\n")
stopifnot(cl$ncell100 == prepped$constants_list$ncell100)
stopifnot(identical(cl$adj, prepped$constants_list$adj))
stopifnot(identical(cl$num, prepped$constants_list$num))
stopifnot(length(prepped$ecoregion_of_cell100) == cl$ncell100)
cat("PASS: adj/num byte-identical -- ecoregion_of_cell100 is safe to reuse as-is.\n\n")

cat("per-ecoregion cell counts (K =", prepped$nregion, "):\n")
print(prepped$region_table)

# ------------------------------ fork 1: + ecoregion fields -------------------
cl_eco <- cl
cl_eco$ecoregion_of_cell100 <- prepped$ecoregion_of_cell100
cl_eco$nregion <- prepped$nregion

eco_dir <- paste0(PROJ, "/HPC/", species_token, "_", TAG, "_ecoregion")
dir.create(eco_dir, showWarnings = FALSE, recursive = TRUE)

real_eco <- real
real_eco$constants_list <- cl_eco
saveRDS(real_eco, file.path(eco_dir, paste0("input_data_", species_token, "_", TAG, "_ecoregion.RDS")))
cat("\nwrote", file.path(eco_dir, paste0("input_data_", species_token, "_", TAG, "_ecoregion.RDS")), "\n")

# ------------------------------ fork 2: unmodified, v2 national_scalar --------
ns_dir <- paste0(PROJ, "/HPC/", species_token, "_", TAG, "_national_scalar")
dir.create(ns_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(real, file.path(ns_dir, paste0("input_data_", species_token, "_", TAG, "_national_scalar.RDS")))
cat("wrote", file.path(ns_dir, paste0("input_data_", species_token, "_", TAG, "_national_scalar.RDS")), "\n")

cat("\ndone. HPC/", species_token, "_v2b/ base bundle untouched; originals untouched.\n", sep = "")
