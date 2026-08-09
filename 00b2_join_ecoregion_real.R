#!/usr/bin/env Rscript
# =============================================================================
# 00b2_join_ecoregion_real.R
# Real-data extension of 00b_prep_ecoregion.R. Read-only w.r.t. Arielle's
# original per-species directories (HPC/<species>/) -- writes new, separately
# tokened bundles under HPC/<species>_ecoregion/ and HPC/<species>_national_scalar/,
# never touching HPC/<species>/ itself.
#
# What this does NOT do: a fresh EPA/CEC spatial join. 00b_prep_ecoregion.R
# already built ecoregion_of_cell100 (length 908) against the FULL shared
# continental cell100 grid (grid100.tif via 00_prep_sim_inputs.R) -- not a
# simulation-only reduced set, despite that script's header comment. Verified
# 2026-07-18: WTD's and moose's own constants_list$adj/$num (the CAR
# adjacency graph that defines cell100 identity) are byte-identical to
# bobcat's, which is what prepped_sim_inputs.RDS's ecoregion_of_cell100 was
# built against. So this script's real job is (a) proving that identity
# holds for the species bundle at hand via a hard adj/num check, then
# (b) appending the already-computed ecoregion_of_cell100/nregion, and
# (c) forking the bundle into two new-token copies (one plus ecoregion
# fields, one an untouched copy for the national_scalar baseline).
#
# RUN (login node; no internet/spatial libs needed here -- nimble_env is fine):
#   Rscript 00b2_join_ecoregion_real.R <species_dir> <new_token>
# Example:
#   Rscript 00b2_join_ecoregion_real.R white-tailed_deer wtd
# Writes:
#   HPC/<new_token>_ecoregion/input_data_<new_token>_ecoregion.RDS
#   HPC/<new_token>_national_scalar/input_data_<new_token>_national_scalar.RDS
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("usage: Rscript 00b2_join_ecoregion_real.R <species_dir> <new_token>")
}
species_dir <- args[1]   # e.g. "white-tailed_deer" -- matches HPC/<species_dir>/input_data_<species_dir>.RDS
new_token   <- args[2]   # e.g. "wtd" -- new bundles named input_data_<new_token>_{ecoregion,national_scalar}.RDS

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")

# ------------------------------ load inputs ---------------------------------
real_path <- paste0(PROJ, "/HPC/", species_dir, "/input_data_", species_dir, ".RDS")
stopifnot(file.exists(real_path))
real <- readRDS(real_path)
cl <- real$constants_list

prepped <- readRDS(paste0(PROTO_DIR, "/prepped_sim_inputs.RDS"))
stopifnot(!is.null(prepped$ecoregion_of_cell100), !is.null(prepped$nregion))

cat("=== identity check: real bundle's cell100 CAR graph vs prepped_sim_inputs' ===\n")
cat("real ncell100:", cl$ncell100, " prepped ncell100:", prepped$constants_list$ncell100, "\n")
stopifnot(cl$ncell100 == prepped$constants_list$ncell100)
stopifnot(identical(cl$adj, prepped$constants_list$adj))
stopifnot(identical(cl$num, prepped$constants_list$num))
stopifnot(length(prepped$ecoregion_of_cell100) == cl$ncell100)
cat("PASS: adj/num byte-identical -- ecoregion_of_cell100 is safe to reuse as-is,\n")
cat("      no fresh spatial join needed for this bundle.\n\n")

cat("per-ecoregion cell counts (K =", prepped$nregion, "):\n")
print(prepped$region_table)

# ------------------------------ fork 1: + ecoregion fields -------------------
cl_eco <- cl
cl_eco$ecoregion_of_cell100 <- prepped$ecoregion_of_cell100
cl_eco$nregion <- prepped$nregion

eco_dir <- paste0(PROJ, "/HPC/", new_token, "_ecoregion")
dir.create(eco_dir, showWarnings = FALSE, recursive = TRUE)

real_eco <- real
real_eco$constants_list <- cl_eco
saveRDS(real_eco, file.path(eco_dir, paste0("input_data_", new_token, "_ecoregion.RDS")))
cat("\nwrote", file.path(eco_dir, paste0("input_data_", new_token, "_ecoregion.RDS")), "\n")

# ------------------------------ fork 2: unmodified, new token -----------------
# for the national_scalar baseline fit -- identical constants_list, just
# relocated under its own token so it never collides with HPC/<species_dir>/'s
# pre-existing (partial, differently-named) fit artifacts.
ns_dir <- paste0(PROJ, "/HPC/", new_token, "_national_scalar")
dir.create(ns_dir, showWarnings = FALSE, recursive = TRUE)

saveRDS(real, file.path(ns_dir, paste0("input_data_", new_token, "_national_scalar.RDS")))
cat("wrote", file.path(ns_dir, paste0("input_data_", new_token, "_national_scalar.RDS")), "\n")

cat("\ndone. Original HPC/", species_dir, "/ untouched.\n", sep = "")
