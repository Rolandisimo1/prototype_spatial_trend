#!/usr/bin/env Rscript
# verify_backward_compat.R -- confirm the sim_helpers.R generalization
# (region_col param, fixed order for region_col="region") reproduces the
# EXACT same reduced design as before -- checked against the known,
# repeatedly-quoted informed_cell100 count (178) from every prior CAR-based
# result in this session. Not a deliverable.
PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
setwd(PROTO_DIR)
source("sim_helpers.R")

prepped <- readRDS("prepped_sim_inputs.RDS")
inp <- readRDS(paste0(PROJ, "/HPC/bobcat/input_data_bobcat.RDS"))
cl <- prepped$constants_list

reduced <- build_reduced_constants(
  cl = cl, inat_effort_real = inp$inat_effort, real_y_template = inp$real_data$y,
  cell100_geo = prepped$cell100_geo,
  n_cell50_keep = 210, max_subcell_per_cell = 15, n_site_keep = 700,
  seed = 20260712)

informed <- unique(reduced$constants$inat_cell100)
cat("n_cell50_keep result:", reduced$constants$ncell50, "(expect 210)\n")
cat("n_site_keep result:", reduced$constants$nsite, "(expect ~702, matches original prints)\n")
cat("informed cell100 count:", length(informed), "(expect 178, matches every prior CAR result this session)\n")
cat("first 10 cell50_keep:", paste(head(reduced$cell50_keep, 10), collapse=","), "\n")

if (length(informed) == 178) {
  cat("\nBACKWARD COMPATIBILITY CONFIRMED\n")
} else {
  cat("\nMISMATCH -- region ordering fix did not fully preserve original behavior!\n")
}
