#!/usr/bin/env Rscript
# =============================================================================
# run_data_prep_v2.R  --  data-prep entrypoint with BOTH fixes wired in
# -----------------------------------------------------------------------------
# Combines Fix 1 (iNat effort/count column-order bug) and Fix 2 (range-mask
# replacement) into one entrypoint, sourced in the correct order. This is
# the file to run/adapt on Hazel for a real re-fit -- do not source
# integration_helper.R alone and expect either fix to be present.
#
# Neither fix has been run against real data as of this commit. This script
# is the wiring/entrypoint; running it is the next step, on Hazel, where the
# real iNat pull and grid files live.
#
# USAGE (on Hazel, from the species prep directory):
#   Rscript run_data_prep_v2.R <species> <taxon_key_path>
#
# WHAT THIS DOES, IN ORDER:
#   1. Sources Arielle's integration_helper.R UNCHANGED (defines everything
#      except the two functions patched below).
#   2. Sources integration_helper_v2.R, which OVERWRITES
#      make_inat_cell_year_matrix() and make_inat_effort_matrix() with the
#      names_sort = TRUE fix (Fix 1). Also defines
#      assert_inat_matrices_aligned(), called after both matrices are built.
#   3. Sources prep_inat_data_grid_v2.R, which defines
#      prep_inat_data_grid_v2() -- Fix 2, the range-mask replacement -- as a
#      NEW function name (does not overwrite the original
#      prep_inat_data_grid(), so both remain callable for comparison if
#      needed).
#   4. Runs prep_inat_data_grid_v2() for the requested species, then builds
#      inat_y / inat_effort with the patched functions and asserts they're
#      aligned before writing anything out.
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop("Usage: Rscript run_data_prep_v2.R <species> <taxon_key_path>")
species <- args[1]
taxon_key_path <- args[2]

# --- source in the documented order ---
source("integration_helper.R")       # Arielle's original, unchanged
source("integration_helper_v2.R")    # Fix 1: overwrites the two matrix builders
source("prep_inat_data_grid_v2.R")   # Fix 2: new function, prep_inat_data_grid_v2()

cat("Loaded: integration_helper.R + Fix 1 (effort/count sort) + Fix 2 (range mask replacement)\n")
cat("Species:", species, "\n")

taxon_key <- read_csv(taxon_key_path)

# Fix 2: range-mask-free iNat grid (captive + spatial-isolation filtering instead)
inat_cell_summary <- prep_inat_data_grid_v2(taxon_key, species, redo = TRUE)
cat("prep_inat_data_grid_v2() complete --", nrow(inat_cell_summary), "cell x year rows\n")

# Fix 1: chronologically-correct matrices, with the guardrail assertion
inat_y      <- make_inat_cell_year_matrix(inat_cell_summary, species)
inat_effort <- make_inat_effort_matrix(inat_cell_summary)
assert_inat_matrices_aligned(inat_y, inat_effort)
cat("make_inat_cell_year_matrix() / make_inat_effort_matrix() aligned and verified --",
    ncol(inat_y), "years x", nrow(inat_y), "cells\n")

cat("\nDONE. inat_y / inat_effort / inat_cell_summary are in memory --",
    "wire these into the existing bundle-forking step",
    "(00b_prep_ecoregion.R / 00b2_join_ecoregion_real.R equivalent)",
    "to produce a real input_data_<species>_v2.RDS.\n")
