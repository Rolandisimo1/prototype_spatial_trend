#!/usr/bin/env Rscript
# =============================================================================
# integration_helper_v2.R  --  Fix 1: iNat effort/count column-order bug
# -----------------------------------------------------------------------------
# Forked from Arielle's `integration_helper.R` (never touches her original).
# Only TWO functions are changed here -- `make_inat_cell_year_matrix()` and
# `make_inat_effort_matrix()` -- both copied byte-for-byte from her file
# except for the single documented fix in each. Every other function in
# `integration_helper.R` (including `prep_inat_data_grid()`, which Fix 2
# replaces separately in `prep_inat_data_grid_v2.R`) is UNCHANGED and must
# still be sourced from her original file.
#
# STATUS: this fix (`names_sort = TRUE`) was already applied by hand,
# directly on Hazel, in forked copies of integration_helper.R used by the
# six currently-running/converged v1fix jobs (bobcat/WTD/moose x
# national_scalar/ecoregion, job IDs 499944-499949, launched ~2026-08-03).
# THIS FILE did not exist yet at that time -- it captures that same fix as a
# real, versioned, committed artifact so it lives in this repo instead of
# only as a hand-edited file on the cluster with no local record. Verify
# this file's two functions match whatever is live on Hazel before treating
# them as interchangeable; they were re-derived from the diagnosis
# (`inat_effort_matrix_fix.diff`) and Arielle's original, not copied off
# the cluster.
#
# ROOT CAUSE (full detail in inat_effort_matrix_fix.diff):
#   make_inat_cell_year_matrix() (builds inat_y) and make_inat_effort_matrix()
#   (builds inat_effort) both pivot the same underlying data with
#   pivot_wider(names_from = Year, ...), but neither set names_sort = TRUE.
#   Without it, column order follows first-appearance-by-row, not calendar
#   order. make_inat_cell_year_matrix happened to come out right BY ACCIDENT
#   of its input already being Year-ascending; make_inat_effort_matrix groups
#   by cell50 first, which reorders rows and silently displaces any year
#   missing from the smallest-ID cell50 to the END of the matrix. Confirmed
#   to corrupt moose (2008/2010/2011 shifted to the tail) and white-tailed
#   deer (2010 shifted to the tail) -- both read by every downstream
#   positional index (year_vals[t], n_cells_year[t], inat_cells_by_year[,t])
#   as if they were the LAST years in the series instead of the correct
#   early ones. Bobcat was unaffected only because its smallest in-range
#   cell50 happens to have full 18-year coverage, not because the original
#   code was correct.
#
# USAGE:
#   source("integration_helper.R")     # defines everything else, UNCHANGED
#   source("integration_helper_v2.R")  # OVERWRITES make_inat_cell_year_matrix
#                                       # and make_inat_effort_matrix with the
#                                       # patched versions below (source this
#                                       # SECOND, after the original, so the
#                                       # patched definitions win)
# =============================================================================

# Make the cell x year matrix for model input -- FIXED: names_sort = TRUE
make_inat_cell_year_matrix <- function(df, species) {

  mat_df <- df %>%
    select(cell50, Year, all_of(species)) %>%
    pivot_wider(
      names_from = Year,
      values_from = all_of(species),
      values_fill = NA,
      names_sort = TRUE   # FIX: enforce chronological column order
                          # explicitly, rather than relying on row order
                          # being Year-major by construction (true for this
                          # function's current input by accident, not
                          # guaranteed by anything that would catch a future
                          # change to the upstream data prep).
    ) %>%
    arrange(cell50)

  mat <- as.matrix(mat_df[,-1])
  rownames(mat) <- mat_df$cell50

  return(mat)
}

# Make the effort matrix for model input -- FIXED: names_sort = TRUE
# (this is the function that actually corrupted moose and WTD)
make_inat_effort_matrix <- function(df) {
  # df must have columns: cell50, Year, effort

  effort_df <- df %>%
    group_by(cell50, Year) %>%
    summarise(effort_sum = sum(effort, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = Year,
      values_from = effort_sum,
      values_fill = 0,  # fill zeros where no effort
      names_sort = TRUE  # FIX: the actual bug fix. Without this, column
                         # order follows first-appearance-by-cell50 (this
                         # function groups by cell50 BEFORE pivoting, which
                         # reorders rows), silently displacing any Year
                         # missing from the lowest-ID cell50 to the end of
                         # the matrix -- confirmed to corrupt moose
                         # (2008/2010/2011 shifted to the tail, read as if
                         # they were the LAST 3 years by every downstream
                         # positional index) and white-tailed deer (2010
                         # shifted to the tail, read as the last year).
    ) %>%
    arrange(cell50)

  # convert to matrix
  mat <- as.matrix(effort_df[,-1])
  rownames(mat) <- effort_df$cell50

  return(mat)
}

# GUARDRAIL (recommended addition from inat_effort_matrix_fix.diff, cheap
# insurance against this class of bug recurring under a future refactor of
# either function above). Call this immediately after building both
# matrices in the data-prep script, e.g.:
#   inat_y <- make_inat_cell_year_matrix(inat_df, species)
#   inat_effort <- make_inat_effort_matrix(inat_df)
#   assert_inat_matrices_aligned(inat_y, inat_effort)
assert_inat_matrices_aligned <- function(inat_y, inat_effort) {
  stopifnot(
    "inat_y and inat_effort column (year) order diverged" =
      identical(colnames(inat_y), colnames(inat_effort)),
    "inat_y and inat_effort row (cell50) order diverged" =
      identical(rownames(inat_y), rownames(inat_effort))
  )
  invisible(TRUE)
}
