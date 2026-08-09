# Two iNaturalist data-pipeline issues found during the ecoregion-trend prototype

Found while validating the spatially-varying trend extension against real
white-tailed deer and moose data. Both are in Arielle's original
`integration_helper.R` (unmodified by the ecoregion-trend prototype — the
model code itself was audited separately and confirmed byte-identical to the
original where unextended). **Neither fix has been applied anywhere.** Both
are review artifacts for the pipeline owner to evaluate and apply
deliberately. The WTD chains that were running on the affected code (Issue 1)
have been killed rather than allowed to finish on data known to be corrupted.

---

## Issue 1: iNaturalist effort data gets shuffled out of calendar order

**Impact:** silently corrupts the fitted temporal trend for any species whose
smallest in-range grid cell lacks full 18-year coverage. Confirmed affected:
**moose** (2008, 2010, 2011 displaced) and **white-tailed deer** (2010
displaced). Confirmed NOT triggered for bobcat, but the same code is running
underneath its fit too — it simply got lucky.

### What it does

`inat_effort` (the all-species reporting-effort denominator, used to
calibrate how many observations of a species we'd expect given search
activity) and `inat_y` (the species' actual observation counts) are both
built as cell-by-year matrices via `pivot_wider()`. Every downstream step
(`year_vals[t]`, `n_cells_year[t]`, `inat_cells_by_year[,t]`, and the model
likelihood itself) reads these matrices by **column position** — column 3 is
assumed to be year 3 (2010), with no name-based check anywhere.

`pivot_wider()` does not sort its output columns chronologically unless
explicitly told to. `make_inat_cell_year_matrix()` (builds `inat_y`) happens
to scan the data in year-ascending order already, so its columns come out
correct **by accident**. `make_inat_effort_matrix()` (builds `inat_effort`)
groups by `cell50` first, which reorders the underlying rows — so its column
order ends up determined by whichever years the *first* (lowest-ID) grid cell
happens to have data for. If that cell is missing an early year (likely for
a species with a small, patchy range), that year gets pushed to the **end**
of the effort matrix instead of its correct position — misaligning it against
`inat_y` and every other year-indexed input.

### Original code (`integration_helper.R`, unmodified, current on Hazel)

```r
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
```

### Proposed fix

```r
make_inat_cell_year_matrix <- function(df, species) {

  mat_df <- df %>%
    select(cell50, Year, all_of(species)) %>%
    pivot_wider(
      names_from = Year,
      values_from = all_of(species),
      values_fill = NA,
      names_sort = TRUE   # FIX: enforce chronological column order explicitly,
                          # rather than relying on row order being Year-major
                          # by construction (true today, but not guaranteed by
                          # anything the type system or a future refactor
                          # would catch).
    ) %>%
    arrange(cell50)

  mat <- as.matrix(mat_df[,-1])
  rownames(mat) <- mat_df$cell50

  return(mat)
}

make_inat_effort_matrix <- function(df) {
  # df must have columns: cell50, Year, effort

  effort_df <- df %>%
    group_by(cell50, Year) %>%
    summarise(effort_sum = sum(effort, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(
      names_from = Year,
      values_from = effort_sum,
      values_fill = 0,  # fill zeros where no effort
      names_sort = TRUE  # FIX: this is the actual bug fix. Without this,
                         # column order follows first-appearance-by-cell50,
                         # which silently displaces any Year missing from
                         # the lowest-ID cell50 to the end of the matrix --
                         # confirmed to corrupt moose (2008/2010/2011 shifted
                         # to the tail) and WTD (2010 shifted to the tail).
    ) %>%
    arrange(cell50)

  # convert to matrix
  mat <- as.matrix(effort_df[,-1])
  rownames(mat) <- effort_df$cell50

  return(mat)
}
```

### Recommended additional guardrail

In `data_prep()` (or the species-specific prep runner), immediately after
both matrices are built:

```r
inat_y <- make_inat_cell_year_matrix(inat_df, species)
inat_effort <- make_inat_effort_matrix(inat_df)

# GUARDRAIL: inat_y and inat_effort must be indexed identically -- every
# downstream use (year_vals[t], n_cells_year[t], inat_cells_by_year[,t])
# treats column position as calendar year with no name-based check. A
# silent divergence here previously baked a multi-year column misalignment
# into two converged bundles (moose, WTD) without any error or warning.
stopifnot(
  "inat_y and inat_effort column (year) order diverged" =
    identical(colnames(inat_y), colnames(inat_effort)),
  "inat_y and inat_effort row (cell50) order diverged" =
    identical(rownames(inat_y), rownames(inat_effort))
)
```

### What this fix does not cover

Does not retroactively repair `input_data_moose_ecoregion.RDS`,
`input_data_moose_national_scalar.RDS`, `input_data_wtd_ecoregion.RDS`, or
`input_data_wtd_national_scalar.RDS` — all four bundles were built from the
buggy code path. Repairing them requires re-running data prep for moose and
WTD with the fix applied, then re-forking the ecoregion/national_scalar
bundles and refitting — a separate, larger decision, not made here. The WTD
chains have already been killed; moose's converged posteriors
(`moose_ecoregion_posterior.csv`, etc.) are being held as diagnosis reference
only, not as trustworthy results.

---

## Issue 2: the range-polygon mask discards real observations outside a species' mapped range

**Impact:** for species with small or patchy ranges, this can discard a large
fraction of real, non-captive observations — including genuine range
expansions, which are exactly the kind of signal a temporal-trend model
should be able to detect. Quantified for moose below.

### What it does

`prep_inat_data_grid()` rasterizes each species' IUCN range polygon onto the
model's 50km grid, then sets that species' observation count to `NA` (not
zero) in every grid cell falling outside the polygon. An `NA` cell is dropped
from the likelihood entirely — treated as "no information," not "confirmed
absence." Note only the **target species'** column is masked this way; the
`effort` column (built from the same file, used for all species) is left
unmasked.

### Quantified for moose

Using the real observation data (2008–2025, all years combined), before any
range masking is applied:

| | in current range polygon | outside current range polygon |
|---|---|---|
| observations | 7,292 (63%) | **4,300 (37%)** |
| occupied grid cells | 278 | 149 |
| per-cell count (median / mean / max) | 6 / 26.2 / 1,290 | 3 / 28.9 / 751 |

The out-of-range observations are not scattered noise — their per-cell count
distribution closely resembles the in-range distribution, and they form a
dense, spatially coherent cluster in Colorado, consistent with Colorado's
known introduced/expanding moose population (established there since the
1970s). All 4,300 of these real observations are currently discarded by the
mask. See `fig15_moose_range_mask_review.png` for the map and histogram.

### Original code (`integration_helper.R`, `prep_inat_data_grid()`, unmodified)

```r
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
```

### Proposed fix (forked as a new function; not applied)

Drop the range-polygon mask; keep a captive-record filter (already available
in the raw iNat download but not currently used) and add a spatial-isolation
outlier filter as the mask's replacement -- this targets true vagrants/
mis-IDs without assuming the historical range polygon is current or
complete, so genuine range expansions (like Colorado's moose) are retained.

```r
for (i in 1:nrow(taxon_key)) {
  this_spec_dat <- inat_cells %>%
    filter(sciname == taxon_key$sci_name[i],
           captive_cultivated == FALSE) %>%          # NEW: exclude known-captive records
    group_by(Year) %>%
    count(cell50, cell100)
  colnames(this_spec_dat)[4] <- taxon_key$common_name_clean[i]

  inat_cell_summary <- left_join(inat_cell_summary, this_spec_dat, by = c("Year","cell50", "cell100"))
  inat_cell_summary[[taxon_key$common_name_clean[i]]][is.na(inat_cell_summary[[taxon_key$common_name_clean[i]]])] <- 0

  # RANGE-POLYGON MASK REMOVED. Proposed replacement: a spatial-isolation
  # outlier filter -- exclude only observations with zero same-species
  # records within 100km and within a 3-year window (flags likely
  # vagrants/mis-IDs without discarding real range-edge/expansion activity).
  # Not yet implemented as code; parameters (100km / 3yr) are a starting
  # proposal, open to adjustment.
}
```

### Open question: minimum observations per grid cell

Considered as a substitute noise filter, but a flat per-cell threshold
penalizes genuinely narrow-but-real detections (e.g. a newly colonizing
population's first 1-2 sightings in a cell) the same way it penalizes actual
noise. The spatial-isolation filter above is a better-targeted primary
filter. A minimum-count threshold may still be useful as a secondary,
reporting-stage flag (e.g., grey out single-detection cells on a map) --
open for discussion.

### Status

Not applied anywhere. This is a genuine pipeline design change (not a
bugfix), so it needs deliberate review and a decision before touching any
bundle. Quantifying this for a few more narrow-range fleet species (elk,
pronghorn, wolf, fisher, kit fox) before deciding is a reasonable next step,
since narrow-range species are the ones most sensitive to range-polygon
accuracy.

---

## Addendum (2026-08-09): range-mask scope confirmed across all 3 species

Launch record for the current v1fix real fits (from the original relaunch
prompt, ~2026-08-03) confirms **Fix 1 only** was applied before all 6 jobs
went out; Fix 2 (this range-mask replacement) was explicitly scoped out of
that launch and deferred pending team review:

| Species | Model | Job ID |
|---|---|---|
| Bobcat | national_scalar | 499944 |
| Bobcat | ecoregion | 499945 |
| WTD | national_scalar | 499946 |
| WTD | ecoregion | 499947 |
| Moose | national_scalar | 499948 |
| Moose | ecoregion | 499949 |

So the range-mask gap quantified for moose in Issue 2 above applies
identically to bobcat and WTD's currently-running/converged v1fix fits —
this is not moose-specific. The 12 additional windowed-refit bundles (5yr/
10yr x 3 species x 2 models, built 2026-08-09, not yet launched) are forked
from the same 18-year source data and inherit the same gap.

## Addendum 2 (2026-08-09): both fixes now committed as real code

`inat_effort_matrix_fix.diff` (Issue 1) was a diagnosis-only artifact -- its
`names_sort = TRUE` patch was applied by hand directly on Hazel for the six
v1fix launches, but never existed as a versioned file until now.
[integration_helper_fix1.R](integration_helper_fix1.R) forks the two affected
functions from Arielle's original with the fix applied and adds the
recommended `assert_inat_matrices_aligned()` guardrail. [Verify this file's
functions match whatever is live on Hazel before assuming they're
interchangeable -- it was re-derived from the diff and the original, not
pulled from the cluster.]

[prep_inat_data_grid_v2.R](prep_inat_data_grid_v2.R) (Issue 2's fix, drafted
earlier this session) and `integration_helper_fix1.R` are wired together in
[run_data_prep_v2.R](run_data_prep_v2.R), the single entrypoint to run on
Hazel for a real re-prep. **Neither fix has been run against real data as
of this commit** -- both are logic-tested (synthetic data only) but not yet
exercised on Hazel's actual iNat pull, and no bundle has been re-forked or
re-fit with either.
