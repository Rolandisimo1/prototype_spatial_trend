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
Hazel for a real re-prep.

## Addendum 3 (2026-08-09): dropping the mask entirely was untractable --
## replaced with a data-driven presence mask (Fix 2b)

The first version of Fix 2 (drop the IUCN mask, keep every cell as a real
0/positive count, rely only on the captive + isolation filters) WAS run
against real data on Hazel this session, isolated from the live fleet. It
reproduced July's numbers exactly: ncell50 goes 382 -> 3,322 for moose (the
full national grid), Colorado goes from 0 moose records (100% masked) to
3,145 recovered records across 43 cells, guardrail passed, chronological
columns confirmed. But only 422 of those 3,322 cells actually have ANY
moose records -- the other ~2,900 are genuine, real zeros (cells with no
moose activity nationwide, not excluded data). Converting those cells from
NA (dropped from the likelihood entirely, as under the old mask) to real
zeros means NIMBLE has to build a live `mu[g,t]` node for each of them at
every one of 18 years, which is what caused `nimbleModel()`'s graph-build
to hang for all six species/model combinations back in July -- confirmed
species-independent, since moose (tiny camera dataset) hung just as hard as
bobcat and WTD.

**Fix 2b, replacing the "drop the mask entirely" approach:** instead of
choosing between the static IUCN polygon and the full 3,322-cell national
grid, build the mask from the species' OWN iNat presence records. A cell50
counts as "in range" only if it accumulates >= `PRESENCE_MIN_RECORDS` total
records of that species across all 18 years, counted AFTER the
spatial-isolation filter has already removed likely vagrants/mis-IDs.
Cells below that threshold go back to NA -- same behavior as the old mask,
but the boundary is drawn from real observations instead of a stale range
polygon.

**At the time this addendum was first written (2026-08-09), the threshold
was set to 3 and only logic-tested on synthetic data.** See Addendum 4
below for the real-data numbers and the resulting threshold revision.

## Addendum 4 (2026-08-10): PRESENCE_MIN_RECORDS lowered from 3 to 1 --
## real Hazel numbers showed the stricter cut wasn't buying tractability

Fix 2b was run against real moose data on Hazel (still isolated from the
live fleet). Measured `ncell50` at both candidate thresholds:

| Threshold | ncell50 | vs. the 3,322-cell wall that broke `nimbleModel()` |
|---|---|---|
| >= 3 records | ~290 | far below |
| >= 1 record (post-isolation-filter) | ~422 | far below |

Both are enormously smaller than 3,322 -- the gap between 290 and 422 is
unlikely to be the difference between "the model build succeeds" and
"it doesn't." That means the >=3 cut's *stated* purpose (tractability) was
not actually being served by the stricter threshold; what it WAS doing was
re-excluding real, thin edge-presence cells -- concretely, Utah's tail of
1-2-record cells, the same *kind* of signal (a real, sparsely-documented
population near a range edge) that motivated dropping the IUCN mask for
Colorado in the first place.

**Decision (2026-08-10): `PRESENCE_MIN_RECORDS` lowered to 1.** The spatial-
isolation filter (100km/3yr, applied BEFORE this threshold) remains the
layer responsible for screening out true vagrants/mis-IDs. `>=1` here means
"at least one record survived that screening," not "no additional
screening at all" -- it is not equivalent to dropping the mask entirely
(that was the ncell50=3,322 version, already ruled out above). The
nimbleModel() build test at ~422 cells (queued as of this addendum) is the
concrete next check on whether this actually resolves tractability for
moose, and by extension bobcat/WTD.

**Also confirmed on Hazel, real data, both thresholds:** the captive-record
filter is a no-op across the ~3M-row master file (`captive_cultivated` is
all-FALSE) -- consistent with the raw pull already being restricted to
iNat's "Research Grade" quality tier, which requires community consensus
the organism is wild. The isolation filter itself flags 0.3-10% of records
depending on species (0.41% for moose specifically).

---

## Issue 3 (found 2026-08-14): soil_clay/soil_silt/soil_sand are compositional
## -- including all three in occ_beta builds an unidentified ridge into the
## design matrix

**This is a separate class of problem from Issues 1-2 above** -- not an
iNat data-pipeline bug, but a latent flaw in the occupancy-covariate design
matrix itself, present in the ORIGINAL 10-covariate KEEP list
(`make_reduced_input.R`, "final" decision) since before any of the fixes
above existed.

### What it does

Soil clay, silt, and sand fractions are compositional data: by construction
they sum to ~100% of a soil sample, so any two of the three fully determine
the third. `make_reduced_input.R`'s KEEP list retains all three
(`soil_clay`, `soil_silt`, `soil_sand`) as separate linear terms in
`occ_beta`. Because the model code indexes `occ_beta` purely by position
(`occ_beta[1:numOccCovars]` in `HPC_<species>_<model>_chain1.R`, no
name-based lookup anywhere), nothing catches that two of these three
coefficients can trade weight against each other along the "sum" direction
with zero change in the model's fit -- an exactly-unidentified ridge.

### How it was found

`moose_v2b_national_scalar` and `moose_v2b_ecoregion` (both built on the
Fix-2b presence-mask data) plateaued at the 50k-iteration cap with
`occ_beta[8]` and `occ_beta[10]` -- confirmed against the bundle's own
column names to be `soil_clay` and `soil_sand` -- stuck at R-hats matching
to four decimal places. Matching R-hats across two parameters is the
classic ridge signature: the chains are drifting along a shared direction,
not failing to mix independently.

Direct correlation check, read off each species' actual camera-site
covariate data:

| Species | n camera sites | cor(soil_clay, soil_sand) |
|---|---|---|
| Moose  | 1,654  | -0.9966 |
| Bobcat | 20,531 | -0.5662 |
| WTD    | 21,559 | -0.5664 |

The ridge is present in **all three species' bundles** at essentially the
same correlation (it depends only on camera-site soil composition, which
the presence-mask fix never touches -- Fix 2b changes iNat cell membership,
not camera sites). It was not *binding* for bobcat/WTD because their camera
sites span enough geography for clay/sand to vary semi-independently;
moose's narrow northern camera-site band pushed the correlation to -0.997
and made the ridge the actual constraint on convergence. **This means the
same ridge is latent in every v1fix bobcat/WTD fit too** (unconfirmed
whether it was ever binding there, since v1fix ecoregion/national_scalar
non-convergence has other causes -- see the WTD ecoregion funnel discussion
elsewhere -- but the redundancy in the design matrix is there regardless of
whether it happens to be the tightest constraint in a given fit).

### Fix

Drop `soil_silt` from the covariate list -- keep `soil_clay` and
`soil_sand` as the two explicit texture-axis effects; silt becomes the
implicit reference level. This is a real, permanent correction to the
covariate design, applicable to any species/model, not a moose-specific
patch. Implemented in
[make_reduced_input_v2.R](make_reduced_input_v2.R), forked from
`make_reduced_input.R`, taking an already-10-covariate bundle down to 9.

### Status

Not yet applied to any running or completed fit as of this addendum.
`moose_v2b_national_scalar`'s `CONVERGED.flag` (written 2026-08-12 off a
20k/20k/20k read that squeaked under threshold) is a false positive --
the very next equal-window check came back 1.1394 with `occ_beta[8]`,
`occ_beta[10]`, `theta0`, and `theta1` all over 1.1. Both moose v2b fits
are being refit from scratch on the corrected 9-covariate design rather
than continued. bobcat/WTD v2b builds should apply this fix before their
first launch, not after a plateau.

## Addendum to Issue 3 (2026-08-14, same day): both the mechanism and the
## drop choice above were wrong -- corrected via a direct VIF check

The "Fix" and framing above were written from the "compositional data, sums
to ~100%" argument alone, without checking the actual retained-pair
correlation or running this against bobcat/WTD data first. A direct
VIF/correlation check on all three species' real bundles caught two errors
before either mattered in practice (no fit was launched on the wrong
version):

**1. Not a universal compositional identity.** R² of `sand ~ clay + silt`
is 0.9995 for moose (nsite=1,654, narrow northern camera-site footprint --
a genuine near-collinearity specific to that footprint) but only 0.606
(bobcat, nsite=20,531) and 0.617 (WTD, nsite=21,559). Max VIF across all 10
original covariates: moose 2,267 (severe), bobcat 3.0, WTD 3.1 (both
completely benign). So this is **not** "exactly unidentified for every
species by construction" -- it's a moose-specific near-collinearity that
happens to be severe enough to bind. The claim above that "the ridge is
latent in every v1fix bobcat/WTD fit too" is not supported by their data.

**2. Dropping `soil_silt` was the wrong fraction.** It leaves `soil_clay` +
`soil_sand` retained together -- exactly the -0.9966-correlated pair that
was `occ_beta[8]`/`occ_beta[10]`, the two parameters that actually
plateaued both moose fits. Max VIF over the retained pair, moose, by drop
choice:

| Drop | max VIF (moose) | retained pair correlation |
|---|---|---|
| none (original 10-cov) | 2,267 | -- |
| soil_silt (this addendum's original fix) | 152 | cor(clay,sand) = -0.997 |
| soil_sand | 4.5 | cor(clay,silt) = -0.297 |
| soil_clay | 4.5 | cor(silt,sand) = +0.220 |

Dropping silt takes moose's max VIF from 2,267 to 152 -- it removes the
exact singularity but leaves a severe ridge that would very likely plateau
again after another multi-day fit.

**Corrected fix: drop `soil_sand`, keep `soil_clay` + `soil_silt`.** Best
conditioning across all three species (moose max VIF 4.5; bobcat/WTD both
1.9). `make_reduced_input_v2.R` has been updated to drop `soil_sand`
accordingly -- the file now matches this corrected fix, not the original
text above.

One conclusion from the original write-up still holds: applying one
uniform 9-covariate design across all five v2b builds remains the right
call for report comparability, and costs bobcat/WTD nothing, since their
soil covariates were never problematically collinear to begin with.

## Addendum, part 2 (2026-08-14, same day): fleet-wide validation --
## mechanism confirmed, one uniform design works for all 8 species

Two follow-up checks, run against all 8 fleet species (moose, bobcat, WTD,
grey wolf, fisher, elk, pronghorn, kit fox) rather than just the three with
launched v2b builds:

**1. Why moose's clay+silt+sand sums to ~100% but bobcat/WTD's doesn't.**
Not an algebraic identity (a hard identity would give an exact, constant
sum for every species). Recovering the raw 0-1 fractions from each
species' cached input confirms the sum is bimodal: a spike at exactly 1.0
plus a tail running down toward 0. Moose and grey wolf sit at 1.0 in ~98%
of camera sites; bobcat/WTD only ~46-51%, with ~30% of sites below 0.90.
This is the compositional identity holding almost exactly wherever ground
cover is genuinely just mineral soil (moose/wolf's narrow, forested,
northern ranges), and breaking down wherever a meaningful fraction of the
ground is something else -- rock, wetland/organic soil, urban fill --
which shows up more as a species' site set gets geographically larger and
more diverse (bobcat, WTD).

**2. The ridge is not moose-specific -- it's a narrow-range-species
problem, and it would have hit at least two more fleet species before
their v2b builds even started.** Max VIF on the original 10-covariate
design, by species:

| Species | n sites | max VIF (10-cov) |
|---|---|---|
| Moose | 1,654 | 2,267.4 |
| Grey wolf | 1,391 | 1,732.7 |
| Fisher | 3,142 | 242.9 |
| Elk | 2,296 | 23.3 |
| Kit fox | 964 | 4.5 |
| Bobcat | 20,531 | 3.0 |
| WTD | 21,559 | 3.1 |
| Pronghorn | 1,244 | 2.7 |

Grey wolf and fisher are in the same severe-collinearity class as moose;
elk is over the conventional VIF>10 danger threshold. All three would have
plateaued in the same way moose did, had their v2b builds gone out on the
original 10-covariate design -- caught here before that cost any cluster
time.

**3. One uniform 9-covariate design (drop `soil_sand`) works fleet-wide --
no need for a species-specific or narrow-range-specific reduced covariate
set.** Max VIF after dropping `soil_sand`, across all 8 species: 4.5
(moose, worst case), all others <=3.2. Dropping `soil_sand` ties dropping
`soil_clay` on conditioning everywhere, and was chosen because sand is the
largest mean fraction (most natural implicit reference level).
`make_reduced_input_v2.R` already implements this (drop `soil_sand`, keep
`soil_clay` + `soil_silt`) -- confirmed sufficient for the whole fleet, not
just the three species with builds already in flight.

**Caveat on the "implicit reference level" framing:** this is only
strictly literal for the narrow-range species (moose, wolf, fisher, elk),
where clay+silt+sand really do sum to ~1 and dropping sand means "clay and
silt are deviations from a fixed whole." For bobcat/WTD/pronghorn/kit_fox,
where the sum isn't constant, the retained `soil_clay`/`soil_silt`
coefficients are just two texture-axis effects, not deviations from a
literal 100% reference -- a minor interpretive nuance worth keeping in mind
when discussing coefficients across species, not a statistical problem.

---

## Issue 4 (found 2026-08-20): resumed chain chunks skip re-adaptation
## burn-in, contaminating ~7% of every round's draws with a sampler-restart
## transient

**Fleet-wide, all species/models, every round boundary** -- distinct from
Issues 1-3, this is a bug in the MCMC chunk-and-resume mechanism itself
(`run_chain_chunk()` / `save_chain_state()` / `restore_chain_state()` in
`HPC_<species>_<model>_chain*.R`), not the data pipeline or the covariate
design.

> ### CORRECTION (2026-08-21): this issue is materially worse than described below
>
> As originally written, Issue 4 says the resume mechanism carries the chain's
> node VALUES across a boundary and loses only the samplers' ADAPTATION. **That
> is wrong. Neither is carried across.** At every resume boundary the chain
> restarts from the model's initial values.
>
> The claim that the restore worked came from a probe that read its verdict off
> a **conjugate-sampled node** (`mu` in a toy model). A conjugate sampler draws
> from its full conditional and lands at the posterior in one step *from any
> starting value*, so the probe could not distinguish "state restored" from
> "state reset to init, then immediately jumped" -- the two hypotheses predict
> identical output. It was never evidence for restoration.
>
> Measured on the real production chains instead: `theta0` snaps to exactly
> -5, `theta1` to 1, `overdisp_inat` to 0.1, `year_var` to 0 at each boundary
> and holds while samplers reject. Audited across 164 `chain_*.RDS` in 58 fit
> dirs (`resume_boundary_audit.csv`): **1,132 of 1,234 boundaries automatically
> confirmed restarting from inits**; the other 102 are indeterminate only
> because no inits file could be loaded for them, and their post-boundary
> values sit exactly on those same init constants. **Zero boundaries were
> confirmed clean.**
>
> Likely mechanism (leading hypothesis, not proven): `Cmcmc$mvSaved` is a
> compiled `CmodelValues` behind an external pointer; `saveRDS` writes only the
> R-level shell, the pointer is dead on read-back, and the assignment silently
> no-ops.
>
> **What this invalidates:** Fix part 2 below (the resume burn-in, drafted in
> `run_chain_chunk_fix_resume_burnin.R`) addresses only the lost adaptation. It
> cannot recover the lost position, and by suppressing the visible transient it
> would *hide* the defect while every round still restarts from inits. It is
> not a fix. Fix part 1 (dropping transient draws at extraction) remains valid
> as far as it goes, but it is mitigation of a symptom, not a repair.
>
> Read the rest of this section with that correction applied.

### What it does

`save_chain_state()` persists only `Cmcmc$mvSaved` (current node values) --
never the adaptive samplers' internal proposal-scale/covariance state.
Every resumed chunk rebuilds `nimbleModel -> configureMCMC -> buildMCMC ->
compileNimble` from scratch, so the new `Cmcmc`'s samplers always restart
at their default proposal scale, regardless of what was learned in prior
rounds. The burn-in call (`Cmcmc$run(burnin_once)`) exists only in the
fresh-start branch; the restore branch jumps straight to
`Cmcmc$run(chunk_iter)` with zero re-adaptation burn-in, so roughly the
first 1,000 draws after every resume are the samplers quietly re-learning
their proposal scale, recorded as if they were converged posterior draws.

### How it was found and confirmed

Investigating an apparent mode-split in `bobcat_v2b_national_scalar`'s
first R-hat read (chain 1 alone permanently displaced) led to a block-mean
trace of `theta0` across all three chains: a shared disturbance at block 11
(the first 1,000 draws after the round-1->round-2 resume) in EVERY chain,
not just chain 1 -- chains 2/3 recovered, chain 1 did not. The same
signature (spike at every resume-point block, all chains) was confirmed on
`bobcat_v2b_ecoregion` and `moose_v2b_national_scalar` via `overdisp_inat`
block means -- moose shows it at both of its round boundaries, in all
three chains, perfectly reproducibly.

### Impact

Point estimates are essentially unaffected -- confirmed on both moose v2b
models, the headline trend numbers barely move once the transient draws
are dropped. Credible intervals were inflated on one tail by roughly 30%
(e.g. moose_v2b_national_scalar theta0: as-reported [-5.865, -3.646] ->
cleaned [-4.497, -3.640]). R-hat is distorted in a DIRECTION-DEPENDENT way:
a transient shared by all chains that all chains recover from deflates
R-hat (makes convergence look BETTER than it is); a transient that
permanently displaces one chain (bobcat_v2b_national_scalar chain 1)
inflates it (makes convergence look WORSE than it is, and can look like a
mode-split when it's actually a stuck sampler).

### Fix (two parts)

1. **Immediate, no refit needed:** regenerate the posterior CSVs dropping
   the first ~1,000 draws of every resumed chunk. Pure re-summarization of
   already-collected chain files -- fixes every already-run model's
   intervals without spending any cluster time.
2. **Going forward:** give resumed chunks their own discarded
   re-adaptation burn-in, drafted in
   [run_chain_chunk_fix_resume_burnin.R](run_chain_chunk_fix_resume_burnin.R).
   **Not yet validated** -- the exact NIMBLE `reset`/`resetMV` semantics in
   the fix need confirming against whatever NIMBLE version is loaded on
   Hazel before fleet rollout (see the file's own validation-plan comment):
   apply to one isolated/throwaway resumed chunk first, re-run the same
   block-mean diagnostic that found this bug, and confirm both that the
   transient is gone AND that the restored `mvSaved` state is genuinely
   still being used (not silently reset to initial values). Sent to
   Arielle for review in parallel; not blocking the isolated-chain
   validation test.

### Status

Not yet applied to any file as of this addendum. `bobcat_v2b_national_scalar`
chain 1 -- the one chain that was permanently displaced rather than
recovering -- should be restarted once the fix is validated (still free:
its remaining chain jobs are unstarted in the queue as of 8/20).

---

## Issue 5 (added 2026-09-01): a third range-mask option -- IUCN polygon
## expanded by contiguous iNat presence, not a replacement for Fix 2b

**Not a bug -- a new candidate mask, requested to compare against Fix 2b's
pure-presence approach.** Per a convention used in a different project: the
IUCN polygon serves as a base "in range" set, then expands outward through a
**contiguous chain** of grid100 (100 km) cells that individually clear a
detection threshold (**>2 detections**, after the same isolation/captive
filtering Fix 2b already applies). A cell joins the mask only if reachable
from the polygon through an unbroken chain of qualifying cells -- a
genuinely isolated record far from both the polygon and any other detection
cluster does not pull its own cell in, even if that one cell individually
clears the threshold. 4-connectivity (edge-adjacent neighbors only), per
2026-09-01 decision.

**Why a third option, not a replacement:** Fix 2b (pure presence threshold,
no polygon at all) is already validated on real fits -- `moose_v2b`
converged and recovered the Colorado/Utah cells the original IUCN mask
excluded. This candidate has not been fit to anything yet. The point of
building it is to compare mask choice against real fleet results, not to
assume one is better going in.

**Implementation:** [build_iucn_expanded_mask.R](build_iucn_expanded_mask.R)
-- the core flood-fill function (`flood_fill_expand_mask()`) is logic-tested
against synthetic data (4 cases: chain propagation through qualifying cells,
threshold correctly blocking propagation, an isolated-but-eligible cell
correctly excluded for lack of a connected path, base cells always
retained -- all passed). **The driver is NOT yet wired to real data** -- it
deliberately `stop()`s rather than guess two things that need confirming on
Hazel directly:
1. The exact per-species file naming convention in the IUCN range directory
   (`RANGE_DIR`) -- referenced by `build_moose_unmasked_review.R` but never
   actually read from there in this repo, so the naming pattern is unverified.
2. Where Fix 2b's isolation/captive-**filtered** cell100 detection counts are
   cached per species. The `>2 detections` threshold must be evaluated on the
   SAME filtered counts Fix 2b uses to define "presence" -- recomputing from
   raw, unfiltered records would silently disagree with what the rest of the
   pipeline already calls a valid detection for the same species.

**Once wired:** run for moose, bobcat, and WTD, and compare resulting
cell50/cell100 counts against both the original IUCN-only mask and Fix 2b's
presence-only mask (the moose/bobcat/WTD range-mask-exclusion table already
in Issue 3's addendum gives the baseline to compare against). Do not wire
this into any HPC run script until that three-way comparison has been made
and reviewed -- this is a candidate under evaluation, not a settled fix.

---

## Issue 5 (2026-09-05): the trend term bypasses Goldstein's gold-standard mechanism

**This is a model-structure finding, not a bug in our code.** Our
`model_code_national_scalar.R` is a byte-for-byte fork of Arielle's production
nimbleCode (sha256 recorded in its header), so nothing on our side changed it.
The issue is in how the temporal trend was originally grafted onto Ben's
integrated model.

### What Goldstein et al. actually do

The preprint (bioRxiv 2025.01.17.633640, "Mammal niches are not conserved over
continental scales") is explicit that cameras are the gold standard, and it
names the mechanism. iNaturalist counts are modelled as

    N_c ~ NegBinomial(mu_c * E_c, phi)
    log(mu_c) = theta0 + theta1 * log( sum_{s in c} lambda_s )

and the latent intensity is

    log(lambda_s) = beta0 + x_s' beta_g(s) + eps_g(s),   eps_g ~ CAR(0, sigma^2)

The paper states this formulation "treats the camera trap data as the 'gold
standard' dataset directly informing the latent intensity process, allowing the
iNaturalist data to more weakly influence its value if the two datasets
diverge." The down-weighting is not a prior or a likelihood weight -- it is
architectural. **iNaturalist never writes into `lambda` directly.** It sees
`lambda` only through the `theta0`/`theta1` link, so when the two streams
disagree, `theta0`/`theta1` absorb the disagreement and `lambda` (camera-defined)
holds. Our own fits estimate theta1 at 0.41-0.61 with every CI excluding 1.0,
i.e. the attenuation is active and substantial.

Critically: **Ben's `lambda_s` has no year term.** A full-text search of the
preprint returns zero matches for "temporal trend", "population trend", "trend
over time", "year effect", or "annual trend". The model is purely spatial.

### Where our structure diverges

The temporal trend is an extension added for this project, and it was inserted
*inside* `lambda`, not through the theta link:

    log(lambda_s) = intensity_intercept + x*beta + MWMT + MCMT
                    + total_var_beta * year          <-- trend lives in lambda
    total_var_beta <- year_beta + year_var

The camera occupancy submodel contains `year_beta` only:

    cloglog(psi_i) = link_occ_intercept + ... + year_beta * year_occ_i

`year_var` therefore appears in **no** camera likelihood term. It is informed
solely by iNaturalist -- and because it sits inside `lambda`, it modifies the
gold-standard latent field directly rather than being attenuated by theta1.
This is the one place in the model where iNaturalist has unmediated write access
to the camera-defined intensity process, and it is precisely the parameter that
carries our headline trend.

The priors give no compensating protection; they are symmetric and vague:

    year_beta ~ dnorm(0, sd = sigma_year_beta);  sigma_year_beta ~ dunif(0, 2)
    year_var  ~ dnorm(0, sd = sigma_year_var);   sigma_year_var  ~ dunif(0, 2)

### What it costs us, measured

Decomposition of the reported trend (single-shot v2b national fits, see
`table_trend_decomposition.csv`):

| Species | camera-anchored `year_beta` | iNat-only `year_var` | reported (sum) |
|---|---|---|---|
| Bobcat | -0.0001 (-0.061, 0.061) | -0.242 (-0.336, -0.150) | -0.242 |
| White-tailed deer | +0.155 (0.116, 0.192) | +0.022 (-0.030, 0.080) | +0.177 |
| Moose | -0.121 (-0.313, 0.042) | -0.001 (-0.195, 0.206) | -0.121 |

Bobcat's entire reported decline is the iNaturalist-only term; its cameras
estimate no trend at all, with a tight interval around zero. WTD is the
opposite and is genuinely camera-anchored. Moose's two components each span
zero and are negatively correlated, so only their sum resolves.

`trend_robust_indicator` was intended to flag exactly this, but it is
`step(year_beta/year_var - 1)` -- a ratio whose denominator sits near zero, so
its posterior mean is NaN for all three species. It cannot do the job; the
decomposition table can.

### The tension in fixing it

`year_var` is not purely an artifact channel. It is also the only way the model
can express "abundance changed while occupancy did not" -- and our own
array-level simulation study found occupancy-based trend estimators badly biased
at high abundance for exactly that reason (psi saturation; occupancy overstated
the decline 3-6x at deer-like abundance). Deleting `year_var` would restore
Ben's architecture but re-import saturation bias for abundant species.

The problem is not that `year_var` exists; it is that `year_var` is identified
by one stream only, so real abundance-only change and residual iNaturalist
effort drift are not separable within it.

### Options (needs a team decision -- all three require refits)

1. **Drop `year_var`**: intensity trend = `year_beta`. Exactly Ben's
   architecture; iNaturalist influences the trend only through theta. Cost:
   saturation bias for abundant species, per our own simulation.
2. **Shrink `year_var` hard toward zero**: keep the term, replace
   `sigma_year_var ~ dunif(0, 2)` with a tight prior (e.g. half-normal with sd
   ~0.05) so the trend defaults to the camera value unless iNaturalist evidence
   is strong. Makes cameras primary in the statistical sense while preserving
   the ability to detect abundance-only change. Requires a prior-sensitivity
   check, since the answer will depend on the chosen scale.
3. **Keep as-is, report the decomposition.** No refit; the report states the
   camera-anchored and iNaturalist-only components separately for every
   species and never quotes the sum alone.

Option 3 is already implemented in the report and is not exclusive of 1 or 2 --
it should be kept regardless of which refit we choose, because the split is
informative in its own right.

### Issue 5b: what this means for the simulation study

**The estimator sweep is structurally blind to Issue 5, by construction.**
`01i_run_estimator_sweep.R:240` generates each replicate with
`simulate_replicate_data(model_code_national_scalar, constants, ...)` -- the
data-generating model *is* the fitted model. The simulated iNaturalist counts
are drawn conditional on the real effort matrix exactly as the likelihood
specifies, so there is no unmodelled iNaturalist effort drift anywhere in the
simulated data. `year_var` is therefore a legitimately identified parameter in
simulation, and the fits recover the total well (median bias 0.002 for
`array_rn`, 0.030 for `camera_rn`).

Consequences, in order of importance:

1. **The sweep's estimator comparison is unaffected.** All four arms faced the
   same correctly-specified truth, so the bias / precision / coverage /
   false-positive results stand as measured. The RN-versus-occupancy finding
   does not depend on Issue 5 in any way.

2. **The sweep's power and coverage are optimistic for real data.** They
   characterise the estimators under conditions kinder than reality: zero
   iNaturalist effort drift. Real-data intervals on `total_var_beta` are
   consequently narrower than warranted, and the real rate of declaring a
   trend when none exists is higher than the 1.1-7.8% measured here, by an
   amount nobody has quantified.

3. **The sweep never scores the split.** Only `tvb_true` is recorded; there is
   no `year_beta_true` or `year_var_true` column, so no replicate was ever
   checked for whether the camera-anchored and iNaturalist-only components are
   individually recovered -- only their sum.

4. **Even in the correctly-specified simulation, `year_var` is the noisier
   component.** Across replicates in the varying scenario:

   | arm | sd `year_beta` | sd `year_var` | sd of sum | corr |
   |---|---|---|---|---|
   | camera_rn | 0.104 | 0.173 | 0.184 | -0.19 |
   | camera_occ | 0.275 | 0.457 | 0.465 | -0.27 |
   | array_rn | 0.150 | 0.256 | 0.236 | -0.42 |
   | array_occ | 0.426 | 0.496 | 0.699 | +0.15 |

   `sd(year_var) > sd(year_beta)` in every arm, and the two are negatively
   correlated (they trade off). So the reported trend's uncertainty is
   dominated by the one component that a single data stream identifies -- and
   the sum is no better determined than `year_var` alone. This is the same
   trade-off seen in the real moose fit, where neither component's CI excludes
   zero but their sum's does.

### Proposed simulation arm (cheaper than the real-data refits)

The question "how much does Issue 5 actually cost us" is answerable in
simulation without touching the real fits, and the sweep infrastructure already
supports it. Generate replicates with `truth$year_var = 0` -- so the entire
trend is camera-anchored and the iNaturalist stream carries no independent
temporal signal -- then inject an unmodelled temporal drift into the simulated
iNaturalist counts (a multiplicative year effect NOT present in the effort
matrix the model conditions on). Score how much of that drift the fit
attributes to `year_var`, and how often it reports a population trend that is
purely injected observer drift.

That is a direct false-trend test aimed at the actual risk, it costs a fraction
of a real-data refit, and its result is what should inform the choice between
Issue 5's options 1, 2 and 3. Recommend running it before committing cluster
time to the structural refits, and sending its result to Arielle and Ben
alongside Issue 5.
