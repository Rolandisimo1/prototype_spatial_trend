# Effort sweep — how much sampling to resolve a per-ecoregion trend

## Question
The abundance sweep varied the ANIMAL (fixed sampling). This sweep varies
SAMPLING at fixed (bobcat) abundance, to answer: how many cameras / how many
iNat observations does a bobcat-density species need for usable per-ecoregion
trend resolution — and does the approach transfer to regions/countries with
different effort?

## Levers
- CAMERAS = controllable effort: number of camera sites (`n_site_keep` in
  build_reduced_constants). Detections per site are NOT set directly — they
  emerge from the occupancy submodel at local intensity, so a rarer region
  yields fewer detections automatically (the point of the exercise).
- iNAT = take-what-you-get: scale the inat_effort data matrix by a global
  factor, and (secondary) by a per-region factor so one ecoregion can be
  "low-effort" like an under-sampled country.

## Design (fixed bobcat abundance: occ_shift=0, count_log_mult=0)
- Camera axis: n_site_keep in {0.5x, 1x, 2x, 4x} of the REAL site count.
  Below real = subsample without replacement. Above real = resample WITH
  replacement (covariates reused + optional jitter) -- flagged as
  extrapolation.
- iNat axis: global effort multiplier in {0.25x, 1x, 4x, 16x}.
- Full cross 4x4 = 16 cells x reps (30, or 20 if queue-limited) = 480 (320) fits.
- Scenario = varying (real regional trend) so per-region RMSE is meaningful;
  optionally add null to re-check false positives under low effort.

## Metric
Fit the ecoregion model; compute_ecoregion_metrics per region:
per-region trend RMSE, coverage, sign recovery. Record realized cam detections
and iNat counts PER REGION so effort -> information -> resolution can be mapped.

## Secondary: per-region effort scenario
Hold cameras+global iNat at 1x; set ONE focal ecoregion's iNat effort to 0.1x
and another to 1x; compare their per-region RMSE. Shows how a low-effort region
(other country) fares vs a well-sampled one at the same abundance.

## Optional phase 2 (feasibility-check first): clustered camera arrays
Snapshot-USA-style clustered arrays instead of independent sites. Changes
spatial independence (clustered cameras cover less independent area per unit).
Claude Code should SCOPE and report feasibility before building, not attempt
blindly.

## Guardrails carried from prior runs
- Per-rep reseed: set.seed(EFFORT_SEED + row_id) AFTER build_reduced_constants,
  BEFORE run_one_replicate (the degeneracy bug fixed twice already).
- Collector anti-degeneracy guard: stop() if per-cell RMSE all identical.
- Never modify Arielle's originals; fork/copy only.

## Deliverable
effort_sweep_summary.csv + figures: per-region RMSE vs #cameras (faceted by
iNat level), per-region RMSE vs iNat counts (faceted by cameras), and an
inverted "effort needed to reach target RMSE" table. Stop before any real-data
refit or fleet action.
