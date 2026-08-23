# Array-level occupancy test — specification

**Goal:** does moving the camera submodel from camera-level to array-level
(a) meaningfully reduce build time, and (b) preserve our ability to detect a
trend? Both halves matter; (a) alone is not a reason to adopt it.

**Status:** test only. Not a decision to switch the production analysis.

## The estimator — read this before coding

NOT "array is occupied if any camera detected." That collapses a 2-of-20
array and a 19-of-20 array to the same value, discards the graded signal, and
pushes psi toward 1 where cloglog(psi)=log(lambda) stops discriminating
abundance (the saturation problem: camera information about lambda scales as
1-psi).

INSTEAD: **cameras are replicates within an array.** For array `a` with
`n_a` cameras, the observation is *how many of those cameras ever detected the
species*, with a per-camera detection probability that accounts for how long
that camera ran.

    z_a ~ Bernoulli(psi_a)                      # array truly used
    w_{a,i} | z_a=1 ~ Bernoulli(q_{a,i})        # camera i in array a detected
    w_{a,i} | z_a=0 = 0

This is exactly the standard occupancy structure with the replicate dimension
changed from 10-day windows to cameras. `dOcc_v` handles it unchanged --
`probOcc = psi_a`, `probDetect = q[a, 1:n_a]`, `len = n_a`.

Interpretation shift to state plainly in any writeup: `q` now absorbs
within-array habitat patchiness as well as detection. `psi_a` means "the array
was used," not "each camera's site was used." That is a real change in what
the parameter means, not just a change of scale.

## Detection covariates for q (the "adjusting for detectability" part)

Per-camera, since cameras are now replicates:

- **effort** -- `log(J[i])` or log camera-days. A camera running 3x longer is
  more likely to have ever detected. THIS IS THE CRITICAL ONE; without it,
  variable deployment length is confounded with occupancy.
- `canopy_height[i]`, `log_roaddist[i]` -- already per-camera in the current
  model, carry over directly.
- `yday` -- currently per camera-window. At array level use the camera's mean
  yday (and its square), or drop it and note the loss. Do not silently pick
  one; report which.

## Node accounting (bobcat, 18-year, nsite=20,531)

| block | camera-level now | array-level |
|---|---|---|
| `p` / `q` | ~66,157 (per camera-window) | 20,531 (per camera) |
| `psi` | 20,531 | n_array (~1,000-2,000) |
| `y` | 20,531 (dOcc_v over windows) | n_array (dOcc_v over cameras) |
| **camera side total** | **~107,200** | **~23,500** |

~78% off the camera side. Combined with Design A's `mu`/`y_inat` work this is
a much larger cut than either alone. Note `q` is still per-camera and remains
the largest remaining block -- folding it into a custom distribution (Design
A option b) would take it to ~0, but that is a separate change and should not
be bundled into this test.

## Phase 1 -- discovery. Report before writing any model code.

Do not assume the array structure. Establish and report:

1. **Is there an explicit array/subproject field** in the camera data? Name
   it, report how many sites have it populated vs missing.
2. **Distribution of cameras per array** -- min / median / max / n_arrays,
   and how many arrays have only 1 camera (those carry no within-array
   replication; decide and report whether they are kept as n_a=1 or dropped).
3. **Are arrays within-year?** Our trend model needs array-year as the unit,
   since the same physical array in 2019 and 2024 must be separate
   observations. Confirm whether the field is array or array-year, and
   construct array-year if needed. (Note: the array grouping in the
   activity-patterns project combined across years -- we cannot.)
4. **For sites with no array field**, report how many, and propose a spatial
   grouping rule (e.g. cluster within X km within a year) WITHOUT implementing
   it yet -- report the proposed radius and resulting array-size distribution
   for review first.
5. **Effort per camera** -- confirm `J[i]` (number of windows) is available
   per camera, and report its spread. This is required for the detection
   covariate.

**Stop here and report.** Steps 2-4 involve judgment calls that change the
science; they should not be made silently.

## Phase 2 -- build one species, both ways

Use **moose** (smallest: ncell50=382, fastest to iterate) for the mechanics,
then **bobcat** for the real comparison (it is the validation species and has
the full-size camera side where the speed question actually matters).

Fork, never modify originals:
- `model_code_array_level.R` -- from `model_code_national_scalar.R`, camera
  block replaced per the estimator above. **The iNat side, the CAR fields,
  the climate SVCs, the trend parameters, and `total_var_beta` all stay
  byte-identical.** The only change is the camera observation model.
- `00_build_array_data.R` -- aggregation, from the phase-1 findings.

## Phase 3 -- the two measurements

Run both scales on the same species, same data, same iterations. **Single-job
runs, no chunk-resume** (the resume defect is unfixed; a chunked run here
would be uninterpretable).

**(a) Build time.** `nimbleModel()` timed alone, both scales. Report the
node count from each so the reduction is measured, not projected.

**(b) Trend detectability -- the part that decides it.** Compare:

| quantity | why |
|---|---|
| `total_var_beta` posterior mean | did the trend estimate move? |
| its 95% CI width | how much precision was lost? |
| `year_beta` and `year_var` separately | which stream lost information? |
| `trend_robust_indicator` | **the decisive one** -- see below |
| per-region trends (ecoregion variant, if run) | did regional resolution survive? |

**Why `trend_robust_indicator` is decisive:** it is `step(snr-1)` where
`snr = year_beta/year_var`, i.e. it fires only when the camera-anchored trend
component matches the iNat-specific one in sign and exceeds it in magnitude.
It is the model's own check that the trend is corroborated by cameras rather
than driven by iNat effort. Moose already sits at 0.478 (national) / 0.432
(ecoregion) -- BELOW the 0.5 threshold. If array aggregation pushes it lower,
the speedup is buying a weaker version of the only guard against reporting an
iNat reporting artifact as a population trend. **A large build-time win with a
degraded indicator is a fail, not a pass.**

## Prior evidence to weigh against

This lab has already compared these two scales, in the Brazil workshop
analysis (`camera_vs_array_slopes.png`, `m2_correlation_camera_vs_array_final.png`):
20 species pairs, camera-level N~1,090 sites vs array-level N=60 arrays. Most
pairs were **significant at camera level only** -- the effect vanished on
aggregation. Only 3 of 20 were significant at both.

Two reasons that is a caution rather than a prediction: it was a different
model (community abundance interactions, not a temporal trend), and it used
plain aggregation at an ~18x reduction, whereas this test uses the graded
proportion estimator at a ~10-20x reduction. But it is direct evidence from
this lab's own data that aggregation can cost real detectability, and it is
the reason (b) is not optional.

## Guardrails

- Fork/copy only; never modify Arielle's originals.
- No chunk-resume in any run here.
- Commit locally; do not push.
- If a phase-1 finding contradicts this spec, stop and say so rather than
  substituting an approach.
