# Estimator simulation — feasibility audit (step 1)

**Question:** can we simulate array-level occupancy, Royle–Nichols, and
N-mixture camera estimators to test their power to detect a trend when
integrated with iNat?

**Answer:** two of three, yes. Classic N-mixture is not simulable from the
data this pipeline retains. Details below.

## Why the simulation framework is the right place to work right now

`fit_replicate()` (sim_helpers.R:414) calls `Cmcmc$run(n_burnin + n_iter)` in
a single call, and `fit_replicate_switch()` (line 739) uses `runMCMC()`.
**Neither uses chunk-resume.** The checkpoint defect that invalidates every
real fit in this project does not touch the simulation results. While real
data is untrustworthy, simulation output is sound — so this is the one place
that can currently produce a defensible answer.

## Estimator feasibility, one at a time

### Array-level occupancy — FEASIBLE

Cameras become replicates within an array-year:

    z_a ~ Bernoulli(psi_a)
    w_{a,i} | z_a ~ Bernoulli(q_{a,i})

`dOcc_v` handles this unchanged — `probOcc = psi_a`,
`probDetect = q[a, 1:n_a]`, `len = n_a`. Only the replicate dimension changes
(10-day windows -> cameras). No new distribution needed.

**Effort caveat, important.** At camera level, effort is handled
*automatically* by the structure: `J[i]` is the number of 10-day windows, so
an 8-week camera contributes ~6 Bernoulli trials against a 3-week camera's ~2,
and `dOcc_v(..., len = J[i])` accounts for it exactly. Aggregating to arrays
collapses the window dimension, so **effort stops being automatic and must be
added back as a per-camera detection covariate** (`log(J[i])` or log
camera-days). Without it, deployment length is confounded with occupancy.

### Royle–Nichols — FEASIBLE, and the natural fit at array level

    N_a ~ Poisson(lambda_a)
    per-camera detection = 1 - (1 - r)^{N_a}
    w_{a,i} ~ Bernoulli(that)

RN is defined on detection/non-detection data, which is exactly what this
pipeline retains — no counts required. And at array level the observation is
"how many of n_a cameras detected," a count out of a known total, which is
RN's native data structure.

**Why it may matter more than a speed fix.** Camera information about the
shared intensity scales as `dpsi/dlambda = 1 - psi`. For abundant species most
cameras in an array detect, `psi_a -> 1`, and the occupancy link goes flat
exactly where you want resolution. RN's expected proportion of detecting
cameras keeps changing as N_a grows, so it does not saturate the same way.
This is the mechanism behind the prediction that array-level occupancy will be
"very high for deer and squirrel."

**Compatibility with the shared intensity field.** The integration rests on
`cloglog(psi) = log(lambda)`, which follows from `psi = 1 - exp(-lambda)` under
a Poisson process of rate lambda. RN's `N_a ~ Poisson(lambda_a)` is the *same*
latent quantity, observed less lossily. So RN is arguably a more faithful
observation model for the existing latent field, not a departure from it.

**Implementation note:** check nimbleEcology for an existing RN distribution
before hand-writing one. It ships occupancy and N-mixture families
(`dOcc_*`, `dNmixture_*`); whether an RN form is included must be verified,
not assumed.

**Known weakness to carry into the writeup:** RN assumes cameras within an
array share detection probability given `N_a`. Real arrays have habitat
variation, and unmodelled detection heterogeneity is RN's documented failure
mode — it inflates abundance estimates. Per-camera covariates (canopy,
road distance, effort) mitigate but do not eliminate this.

### Classic N-mixture — NOT FEASIBLE from retained data

N-mixture needs *counts per visit* (how many individuals seen). This pipeline
does not retain them:

- Real data: camera detections are discretized to **binary 0/1 per 10-day
  window** before the model sees them. Two detections in one window are
  indistinguishable from twenty.
- Simulated data: `simulate_replicate_data()` (line 344) draws `y` from
  `dOcc_v` — binary by construction. Simulating counts would mean simulating
  from a model we are not fitting.

Recovering counts means returning to raw trigger data upstream of this
pipeline. Out of scope here.

**However:** the array-level observation "k of n_a cameras detected" *is* a
binomial count, and modelling it with abundance driving detection **is**
Royle–Nichols. So the useful part of the N-mixture idea is already covered by
RN; there is no separate third estimator to build. Reporting a distinct
"N-mixture" arm would be labelling the same model twice.

## Extension points in sim_helpers.R

| function | line | what to fork / reuse |
|---|---|---|
| `build_reduced_constants` | 112 | reuse as-is; add array assignment downstream of it |
| `simulate_replicate_data` | 344 | fork per estimator — simulates from `dOcc_v`/`dnbinom` via nimble `simulate()` |
| `fit_replicate` | 388 | fork per estimator; single-run MCMC, no resume |
| `compute_replicate_metrics` | 506 | extend, don't replace — add discrimination + power metrics |
| `run_one_replicate` | 570 | fork as the per-row driver |
| `true_param_list` | 251 | reuse; truth stays on the intensity field, unchanged by estimator |
| `scale_truth_abundance` | (sim_helpers_abundance.R) | reuse — the 3-level ladder is the axis RN's advantage should appear along |

Design invariant across all arms: **the iNat side, CAR fields, climate SVCs,
and every trend parameter (`year_beta`, `year_var`, `total_var_beta`, `snr`,
`trend_robust_indicator`) stay byte-identical.** Only the camera observation
block differs. Otherwise a difference in trend recovery cannot be attributed
to the estimator.

## Compute placement (measured, not estimated)

From the completed 180-task abundance sweep: **8 GB per fit**, elapsed
min 24.4 / median 42.3 / max 46.6 minutes, on the reduced design
(`n_cell50_keep = 210`, `n_site_keep = 700`).

Local machine: 36 GiB RAM, 14 cores. **Memory binds before cores** — at 8 GB
per fit, 3–4 concurrent, not 14.

| stage | tasks | laptop (3–4 parallel) | Hazel |
|---|---|---|---|
| pilot: 3 estimators x 1 abundance x 2 scenarios x 3 reps | 18 | ~4 h | queue wait ~= runtime |
| full: 3 x 3 abundance x 2 x 20 reps | 360 | ~65–85 h | ~4 h wall (`%10` throttle) |

**Recommendation: develop and pilot locally, run the full sweep on Hazel.**
The pilot's job is "does it compile, run, and produce sane metrics" — that
needs 3 replicates, not 30 — and iterating locally avoids the SSH/2FA
interruptions that have repeatedly cost this project time.

## Conventions to follow (per the ben-golstein skill)

- `design_df`: one row per (estimator x abundance x scenario x rep).
- Wrap each run in `tryCatch`; a failure returns an "Errored" summary row
  rather than killing the batch.
- Fixed metric set returned as tidy data frames: absolute error, 95% CI
  coverage, in- and out-of-sample RMSE, plus a GOF check.
- Roxygen block (`@name`, `@description`, `@param`) on every function; banner
  header on every file.
- Per-replicate reseed after `build_reduced_constants()` — this project has
  had to fix that same degeneracy bug twice.
- Simulate from the real dataset's covariate and effort patterns rather than
  idealized ones (already how `build_reduced_constants` works).

## Revised arm count

**Three arms, not four:**

1. camera-level occupancy (baseline — the current production estimator)
2. array-level occupancy
3. array-level Royle–Nichols

Classic N-mixture is dropped for lack of count data, and its usable content is
subsumed by arm 3.
