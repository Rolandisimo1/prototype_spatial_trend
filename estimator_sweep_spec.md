# Estimator-comparison sweep — run specification

**Question.** Does moving the camera submodel to array level, or to
Royle–Nichols at array level, change our ability to detect a **trend** when
integrated with iNaturalist data?

**Why now.** Every real fit in this project is compromised by the chunk-resume
defect. The simulation framework is not: `fit_replicate()` uses a single
`Cmcmc$run()` call (sim_helpers.R:414). This is currently the only part of the
project that can produce a defensible answer.

**Status.** Code written and syntax-checked. **Blocked on two prep additions
and one sign-off.** Do not launch the sweep until all three clear.

## The three arms

| arm | camera observation | replicates |
|---|---|---|
| `camera_occ` | baseline, current production estimator | 10-day windows within a camera |
| `array_occ` | array-level occupancy | cameras within an array-year |
| `array_rn` | Royle–Nichols at array level | cameras, with `N_a` driving detection |

Classic **N-mixture was evaluated and dropped**: it needs counts per visit, and
this pipeline binarizes camera data to 0/1 per 10-day window before the model
sees it. Its usable content is already arm 3 — "k of n_a cameras detected" with
abundance driving detection *is* Royle–Nichols. See
`sim_estimator_feasibility.md`.

**Design invariant (verified, not asserted):** the trend block is byte-identical
across all three model files —
`md5(normalized trend block) = 69213000bc266620f68b21c71316358b` for
`model_code_national_scalar.R`, `model_code_array_occ.R`,
`model_code_array_rn.R`. The iNat likelihood, `calcIntensity_SVC`, the CAR
fields, and the climate SVCs are likewise unchanged. Only the camera
observation block differs, which is what makes a difference in trend recovery
attributable to the estimator. **Re-run that check after any edit.**

## Design matrix

3 estimators x 3 abundance levels x 2 scenarios x `N_REP`.

- **abundance** (`scale_truth_abundance`, existing ladder):
  `bobcat_baseline` / `moderate` / `common_deerlike`.
  **Crossed in, not fixed** — the saturation problem RN addresses *is*
  abundance-dependent (camera information about lambda scales as `1 - psi`), so
  a single-abundance sweep would answer the wrong question and could show RN
  offering nothing or everything purely as an artifact of the level chosen.
- **scenario**: `varying` (true trend present) and `null` (no trend).
  Under `varying`, `detect_rate` is **power**. Under `null`, the same column is
  the **false-positive rate**. Reporting power without the null rate is how a
  method that "detects" trends everywhere looks good.

At `N_REP = 20`: 360 tasks. At `N_REP = 3` (pilot): 54 tasks.

## Metrics — two families, because one would mislead

The Brazil workshop comparison (camera-level N~1,090 vs array-level N=60, 20
species pairs) found aggregation **generally reduced power but sometimes
increased AUC**. Those are different questions, so both are recorded:

| family | question | columns |
|---|---|---|
| power | can we detect *that* a trend exists? | `detect_rate`, `bias`, `ci_width`, `coverage` |
| discrimination | can we say *where* intensity is high? | `auc`, `spearman`, `field_rmse` |
| corroboration | is the trend camera-backed or iNat-driven? | `tri_rate` |

**`tri_rate` may be the most important column.** `trend_robust_indicator` is
`step(snr - 1)` with `snr = year_beta/year_var` — it fires only when the
camera-anchored trend component matches the iNat-specific one in sign and
exceeds it in magnitude. Real moose currently sits at **0.478 (national) /
0.432 (ecoregion)**, i.e. *below* the 0.5 threshold: the fitted trend is
iNat-driven, not camera-corroborated. **An estimator that raises this is more
valuable than one that is merely faster**, because this is the model's own
guard against reporting an iNat effort artifact as a population trend.

## Blockers — all three must clear before launch

### 1. `inputs$site_array` (prep addition)

`00_prep_sim_inputs.R` does not carry array labels. Add the real
`subproject_name` field from `umflist.RDS` `siteCovs` (943 groups, 99.9%
populated, species-independent — one master deployment table serves all 36
species), **spatially split per the agreed rule**, as `inputs$site_array`.

The array arms **refuse to run** without it rather than falling back to
`assign_arrays_spatial()`. That spatial approximation exists for local
development only; given this project's history with silent fallbacks, an
explicit stop is worth more than convenience.

### 2. `constants$site_keep` (helper addition)

`build_reduced_constants()` computes the retained-site index internally but
does not return it. It must, so `site_array` can be aligned to the reduced
design. Additive change; nothing else reads it.

### 3. Array construction sign-off — PI review required

`subproject_name x year` is **necessary but not sufficient**:

- **Time:** 255 of 943 subprojects (27%) span multiple years, up to 7 — hence
  the `x year`, giving 1,483 array-years. Independently corroborated:
  `deployment_id` has 246 duplicated rows and all span multiple years, so the
  table is already at camera-year granularity.
- **Space, still unfixed:** the largest array-year has **163 cameras spanning
  162 km**. That is an administrative aggregate, not an array. Averaging
  occupancy covariates across 162 km is meaningless.

Within-array-year nearest-neighbour spacing is median 0.28 km / 95th pct
2.42 km, so genuine arrays are compact and oversized units should be split
spatially. Report, **for review before finalizing**, at candidate radii
(2.5 / 5 / 10 km):

- resulting unit count; cameras per unit min/median/mean/max
- maximum spatial extent (km) per unit; how many exceed ~10 km
- how many original array-years were split, and into how many pieces
- graded-proportion and informative-unit counts recomputed **post-split and
  post-range-mask**, for **both moose and bobcat**

This gate is deliberate: the PI flagged that array construction went wrong in
the activity-patterns analysis, and does not want units that are too large or
that overlap too much time.

## Settled: 5 km radius, and NO array random effect

**Radius = 5 km, by protocol not by optimization.** Snapshot USA deploys
cameras on a nominal 5 km spacing, so 5 km is the design's own unit and the
split should reproduce it rather than maximize a statistic. A sample-size
argument favoured 2.5 km (bobcat 1,429 vs 1,151 graded informative units;
moose 129 vs 106), but matching the sampling protocol is the stronger
consideration: the units then correspond to what the field crews actually
laid out. Retention at 5 km is 97.0% of cameras, 4,002 array-years,
median 4 cameras per unit, max extent 5 km guaranteed by complete linkage.

**No array random effect in the array-level arms.** Once the analysis unit
IS the array, cameras within it are replicates of that array's occupancy
state, not separate sites that happen to be near each other. The clustering
a random effect would absorb is handled structurally by the aggregation
itself. Adding `array_effect[a] ~ dnorm(0, sigma_array)` on top of a
one-row-per-array-year likelihood would be modelling a grouping that no
longer has anything nested inside it.

Two supporting facts, both measured:

- **The CAR field already carries place.** 89.5% of 5 km arrays nest inside
  a single cell100, and 70.3% of revisited arrays have all their year-visits
  in one cell100 — so two array-years at the same location already share
  `link_occ_intercept[cell]`. Location correlation is modelled; it is not
  missing.
- **There is almost no within-array trend information to exploit.**
  Within-array year variance is 0.49 against a total of 14.13, i.e. **3.5%**.
  Median revisit span is 2 years, only 15 arrays span >= 5 years, one spans
  >= 10. A model leaning on within-array year contrasts would lose precision,
  not gain it. The revisited arrays are valuable because they remove
  site-composition differences from the comparison, which the CAR field
  already delivers — not because they supply long per-array time series.

### Consequence for reading the sweep — important

The **camera-level baseline** arm treats cameras as independent sites and has
no array term, so it carries a pseudo-replication problem the array arms do
not. Variance-weighted mean cluster size at 5 km is 16.5 cameras; at a
plausible ICC of 0.15 the implied design effect is 3.33, i.e. genuine
intervals on `year_beta` about **82% wider** than the camera-level model
reports.

Therefore: **if the array arms return wider credible intervals than the
camera-level arm, that is NOT automatically evidence that aggregation cost
precision.** It may be the array arms reporting honest uncertainty while the
camera-level arm is overconfident. The two are not on the same footing for CI
width, and the write-up must say so. Coverage — not width — is the metric that
separates these cases, which is why `tvb_covered` is recorded per replicate:
an arm whose 95% interval covers truth ~95% of the time is calibrated, and a
narrower interval with 70% coverage is simply wrong.

This makes calibration a first-class result of the sweep rather than a
footnote, and it is a point in favour of array-level aggregation that has
nothing to do with build time.

## Agreed construction decisions

| decision | choice | rationale |
|---|---|---|
| singletons (`n_a = 1`) | **drop** | one Bernoulli trial cannot separate `psi_a` from `q`; uninformative under occupancy, degenerate under RN; 0.2% of cameras |
| `yday` | array mean **+ square** | detection is seasonal; cheap. Within-array seasonal spread is lost — state this |
| 29 unassigned sites | **drop** | 0.12%, one project, one year |
| effort | **required** `q_beta[5] * log_effort` | see below |

**Effort is not optional.** At camera level it enters *automatically*: `J[i]`
is the window count, so an 8-week camera contributes ~6 Bernoulli trials
against a 3-week camera's ~2, and `dOcc_v(..., len = J[i])` accounts for it
exactly. Aggregating to arrays **collapses the window dimension**, so effort
must re-enter as a per-camera detection covariate. Measured `J` spread is
1–10 occasions (median 3, mean 3.22) — a 10x range. Omitting it confounds
deployment length with occupancy.

## Species

- **moose** — mechanics only. Data-poor: ~9 informative array-years per year
  to estimate `year_beta`. **A moose failure is a sample-size result, not a
  verdict on the method.** Both the PI and Claude Code reached this
  independently; it is written into the design so a moose result cannot stand
  as the headline.
- **bobcat** — **the deciding arm.** Validation species, 762 informative
  array-years, and the abundance level the ladder is anchored to.

## Guardrails

- **Per-replicate reseed** — `set.seed(DESIGN_SEED + row_id)` **after**
  `build_reduced_constants()` (which calls `set.seed()` internally at
  sim_helpers.R:97) and **before** the fit. This project has shipped the
  missing-reseed bug **twice** (job 436723; then again in
  `01e_run_abundance_sweep.R`), each time producing an n=1 result that looked
  like n=30. Do not move this line.
- **Degeneracy guard** — `summarize_estimator_sweep()` stops if all
  replicates in a cell share an identical `tvb_mean`. Tested in both
  directions (fires on identical, passes on varied). The collector must call
  it; do not summarize around it.
- **No chunk-resume.** Single `Cmcmc$run()` only. The project's
  checkpoint/restore is known broken (restarts from inits at every boundary),
  so a chunked run here would be uninterpretable.
- **`tryCatch` per row** returning an `Errored` row rather than killing the
  batch (ben-golstein convention).
- **Fork/copy only** — never modify Arielle's originals.
- **Commit locally; do not push.** The PI pushes manually.

## Compute placement (measured)

From the completed 180-task abundance sweep: **8 GB per fit**, elapsed
min 24.4 / median 42.3 / max 46.6 min, on the reduced design
(`n_cell50_keep = 210`, `n_site_keep = 700`).

| stage | tasks | laptop (36 GiB -> 3–4 concurrent) | Hazel |
|---|---|---|---|
| pilot, `N_REP = 3` | 54 | ~10 h | ~1 h at `%10` |
| full, `N_REP = 20` | 360 | ~65–85 h | ~4 h at `%10` |

**Memory binds before cores** — 8 GB per fit means 3–4 concurrent on the
laptop, not 14. Recommendation: **pilot locally** (iterating without SSH/2FA
interruptions is worth real time), **full sweep on Hazel**. Note the pilot
still needs `sim_inputs.RDS` fetched from Hazel once.

RN caveat for run sizing: `N_a` is a latent discrete parameter, so NIMBLE
assigns it a slice/RW sampler and the RN arm may mix more slowly than the
occupancy arms. Check the `N_a` posterior is not piling up at `N_max` — that
would mean the cap is binding and must be raised.

## Known limitation to state in any writeup

These are simulations from a **correctly specified** model: the data are
generated by the same structure being fitted. Real data will be worse. In
particular the RN arm's known weakness — unmodelled within-array detection
heterogeneity inflating abundance — is *absent by construction* here, so the
RN arm's simulated performance is optimistic relative to real data.
