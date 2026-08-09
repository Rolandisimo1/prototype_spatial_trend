# Moose model results — full portfolio report

**Model:** `moose_v1fix_national_scalar` and `moose_v1fix_ecoregion`, converged 2026-08-05.
**Fixes applied:** iNat effort-matrix sorting fix (Fix 1), 10-covariate KEEP list
(matches original reviewed design), range mask as currently configured (see
**important caveat** below — this is a real, unresolved data gap, not just a
map-rendering choice).
**Status:** Sections 1–7 are built from data already on disk. Three more
items — theta0/theta1's posterior, occ_beta's full posterior, and overall
model-fit metrics (WAIC, occupancy AUC) — need raw chain `.RDS` files on
Hazel; extraction scripts are ready but this session's Hazel SSH target was
unreachable ("Permission denied"). **The range mask problem (see below) is
confirmed across all three species' current v1fix fits, not moose-specific**
— a real fix has been drafted as code (not yet run) and should be applied
before any further Hazel launches, including the 12 windowed-refit bundles
built 2026-08-09.

**Version-control note:** a GitHub push for this pipeline was requested
2026-07-12 but never followed through — the next work session pivoted
straight into the Fix-1-only real-fit launch decision and no `git` repo was
ever initialized for `prototype_spatial_trend/`. This report and its
companion scripts are being prepared now for that initialization.

This report is intended as the **template** for the equivalent bobcat and
white-tailed deer reports once their fits converge — section structure and
extraction script should carry over with only the species name changed.

---

## ⚠ Important caveat discovered this session: the range mask problem is bigger than previously stated

Earlier sessions established that the IUCN range polygon (`moose_range.geojson`)
is unreliable and should not be used to CLIP MAPS — the model's own 382-cell
fitted footprint was treated as the safe, mask-free alternative. **That
footprint is not actually mask-free.** Checking it directly against
`moose_cell_counts_unmasked.csv` (real per-cell iNat counts with no filtering)
shows:

- **Every one of the model's 382 fitted cells has `in_range == TRUE` under the
  old IUCN mask, and zero fitted cells have `in_range == FALSE`.** The mask was
  applied at **data-prep time**, before fitting — not just at plotting time.
- Of 427 grid cells nationwide with real, non-zero iNat moose records, **149
  were excluded from the model fit entirely** — the model was never shown
  this data and cannot produce `mu`/`psi` for these locations.
- These 149 excluded-but-documented cells break down by state as: **Colorado
  (43), Wyoming (25), North Dakota (20), Montana (18), Utah (8)**, plus smaller
  counts in NV, WA, SD, CT, MA, ID, NM, MI, MN, WI, NE.
- This directly explains the Colorado/Utah gap flagged this session: the
  documented Colorado moose range expansion (and a real Utah population) are
  both concentrated in cells the model literally never fit on.

**Consequence:** every map in this report showing "the model's real
382-cell data footprint" is still showing the outcome of the flawed mask, one
step removed — it accurately reflects what the *converged fit* covers, but
that coverage itself was set by the same bad polygon. All spatial figures
below now overlay the 149 excluded-but-documented cells explicitly so this
gap is visible rather than silently absent. **A re-fit including these cells
is required before Colorado/Utah moose intensity, or a Colorado-inclusive
national/regional trend, can be reported** — this is now the top blocking
item for a fully trustworthy moose portfolio, ahead of the windowed 5yr/10yr
refits.

[moose_excluded_cells_with_inat_activity.csv]({{artifact:0ab80c19-bdfb-4f9b-b45e-30ac6d8fc33b}})
— the 149-cell list with coordinates, state, and raw iNat count.

**Scope confirmed beyond moose (2026-08-09):** the launch record for the
current v1fix real fits shows all 6 jobs (bobcat + WTD + moose, both models
each — job IDs 499944–499949, launched ~2026-08-03) went out with **Fix 1
only**; Fix 2 (this range mask) was explicitly deferred pending team review
and never applied. So bobcat's and WTD's currently-running/converged v1fix
fits carry the identical gap — this is a project-wide data-prep issue, not
a moose peculiarity. The 12 windowed-refit bundles (5yr/10yr × 3 species ×
2 models, built 2026-08-09 from the same 18-year source, not yet launched)
inherit it too.

**Fix implemented as real code (not yet run):**
[prep_inat_data_grid_v2.R]({{artifact:94fdb205-b92f-46a8-867a-4b50062c96d3}})
forks `prep_inat_data_grid()` from `integration_helper.R` (Arielle's
original is untouched) and replaces the range-polygon mask with (1) a
captive-record exclusion filter and (2) a spatial-isolation outlier filter
— excluding an observation only if it has zero same-species records within
100 km and 3 years, rather than by historical range boundary. The 100km/3yr
thresholds are the team memo's starting proposal, not yet independently
reviewed. **Recommended next step: implement this fix, re-fork all six
national_scalar/ecoregion input bundles, re-verify the 10-covariate KEEP
list, and re-launch — before running the 12 windowed jobs or trusting any
further bobcat/WTD convergence output.**

---

## 1. Convergence diagnostics

| Model | Status | N params monitored | Max R-hat | N over 1.1 | Mean R-hat | Chains × iters |
|---|---|---|---|---|---|---|
| national_scalar | CONVERGED (2026-08-05) | 2751 | 1.052 | 0 | 1.004 | 3 × 50,000 |
| ecoregion | Stopped at cap, 1 param over | 2760 | 1.1013 | 1 (`link_occ_intercept[7]`) | 1.010 | 3 × 50,000 |

**Ecoregion caveat, accepted as converged:** the single miss is 0.0013 over
threshold on one CAR-field spatial-intercept node, out of 2,759 parameters —
every trend, detection, and iNat parameter cleared 1.1. Treated as a marginal,
practically-converged result rather than a real non-convergence signal.

### Per-parameter-family R-hat breakdown

![Moose R-hat by family]({{artifact:40784d2d-eba1-4bc0-a82a-db73de01c112}})

| Model | Family | N params | Max R-hat | Mean R-hat | N over 1.1 |
|---|---|---|---|---|---|
| national_scalar | CAR_fields | 2724 | 1.042 | 1.004 | 0 |
| national_scalar | other | 15 | 1.052 | 1.013 | 0 |
| national_scalar | iNat | 3 | 1.021 | 1.013 | 0 |
| national_scalar | trend | 4 | 1.006 | 1.003 | 0 |
| national_scalar | detection | 5 | 1.003 | 1.002 | 0 |
| ecoregion | **CAR_fields** | 2724 | **1.101** | 1.010 | **1** (`link_occ_intercept[7]`) |
| ecoregion | other | 15 | 1.016 | 1.006 | 0 |
| ecoregion | iNat | 3 | 1.015 | 1.009 | 0 |
| ecoregion | ecoregion | 9 | 1.004 | 1.001 | 0 |
| ecoregion | trend | 4 | 1.004 | 1.003 | 0 |
| ecoregion | detection | 5 | 1.002 | 1.001 | 0 |

The ecoregion model's single over-threshold parameter is isolated to the
CAR-field family (spatial intercepts on the cell100 grid) — every
scientifically load-bearing family (trend, iNat, detection, ecoregion)
clears 1.1 comfortably in both models. This strengthens the accept-as-
converged call: the marginal miss sits in a nuisance spatial-smoothing term,
not in any parameter this report draws conclusions from.

---

## 2. National trend

| parameter | mean | median | sd | 95% CI | P(<0) | R-hat |
|---|---|---|---|---|---|---|
| total_var_beta | -0.299 | -0.300 | 0.085 | [-0.463, -0.126] | 0.999 | 1.006 |
| year_beta (camera-anchored) | -0.242 | -0.245 | 0.116 | [-0.464, -0.012] | 0.982 | 1.003 |
| year_var (iNat-specific) | -0.057 | -0.046 | 0.124 | [-0.314, 0.181] | 0.678 | 1.002 |
| trend_robust_indicator | 0.478 | 0 | 0.500 | [0, 1] | — | 1.000 |

![Moose trend results]({{artifact:2732181d-bfe9-4730-a586-404dc7807278}})

**Reading this:** `total_var_beta` is the trend actually used in the intensity
equation (`year_beta + year_var`); it's the number to report as "the national
trend." Both models agree closely (-0.299 national-scalar vs. -0.287 ecoregion
global). `year_beta` alone (camera-anchored) is directionally consistent but
noisier; `year_var`'s 95% CI spans zero.

**Caveat — data density:** camera trap data for moose exists only 2019–2025;
2008–2018 in this model is carried by iNat data alone. The "18-year trend"
label is not uniformly supported across its own window. **A genuine 5-year
(2021–2025) and 10-year (2016–2025) refit is in progress separately** (see
`real_fit_wtd_moose_plan.md` and the windowed-refit thread) — those results
will supersede the 18-year headline number once available.

---

## 2b. Overall model-fit metrics (requested addition — AUC/WAIC-style)

R-hat (Section 1) tells you the chains **converged** — it says nothing about
whether the converged model actually **fits the data well**. Two standard
metrics fill that gap, and **neither could be computed this session**:

- **WAIC** — nimble was configured with `enableWAIC = TRUE` at build time
  (`HPC_run_model_chunks_chain1.R` L220, inherited unchanged by all moose
  chain scripts), but no step in the checkpointed chunking workflow ever
  called `Cmcmc$getWAIC()` or saved a WAIC value — the `chain_*.RDS` files
  contain `obj$samples` only. Computing it now requires `calculateWAIC()`
  against the live compiled nimble model plus the full pooled posterior,
  which needs Hazel (cannot be reconstructed from CSVs off-cluster).
- **Occupancy AUC** (camera side only — how well fitted `psi` separates
  sites with at least one detection from sites with none) requires the raw
  `y[i,j]` detection-history matrix from `input_data_moose_*.RDS`, which was
  never extracted locally either.

[extract_moose_v1fix_fit_metrics.R]({{artifact:4aec188a-ad46-4d6a-8c67-d442c0611344}})
— ready-to-run on Hazel: computes WAIC via `nimble::calculateWAIC()` for
both models, and camera-side occupancy AUC via `pROC::auc()` (posterior
mean + 95% CI across a draw subsample). No iNat-side AUC equivalent exists
— `y_inat` is a count, not a binary detection, so a posterior-predictive
count check would be the analogous iNat-side diagnostic; not attempted here.

---

## 3. Regional trend (ecoregion model)

| Region | Deviation (year_region) | 95% CI | P(<0) | Absolute trend | 95% CI | Camera sites | iNat cells |
|---|---|---|---|---|---|---|---|
| Eastern Temperate Forests | -0.058 | [-0.416, 0.168] | 0.654 | -0.345 | [-0.696, -0.078] | 78 | 42 |
| Great Plains | 0.072 | [-0.162, 0.482] | 0.332 | -0.216 | [-0.514, 0.244] | 0 | 19 |
| Marine West Coast Forest | 0.001 | [-0.380, 0.387] | 0.498 | -0.286 | [-0.705, 0.172] | 0 | 0 |
| Mediterranean California | 0.000 | [-0.383, 0.385] | 0.500 | -0.287 | [-0.714, 0.176] | 0 | 0 |
| North American Deserts | 0.048 | [-0.191, 0.411] | 0.379 | -0.239 | [-0.522, 0.137] | 43 | 38 |
| Northern Forests | -0.039 | [-0.298, 0.150] | 0.646 | -0.326 | [-0.545, -0.132] | 925 | 118 |
| Northwestern Forested Mtns | -0.010 | [-0.234, 0.199] | 0.537 | -0.297 | [-0.494, -0.101] | 602 | 141 |
| Water | -0.007 | [-0.308, 0.281] | 0.522 | -0.295 | [-0.618, 0.038] | 6 | 24 |

sigma_region posterior mean = 0.136 (95% CI [0.007, 0.465]).

**Reading this:** no single region's deviation CI excludes zero — the
regional structure is directionally informative (see the model-vs-agency
comparison below) but individually uncertain. Marine West Coast Forest and
Mediterranean California have zero camera and zero iNat cells — their values
are essentially the prior/spatial-smoothing default, not real information.

---

## 4. Model vs. state agency comparison

### 4a. Raw agency data, state by state (restored — this is the reference view)

![Moose agency trend by state, raw]({{artifact:124ecd03-a70b-475f-abaf-56615733130f}})

This is the original per-state map from `deer_moose_trends_master.csv`,
restored exactly as first built (2026-07-20) — one color per **state**
(not ecoregion), directly from each state agency's own reported 10-year
direction, with **no model overlay, no ecoregion aggregation, and no
composite/derived score**. This is the ground truth the model is being
checked against, and it is the primary reference for that comparison —
the ecoregion-level ecoregion comparison below is a secondary, model-side
re-aggregation of the same underlying facts, not a replacement for this
view.

### 4b. Ecoregion-level categorical comparison (model-side aggregation)

![Moose model vs agency, rebuilt]({{artifact:de26095a-d8e5-46c5-84a5-f8b2c545e4b0}})

Rebuilt per the 2026-08-08 decision: the old `agency_score` continuous number
(categorical agency labels area-weighted into a single figure — an
undisclosed data translation) is **dropped entirely**. Comparison is now
strictly categorical: the model's regional posterior is thresholded into
`decreasing` (P(decline) ≥ 0.90) / `increasing` (P(decline) ≤ 0.10) /
`uncertain-stable` (in between), and compared region-by-region against the
agency's own plain `dominant_dir` label (`increasing`/`stable`/`decreasing`/
`mixed`) — no correlation coefficient, no ordering statistic. This view
exists only because the model itself reports at ecoregion, not state,
resolution — the state-by-state map in 4a above is the one to read first.

| Region | Model direction | Agency direction | Sign agreement |
|---|---|---|---|
| Eastern Temperate Forests | decreasing | stable | **disagree** |
| Great Plains | uncertain/stable | stable | agree |
| Marine West Coast Forest | no data (prior only) | stable | not comparable (0 camera, 0 iNat cells) |
| Mediterranean California | no data (prior only) | no agency data | not comparable |
| North American Deserts | decreasing | increasing | **disagree** |
| Northern Forests | decreasing | stable | **disagree** |
| Northwestern Forested Mountains | decreasing | stable | **disagree** |
| Water | decreasing | mixed | not comparable |

[moose_agency_sign_agreement_table.csv]({{artifact:27cd664f-b0b7-4d5f-8a66-3b4c3c63a925}})
— full table with posterior P(decline), n camera sites, n iNat cells per region.

**Reading this — and the range-mask finding changes the interpretation:**
Of the 8 regions in the table, 3 are marked not comparable (zero/mismatched
data on one side), leaving 5 comparable regions — of those, 4 disagree in
sign and only 1 (Great Plains) agrees. That is a starker divergence than the
"5 of 6 agree" read from the earlier (flawed, ordering-only) scatter
comparison. But the
North American Deserts disagreement is now explainable rather than a genuine
model failure: **43 of the 149 cells with real, undocumented-by-the-model
iNat moose activity are in Colorado**, which sits in this exact ecoregion —
the model's "decreasing" call for North American Deserts is based on the
cells it *did* fit (mostly Utah/Nevada margins), while the agency's
"increasing" call is driven specifically by the Colorado reintroduction
population the model never saw. This is not evidence the model is wrong
about the cells it covers — it is evidence the model's covered footprint
and the agency's basis for comparison are answering about different
geography. The other three disagreements (Eastern Temperate Forests,
Northern Forests, Northwestern Forested Mountains — all "decreasing" vs.
agency "stable") are not explained by the excluded-cell gap and stand as
real, unresolved model-vs-agency divergence.

---

## 5. Spatial intensity over time — snapshot maps

![Moose intensity snapshots]({{artifact:06d506c7-e2a6-4dde-ba50-4f707ce5313c}})

`mu[g,t]` was not directly monitored in either model; reconstructed exactly
(not approximated) from monitored parameters (`theta0`, `theta1`,
`total_var_beta` + covariate data) for snapshot years 2008, 2013, 2019, 2025
(t = 1, 6, 12, 19). Shown at the model's fitted 382-cell footprint, zoomed to
its own bounding box — no IUCN range polygon applied to the map itself.

**Red-outlined cells mark the 149 real, iNat-documented cells the model was
never fit on** (see the range-mask caveat above) — most visibly a dense
cluster in Colorado/Utah that is entirely absent from the colored (fitted)
cells. The 2008→2025 decline visible in the colored cells is national-scalar
and spatially near-uniform by construction (consistent with
`total_var_beta = -0.299`); it says nothing about the excluded cells, which
have no model output at any year.

A dedicated 2008→2025 percent-change map was built and then dropped this
session per user direction — it spans a period moose camera data does not
cover (cameras start 2019), which contradicts the report's own data-density
caveat (Section 8, item 1). The 5-yr (2021–2025) and 10-yr (2016–2025)
windowed refits, once they land from the parallel Hazel thread, will provide
the change maps that respect this caveat.

---

## 6. Climate SVC response maps (MWMT / MCMT)

![Moose climate SVC maps]({{artifact:5ca090d8-fca7-43c7-af97-68e8cf287842}})

`MWMT_effect[1:908]` and `MCMT_effect[1:908]` posterior means, cell100 grid
(100 km), both models — no range mask, zoomed to the data footprint's
bounding box. Black cell outlines mark where the 95% CI excludes zero.

**MWMT (max summer temp):** only 9/908 cells have a CI excluding zero in
either model — the climate-heat-stress signal is directionally present
(the strongest positive/blue cluster sits in the northern Rockies, and the
southern-edge cells trend more negative) but not sharp enough at most
individual cells to call significant. This is a weaker signal than the
sanity-check expectation of a clear southern-edge heat-stress signature —
worth flagging as a modest, not strong, climate effect.

**MCMT (min winter temp):** a stronger, more spatially coherent signal —
31–33/908 cells (national-scalar/ecoregion) have CIs excluding zero, mostly
in the Rockies/Northwestern Forested Mountains and clustered around the
Great Lakes/Northeast. Consistently negative in the significant cells,
meaning colder winters associate with *higher* CAR-field values there —
plausible given moose's cold-adapted physiology and heat intolerance (i.e.
milder winters, not just hot summers, may be the more informative climate
covariate for this species).

Green diamonds mark the same 149 excluded-but-documented cells from the
range-mask caveat above. The CAR fields are fit on the cell100 grid, one
level coarser than the cell50 grid the mask cut at, so climate response
*surfaces* still exist under those markers (interpolated from neighboring
cells) — but the local intensity/occupancy signal at those specific
locations was never incorporated, so treat the climate response there as
extrapolated, not locally informed.

---

## 6b. Other occupancy covariate effects (requested addition — beyond MWMT/MCMT)

MWMT and MCMT (Section 6) are the two **spatially-varying** climate
coefficients — the only ones with their own CAR field per cell. The model
also carries **10 ordinary (non-spatial, single global coefficient each)
occupancy covariates** in `occ_beta[1:10]`, decided in `make_reduced_input.R`
(dropping 7 collinear/NaN-prone candidates — CWD, NDVI_sd, Impervious, PDSI,
PPT, elevation, Temp — while keeping the SVC climate effects on):

**Human_pop, NDVI_mean, Ag, Deciduous, Evergreen, Mixed, terrain_ruggedness,
soil_clay, soil_silt, soil_sand**

These enter the camera occupancy submodel as
`inprod(occ_beta[1:10], occ_covars[i, 1:10])` and the iNat intensity
pathway identically (same global `beta`, per `calcIntensity_SVC`) — one
coefficient per covariate, shared across all cells and years (unlike
MWMT/MCMT, which get a full per-cell CAR field each).

**Not yet extracted for moose.** `occ_beta` is a monitored node (confirmed
present in the chain files, since its `[1]` element shows up as the
national-scalar model's `worst_param` in the R-hat-by-family table — R-hat
1.052, well converged) but no per-covariate posterior summary was ever
pulled, only that single R-hat number. The bobcat project ran this exact
extraction already (`bobcat_env_beta_forest.png`, a 10-covariate forest
plot with posterior means and 95% CIs) — moose has no equivalent yet.

[extract_moose_v1fix_occbeta.R]({{artifact:3f9b5120-7177-499b-bab2-b290cda6872c}})
— ready-to-run on Hazel, matching the bobcat forest-plot pipeline, both
models. Includes a built-in check that `occ_covars` column order actually
matches the 10-name list above before labeling anything — it fails loudly
rather than silently mislabeling if the order doesn't match.

---

## 7. Camera vs. iNat congruence

### 7a. Temporal (parameter-level)

![Moose temporal congruence]({{artifact:46d23e11-1ec0-4cac-9e48-6f71a92c1fe7}})

`trend_robust_indicator = step(snr-1)`, `snr = year_beta/year_var`, evaluated
per posterior draw — 1 only when `year_beta` and `year_var` share sign AND
`|year_beta| >= |year_var|`. Posterior mean 0.478 (national-scalar) / 0.432
(ecoregion) — below 0.5 in both models, meaning the camera-anchored signal
does not dominate the iNat-specific signal by the model's own criterion. Not
camera-corroborated.

**[STILL PENDING HAZEL]** `theta0`/`theta1` full posterior summary
(mean/median/CI/R-hat) is not yet pulled — this session's Hazel SSH access
was unreachable. Per Goldstein et al. (bioRxiv 2025.01.17.633640), `theta1`
is the paper's own primary camera/iNat congruence metric (they report
inspecting the posterior estimates of theta_1 to characterize how often the
two datasets corresponded) —
our model's `theta1` plays the identical structural role
(`log(mu) = theta0 + theta1*log(sum(lambda))`). `theta1`'s R-hat is already
available by family (1.021 national-scalar, 1.015 ecoregion — both well
converged, see Section 1) but its posterior mean/median/CI is not. The
extraction script (`extract_moose_v1fix_theta.R`, ready to run) pulls both
nodes directly from the already-loaded samples matrices — no reconstruction
needed, just a Hazel connection.

[extract_moose_v1fix_theta.R]({{artifact:753619bd-cf2f-4b10-8da4-1f617dd6257b}})

### 7b. Spatial (per-cell/site)

![Moose spatial congruence]({{artifact:ed56d674-7816-42c0-bcfa-acb932044d3f}})

Shown as adjacent panels rather than a forced bivariate map: `mu` is a
cell50-level intensity (iNat-informed), `psi` is per-site occupancy
probability (camera-informed) — different spatial units, not directly
overlaid. Camera data for moose exists only 2019 and 2025 (179 and 231
sites respectively; zero sites in 2008/2013, not shown or approximated).

Both panels are zoomed to the model's fitted footprint and show the same
three broad clusters (Rockies/PNW, upper Midwest, Northeast). `psi`'s
per-site granularity reveals sub-cluster heterogeneity (some sites near-zero
occupancy adjacent to high-occupancy neighbors) that the coarser cell50 `mu`
field smooths over — consistent with occupancy saturating at high
underlying intensity (the `cloglog(psi) = log(lambda)` link) while `mu`
continues to track continuous variation via the iNat count model. This
qualitative visual congruence is the spatial complement to the
`trend_robust_indicator` reading in 7a; a full quantitative site-to-cell
join (e.g., aggregating `psi` to cell50 by mean and computing a rank
correlation against `mu`) was not attempted this session and would be a
reasonable small follow-up.

---

## 8. Known open issues / caveats carried into this report

1. **18-year trend spans a camera-data-sparse period** (moose cameras only
   2019–2025) — windowed 5yr/10yr refits in progress separately on Hazel
   (bobcat-scale checkpoint test was the last reported status).
2. **trend_robust_indicator < 0.5 in both models** — trend is iNat-driven, not
   camera-corroborated (this is the SAME caveat flagged for bobcat in the
   original provenance memo — a project-wide pattern, not moose-specific).
3. **Range mask excludes 149 real, iNat-documented cells from the model FIT
   itself, not just from map rendering** (Colorado 43, Wyoming 25, North
   Dakota 20, Montana 18, Utah 8, plus smaller counts elsewhere) — see the
   top-of-report caveat. **Confirmed 2026-08-09 to affect bobcat and WTD's
   current v1fix fits identically** (launch record: all 6 jobs, IDs
   499944–499949, went out with Fix 1 only; Fix 2 was deferred) — this is a
   project-wide data-prep issue, not moose-specific, and the 12 windowed-
   refit bundles built the same day inherit it too. It directly explains 1 of
   4 "disagree" calls in the rebuilt agency comparison (North American
   Deserts) and means Colorado/Utah moose intensity cannot be reported at all
   from the current fit. **The fix is now implemented as real, ready-to-run
   code** — [prep_inat_data_grid_v2.R]({{artifact:94fdb205-b92f-46a8-867a-4b50062c96d3}})
   forks `prep_inat_data_grid()`, dropping the IUCN polygon mask in favor of
   a captive-record filter plus a spatial-isolation outlier filter (100km/
   3yr, not yet independently reviewed) — but it has not yet been run, and no
   bundle has been re-forked or re-fit against it. **Recommended: apply this
   before launching the 12 windowed jobs or trusting further bobcat/WTD
   convergence output.**
4. **Regional deviations individually uncertain** — no single region's CI
   excludes zero; only the aggregate ordering (vs. agency direction) is
   informative.
5. **Agency comparison, categorical version, shows more disagreement than the
   earlier flawed ordering-only comparison suggested** — of 5 comparable
   regions, 4 disagree in sign (only 1 of those 4 is explained by the excluded-
   cell gap above). Eastern Temperate Forests, Northern Forests, and
   Northwestern Forested Mountains all show the model calling "decreasing"
   against an agency "stable" call, with no known data-gap explanation —
   flagged as real, unresolved divergence rather than resolved by this
   session's work.
6. **theta0/theta1 full posterior summary still pending Hazel** — R-hat is
   available (both converged) but mean/median/CI is not; extraction script
   is ready (`extract_moose_v1fix_theta.R`).
7. **Spatial mu-vs-psi congruence is qualitative only** (adjacent panels,
   visual read) — a quantitative site-to-cell join was not attempted this
   session.
8. **No overall model-fit metric (WAIC, occupancy AUC) yet** — R-hat confirms
   convergence, not fit quality. Both require the live nimble model / raw
   detection-history matrix on Hazel; extraction script is ready
   (`extract_moose_v1fix_fit_metrics.R`).
9. **Full occ_beta posterior (10 covariates) not yet extracted** — only its
   R-hat summary exists. Extraction script ready
   (`extract_moose_v1fix_occbeta.R`), matching the existing bobcat forest-plot
   pipeline.
10. **No version control on this pipeline** — a GitHub push was requested
    2026-07-12 and never completed; no git repo exists for
    `prototype_spatial_trend/` as of this report. All code changes since then
    (Fix 1, the v1fix relaunch, the windowed-refit builder, this session's
    Fix 2 draft) have no commit history. Recommended: initialize the repo and
    commit going forward, starting with this report and its companion
    scripts.

---

## Reusable template notes (for bobcat / WTD / future species)

- Sections 1–2 (convergence, national trend) need only the model's own
  extracted posterior CSVs — no new Hazel work beyond what any converged fit
  already produces via the standard convergence-check pipeline.
- Section 3 (regional) only applies to species with an ecoregion fit.
- Sections 5–7b require the `mu`/`psi` reconstruction extraction — **for
  future species, add `mu` and `psi` to the monitor list before launching**
  (confirmed cheap to do for the not-yet-launched windowed refits) so this
  reconstruction workaround is not needed again.
- Section 4 (agency comparison) only applies to species with compiled agency
  data (`deer_moose_trends_master.csv` currently covers WTD and moose only).
- cell100/cell50 coordinate crosswalks and range polygons are reusable across
  species once built — check `cell50_coordinates.csv` and species-specific
  range GeoJSONs before re-fetching.
- **New lesson from this report, check for every future species before
  calling any spatial figure "unmasked":** don't assume the model's fitted
  cell footprint equals "no mask" just because no range polygon is applied
  at *plot* time. Cross-check the fitted cell list against an unfiltered raw
  count table (the moose check used `moose_cell_counts_unmasked.csv`) to
  confirm the mask wasn't already baked in at *data-prep* time. For moose
  this check reversed the earlier assumption entirely — the 382-cell
  footprint looked mask-free but was in fact 100% coincident with the old
  IUCN mask's inclusion set. Do this check once per species before writing
  any "zoomed to the model's real data footprint, no mask applied" caption.
