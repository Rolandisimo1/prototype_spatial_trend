# Camera-anchored regional trends: design specification

This file is the specification referred to in the 23 November handoff. It is
reproduced below unchanged, with two sections added at the end: an independent
confirmation of the diagnosis against the model source, and runtime evidence
for sizing the job array.

Nothing has been built. No fork exists, no job has been submitted.

---

Camera-anchored regional trends: new model fork and validation simulation
Roland's decision: cameras are the gold standard for every reported trend, national AND regional. The current regional model does not meet that standard, so the regional simulations and regional maps are to be redone on a camera-anchored model.

The problem, verified in the code
model_code_ecoregion_trend.R:

camera occupancy: cloglog(psi[i]) <- ... + year_beta * year_occ[i] (lines ~117-126), no regional term;
iNat intensity: total_var_beta + year_effect[inat_cell100[g]] (lines ~57-71), with year_effect[c] <- year_region[ecoregion_of_cell100[c]] (lines ~183-188).
So year_region is identified by iNaturalist alone. Every simulation that scored it (abundance sweep 01e, indicator test 01f, effort sweeps 01g/01h) generated data with simulate() from this same model, so the simulated CAMERA data carried no regional trend at all. Those results describe an iNat-only quantity, and the effort sweep's "iNat matters 3x more than cameras" follows from the model structure. REVISED_camera_anchored_trends.md §3 already recorded this for the real fits.

Step 1. New model fork (do not edit any existing model file)
Create model_code_ecoregion_camanchored.R, forked from model_code_ecoregion_trend.R, banner header + roxygen, lab style. Mirror the national shared-parameter + sink structure at the regional level:



# camera occupancy: regional deviation is SHARED, so camera-identified
cloglog(psi[i]) <- ... + (year_beta + year_region[region_of_site[i]]) * year_occ[i]

# iNat intensity: shared regional deviation plus its own iNat-only regional sink
total_var_beta_cell[c] <- year_beta + year_var + year_region[ecoregion_of_cell100[c]]
                                              + year_region_var[ecoregion_of_cell100[c]]

year_region[r]     ~ dnorm(0, sd = sigma_region)       # reportable, camera-anchored
year_region_var[r] ~ dnorm(0, sd = sigma_region_var)   # iNat-only sink, never reported
The reportable regional trend is year_beta + year_region[r], converted to %/yr with the CAMERA year scaling.
year_region_var is the regional analogue of year_var: it absorbs iNat-specific regional change (observer growth, campaigns) so it cannot leak into the camera-anchored term. It is not optional.
Camera sites need a region_of_site constant (ecoregion of the site's cell100).
Keep EPA/CEC Level I ecoregions as the regional unit, so results line up with the earlier runs.
Step 2. Validation simulation
Generate data from the NEW model so the regional trend reaches both streams. Three truth scenarios, each defined by zeroing or setting trend parameters (not the spatial field):

Scenario	year_region	year_region_var	Tests
NULL	0	0	false positive rate of the regional term
REAL	sd 0.2 across regions	0	power, bias, coverage of camera-anchored regional trends
SPURIOUS	0	sd 0.2	does iNat-only regional change leak into year_region?
REAL+SPURIOUS	sd 0.2	sd 0.2, independent	same, with a real signal present
Cross with data volume at REAL fleet scale, using distinct real camera sites (no resampling with replacement; the earlier 2x/4x camera arms duplicated sites and added no new places). Build on the Tier 1 design infrastructure (build_tier1_constants()), which already handles 1,000-25,000 sites:

camera sites: 1,000 / 3,000 / 9,000 (9,000 is the Tier 1 mid level, near the real fleet)
iNat records: 5,000 / 20,000 / 100,000
one abundance level (bobcat-like) and one range (300 cell100s), window_frac 0.4
4 scenarios x 3 x 3 x 30 reps = 1,080 fits
Smoke-test first: confirm each axis is live and that a camera-only dataset (iNat likelihood removed) recovers year_region in the REAL scenario. Per-rep reseed after building constants.

What to record per replicate and region
For year_beta + year_region[r] (reportable) and, separately, year_region_var[r]: true value, posterior mean, posterior SD, 95% interval, interval excludes 0, correct sign, R-hat, ESS. Also realized camera detecting deployments and iNat records PER REGION, so regional precision can be related to regional data volume.

Report back
Per scenario x volume: regional sign recovery, power (interval excludes 0 when true |deviation| > 0.1), false positive rate (NULL), bias, coverage, median posterior SD.
Leakage: in SPURIOUS, the share of regions whose year_region interval excludes 0, and the correlation between true year_region_var and estimated year_region. Near zero means the sink works.
Does camera volume or iNat volume buy regional precision now? Report the log-log slope of regional posterior SD on regional camera detections and on regional iNat records.
Convergence failures by cell, separately.
Only after this validates: refit bobcat, white-tailed deer and moose on real data with the new model, plus a camera-only regional fit for each as the control.

Do not
Edit model_code_ecoregion_trend.R, model_code_national_scalar.R or any existing model.
Bundle this with the Tier 1 launch: separate output directory, separate job array.
Report year_region_var or any iNat-side regional quantity as a population trend.

---

## Added: independent confirmation of the diagnosis

Checked against `prototype_spatial_trend/model_code_ecoregion_trend.R` on the
laptop, 23 November. The diagnosis holds exactly as stated.

- Camera occupancy, lines 117-126: `cloglog(psi[i]) <- link_occ_intercept[cell[i]]
  + ... + year_beta * year_occ[i]`. No regional term. The file's own header
  comment at line 36 states that this expression is untouched from the
  pre-ecoregion model.
- iNaturalist intensity, line 71: `total_var_beta = total_var_beta +
  year_effect[inat_cell100[g]]`.
- Line 188: `year_effect[c] <- year_region[ecoregion_of_cell100[c]]`.

So `year_region` enters one likelihood and it is not the camera one.

The four sweeps named in the specification all source this model file and build
truth through `make_true_year_region` / `true_param_list_ecoregion` before
simulating: `01e_run_abundance_sweep.R` (line 12), `01e_run_ecoregion_sim.R`
(line 41), `01f_run_indicator_test.R` (line 52), `01g_run_effort_sweep.R`
(line 51), `01h_run_effort_secondary.R` (line 40).

The consequence for the reported result is therefore structural rather than a
question of power or scale: in the data those sweeps generated, the camera
stream carried no regional trend, so no camera arm could have moved regional
recovery whatever its volume.

## Added: runtime evidence for sizing the array

From `sim_results/effort_sweep_summary.csv`, 480 successful replicates:

- median 23 minutes per replicate, maximum 28 minutes
- largest camera arm in that sweep: 2,800 sites

This design goes to 9,000 camera sites, 3.2 times that arm, across
1,080 fits:

| Assumption | Compute-hours | At 80 concurrent slots | At 11 slots |
|---|---|---|---|
| Runtime flat in site count | 419 | 5 h | 1.6 days |
| Runtime linear in site count | 1,345 | 17 h | 5.1 days |

Two requests follow from this, and neither changes the design:

1. Report the smoke-test result on its own as soon as it lands, before the full
   grid finishes. The check that decides whether the fork works is whether a
   camera-only dataset recovers `year_region` under REAL. That is worth having
   days before the grid.
2. Measure per-fit runtime in the smoke test and size the array against it
   before submitting all 1,080. If something has to give, cut replicates
   before cutting scenarios: SPURIOUS is the scenario that earns the fork.
