# Seed document — array-level iSDM trend analysis

**Purpose of this file.** Self-contained brief for a new conversation. Assume the
reader has no prior context. Everything asserted here was measured in-session;
where a number is quoted, the file it came from is named.

**PI:** Roland Kays, Pacifici lab, NC State. **Cluster:** NCSU Hazel (Slurm).
**Working folder:** `/Users/rwkays/claude_code/data_integration_arielle/prototype_spatial_trend`

---

## 1. What this project is

Fuse camera-trap and iNaturalist data in a NIMBLE hierarchical model to estimate
**temporal trends** in North American mammal abundance/occupancy, with spatially
varying coefficients. The camera and iNat streams inform one shared latent
intensity field: cameras through an occupancy link, iNat through an
effort-offset negative-binomial count likelihood at 50 km grid cells.

**This project is the array-level variant.** The parent project models each
camera as an independent site. Here the camera observation unit is the
**array-year** — a spatially coherent group of cameras surveyed in one year —
with cameras acting as replicates within it. Everything else in the model
(iNat likelihood, CAR spatial fields, climate SVCs, the trend block) is held
identical so any difference is attributable to the camera estimator alone.

**First task: a simulation study.** Do array-level estimators, integrated with
iNat, have the statistical power to detect trends — compared with the
camera-level simulations already completed? Nothing about real-data fitting is
in scope until the simulation answers this.

---

## 2. Why simulate rather than fit real data

Two independent reasons, both established:

1. **A trend signal cannot be validated against itself.** Only a simulation
   with known truth can measure power, bias, and coverage.
2. **The real-data fitting path is currently broken.** The inherited chunked-MCMC
   checkpoint code does not restore chain state: the saved object is a compiled
   structure whose external pointer dies on deserialization, so the assignment
   silently does nothing and each resumed chunk restarts from the model's
   initial values. An audit found restart at **every** resume boundary across
   the whole fleet. No fit in the parent project currently carries a
   trustworthy convergence assessment.

   **The simulation framework is unaffected** — it uses single-call MCMC
   (`Cmcmc$run(n_burnin + n_iter)` and `runMCMC`), with no chunk-resume
   anywhere. It is the one part of this codebase currently able to produce
   defensible results. Do not introduce chunk-resume into simulation code.

---

## 3. The camera-level baseline to compare against

These are the numbers the array-level arms must be judged against. All are from
the **corrected** re-run of the camera-level ecoregion sweep (an earlier run was
retracted — see §7).

**Regional trend sign recovery** (fraction of ecoregions whose estimated trend
has the correct sign, `varying` scenario):

| abundance level | sign recovery | vs 50% chance |
|---|---|---|
| bobcat-like (low) | **66%** | p = 6e-06 |
| ~3x bobcat | 73% | — |
| deer-like | 85% | — |

Genuinely above chance at every level, improving smoothly with information.

**Structure-detection gate.** WAIC could NOT discriminate a real regional trend
from none at any abundance level (p = 0.80 / 0.91 / 0.61) — **WAIC was dropped
as a gate.** It was replaced by a Bernoulli inclusion indicator on the regional
trend deviation, fit by reversible-jump MCMC, following Goldstein et al.
(bioRxiv 2025.01.17.633640). Results at the 0.5 support threshold:

- False-negative rate at bobcat abundance: **70%** — the test misses a real
  regional trend more often than it catches it.
- At the 0.9 threshold, false positives essentially vanish (0-3%) but
  false negatives reach **100%** at bobcat and moderate abundance.
- Published benchmark for comparison: FP 9.1% / FN 14.9% at 0.5;
  FP 0.7% / FN 38.1% at 0.9.

**The headline limitation the array-level work is trying to move:** at bobcat's
real density, neither WAIC nor the indicator can decide whether regional trend
structure exists. The camera stream is the constraint. That is the bar.

---

## 4. Settled design decisions — do not re-litigate without new evidence

**Analysis unit: array-year at a 5 km radius.** Groups are built from the raw
data's array field (`camera_trap_array` where populated, falling back to
`subproject_name`), crossed with year, then split spatially so no unit exceeds
5 km diameter.

- **5 km because Snapshot USA deploys on a nominal 5 km spacing** — the protocol
  is the design, so the unit should reproduce it. A sample-size argument favours
  2.5 km (bobcat 1,429 vs 1,151 graded informative units; moose 129 vs 106) and
  was deliberately declined: matching the field protocol is the stronger
  consideration.
- **Complete linkage, not single linkage.** Complete linkage *guarantees* the
  diameter bound. Single linkage chains and can produce a unit far larger than
  the threshold.
- Array-year, not array: the array field pools across years for ~27% of groups,
  so the year crossing is mandatory.
- Singleton units (1 camera) are dropped.

**At 5 km:** 4,002 array-years, 26,748 cameras, 963,494 trap-nights, median 4
cameras/unit, 97.0% of cameras retained in units of >= 2.

**NO array-level random effect.** Once the unit *is* the array, cameras within
it are replicates of that array's state, not separate nearby sites — a random
effect would be a group with one observation. Two supporting measurements:

- 89.5% of 5 km arrays nest inside a single 100 km CAR cell, and 70.3% of
  revisited arrays have all their year-visits in one cell — so location
  correlation is already carried by `link_occ_intercept[cell]`.
- Within-array year variance is 0.49 against a total of 14.13, i.e. **3.5%**.
  Median revisit span is 2 years; only 15 arrays span >= 5 years. There is
  almost no within-array trend information to exploit.

**NO seasonality term added.** Tested and declined. The camera side already has
a quadratic day-of-year term on *detection*
(`p_beta[1]*yday + p_beta[2]*yday^2`); the iNat side has none. iNat seasonality
is real and species-specific (moose 6.7x amplitude, peak June) and is **not**
absorbed by the effort offset — for bobcat the offset actually raises amplitude
from 2.1x to 4.6x, because bobcat's seasonality runs opposite to the all-mammal
average. But for a *trend* estimand what matters is whether the seasonal shape
**drifts**, and it is stable:

| species | mean day-of-year drift | as % of the count trend it would bias |
|---|---|---|
| Bobcat | -0.09 d/yr (p=0.89) | 0.2% |
| Moose | -0.53 d/yr (p=0.043) | **2.4%** |
| WT deer | -0.14 d/yr (p=0.60) | 0.6% |

Only moose drifts significantly, and at 2.4% of the trend signal. Record as a
caveat; do not add the parameter.

**One asymmetry that must be stated when results are read.** The camera-level
baseline treats cameras as independent and has no array term, so it carries
pseudo-replication the array arms do not. Variance-weighted cluster size at
5 km is 16.5 cameras; at ICC 0.15 the implied design effect is 3.33, i.e. true
intervals on `year_beta` roughly **82% wider** than the camera-level model
reports. **Therefore wider credible intervals from the array arms are not
automatically evidence that aggregation cost precision** — they may be honest
where the camera-level arm is overconfident. Judge on **coverage**, not width:
`tvb_covered` is recorded per replicate for exactly this reason.

---

## 5. The estimator arms — and an unresolved question about N-mixture

The requested comparison is array occupancy, Royle-Nichols, and N-mixture.
Occupancy and RN are settled and their model code is written and compile-tested.
**N-mixture needs a decision in the new conversation**, because its feasibility
changed mid-analysis and the arms may not be as distinct as they look.

**Arm 1 — array occupancy.** Response is *k of n_a cameras detected*, so
cameras are graded replicates. This was the PI's own specification: the
proportion of cameras detecting, adjusted for detectability. Confirmed sound by
the data — 65 of 66 informative moose array-years and 751 of 762 bobcat
array-years are **graded** (0 < proportion < 1), median proportion 0.17-0.18,
only ~1.4% saturating at 1.0. A plain any-detection aggregation would have
collapsed nearly all of them to a single value. Model:
`model_code_array_occ.R`.

**Arm 2 — array Royle-Nichols.** Latent abundance `N_a` drives detection
probability via `1 - (1 - r)^N_a`. Runs on detection/non-detection; it never
needed counts. `nimbleEcology` ships **no** RN distribution (verified, not
assumed) — the hand-written formulation in `model_code_array_rn.R` is
deliberate. It compiles and returns a valid log-probability, including the
construct with a latent discrete exponent, which was the main uncertainty.

**Arm 3 — N-mixture: status changed twice; resolve before building.**

- Originally dropped: the prepared pipeline binarizes camera data to 0/1 per
  10-day window, and `simulate_replicate_data()` draws from `dOcc_v`, so no
  counts existed to fit.
- Raw camera data was then supplied, and it *does* carry counts — but
  **`group_size` is unusable**: mean 1.08 for moose and 1.04 for bobcat (93%
  and 97% of detections are a single animal), and it is 100% missing for 2 of
  27 project groups — structured whole-project absence, so dropping NAs would
  drop two projects and their geography.
- What the raw data *does* give is **independent detection frequency per
  camera-window**, which is graded and fully covered. See §6 for the numbers.

**The open question:** an N-mixture fit to detection-*frequency* counts is not
classic N-mixture (which counts individuals per visit). Its observation process
is repeat visits by possibly the same animals, which is exactly what RN already
models. So arm 3 risks being a reparameterization of arm 2 rather than a third
estimator. Two defensible resolutions — pick one explicitly:

1. **Binomial-N-mixture on cameras-detecting:** `k_a ~ dbin(p_eff, n_a)` with
   `N_a` latent. Genuinely distinct from RN in its observation model, and uses
   only data we are certain of.
2. **Poisson/NB-N-mixture on per-window detection counts:** closer to what the
   PI likely means, but must be reported as *detection-frequency* N-mixture,
   with the non-independence of repeat visits stated as a known bias source.

Do not build a third arm until this is settled — and say plainly if the
conclusion is that RN and N-mixture cannot be separated on this data.

---

## 6. Data on hand, and how it was verified

Raw files at `/Users/rwkays/claude_code/data_integration_arielle/raw_cam_data/`:

| file | size | contents |
|---|---|---|
| `combined_deployments_all.csv` | 4.2 MB | 26,798 deployment rows |
| `combined_sequences_all.csv` | 574 MB | 2,517,928 sequence-class rows |
| `inat_combo_nam_mams.csv` | 2.2 GB | 3,014,977 iNat records, full `observed_on` date |

**Verified as a clean superset of what the models use** (`umflist.RDS`
`siteCovs`): 26,798 vs 24,869 deployment rows, 60 vs 58 projects, 962 vs 943
subprojects, identical 2008-2025 year range; all 26,517 deployment IDs appear in
both raw files with zero orphans either direction. Independent effort check:
`survey_nights` median 31 days -> 3.1 ten-day windows, against the prep
pipeline's measured mean `J` of 3.22.

**Three QC steps are mandatory and each is load-bearing** (see the
`camera-trap-qc` skill, which was built on this exact dataset):

1. **`sequence_id` is NOT a unique key.** The table carries one row per
   (species x age x sex) class, so a doe with two fawns is legitimately two
   rows on one sequence. 2,517,928 rows -> 2,248,930 distinct sequence IDs;
   122,369 multi-row sequences. Counting raw rows inflates per-deployment
   detections — moose max was 146 by raw rows, **73** after collapsing. Collapse
   to one row per (deployment, sequence, species) for a frequency response.
2. **30-minute independence filter.** Sequences are algorithmically grouped
   triggers, so one lingering animal produces several; RN reads that clustering
   as more animals. Thinning to >= 30 min between consecutive same-species
   detections at a deployment removes **22.6%** of moose events, 7.7% of
   bobcat, **30.0%** of deer.
3. **Effort bounds (0, 365] nights.** `survey_nights` has 87 zero-night
   deployments, 19 NA-year rows, and a max of 2,688 days (~7.4 years — a
   deployment never closed). Bounding drops 115 of 26,798 (0.43%). Effort now
   enters the array arms as a covariate, so an outlier would propagate.

**Measured abundance ladder** — independent detections per occupied 10-day
window, after all three QC steps (`rn_count_calibration.csv`):

| rung | species anchor | mean count | occupied windows | % > 1 | multiplier |
|---|---|---|---|---|---|
| low | bobcat | **1.555** | 6,079 | 28.3% | 1.00x |
| intermediate | (geometric mean) | **3.156** | — | — | 2.03x |
| high | white-tailed deer | **6.405** | 60,125 | 78.6% | 4.12x |

Moose measures 1.757 — nearly identical to bobcat. **Moose is data-POOR, not
low-abundance-per-window**, so it cannot stand in for the low rung; bobcat
anchors it. The intermediate rung is the geometric mean, giving equal log
spacing, correct because the count enters the linear predictor on the log scale.

These multipliers (1.00 / 2.03 / 4.12) **replace** an earlier placeholder ladder
of 1x / 3x / 8x whose own docstring admitted the values were not fitted. The
real span from a widespread-rare carnivore to the most abundant ungulate is
~4x, not ~8x; the placeholder overstated the gradient.

**Sampling season is drifting on the camera side** — mean deployment
day-of-year +3.6 days/yr (p=0.027, 51-day drift over 2011-2025), survey window
narrowing 15 days/yr (p=0.006, 196 -> 51 days wide), September rising from 17%
to 51% of deployments. The quadratic `yday` detection term is what keeps this
out of `year_beta`, so it is load-bearing and **must survive aggregation**. It
is carried as the array mean plus its square, which loses within-array spread —
state that in methods.

---

## 7. Traps this project has already fallen into

**The missing per-replicate reseed — shipped twice.** `build_reduced_constants()`
calls `set.seed()` internally with a fixed seed, and `simulate_replicate_data()`
does not reseed. A driver that omits its own reseed leaves the RNG in an
identical state for every array task, so all "replicates" simulate and fit the
**same dataset**. This produced a sweep where every scientific quantity was
byte-identical across all 30 replicates per cell and only `elapsed_sec` varied —
effectively n=1. Every rate, coverage, and information-criterion claim from that
run was retracted.

Two defences, both already in the code — keep them:
- `set.seed(DESIGN_SEED + row_id)` **after** `build_reduced_constants()` and
  before the fit.
- A **degeneracy guard** in the collector that asserts per-cell `rmse_all` /
  `waic_primary` are not all identical before writing output. It has been tested
  in both directions: it fires on identical replicates and passes on varied ones.

**Two data-pipeline bugs found in the parent project**, both fixed but worth
knowing: an iNat effort matrix built with unsorted column order, so the effort
covariate and count matrix could disagree on year alignment (caught because one
species' trend was implausibly large); and a static range mask that set counts
to missing outside a polygon while leaving the effort denominator unmasked —
dropped, since an observation proves occurrence and the estimand is change.

**Do not trust a `PASS` written in a file.** A prior diagnosis of the checkpoint
defect was confidently worded and its central claim was false: the probe read
its verdict off a node whose sampler draws from its full conditional in one
step, and so cannot distinguish a restored state from a fresh start.

---

## 8. Files — all in `prototype_spatial_trend/`

**Written, parse-clean, compile-tested where applicable:**

| file | role |
|---|---|
| `sim_helpers_array.R` | array aggregation; carries `log_effort` forward since collapsing windows removes automatic effort accounting; `summarize_array_structure()` reports `mean_prop_detect` / `frac_saturated` |
| `model_code_array_occ.R` | arm 1; builds under `dOcc_v` with ragged `len = n_a[a]` |
| `model_code_array_rn.R` | arm 2; hand-written RN, compiles, logProb verified |
| `01i_run_estimator_sweep.R` | design matrix + driver; per-replicate reseed; `tryCatch` returns an "Errored" row rather than killing the batch; radius default 5.0 |
| `sim_helpers_estimator_metrics.R` | both metric families (see §9); degeneracy guard |
| `sim_helpers_abundance.R` | `abundance_levels_measured()` — the §6 ladder |
| `00d_prep_rn_counts.R` | per-window independent-detection counts with all three QC steps, plus `verify_counts_against_binary()` which **stops** on inconsistency rather than patching |
| `estimator_sweep_spec.md` | full run spec, decisions, and their rationale |
| `rn_count_calibration.csv` | measured anchors |
| `array_summary_by_year.csv`, `array_units_r5km.csv`, `array_radius_retention.csv` | array structure tables |
| `raw_data_verification.md`, `sim_estimator_feasibility.md` | verification and feasibility records |

**Design invariant, verified not asserted:** the trend block is byte-identical
across `model_code_national_scalar.R`, `model_code_array_occ.R`, and
`model_code_array_rn.R` (md5 of the normalized block:
`69213000bc266620f68b21c71316358b`). An earlier version of the array models had
the same statements in a different order — NIMBLE is declarative so it compiled
identically, but the stated invariant was broken and was fixed. **Re-check this
hash after any edit to a model file.**

**Three blockers before the sweep can run:**

1. `inputs$site_array` — the 5 km array-year label per camera, added in
   `00_prep_sim_inputs.R` and carried in `sim_inputs.RDS`.
2. `constants$site_keep` — `build_reduced_constants()` computes the retained
   site index internally but does not return it; needed to align `site_array`
   with the reduced design.
3. `sim_inputs.RDS` is built on Hazel and is not local, so the pilot needs it
   fetched.

The driver **stops** if `site_array` is absent rather than silently falling back
to the development-only spatial approximation. Keep that behaviour.

---

## 9. Metrics — record both families

The lab's own prior camera-vs-array comparison (a different model, plain
any-detection aggregation, N ~ 1,090 cameras -> 60 arrays) found aggregation
**reduced statistical power** while sometimes **improving AUC**. Those are
different quantities — detecting *that* there is a trend versus discriminating
*where* — and a study measuring only one would give a misleading verdict.

- **Power family:** `detect_rate` (power under the `varying` scenario, and the
  false-positive rate under `null` — same column, different scenario), bias,
  CI width, and `tvb_covered` for coverage.
- **Discrimination family:** AUC (via the Mann-Whitney identity), Spearman
  correlation, and field RMSE.
- Plus `trend_robust_indicator` = `step(snr - 1)` where `snr = year_beta /
  year_var` — whether the camera stream corroborates the iNat-driven trend.
  Moose currently fails this. **If an array estimator raises it, that matters
  more than any speed gain**, because it is the difference between a trend the
  camera data supports and one it does not.

---

## 10. Compute

Measured from a completed sweep: **8 GB and a median 42 min per replicate fit**
(min 24, max ~90). The laptop (36 GB, 14 cores) fits 3-4 concurrently — memory
binds before cores. A 3-replicate pilot is ~4 h locally; the full design
(3 estimators x 3 abundance x 2 scenarios x 20-30 reps) is ~3 days locally
versus ~4 h as a Hazel array job.

**So: develop and pilot locally, sweep on Hazel.** The local half matters
independently — cluster access needs Duo push authentication that only works
from the PI's own machine, so iterating locally avoids a round-trip per attempt.

Local R environment `isdm-sim` has **NIMBLE 1.4.3 and nimbleEcology installed**
in a workspace-local library (`rlibs`; set `R_LIBS_USER` to it). NIMBLE is
CRAN-only, not available through conda — which is also part of why model builds
are slow: it compiles C++ at runtime.

**Hazel notes:** login nodes have internet, compute nodes do not — install
packages on login. **Never SSH directly to a compute node** — this drew two
acceptable-use warnings; use `srun --jobid=<id> --overlap --nodelist=<node>
--pty /bin/bash`, or prefer `sacct` / `squeue` / `sstat` from the login node.

---

## 11. Working rules

- **Never modify the original analyst's (Arielle / ahwaldst) files or
  environments.** Fork or copy only. This is invariant.
- **Write code to the PI's hard drive and commit locally; never push.** The PI
  pushes manually so he knows when and why code becomes public. Relatedly:
  never claim in a commit message, README, or file header that something is
  validated when it is not — that overclaiming is a mistake this project has
  already had to correct.
- **Claude Science drafts specs, code, analysis, and figures; Claude Code
  executes on Hazel.** Prompts are relayed by hand because of the Duo
  requirement.
- **State the limits of the simulation.** These are correctly-specified-model
  simulations: data are generated by the structure being fitted, so real data
  will be worse. The RN arm's known weakness — unmodelled within-array
  detection heterogeneity inflating abundance — is absent by construction, so
  its simulated performance is optimistic.
- **Moose is a mechanics test, not a verdict.** ~9 informative array-years per
  year for `year_beta`; a moose failure is a sample-size result. **Bobcat is
  the deciding arm.** Report both and do not let moose stand as the headline.
- If a finding contradicts this document, say so rather than substituting an
  approach quietly.

---

## 12. Suggested first steps

1. Resolve the **arm-3 / N-mixture** question in §5 explicitly, and record the
   reasoning. Be willing to conclude there are two arms, not three.
2. Clear the three blockers in §8 — the two prep additions and fetching
   `sim_inputs.RDS`.
3. Run a **3-replicate local pilot** of the existing arms: confirm the driver
   runs end-to-end, the metrics compute, and the degeneracy guard passes on
   genuinely varied replicates.
4. Verify the array-structure diagnostics on simulated data —
   `mean_prop_detect` and `frac_saturated` — so the graded-replicate assumption
   is checked rather than assumed.
5. Only then launch the full sweep on Hazel.
6. Compare against the §3 camera-level baseline. The specific question:
   **does any array-level estimator improve on 66% sign recovery and a 70%
   false-negative rate at bobcat-like abundance?** State plainly if the answer
   is that none of them helps.
