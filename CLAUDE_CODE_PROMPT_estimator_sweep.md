## STATUS UPDATE (read this first — supersedes the original prompt below)

Since the original prompt below was written, three more rounds of fixes
happened (see `git log --oneline` — commits `c69ec6e`, `1c5d8a0`,
`96b6c9f`, on top of the two the original prompt already reflects). Do
**not** restart from the original prompt's Step 1-3 as written; pick up
from here.

**What's now done and confirmed:**

- `inputs$site_array` is built and validated
  (`00b3_validate_array_split.R`, `00c_build_site_array.R`) — the
  re-derived 5 km complete-linkage split matches the agreed reference
  numbers essentially exactly (4,003 vs 4,002 array-years; the single
  discrepancy is one cluster boundary sitting on the 5 km cut height,
  not a rule mismatch — already root-caused and reported, not silently
  accepted). `prepped_sim_inputs_with_array.RDS` is the input file to
  use; the original `prepped_sim_inputs.RDS` is untouched.
- Two upstream data-parsing bugs were found and fixed during that work
  (a CSV quoting bug that shifted columns, and an integer-factor-code
  vs. label mismatch in `occ.covs` — see commit `c69ec6e` for detail).
- Five more bugs in `01i_run_estimator_sweep.R` beyond the original
  four were found and fixed: wrong input filename, wrong field names
  (`inputs$cl` etc. instead of the real `inputs$constants_list` etc.),
  a missing `source()` for `integration_helper.R`, and (this session,
  locally) a stale `inputs$y_template` reference that should be
  `reduced$y_ncol`.
- **The pilot's site-reduction design has changed.** Independent
  random-site sampling (the original `build_reduced_constants()`
  default) was found to destroy the array structure being tested —
  see `sample_sites_by_array()` in `sim_helpers_array.R` (commit
  `96b6c9f`) for the diagnostic numbers and the fix. `01i` now samples
  **whole array-year units** (up to a ~700-camera budget) instead of
  independent sites, via `build_reduced_constants(...,
  site_keep_override = sample_sites_by_array(...)$site_keep)`. This
  changes which sites `camera_occ` fits too (not just the array arms)
  — that's intentional; all three arms should be compared on the same
  retained-site set either way.
- A latent bug in `assign_arrays_from_field()` was also fixed: sites
  with no real array label (`NA`) were being spuriously pooled into a
  shared fake array with other same-year `NA` sites. Now excluded
  correctly.
- **The current blocker**, per your last report: all three arms
  (including `camera_occ`) got past model *definition* and into C++
  *compilation*, then failed identically with:

  ```
  Dimension of 'y[i, 1:J[i]]' does not match required dimension for the
  distribution 'dOcc_v'. Necessary dimension is 1.
  ```

  Root cause (diagnosed locally, fix applied in commit `1c5d8a0`, but
  **not yet verified past `nimbleModel()` construction** since
  `compileNimble()`/`dyn.load()` are blocked in this sandbox): this is
  a ragged vectorized declaration (`y[i, 1:J[i]] ~ dOcc_v(...)`), and
  whenever some `J[i] == 1` (a real, legitimate single-occasion camera
  site — not a droppable degenerate case), `1:1` collapses to a scalar
  rather than a length-1 vector, which `dOcc_v` — a vector-valued
  distribution — rejects. Fixed by passing an explicit `dimensions=`
  argument to both `nimbleModel()` calls (`simulate_replicate_data()`
  and `fit_replicate()` in `sim_helpers.R`), pinning `y`'s (or `w`'s)
  shape to the actual supplied matrix's dimensions regardless of any
  individual row's real length.

## What to do next

1. `git pull`/sync so you have commits `c69ec6e` through `96b6c9f`.
2. Re-run the smoke test (one rep of `camera_occ`, `array_occ`,
   `array_rn` at bobcat abundance, `varying` scenario). Use the same
   `--time=02:00:00` budget as below — the earlier `dOcc_v` failure
   happened well before the MCMC budget mattered, but the corrected
   run needs real wall-clock room.
3. **If the `dOcc_v` compile error is gone**: check the resulting
   `array_diag$mean_prop_detect` / camera-count numbers against the
   whole-array-sampling design's expectations (should look like a
   genuine, if smaller, version of the real array structure — not the
   near-singleton distribution the old random-site design produced),
   then proceed to the full N_REP=3 pilot as originally planned.
4. **If the `dOcc_v` error persists** (i.e. the `dimensions=` fix
   wasn't sufficient): this means the fix needs to go further than
   pinning dimensions — likely reformulating the ragged declaration
   itself (e.g., splitting out `J[i]==1` sites into an explicit scalar
   `dbern()` likelihood rather than routing them through `dOcc_v` at
   all). Report back rather than attempting that rewrite blind — it
   touches the byte-identical camera-observation block shared
   assumption and needs to be checked against all three arms
   consistently.
5. Re-confirm the reseed/degeneracy guard (SEED §7) once you have
   real `status == "OK"` rows — this has shipped broken twice before
   and neither previous round of fixes touched that logic, so it's
   still an open risk once fits start succeeding.

Everything below this line is the original prompt, kept for the
smoke-test/pilot/reporting steps (3 onward) that are still accurate —
just skip its Step 1 (site_array prep, now done) and don't reintroduce
independent random-site sampling.

---
# Prompt for Claude Code (Hazel) — array-level estimator sweep, bobcat only

Paste everything below this line into Claude Code on a machine with working
Hazel SSH/Duo access.

---

## Context

Working directory: `/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/prototype_spatial_trend`
on Hazel, mirrored locally at
`/Users/rwkays/claude_code/data_integration_arielle/prototype_spatial_trend`
(already has git history — this is a git repo; **commit locally, never push**,
per the project's own working rule).

We are running a simulation study to see whether an array-level camera
estimator (cameras-as-replicates-within-an-array, instead of the current
per-camera occupancy model) preserves the power to detect a temporal trend
when integrated with iNaturalist data — before considering a switch away
from the camera-level production model. Read `SEED_array_level_project.md`,
`estimator_sweep_spec.md`, and `sim_estimator_feasibility.md` in this
directory for full background; the short version:

- 3 estimator arms: `camera_occ` (current production baseline),
  `array_occ` (array-level occupancy), `array_rn` (array-level
  Royle-Nichols). All three share a byte-identical trend block
  (`year_beta`, `year_var`, `total_var_beta`, `snr`,
  `trend_robust_indicator`) — only the camera-observation block differs.
- Array construction (5 km radius, complete linkage, array-year unit,
  singletons dropped, no array random effect) is **already decided** —
  see SEED §4. Do not re-litigate it.
- **Bobcat only for this pass.** Moose is explicitly a sample-size
  mechanics check, not a verdict (SEED §11) — skip it for now.
- Local (laptop) execution of this pipeline is **not possible**: I ran
  `nimbleModel()` construction tests locally and confirmed the model code
  builds correctly, but this sandbox blocks `dyn.load()` of freshly
  compiled shared libraries (a security restriction, not a bug), so no
  actual MCMC can run there. All fitting — smoke test, pilot, and full
  sweep — has to happen here on Hazel.

## Bug fixes already made (local session, committed, not pushed)

Before you start, `git log -1 --stat` in this directory to see what
changed. Four bugs were found and fixed via local `nimbleModel()`
construction testing (build succeeds without compiling/running — this
sandbox permits that much, just not `dyn.load`):

1. **`model_code_array_occ.R` / `model_code_array_rn.R`** — the iNat
   intensity block's `calcIntensity_SVC()` call used the wrong argument
   name (`intercept=` instead of `intensity_intercept=`) and dropped the
   year index `t` from `xdat_inat`/`MWMT_inat`/`MCMT_inat`, contradicting
   both files' own header claim that this block is byte-identical to
   `model_code_national_scalar.R`. Fixed to match baseline exactly.
2. **`sim_helpers_array.R` `build_array_constants()`** — never
   constructed `array_yday` or `weights`, both of which the array model
   files read directly (`array_yday[a,i]` in the detection predictor;
   `weights[1:nadj]` in every `dcar_normal()` call). Added both.
3. **`sim_helpers.R` `fit_replicate()`** — hardcoded
   `data = list(y = sim_data$y, ...)`. The array arms pass
   `sim_data$w` (no `$y` at all); NIMBLE silently ignores an unused
   data-list key rather than erroring, so the array models would have
   sampled `w` from its own prior instead of conditioning on simulated
   detections — a silent-failure bug, not a crash. Now builds the data
   list from whichever of `$y`/`$w` is present, and wires
   `sim_data$N_a_init` into inits for the RN arm's latent `N_a` node.
4. **`01i_run_estimator_sweep.R`** — called `fit_replicate()` with its
   *default* `trend_inits`/`extra_monitors` (`year_effect`/`tau_year`),
   which are CAR-ecoregion-model nodes that none of these three sweep
   models have. `configureMCMC()$addMonitors()` would error `"These
   variables are not in model"` for every arm. Fixed to pass
   `trend_inits = list(year_beta = 0, year_var = 0)`,
   `extra_monitors = character(0)`, matching `sim_helpers.R`'s own
   docstring guidance for the national-scalar case.

**These fixes were verified only at the `nimbleModel()` construction
level** (on synthetic toy data, not the real bobcat structure) — build
succeeds, `w`/`N_a` correctly register as data/inits,
`configureMCMC()+addMonitors()` succeeds for all three arms. **They have
never been compiled or run.** Your first job (the smoke test below) is
the first real test of whether they're actually correct end-to-end.

## Step 1 — `site_array` prep (the one real blocker)

`00_prep_sim_inputs.R` currently does not carry an array-label field.
Add it:

1. Load `umflist.RDS` (or wherever `siteCovs` for bobcat currently
   lives on Hazel — check `input_data_bobcat.RDS` / the real fit
   pipeline for the exact path) and pull `subproject_name` /
   `camera_trap_array` per site, in the **same site order** that
   `cl$cell` (from `input_data_bobcat.RDS$constants_list`) uses. This
   alignment is the part that has to happen on Hazel — the ordering is
   internal to a Hazel-resident RDS and can't be reconstructed
   elsewhere.
2. Apply the already-agreed spatial split: 5 km radius, **complete
   linkage** (not single linkage — see SEED §4 for why), crossed with
   year, singletons dropped. If a working implementation of this split
   already exists from the phase-1 discovery work (check
   `array_units_r5km.csv` / whatever produced it — that file's
   `unit_5` column is exactly this label, already computed for the full
   raw dataset), prefer reusing it over re-deriving the rule, and
   confirm the counts still match: **4,002 array-years, 26,748 cameras,
   963,494 trap-nights, median 4 cameras/unit** (SEED §4). If your
   numbers differ, stop and report the discrepancy rather than
   proceeding on a different structure.
3. Add the result as `inputs$site_array` in `00_prep_sim_inputs.R`'s
   output (`prepped_sim_inputs.RDS`, or wherever `sim_inputs.RDS` is
   assembled for `01i_run_estimator_sweep.R` — check whether that file
   is currently produced by a separate step; `01i` reads
   `sim_inputs.RDS` directly and the seed doc says it "is built on
   Hazel and is not local," so confirm which script currently produces
   it and add `site_array` there).
4. Confirm `build_reduced_constants()` (`sim_helpers.R`) already returns
   `site_keep` — I checked this locally and it does
   (`list(constants = constants, ..., site_keep = site_keep)` at the
   end of the function). This blocker from the original spec appears
   already resolved; just confirm it wasn't reverted.
5. Report the realized array-year count, cameras-per-unit distribution,
   and singleton count for bobcat specifically, post-range-mask —
   compare against the deciding-arm numbers already measured
   (762 informative array-years, 751 graded, SEED §5) so a mismatch is
   caught immediately rather than discovered mid-sweep.

## Step 2 — smoke test (3 fits, not the full pilot)

Before committing several hours of cluster time, run **one replicate of
each of the 3 estimators at bobcat abundance, `varying` scenario only**
— i.e. rows where `abundance == "bobcat_baseline"`,
`scenario == "varying"`, `rep_id == 1`, across the three `estimator`
values. Find their `row_id`s from `build_design_df()` in
`01i_run_estimator_sweep.R` (the design is ordered
`estimator x abundance x scenario x rep_id`, so with `N_REP=20` the
three bobcat/varying/rep1 rows are at estimator-block boundaries —
print `design_df` and grep for them rather than computing the index by
hand).

Run these **as three individual Slurm jobs** (not an array, for this
smoke test — you want each one's log immediately, not batched), with
`N_REP` irrelevant here since you're pinning specific `row_id`s:

```bash
cd /rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/prototype_spatial_trend
export SIM_BASE=.
# find the three row_ids first (see above), then e.g.:
sbatch --job-name=smoke_cam --partition=compute_partners --qos=short \
  --mem=8G --time=02:00:00 \
  --wrap="Rscript 01i_run_estimator_sweep.R <row_id_camera_occ>"
sbatch --job-name=smoke_arrocc --partition=compute_partners --qos=short \
  --mem=8G --time=02:00:00 \
  --wrap="Rscript 01i_run_estimator_sweep.R <row_id_array_occ>"
sbatch --job-name=smoke_arrrn --partition=compute_partners --qos=short \
  --mem=8G --time=02:00:00 \
  --wrap="Rscript 01i_run_estimator_sweep.R <row_id_array_rn>"
```

Use the `nimble_env`-on-PATH pattern from `smoke_test_sbatch.sh` /
`01e_run_abundance_sweep_sbatch.sh` (export the conda env's `bin` onto
`PATH` and invoke `Rscript` by full path directly — do **not** rely on
`conda activate` inside a non-interactive sbatch script; this project
already hit and debugged that exact failure mode, see the comments in
those two files).

**What to check when these finish:**
- All three produced `status == "OK"` rows in
  `estimator_sweep_out/row_*.rds`, not `"Errored"`.
- If any errored, read the error message first against the 4 fixes
  above — check whether the fix didn't fully take, or whether this is a
  genuinely new issue the toy-data construction test couldn't have
  caught (real bobcat data has different dimensionality, real
  `array_yday`/effort ranges, etc.).
- For the RN arm specifically: check the `N_a` posterior isn't piling
  up against `N_max` (the RN model's cap) — if it is, the cap needs
  raising (see `array_level_test_spec.md`'s compute-placement notes on
  `N_a` mixing).
- Compare `array_diag$mean_prop_detect` / `frac_saturated` against the
  measured real-data numbers (proportion 0.17-0.18, ~1.4% saturating,
  SEED §5) — a wildly different number in the simulated array structure
  would mean the simulation's array assignment doesn't match reality.
- Elapsed time per fit, so you can size the full pilot's wall-clock
  request from a real number instead of the abundance-sweep's 8GB/42min
  estimate (that number is from the *camera-level* model; array models
  have a smaller camera-side node count so may run faster, per
  `array_level_test_spec.md`'s node-accounting table — ~78% fewer nodes
  on the camera side — but the RN arm's discrete `N_a` sampler may be
  slower to mix, so don't assume it nets out).

**Do not proceed to the full pilot until all three smoke-test rows
report `status == "OK"` with sane diagnostics.** If one arm keeps
failing, report back rather than silently dropping it from the design.

## Step 3 — full N_REP=3 pilot (54 tasks, bobcat only)

Once the smoke test passes, run the full pilot design as specified in
`estimator_sweep_spec.md` and `01i_run_estimator_sweep.R`'s own header:
3 estimators × 3 abundance levels × 2 scenarios × `N_REP=3` = 54 tasks.

```bash
cd /rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/prototype_spatial_trend
sbatch --job-name=estimator_pilot --partition=compute_partners --qos=short \
  --mem=8G --time=02:00:00 --array=1-54%10 \
  --export=ALL,N_REP=3,SIM_BASE=. \
  --wrap='export PATH="'"$PWD"'/../HPC/conda_envs/nimble_env/bin:$PATH"; Rscript 01i_run_estimator_sweep.R $SLURM_ARRAY_TASK_ID'
```

Adjust `--time` upward from the smoke test's measured per-fit elapsed
time with margin (the abundance-sweep baseline was min 24 / median 42 /
max ~90 min per fit on the camera-level model, per SEED §10; budget
accordingly once you know the array arms' actual numbers). Keep the `%10`
concurrency
throttle — this project's convention on the shared partition.

**Reseed and degeneracy guard — do not skip verifying these actually
fired.** This project has shipped the missing-per-replicate-reseed bug
**twice** already (see SEED §7), each time producing an n=1 result that
looked like n=30. After the pilot completes:

1. Confirm the 3 replicates within each (estimator × abundance ×
   scenario) cell are **not byte-identical** — spot check a couple of
   `tvb_mean` values by hand before trusting the degeneracy guard.
2. Run `summarize_estimator_sweep()` (in
   `sim_helpers_estimator_metrics.R`) on the collected rows — it stops
   with an explicit error if it detects the degenerate case, so a clean
   run through it is itself a check.

## Step 4 — report back

Report, per estimator arm, at bobcat abundance:

- `detect_rate` (power, under `varying`) and the false-positive rate
  (same column, under `null`)
- `bias`, `ci_width`, `coverage` (`tvb_covered`)
- `tri_rate` (the camera-corroboration indicator — SEED calls this
  potentially the most important single number: does an array
  estimator raise it above the real-data baseline of 0.478/0.432?)
- `auc`, `spearman`, `field_rmse` (discrimination family)
- median elapsed time per fit, per arm

And the explicit comparison SEED §12 asks for: **does any array-level
estimator improve on the camera-level baseline's 66% sign recovery and
70% false-negative rate at bobcat-like abundance** (SEED §3)? State
plainly if the answer is that none of them helps — that is a valid and
useful result, not a failure to report around.

**Do not push any commits.** Commit locally with a clear message; I'll
push manually.
