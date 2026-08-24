# Prompt for Claude Code (Hazel) — array-level estimator sweep, bobcat only

Paste everything below this line into Claude Code on a machine with working
Hazel SSH/Duo access.

---

## Context

Working directory: `/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/prototype_spatial_trend`
on Hazel, mirrored locally at
`/Users/rwkays/claude_code/data_integration_arielle/prototype_spatial_trend`
— a git repo; **commit locally with a clear message, never push** (the PI
pushes manually).

We are testing whether an array-level camera estimator
(cameras-as-replicates-within-an-array, instead of the current per-camera
occupancy model) preserves the power to detect a temporal trend when
integrated with iNaturalist data, before considering any change to the
production model. Read `SEED_array_level_project.md`,
`estimator_sweep_spec.md`, and `array_level_test_spec.md` for full
background. Short version:

- 3 estimator arms, one driver script (`01i_run_estimator_sweep.R`):
  `camera_occ` (current production baseline), `array_occ` (array-level
  occupancy), `array_rn` (array-level Royle-Nichols). All three share a
  byte-identical trend block (`year_beta`, `year_var`, `total_var_beta`,
  `snr`, `trend_robust_indicator`) — only the camera-observation block
  differs.
- Array construction (5 km radius, complete linkage, array-year unit,
  singletons dropped, no array random effect) is **already decided** —
  SEED §4. Do not re-litigate it.
- **Bobcat only.** Moose is a sample-size mechanics check, not a
  verdict (SEED §11) — out of scope here.
- All fitting has to happen on Hazel. `nimbleModel()` construction can
  be (and has been) validated off-cluster, but `compileNimble()` →
  `dyn.load()` cannot run anywhere but here.

## Where things stand — see `git log --oneline` for the full history

`inputs$site_array` is built and validated end to end
(`00b3_validate_array_split.R`, `00c_build_site_array.R`): the
re-derived 5 km complete-linkage split matches the agreed reference
almost exactly (4,003 vs. 4,002 array-years — one cluster boundary
sitting on the 5 km cut height, already root-caused, not a rule
mismatch). Use `prepped_sim_inputs_with_array.RDS` as the input file
(the original `prepped_sim_inputs.RDS` is untouched).

Since then, a chain of real bugs in `01i_run_estimator_sweep.R` and its
dependencies has been found and fixed, entirely by construction-level
`nimbleModel()` testing off-cluster plus on-cluster runs. In order:

1. `model_code_array_occ.R`/`model_code_array_rn.R`'s iNat block called
   `calcIntensity_SVC()` with the wrong argument name and missing year
   indices, despite claiming to be byte-identical to the baseline.
2. `build_array_constants()` never built `array_yday` or `weights`,
   both read directly by the array model code.
3. `fit_replicate()` hardcoded `data=list(y=...)`, silently dropping
   the array arms' `w` data (NIMBLE ignores unused data-list keys
   rather than erroring) and never wired `N_a` inits for the RN arm.
4. `01i` passed `fit_replicate()`'s CAR-model default
   `trend_inits`/`extra_monitors`, referencing nodes none of these
   three models have.
5. `01i` never unwrapped `build_reduced_constants()`'s wrapper list
   (`list(constants=, inat_effort=, y_ncol=, cell50_keep=,
   site_keep=)`), so `constants$adj`/`$J`/`$year_occ`/`$site_keep` all
   resolved to `NULL` for every arm, and used the wrong (full-size)
   `inat_effort` matrix.
6. Wrong input filename (`sim_inputs.RDS` doesn't exist — the real
   file is `prepped_sim_inputs_with_array.RDS`), wrong field names
   (`inputs$cl` etc. instead of `inputs$constants_list` etc.), and a
   missing `source()` for `integration_helper.R` (which defines
   `calcIntensity_SVC` — every model code file calls it).
7. A stale `inputs$y_template` reference (that field doesn't exist;
   `reduced$y_ncol`, already computed correctly inside
   `build_reduced_constants()`, is used instead).
8. **Site-reduction redesign.** Independent random-site sampling (the
   original `build_reduced_constants()` default, `n_site_keep=700`)
   was found to destroy the array structure being tested: real bobcat
   arrays average 6.17 cameras/array (median 4, max 65); after
   independent-site reduction the surviving arrays averaged 2.36
   (median 2), because random sampling scatters cameras across units
   and singleton-dropping then discards most of what was kept. Fixed
   by adding `sample_sites_by_array()` (in `sim_helpers_array.R`),
   which samples **whole array-year units** up to a ~700-camera
   budget, wired into `01i` via `build_reduced_constants(...,
   site_keep_override = sample_sites_by_array(...)$site_keep)`. This
   changes which sites `camera_occ` fits too — intentional, so all
   three arms are compared on the same retained-site set.
9. While building (8), found `assign_arrays_from_field()` was pooling
   sites with no real array label (`NA`) into a shared FAKE array with
   other same-year `NA` sites (`as.character(NA)` → literal `"NA"`
   string → shared join key). Fixed to exclude them instead.
10. **The last blocker reported**: all three arms got past model
    *definition* and into C++ *compilation*, then failed identically
    with `Dimension of 'y[i, 1:J[i]]' does not match required
    dimension for the distribution 'dOcc_v'. Necessary dimension is
    1.` Root cause: `y[i, 1:J[i]] ~ dOcc_v(...)` is a ragged
    vectorized declaration, and any real (legitimate) `J[i]==1`
    single-occasion site collapses `1:1` to a scalar, which `dOcc_v` —
    a vector-valued distribution — rejects. Fixed by passing an
    explicit `dimensions=` argument to both `nimbleModel()` calls
    (`simulate_replicate_data()` and `fit_replicate()`, in
    `sim_helpers.R`), pinning `y`'s (or `w`'s) shape to the actual
    supplied matrix regardless of any individual row's length.

**Fix #10 is the one thing in this list not yet confirmed past a
compile step.** `compileNimble()`/`dyn.load()` are blocked in every
off-cluster sandbox used to develop these fixes, so #10 is verified
only via `nimbleModel()` construction with a synthetic `J[i]==1` toy
site. It needs a real compile on Hazel to confirm.

## Step 1 — smoke test

Sync to the latest commit, then run one replicate of each of the 3
estimators at bobcat abundance, `varying` scenario:
`abundance == "bobcat_baseline"`, `scenario == "varying"`,
`rep_id == 1`, across the three `estimator` values. Get their
`row_id`s from `build_design_df()` in `01i_run_estimator_sweep.R`
(ordered `estimator x abundance x scenario x rep_id`) rather than
computing the index by hand.

Run as three individual Slurm jobs (not an array — you want each log
immediately):

```bash
cd /rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/prototype_spatial_trend
export SIM_BASE=.
# find the three row_ids first, then e.g.:
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

Use the `nimble_env`-on-PATH pattern already established in this
project (export the conda env's `bin` onto `PATH`, invoke `Rscript` by
full path — do not rely on `conda activate` in a non-interactive
sbatch script; already debugged once, jobs 435927→435928).
`qos=short` caps below 3h; 2h is the budget known to work.

**If a job fails with the `dOcc_v` dimension error again**: fix #10
wasn't sufficient — the ragged declaration itself likely needs
reformulating (e.g., an explicit scalar `dbern()` path for `J[i]==1`
sites rather than routing them through `dOcc_v`). Report back rather
than attempting that rewrite unilaterally — it touches the
byte-identical camera-observation block shared across all three arms
and needs checking against all of them consistently.

**If all three jobs succeed**, check before moving on:
- All three produced `status == "OK"` rows in
  `estimator_sweep_out/row_*.rds`, not `"Errored"`.
- RN arm: the `N_a` posterior isn't piling up against its cap
  (`N_max`) — if it is, the cap needs raising.
- `array_diag$mean_prop_detect`/`frac_saturated` are a sane order of
  magnitude given the whole-array-sampled subset actually retained
  (not necessarily identical to the full-data 0.17-0.18 /
  ~1.4% figures, since this is now a ~700-camera subsample of whole
  arrays, not the full 26,748 — but wildly different would mean the
  array assignment doesn't match the real structure).
- Elapsed time per fit, to size the pilot's wall-clock request from a
  real number rather than the camera-level abundance-sweep's own
  measured min 24 / median 42 / max ~90 min (array arms may run
  faster on the camera side, ~78% fewer nodes per the node-accounting
  table in `array_level_test_spec.md`, but RN's discrete `N_a`
  sampler may mix more slowly — don't assume it nets out).

**Do not proceed to the full pilot until all three smoke-test rows
report `status == "OK"` with sane diagnostics.**

## Step 2 — full N_REP=3 pilot (54 tasks, bobcat only)

3 estimators × 3 abundance levels × 2 scenarios × `N_REP=3` = 54 tasks:

```bash
cd /rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/prototype_spatial_trend
sbatch --job-name=estimator_pilot --partition=compute_partners --qos=short \
  --mem=8G --time=02:00:00 --array=1-54%10 \
  --export=ALL,N_REP=3,SIM_BASE=. \
  --wrap='export PATH="'"$PWD"'/../HPC/conda_envs/nimble_env/bin:$PATH"; Rscript 01i_run_estimator_sweep.R $SLURM_ARRAY_TASK_ID'
```

Adjust `--time` from the smoke test's measured elapsed time, with
margin. Keep the `%10` concurrency throttle.

**Reseed and degeneracy guard — verify these actually fired.** This
project has shipped the missing-per-replicate-reseed bug twice already
(SEED §7), each time producing an n=1 result that looked like n=30.
After the pilot completes:

1. Confirm the 3 replicates within each (estimator × abundance ×
   scenario) cell are not byte-identical — spot check a couple of
   `tvb_mean` values by hand.
2. Run `summarize_estimator_sweep()` (`sim_helpers_estimator_metrics.R`)
   on the collected rows — it stops with an explicit error on the
   degenerate case, so a clean pass through it is itself a check.

## Step 3 — report back

Per estimator arm, at bobcat abundance:

- `detect_rate` (power, `varying`) and false-positive rate (same
  column, `null`)
- `bias`, `ci_width`, `coverage` (`tvb_covered`)
- `tri_rate` (camera-corroboration indicator — does an array estimator
  raise it above the real-data baseline of 0.478/0.432?)
- `auc`, `spearman`, `field_rmse` (discrimination family)
- median elapsed time per fit

And the explicit question this sweep exists to answer (SEED §12):
**does any array-level estimator improve on the camera-level
baseline's 66% sign recovery and 70% false-negative rate at
bobcat-like abundance** (SEED §3)? State plainly if none of them
helps — that's a valid, useful result, not a failure to report around.

**Do not push any commits.** Commit locally with a clear message; the
PI pushes manually.
