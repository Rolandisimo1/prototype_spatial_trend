# iSDM trend simulations — reproducible code base

Self-contained bundle for independently verifying the simulation results in the
camera/iNaturalist integrated trend-modelling project. Everything here runs
without cluster access and without the real data.

## What you can and cannot reproduce

**You can verify the code does what the methods describe.** Run the smoke test
and the synthetic-input pipeline; confirm the estimators are implemented as
specified, that the four models share a byte-identical trend block, that the
design grid is what we claim, and that the metrics and guards behave.

**You cannot reproduce our reported numbers from this bundle alone.** The
production inputs (`sim_inputs.RDS`) are built by `00_prep_sim_inputs.R` from
the real camera/iNaturalist bundle plus a fitted posterior, both held on the
NCSU Hazel cluster behind institutional access. `00z_make_synthetic_inputs.R`
writes a structurally identical stand-in with plausible but invented values, so
the pipeline executes end to end. Results from synthetic inputs will differ from
the reported results and are not comparable to them.

To reproduce the reported numbers you need the real bundle. Contact the project
PI for access; then skip `00z` and use the real `sim_inputs.RDS`.

## Requirements

R >= 4.3 with `nimble` (>= 1.4) and `nimbleEcology`, both from CRAN:

    install.packages(c("nimble", "nimbleEcology"))

NIMBLE compiles C++ at run time, so a working toolchain is required. No other
packages are needed; the collectors use base R only.

## Quick start

    Rscript 00y_smoke_test.R                                   # 21 structural checks
    Rscript 00z_make_synthetic_inputs.R prepped_sim_inputs_with_array.RDS
    N_REP=30 N_BURNIN=50 N_ITER=100 \
      Rscript 01i_run_estimator_sweep.R 541                    # one camera_rn replicate
    Rscript 01j_collect_estimator_sweep.R                      # aggregate

Note the input **filename matters**: the driver reads
`prepped_sim_inputs_with_array.RDS`, not `sim_inputs.RDS`.

Verified on this bundle: row 541 completes with `status: OK` in ~35 s at
`N_BURNIN=50 N_ITER=100`. Production settings are `N_BURNIN=2000 N_ITER=8000`,
a median 42 min per replicate; the full 720-task sweep is designed for a job
array (`01i_camera_rn_sbatch.sh` shows the pattern).

### Row numbering

Row ids map to design cells in sorted order (estimator, abundance, scenario,
rep). With four arms: rows 1-180 `array_occ`, 181-360 `array_rn`, 361-540
`camera_occ`, 541-720 `camera_rn`.

## The estimator sweep

The central experiment. A 2x2 of camera submodel by aggregation level, crossed
with abundance and trend scenario:

|                  | occupancy                     | Royle-Nichols            |
|------------------|-------------------------------|--------------------------|
| **camera level** | `model_code_national_scalar.R`| `model_code_camera_rn.R` |
| **array level**  | `model_code_array_occ.R`      | `model_code_array_rn.R`  |

4 estimators x 3 abundance levels x 2 scenarios x 30 replicates = 720 tasks.

**Design invariant.** All four models share a byte-identical trend block and an
identical iNaturalist likelihood; only the camera observation block differs.
Any difference in trend recovery is therefore attributable to the estimator.
The smoke test verifies this by hashing the block in all four files.

## File map

| File | Role |
|---|---|
| `00y_smoke_test.R` | 21 structural checks; run this first |
| `00z_make_synthetic_inputs.R` | synthetic `sim_inputs.RDS` for standalone runs |
| `00_prep_sim_inputs.R` | real input prep (needs cluster data) |
| `00c_build_site_array.R` | 5 km array construction from the camera array field |
| `00d_prep_rn_counts.R` | per-window detection counts, with QC |
| `sim_helpers.R` | design reduction, simulation, MCMC fitting, metrics |
| `sim_helpers_abundance.R` | measured abundance ladder |
| `sim_helpers_array.R` | array aggregation and array constants |
| `sim_helpers_effort.R` | camera and iNaturalist effort levels |
| `sim_helpers_estimator_metrics.R` | power, calibration, discrimination, guards |
| `model_code_*.R` | the model variants |
| `01e_run_abundance_sweep.R` | regional sign-recovery sweep (ecoregion model) |
| `01f_run_indicator_test.R` | RJMCMC inclusion-indicator test |
| `01g_run_effort_sweep.R` | camera x iNaturalist effort grid |
| `01i_run_estimator_sweep.R` | the 4-arm estimator sweep |
| `02*.R`, `01j_*.R` | collectors |
| `results/` | our result CSVs, for comparison |

## Traps found while verifying this bundle

Each of these produced a failure that parse checks and `nimbleModel()` build
checks did **not** catch. They are recorded because an independent reviewer
regenerating inputs will hit the same ones.

1. **`dcar_normal` and the `weights` argument.** The baseline camera model omits
   `weights` (NIMBLE defaults to unit weights); the array models pass
   `weights[1:nadj]` explicitly, which works only because
   `build_array_constants()` creates that vector. `build_reduced_constants()`
   does **not**. A camera-level model written from the array template therefore
   fails inside `buildMCMC()` at `CAR_normal_processParams()` with "missing
   value where TRUE/FALSE needed" — after the model has already built and
   samplers have been configured.
2. **Build-time branch flags.** `has_SVC`, `hasSVC` (both, separately), and
   `prior_type` are resolved by `nimbleModel()` at graph-construction time.
   `prior_type` is compared to the literal `"Normal"` and is case-sensitive; a
   missing or mis-cased value fails in `codeProcessIfThenElse()` with a bare
   "argument is of length zero" and no line number.
3. **`inat_cell50_start`/`_end` index ROWS of `xdat_inat`**, not 50 km cell ids.
   Conflating the two yields out-of-range indices that surface much later as
   "NA/NaN argument" inside `build_reduced_constants()`.
4. **`n_site_target` is hardcoded at 700** in the driver's
   `sample_sites_by_array()` call. A small toy grid cannot satisfy it and fails
   with "missing value where TRUE/FALSE needed". This is why the synthetic
   inputs are production-scale; reduce `N_BURNIN`/`N_ITER` for speed, not the
   grid.
5. **`build_reduced_constants()` argument order** is
   `(cl, inat_effort_real, real_y_template, cell100_geo, ...)`. Call it by name.

## Reproducibility notes

**Seeding.** Every driver reseeds per replicate after `build_reduced_constants()`,
which calls `set.seed()` internally. Omitting this reseed makes all replicates
in a cell identical — a real defect that reached results once in this project
and was retracted. `summarize_estimator_sweep()` now refuses to aggregate a cell
whose replicates are all identical, and the smoke test verifies the guard fires.

**Reading the result CSVs.** The scenario column contains the literal string
`"null"`. In pandas this is parsed as a missing value by default and silently
empties half the data; read with
`keep_default_na=False, na_values=['NA','']`. Base R's `read.csv` is unaffected.

**No chunked MCMC.** These fits use a single `Cmcmc$run()` call. They are
therefore unaffected by the checkpoint state-restoration defect documented in
`REVISION_NOTE_resume_defect.md`, which affects the project's real-data fits.

## Interpreting the metrics

Two families are recorded because they can diverge, and reporting only one is
misleading:

- **Power** — does the 95% credible interval exclude zero? Note this is a strict
  threshold; a run can identify a trend's direction reliably while failing it.
- **Sign recovery** — does the posterior mean have the correct sign? The lenient
  threshold on the same underlying quantity.
- **Discrimination** — AUC, Spearman, and field RMSE against the true spatial field.

Both scenarios in the estimator sweep carry the same true global trend
(-0.1795); the `null` scenario zeroes the *regional* deviations only. The
detection rate under `null` is therefore **not** a false-positive rate.

## Limitations

These are correctly-specified simulations: data are generated by the same
structure that is then fitted, so real-data performance will be worse. The
Royle-Nichols arms' known weakness — unmodelled within-unit detection
heterogeneity, which inflates abundance — is absent by construction, so those
arms are flattered here.
