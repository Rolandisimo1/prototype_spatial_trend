# Abundance sweep — seeding regression (must re-run before interpreting)

**Status:** the `abundance_sweep_summary.csv` / `abundance_sweep_recovery.png`
results are NOT interpretable at face value. All 30 replicates within each
abundance×scenario cell are byte-identical — the run is effectively n=1 per
cell (6 datasets total), not n=30.

## Evidence
- Within any cell (e.g. `common_deerlike/varying`), `waic_primary`, `rmse_all`,
  `bias_all`, `info_inat_count_total`, `info_cam_detections_total` each have
  **exactly 1 unique value across all 30 reps**. Only `elapsed_sec` varies
  (30 unique) — the tasks ran as separate processes but computed identical
  results.

## Cause
- `build_reduced_constants()` calls `set.seed(seed)` internally with a FIXED
  seed (sim_helpers.R line 97), by design, so every replicate shares the same
  retained-cell/site subsample.
- `simulate_replicate_data()` does not reseed.
- `01e_run_abundance_sweep.R` — the driver that produced this CSV — has NO
  per-replicate `set.seed()`. So R's global RNG is left in an identical
  deterministic state for every array task, and every "replicate" simulates
  and fits the exact same dataset.
- This is the SAME bug documented in README bug #4 ("Fixed-seed RNG leak",
  job 436723). It was fixed in `01e_run_ecoregion_sim.R` (line 94:
  `set.seed(DESIGN_SEED + (scenario=="null")*1000 + rep_id)`) but that fix was
  not carried into `01e_run_abundance_sweep.R`.

## Fix (one line)
In `01e_run_abundance_sweep.R`, immediately BEFORE the `run_one_replicate(...)`
call (must be AFTER `build_reduced_constants`, which resets the seed):

```r
set.seed(20260712 + row_id)   # row_id is unique per array task
```

Then re-submit the 180-task array and re-run `02e_collect_abundance_sweep.R`.
Each task is ~40 min; full array re-runs in ~1 h wall time.

## What the current run can and cannot say
INVALID until re-run (all rest on across-replicate variation / rates):
- "50% sign-recovery = chance" (it's one draw: 4/8 regions correct)
- "0/30 vs 30/30 WAIC preference"; "false-positive control unambiguous" (single
  WAIC comparison ×30; no false-positive RATE is measured)
- coverage percentages (100%, 62.5%)
- the "threshold, not a curve" shape (curve through 3 single-draw points)

Survives only as single-draw anecdote (directional):
- point estimates near-unbiased in every cell (|bias| ~0.002–0.019)
- sign-recovery on the single draws rose 4/8 → 4/8 → 7/8 with information
  (bobcat → moderate → deer-like) — direction sensible, magnitude unquantified.

## README edits needed after re-run
Section "Abundance × trend-scenario sweep": the result table, the
"single cleanest result" claim, the "threshold not a curve" headline, and both
"not fully understood" bullets must be recomputed from the reseeded run.
