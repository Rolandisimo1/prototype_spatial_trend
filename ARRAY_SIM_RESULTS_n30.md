# Array-level estimator sweep (n=30) — results

**Run:** 540 array tasks, all completed, none failed. Aggregated locally on
2026-08-27 from `hazel_pull_20260827/estimator_sweep_n30/rows/` (540 `.rds`).
Supersedes the n_rep=3 pilot.

**Design:** 3 estimators x 3 abundance levels x 2 scenarios x 30 replicates.
Design coverage verified complete: every one of the 18 cells has exactly 30 OK
rows.

---

## READ THIS FIRST — two things that change how the table is interpreted

**1. The "null" scenario is not a no-trend scenario.** The driver builds it as
`make_true_year_effect(amplitude = 0)`, which zeroes the *regional* deviations.
The **global** trend is untouched: `tvb_true = −0.1795 in BOTH scenarios**
(verified — identical single value in all 540 rows).

`tvb_detected` is "the credible interval on `total_var_beta` excludes zero." So
in the null scenario it is testing a global trend that **is genuinely present**.
The collector's own header labels this column the false-positive rate under
null. **It is not.** Under this design both scenarios measure power on the same
global trend; what differs between them is only whether regional structure
exists. There is no false-positive rate in this sweep.

This does not invalidate anything — but any statement of the form "estimator X
has a Y% false-positive rate" from this run would be wrong, and the pilot
summary carries the same mislabel.

**2. Statistical resolution.** 30 draws gives a rate a 95% interval of roughly
±0.18. **Differences between estimators smaller than about 18 percentage points
are not resolvable by this design.** Almost every difference below is smaller
than that.

**Degeneracy guard: passed.** 2–5 of 16 numeric metrics are constant within a
cell, and all are constant by construction (`tvb_true`, `n_arrays`,
`tri_fires`, and `tvb_covered`/`mean_prop_detect` in 6 cells). The reseed is
working; this is not the n=1 failure of the earlier abundance sweep.

---

## Results — varying scenario (regional structure present)

| estimator | abundance | power | bias | CI width | coverage | AUC | Spearman | field RMSE |
|---|---|---|---|---|---|---|---|---|
| camera_occ | bobcat-like | **0.167** | −0.122 | 2.32 | 0.967 | 0.637 | 0.337 | **2.23** |
| array_occ | bobcat-like | **0.000** | **−0.410** | **5.42** | 1.000 | **0.697** | **0.419** | 17.07 |
| array_rn | bobcat-like | **0.167** | +0.109 | **1.05** | 0.867 | 0.631 | 0.302 | 3.16 |
| camera_occ | moderate | 0.067 | −0.001 | 1.68 | 1.000 | 0.661 | 0.366 | 2.03 |
| array_occ | moderate | 0.033 | −0.237 | 5.09 | 1.000 | 0.654 | 0.316 | 17.13 |
| array_rn | moderate | 0.067 | +0.007 | 1.04 | 0.967 | 0.628 | 0.295 | 3.00 |
| camera_occ | deer-like | 0.167 | −0.148 | 1.74 | 1.000 | 0.672 | 0.400 | 1.83 |
| array_occ | deer-like | 0.033 | +0.076 | 4.94 | 0.967 | 0.661 | 0.333 | 16.04 |
| array_rn | deer-like | 0.133 | +0.045 | 0.96 | 0.933 | 0.636 | 0.328 | 3.12 |

## The headline: nobody has power to detect the global trend

**Power never exceeds 0.167 in any arm at any abundance.** The best case is
17% — 5 of 30 replicates. The true global trend is −0.18 on the link scale, and
essentially no estimator's credible interval excludes zero.

This is the dominant result and it swamps every estimator comparison. **The
answer to "do array-level estimators have the statistical power to detect
trends?" is: no — but neither does the camera-level estimator on this same
design.** No arm is usable as a trend detector here.

**Power does not increase with abundance in any arm** — camera_occ goes
0.167 → 0.067 → 0.167, array_rn 0.167 → 0.067 → 0.133. That flatness is itself
informative: it says the binding constraint is not information about the species
but something structural in the design (see "what to check next").

## Estimator differences that are real

Two differences are large enough to exceed the ±0.18 resolution floor or are
consistent across all three abundance levels:

**Array occupancy is badly behaved.** CI width **4.9–5.4** versus 1.0–1.7 for the
other two — roughly 3–5x wider — with bias reaching **−0.41** at bobcat
abundance and field RMSE of **16–17** against 1.8–3.2. Its coverage is 0.97–1.00,
but that is coverage bought with an interval so wide it is uninformative. This
is consistent across all three abundance levels, so it is not noise.

**Array RN is the best-calibrated arm.** Narrowest intervals (0.96–1.05),
smallest bias (+0.007 to +0.109), coverage 0.87–0.97. If any array-level
estimator is worth pursuing it is RN, not occupancy — which inverts the
expectation that occupancy would be the safer choice.

**Camera-level retains the best spatial discrimination.** Field RMSE 1.8–2.2 vs
3.0–3.2 (RN) and 16–17 (array_occ). AUC is statistically indistinguishable
across arms (0.63–0.70, all weak).

**`trend_robust_indicator` fires in 0 of 540 replicates** — zero across every
estimator, abundance, and scenario. No arm produces a camera-corroborated trend
under any condition tested. This directly matches moose's real-data behaviour
(indicator 0.24–0.26, failing the 0.5 threshold), so the simulation reproduces
the real failure mode.

## Comparison to the camera-level baseline

The earlier camera-level sweep asked a **different question** and the numbers
are not directly comparable:

| | camera-level sweep | this sweep |
|---|---|---|
| quantity | per-ecoregion trend **sign recovery** | **global** trend detection (CI excludes 0) |
| result | 66% → 85% with abundance, above chance everywhere | power ≤ 0.167, flat in abundance |
| structure test | RJMCMC indicator, FN 70% at bobcat | not run |

So this sweep does **not** show that array-level aggregation destroyed the 66%
sign recovery — it never measured sign recovery. What it shows is that on the
*global* trend metric, all three arms including camera-level are near-powerless.

**The one genuinely comparable statement:** the camera-level arm in this sweep
(power 0.067–0.167) performs no better than the array arms, so **on this metric
array aggregation costs nothing measurable.** That is a real, if modest, point
in favour of array-level — the Brazil-derived worry that aggregation would
reduce power is not confirmed here.

## What to check before drawing conclusions

The flat, uniformly low power across all arms *and* abundance levels is the
signature of a design problem rather than an estimator ranking. Three candidates,
in order of my suspicion:

1. **Is the global trend detectable in principle here?** With `tvb_true` =
   −0.18 and typical CI widths of 1.0–1.7, an interval would need to be ~10x
   narrower to exclude zero. The reduced simulation design may simply be too
   small for this effect size — which would make "power ≈ 0" a statement about
   the design, not the estimators.
2. **Array occupancy's 5x CI inflation needs diagnosing** before it is reported
   as an estimator property. `mean_prop_detect` is 0.13–0.15, so the graded
   replicate structure is present as intended, and `n_arrays` = 97 is constant
   across all cells — worth confirming 97 arrays is the intended reduced-design
   size and not an unintended truncation.
3. **Sign recovery was never computed.** The metric that made the camera-level
   sweep informative is absent here. Adding it to the collector is cheap and is
   the only way to make the two sweeps genuinely comparable.

**Recommendation: do not use this run to choose an estimator.** Use it for the
three things it does establish — array_occ is poorly behaved, array_rn is
well-calibrated, and array aggregation costs no measurable power versus
camera-level — then fix the metric mismatch and re-run before ranking anything.

---

## Provenance and caveats

- All 540 tasks OK; design coverage complete; degeneracy guard passed.
- Aggregated with `01j_collect_estimator_sweep.R` locally, NIMBLE not required
  for aggregation. Outputs: `estimator_sweep_n30_rows.csv`,
  `estimator_sweep_n30_summary.csv`.
- **Reading the CSVs in pandas requires `keep_default_na=False,
  na_values=['NA','']`** — the literal scenario string `"null"` is otherwise
  parsed as missing, which silently empties half the dataset. This bit us once
  before on the abundance sweep.
- Correctly-specified-model simulation: data generated by the structure being
  fitted, so real data will be worse. The RN arm's known weakness (unmodelled
  within-array detection heterogeneity) is absent by construction, so its
  performance here is optimistic.
