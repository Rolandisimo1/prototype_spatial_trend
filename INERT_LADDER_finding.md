# The abundance ladder was inert — finding, scope, and fix

**Status:** confirmed independently from both code and output. Fix written and
verified; the affected 720 rows are preserved under tag `n30_4arm` and must be
re-run before any abundance-indexed claim is made.

## The defect

`01i_run_estimator_sweep.R` line 249 (pre-fix):

    truth <- scale_truth_abundance(truth, label = cfg$abundance)

`scale_truth_abundance()` does not look a label up. It only *records* one as an
attribute; all of its scaling comes from `occ_shift` and `count_log_mult`, both
of which default to `0`. So the call was a silent no-op, and the driver never
called `abundance_levels_default()` or `abundance_levels_measured()` anywhere.

The sibling `01e_run_abundance_sweep.R` (lines 103-106) does it correctly,
passing `abn$occ_shift` and `abn$count_log_mult`. This is a copy-omission, and
it is the third instance in this project of a defect that reached production in
one driver while a working version sat in a sibling script.

## Confirmation from the output, not only the code

If the ladder were live, `occ_shift` of 0 / 0.75 / 1.5 on the cloglog scale
implies lambda multipliers of 1.00 / 2.12 / 4.48. Observed, across the 540
array-arm rows:

| level | mean_prop_detect | ratio |
|---|---|---|
| bobcat_baseline | 0.1458 | 1.000 |
| moderate | 0.1435 | 0.985 |
| common_deerlike | 0.1444 | 0.991 |

Kruskal-Wallis across the three levels: `mean_prop_detect` p = 0.58,
`tvb_detected` p = 0.51, `tvb_bias` p = 0.52, `tvb_ci_width` p = 0.58,
`disc_auc` p = 0.98. No metric varies with abundance. The three levels are
i.i.d. replicate sets from one data-generating process, differing only by seed.

Claude Code found this via `lambda_true_med` in the `ncam_diag` side file, which
came back as 0.113748 to six decimals in all 180 camera_rn replicates. That
diagnostic was added for an unrelated purpose (checking N_cam runaway) and
caught this instead.

## Why it matters more than a lost axis

Abundance is the axis along which the two estimator families were *expected to
diverge*: occupancy's information about intensity dies as psi saturates, while
Royle-Nichols keeps reading detection frequency. Fisher information about
log(lambda) at J=4, r=0.1 crosses over near lambda≈5 and RN peaks near
lambda≈20, where occupancy's information is ~0.

Both RN arms were therefore never tested where they should win. Any statement
of the form "RN shows no advantage at deer-like abundance" is unsupported —
there is no deer-like cell in the data to show an advantage in.

## A third mismatch, which would have defeated a naive fix

The design levels were `bobcat_baseline` / `moderate` / `common_deerlike`.
These match **neither** ladder: `abundance_levels_default()` uses
`bobcat_baseline` / `moderate` / `common`, and `abundance_levels_measured()`
uses `bobcat_like` / `intermediate` / `deer_like`. Simply passing the label
through to a lookup would have failed on two of three levels.

## The fix

Applied in `01i_run_estimator_sweep.R`:

1. The design levels are now `abundance_levels_measured()$level` verbatim, so
   the lookup cannot silently miss.
2. The driver looks the level up and passes both scaling arguments, with a
   hard `stop()` if the level is absent and a second `stop()` if a
   non-baseline level resolves to a zero shift.
3. `abundance_levels_measured()` is used rather than
   `abundance_levels_default()`, because the measured ladder is anchored on
   real per-window detection counts from the raw camera data (1.555 for
   bobcat, 6.405 for white-tailed deer, geometric midpoint 3.156), whereas the
   default ladder's own docstring flags its multipliers as placeholders.

Verified: lambda multipliers now 1.00 / 2.12 / 4.48 and iNat count multipliers
1.00 / 2.03 / 4.12, matching the ladder exactly.

Six checks were added to `00y_smoke_test.R` — three static (the driver passes
both scaling arguments; the driver builds a ladder; every design level resolves
in it) and three behavioural (scaling moves lambda across levels; the
label-only call is a documented no-op). The smoke test gates the sweep, so a
recurrence stops one task rather than 720.

**Renaming the levels changes the sort order and therefore what each `row_id`
denotes.** Results from the inert-ladder run are not comparable row-by-row and
the sweep must be re-run whole.

## What the existing 720 rows still support

They are a single-abundance experiment at bobcat-like density with **90
replicates per arm** — larger and cleaner than the n=30 the design intended.
Every within-DGP comparison of estimators holds. Only abundance-indexed claims
are retracted.
