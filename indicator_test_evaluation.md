# RJMCMC indicator test — evaluation

The RJMCMC switch model (one Bernoulli indicator gamma on the ecoregion trend,
fixed slab SD = 0.2 matching the true generating scale, fit with configureRJ)
was run on the same 180 simulated datasets as the abundance sweep. Convergence
judged on R-hat of the data log-likelihood (NOT on gamma). 180/180 tasks ran;
per-cell gamma values are non-degenerate (28-30 unique of 30).

## Headline

The indicator test is a usable detector of "does this species need a regional
trend term" — but ONLY at higher abundance. At bobcat's real abundance it
cannot distinguish real regional structure from none, and it also mixes poorly
there. Detection power is governed by information (abundance), exactly like the
sign-recovery result.

## The separation grows with abundance (Fig 4a)

Mean posterior inclusion P(gamma=1), null vs varying scenario:

| abundance | P(gamma) null | P(gamma) varying | separation |
|---|---|---|---|
| bobcat_baseline | 0.494 | 0.501 | **0.007** |
| moderate | 0.447 | 0.592 | 0.145 |
| common_deerlike | 0.414 | 0.837 | **0.423** |

At bobcat abundance the switch sits at ~0.5 (its prior) whether or not a real
trend is present — no information to move it. At deer-like abundance it clearly
lifts toward 1 under a real trend and drifts below 0.5 under the null.

## False-positive / false-negative rates (Fig 4b, threshold 0.5, all reps)

| abundance | false positive (null) | false negative (varying) |
|---|---|---|
| bobcat_baseline | 20.0% (6/30) | 70.0% (21/30) |
| moderate | 10.0% (3/30) | 33.3% (10/30) |
| common_deerlike | 16.7% (5/30) | **6.7% (2/30)** |

Paper benchmark (single effect, threshold 0.5): FP 9.1%, FN 14.9%. Only the
deer-like level reaches or beats the benchmark FN. At bobcat abundance the
false-negative rate is 70% — the test misses a real regional trend more often
than it catches it. At the strong 0.9 threshold, false positives essentially
vanish everywhere (0-3%) but false negatives climb to 100% at bobcat and
moderate — the same threshold/power trade-off the paper reported, shifted
unfavorably by bobcat's low information.

## Convergence caveat (Fig 5)

RJMCMC mixing (R-hat of data log-lik <= 1.1) was reached in only 57% of
bobcat-abundance replicates, vs 83% (moderate) and 92% (deer-like). So at low
abundance the switch is both underpowered AND hard to fit. Filtering to
converged replicates only (summary rows with convergence_filter =
rhat_logLik_le_1.1_only) barely changes the rates, so the conclusion is not an
artifact of including poorly-mixed chains — but the low convergence rate at
bobcat abundance is itself a warning that this test is operating at the edge of
feasibility there.

## Comparison to WAIC (the point of this exercise)

WAIC could not distinguish null from varying at ANY abundance (p>0.6
everywhere). The RJMCMC indicator CAN distinguish them at moderate-to-high
abundance and gives an interpretable, thresholded FP/FN answer comparable to a
published benchmark. So the indicator is a strictly better gate than WAIC — but
it inherits the same fundamental limit: at bobcat's real data density there is
not enough information to reliably decide whether a regional trend term is
warranted, by any method tried.

## Implications for the project

- The indicator test REPLACES WAIC as the "is the ecoregion term needed"
  screen. Report P(gamma=1) with the 0.5/0.9 thresholds, not a WAIC delta.
- For a bobcat-abundance species the honest answer is "cannot determine from
  data" -- the gate should return abstain/undetermined, not a false yes or no.
  This is consistent with every other result: bobcat is the low-information
  stress case.
- Fleet screening should run this test per species and gate on BOTH the
  indicator AND achieved convergence -- a species that can't mix the switch is
  a species where the term can't be evaluated.
- The FN rates are a best case: the slab SD was fixed at the TRUE trend scale,
  giving the model oracle knowledge of the magnitude. A real species whose
  trend scale differs would do somewhat worse. A diffuse-slab sensitivity run
  (sd=0.4) was deferred; worth adding before relying on these rates fleet-wide.

## Figures
- fig4_indicator.png -- (a) P(gamma=1) distributions null vs varying by
  abundance; (b) FP/FN vs abundance at threshold 0.5 with paper benchmark.
- fig5_convergence.png -- RJMCMC convergence rate by abundance.
