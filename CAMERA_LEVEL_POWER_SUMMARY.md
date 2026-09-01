# Camera-level trend-detection power — results summary

**Shareable summary for another conversation.** Self-contained; assumes no prior
context. Every number here is copied from a result file in
`prototype_spatial_trend/` — the source file is named for each.

---

## What was simulated

A NIMBLE integrated model fusing **camera-trap** and **iNaturalist** data to
estimate mammal abundance/occupancy trends across North America. Camera data
enters at the **camera level** (each camera an independent site, binary
detection per 10-day window, via `dOcc_v`); iNaturalist enters as an
effort-offset negative-binomial count at 50 km grid cells. Both streams inform
one shared latent intensity field.

The question: **can this model detect a spatially varying trend — specifically a
per-ecoregion (EPA/CEC Level I) deviation from the national trend — and how does
that depend on species abundance?**

**Design:** 3 abundance levels x 2 trend scenarios x 30 replicates = 180 fits.
Real camera locations, real iNat effort, real covariates; simulated detections
with known truth. The `varying` scenario has genuine per-ecoregion trend
deviations; the `null` scenario has a single national trend and no regional
structure. Abundance was scaled from a bobcat-calibrated baseline
(bobcat = widespread but rare; "deer-like" = most abundant fleet species).

Two structure-detection tests were run on the **same** 180 datasets: a WAIC
comparison (ecoregion vs national-scalar model) and a reversible-jump MCMC
Bernoulli inclusion indicator on the regional trend term.

---

## Result 1 — trend sign recovery works, and improves smoothly with abundance

Fraction of ecoregions whose estimated trend has the correct sign
(`varying` scenario; p-values are one-sample t vs 50% chance over 30 reps):

| abundance | sign recovery | p |
|---|---|---|
| bobcat-like (low) | **66.2%** | 6e-06 |
| ~3x bobcat | 72.5% | 4e-09 |
| deer-like (high) | 85.4% | 4e-17 |

**Significantly above chance at every level, including bobcat's real density.**
This is a gradual improvement with information, not a threshold effect.

**Important refinement:** broken out per region, recovery tracks the *magnitude*
of the true regional trend, not just abundance. Regions with a strong true trend
(|slope| ~0.20-0.29) recover at **80-100%** even at bobcat abundance, while
regions whose true trend is near zero (~0.04) sit at chance everywhere. That is
correct behaviour — you cannot recover the sign of a trend that is essentially
zero — and it means the aggregate 66% *understates* performance for the regions
that actually have something to detect.

Recovery is governed by **(trend magnitude x information)**, so a fleet screen
should be run per species against its own posterior rather than read off this
bobcat-anchored ladder.

*Source: `abundance_sweep_evaluation.md`, `abundance_sweep_summary.csv`,
`fig1_recovery.png`.*

## Result 2 — the estimator is well calibrated

| abundance | scenario | mean bias | 95% CI coverage |
|---|---|---|---|
| bobcat-like | null | -0.007 | 100% |
| ~3x | null | -0.003 | 100% |
| deer-like | null | -0.003 | 100% |
| bobcat-like | varying | -0.005 | 97.1% |
| ~3x | varying | -0.006 | 90.0% |
| deer-like | varying | -0.009 | 97.5% |

Bias is **< 0.01 logit points/year in all six cells**, and coverage is at or
above nominal almost everywhere. So where the model reports a regional trend,
the point estimate and its uncertainty can be trusted. The difficulty is purely
in resolving *sign* for weak-trend regions at low abundance — not in
overconfidence.

*Source: `abundance_sweep_evaluation.md`, `fig3_calibration.png`.*

## Result 3 — WAIC fails as a structure-detection gate

Fitting both models to each dataset and comparing WAIC, the null and varying
distributions overlap almost completely. Difference **non-significant at every
abundance level: p = 0.80 / 0.91 / 0.61.** Under the *null* scenario WAIC still
"prefers" the ecoregion model 12-18 times out of 30 — a coin flip.

**WAIC was dropped as the "does this species need a regional trend term" gate.**

*Source: `abundance_sweep_evaluation.md`, `fig2_waic.png`.*

## Result 4 — the RJMCMC indicator is better, but hits the same wall

WAIC's replacement: one Bernoulli indicator `gamma` on the regional trend
deviation, fit with `configureRJ`, following Goldstein et al. (bioRxiv
2025.01.17.633640). Slab SD fixed at 0.2, matching the true generating scale.
Convergence judged on R-hat of the data log-likelihood, **not** on `gamma`
(a Bernoulli node can sit at 0 or 1 for long stretches, so the standard
between-chain statistic is inappropriate for it).

Mean posterior inclusion probability:

| abundance | P(gamma) null | P(gamma) varying | separation |
|---|---|---|---|
| bobcat-like | 0.494 | 0.501 | **0.007** |
| ~3x | 0.447 | 0.592 | 0.145 |
| deer-like | 0.414 | 0.837 | **0.423** |

False-positive / false-negative rates at the 0.5 support threshold:

| abundance | FP (null) | FN (varying) |
|---|---|---|
| bobcat-like | 20.0% (6/30) | **70.0% (21/30)** |
| ~3x | 10.0% (3/30) | 33.3% (10/30) |
| deer-like | 16.7% (5/30) | **6.7% (2/30)** |
| *published benchmark* | *9.1%* | *14.9%* |

Only the deer-like level reaches or beats the benchmark. **At bobcat abundance
the false-negative rate is 70% — the test misses a real regional trend more
often than it catches it.** At the stricter 0.9 threshold, false positives
essentially vanish (0-3%) but false negatives climb to **100%** at bobcat and
moderate abundance — the same threshold/power trade-off the paper reported,
shifted unfavourably by bobcat's low information.

**It also mixes poorly where it is needed most:** the switch reached R-hat <= 1.1
in only **57%** of bobcat-abundance replicates, vs 83% (moderate) and 92%
(deer-like). A 57% convergence rate is itself a warning that the test is at the
edge of feasibility there.

**One correction to the source memo.** `indicator_test_evaluation.md` states
that filtering to converged replicates "barely changes the rates." That holds
for the two higher abundance levels but **not for bobcat**, where restricting to
converged replicates makes the false-positive rate substantially *worse*
(re-read from `indicator_test_summary.csv`):

| abundance | FP all reps | FP converged only | FN all reps | FN converged only |
|---|---|---|---|---|
| bobcat-like | 20.0% | **37.5%** | 70.0% | 61.1% |
| ~3x | 10.0% | 11.1% | 33.3% | 30.4% |
| deer-like | 16.7% | 17.9% | 6.7% | 3.7% |

So at bobcat abundance the converged subset gives a slightly better FN but a
**nearly doubled FP** — the test is not merely underpowered there, it is
unstable, and the well-mixed replicates are not a clean subsample. Treat the
bobcat row as "unusable" rather than "usable with a wider interval." This
strengthens conclusion 3 below rather than changing it.

*Source: `indicator_test_evaluation.md`, `indicator_test_summary.csv`,
`fig4_indicator.png`, `fig5_convergence.png`.*

---

## The bottom line

1. **Trend estimation is sound and honest.** Near-zero bias, at-or-above-nominal
   coverage at every abundance level.
2. **Regional sign recovery beats chance everywhere** (66% -> 85%), and reaches
   80-100% for regions with a genuinely strong trend.
3. **No method tried can decide whether regional structure exists at bobcat's
   real data density.** WAIC fails at every level; the indicator fails at low
   abundance (70% FN, 57% convergence). For a bobcat-abundance species the
   honest output is **"cannot determine from data"** — abstain, not a false yes
   or no.
4. **The camera stream is the binding constraint.** Everything above improves
   with abundance, which is why alternative camera estimators are worth testing.

**Practical consequence:** report `P(gamma=1)` with 0.5/0.9 thresholds rather
than a WAIC delta; gate fleet screening on **both** the indicator and achieved
convergence, since a species that cannot mix the switch is a species where the
term cannot be evaluated at all.

---

## Caveats that limit how far these numbers travel

- **Correctly-specified-model simulation.** Data were generated by the same
  structure being fitted, so real data will be worse. These are upper bounds.
- **The indicator FN rates are a best case.** The slab SD was fixed at the
  *true* trend scale, giving the model oracle knowledge of effect magnitude. A
  species whose trend scale differs would do worse. A diffuse-slab sensitivity
  run (sd = 0.4) was deferred and is worth doing before relying on these rates
  fleet-wide.
- **An earlier version of this sweep was retracted.** The first run lacked a
  per-replicate reseed, so all 30 "replicates" per cell simulated and fit the
  identical dataset — effectively n=1. Every rate and coverage claim from it was
  withdrawn. **The numbers above are from the corrected re-run** (180/180 ok,
  per-cell values 30/30 unique, collector degeneracy guard passed). If you see
  different figures attributed to this sweep elsewhere, they are from the broken
  run. See `abundance_sweep_seed_diagnosis.md`.
- **These are simulation results only.** They say nothing about the real-data
  fits in the parent project, which are separately affected by a chunked-MCMC
  checkpoint defect. The simulation framework uses single-call MCMC with no
  chunk-resume and is unaffected.
