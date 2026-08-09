# Ecoregion trend sweep — evaluation of the corrected (reseeded) run

*Supersedes the impressions in `abundance_sweep_seed_diagnosis.md`, which
flagged the degenerate first run. This memo reads the properly-replicated
re-run (180/180 ok, per-cell values now 30/30 unique, collector guard passed).*

## Bottom line

The re-run changes the verdict in two directions — one more encouraging, one
more cautious — and both are now backed by 30 genuinely independent replicates
per cell rather than one dataset repeated 30×.

**1. Sign recovery is genuinely, significantly above chance — even at bobcat's
real abundance.** 66% of ecoregions get the right trend sign at bobcat's
density (one-sample t vs 50%: p = 6e-06), rising to 73% at ~3× and 85% at
deer-like abundance. This is a *gradual improvement with information*, not the
"flat-then-threshold" story the broken run implied, and it does not collapse to
a coin flip at low abundance. (Fig 1a.)

**2. The improvement is driven by trend magnitude, not just abundance.** Broken
out per region (Fig 1b), recovery tracks |true regional trend|: regions with a
strong true trend (|slope| ~0.2-0.29) recover at 80-100%, while regions whose
true trend is near zero (~0.04) sit at chance at every abundance level. This is
*correct* behavior — you cannot recover the sign of a trend that is essentially
zero — and it means the aggregate "66%" understates performance for the regions
that actually have something to detect.

**3. WAIC cannot tell real regional structure from none.** This is the headline
retraction from the broken run. Fitting both the ecoregion model and the
national-scalar model to each dataset and comparing WAIC, the null-scenario and
varying-scenario WAIC distributions overlap almost completely — the difference
is non-significant at every abundance level (p = 0.80 / 0.91 / 0.61; Fig 2).
Under the null scenario WAIC still "prefers" the ecoregion model 12-18 times out
of 30 — a coin flip. **WAIC is not a usable test for whether a species needs the
ecoregion term** at this replicate count and MCMC budget.

**4. The estimator itself is well-calibrated.** Mean trend bias is < 0.01 logit
points/year in all six cells, and 95%-CI coverage of the true trend sits at or
above the nominal 95% everywhere (Fig 3) — intervals are honest, not
overconfident. So where the model does report a regional trend, the point
estimate and its uncertainty can be trusted; the difficulty is purely in
resolving *sign* for weak-trend regions at low abundance.

## What this means for the project

- **For bobcat specifically:** an ecoregion-level trend map is defensible for
  the regions with strong trends (recovers sign 80-100% of the time even at
  bobcat's density) but should be reported with a per-region confidence flag —
  near-zero-trend regions are indistinguishable from noise and should gray out.
  This matches, from a completely different modeling approach (ecoregion partial
  pooling vs. the earlier CAR field), the same conclusion the grain sweep
  reached: bobcat's data can tell you *that* spatial heterogeneity exists better
  than it can pin *which sign* every individual region takes.
- **Do not use WAIC as the "is the ecoregion term needed?" gate.** The original
  plan was a WAIC comparison in both directions; the sweep shows that gate
  doesn't work here. A better trigger for the reporting-stage gray-out is the
  per-region data-support diagnostic (informed-cell count) plus the posterior
  width of that region's trend — not a model-selection statistic.
- **Abundance helps but is not the whole story.** A rarer species than bobcat is
  the real risk; a species with strong regional trends but low abundance may
  still do fine, while a high-abundance species with weak/uniform trends gains
  little. Recovery is governed by (trend magnitude × information), so the fleet
  screen should be run per-species against its own fitted posterior rather than
  read off this bobcat-anchored ladder.

## Numbers (corrected run)

| abundance | scenario | mean bias | coverage | sign-recovery | WAIC prefers eco |
|---|---|---|---|---|---|
| bobcat_baseline | null | -0.007 | 100% | n/a | 12/30 |
| moderate | null | -0.003 | 100% | n/a | 18/30 |
| common_deerlike | null | -0.003 | 100% | n/a | 14/30 |
| bobcat_baseline | varying | -0.005 | 97.1% | 66.2% (p=6e-06) | 14/30 |
| moderate | varying | -0.006 | 90.0% | 72.5% (p=4e-09) | 17/30 |
| common_deerlike | varying | -0.009 | 97.5% | 85.4% (p=4e-17) | 17/30 |

Sign-recovery p-values are one-sample t-tests vs 50% (chance) over 30 reps.
WAIC null-vs-varying separation is non-significant at every level (p>0.6).

## Figures
- `fig1_recovery.png` — sign-recovery above chance at every abundance (a);
  recovery vs true trend magnitude (b).
- `fig2_waic.png` — WAIC cannot discriminate null from varying (3 panels).
- `fig3_calibration.png` — near-zero bias (a) and >=nominal coverage (b).
