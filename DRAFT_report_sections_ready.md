## Two model structures, and why we fit both

For each species we fit the same underlying model twice, in two different
spatial structures. Both structures combine the same two data sources —
opportunistic iNaturalist sightings, which cover the whole country but are
biased by where people go and what they choose to photograph, and camera-trap
detections, which are unbiased about presence/absence at a site but limited
to wherever a camera happens to be. Both models estimate one thing we care
about most: whether the species' population trend, over time, is increasing,
decreasing, or flat.

**The national_scalar model** assumes the trend is a single number for the
entire country: one slope, applied everywhere. It is the simpler of the two,
and it borrows the most statistical strength across the full dataset, which
is why it converges fastest and gives the tightest confidence intervals.

**The ecoregion model** relaxes that assumption. It lets each of the eight
EPA Level I ecoregions in the study area — Eastern Temperate Forests, Great
Plains, Northern Forests, and so on — have its own deviation from the
national trend (`year_region[r]`, one value per region), with those
deviations themselves drawn from a shared distribution whose spread is a
single estimated parameter (`sigma_region`). This is a standard technique
called partial pooling: a data-rich region can pull its own number away from
the national average, while a data-poor region's number stays close to the
national average by default (because there is little region-specific
evidence to overrule it). It answers a genuinely different question than
the national model does — not "what is the trend everywhere" but "does the
trend differ meaningfully by region, and if so, where."

We fit both because they serve different purposes. The national_scalar model
is the number to quote if someone asks "is this species increasing or
decreasing across its range." The ecoregion model is the number to quote if
someone asks about a specific region — including, for readers focused on one
state, reading that state's ecoregion(s) out of this model rather than
re-fitting the model on that state alone (which would answer a different,
less-informed question; see the North Carolina discussion earlier in this
project's record for why).
## What "the CAR fields fail the whole-fit gate" means, and why it does not undermine the trend numbers

Every ecoregion model contains two kinds of parameter that behave very
differently, and telling them apart is the key to understanding this
caveat.

**Trend parameters** — `total_var_beta`, `year_beta`, `year_var`,
`sigma_region`, and each region's `year_region[r]` — are single numbers,
each one estimated from data pooled across the entire country or an entire
region (thousands of camera-nights and iNat records). With that much
supporting data, these parameters converge quickly and reliably: running
three independent MCMC chains from different starting points, they land in
essentially the same place (formally, their Gelman-Rubin R-hat statistic,
the standard convergence check, is close to 1.0, meaning the chains agree).

**CAR fields** — `link_occ_intercept` and `MCMT_effect` — are a different
kind of object entirely. Instead of one number, each is a value assigned to
*every one of the roughly 900 grid cells* covering the study area (a "CAR
field," short for Conditional AutoRegressive field: a common spatial
statistics tool where each cell's value is modeled as similar to its
immediate map neighbors, so the whole field is spatially smooth rather than
cell-by-cell independent). Each individual cell's value is estimated mostly
from whatever data fall in that one cell and its neighbors — far less data
per parameter than any trend parameter gets.

**"Fails the whole-fit gate" means:** when we check convergence for every
parameter in the model — all ~2,750 of them — a small, localized patch of
CAR-field cells (roughly a dozen to a few dozen, always the same patch, in
the corner of the map where a data-rich area sits directly next to a
data-empty one) does not converge cleanly, even though every trend
parameter does. The three MCMC chains land in different places for those
specific cells, because there is a real, hard-to-resolve statistical tension
at exactly that spot — but it does not spread to, or contaminate, the
trend estimates, which draw on a much larger, unrelated pool of data.

**What this means in practice for this report:** every national and
regional trend number quoted here is trustworthy on its own terms. What is
*not* trustworthy is the fine-grained spatial intercept map and the winter-
temperature climate-response map (MCMT) at that specific patch of cells —
those should not be read at all in that location, on either species, and
any figure showing them should grey the patch out rather than show a
number that might be an artifact of unresolved chain disagreement.
## Convergence: what converged, what did not, and one number that looks alarming but is not

We ran three independent MCMC chains per model and assessed convergence with
the Gelman-Rubin statistic (R-hat), which compares variance between chains to
variance within them; values near 1.0 indicate the chains have found the same
answer from different starting points.

All four species-tracks converged on every parameter that carries a population
trend. Across the eight fits, the largest R-hat among the trend parameters
(`total_var_beta`, `year_beta`, `year_var`, `sigma_region`, `year_region[1..8]`)
was 1.086, and all 72 occupancy-covariate coefficients across all eight fits
were at R-hat ≤ 1.024.

**One apparent exception, which is not one.** Bobcat's national_scalar fit shows
R-hat = 1.115 on `trend_robust_indicator` — above the conventional 1.1
threshold. This is an artifact of what that quantity is rather than evidence of
a failed fit. `trend_robust_indicator` is a derived probability, not a sampled
parameter, and bobcat's value is 0.0005: pinned against the lower boundary of
[0, 1]. R-hat is poorly behaved for a quantity compressed against a boundary,
because between-chain variance in the last few decimal places is compared
against a within-chain variance that is nearly zero. Every actual sampled trend
parameter in that same fit is at R-hat ≤ 1.0006. The fit is sound; the
indicator's R-hat is not interpretable.

**Separately, and more consequentially:** both ecoregion fits for bobcat and
white-tailed deer fail a whole-model convergence gate on the spatial CAR fields
only — a localized, contiguous patch of grid cells, with 11 of 12 cells shared
between the two species. Because the same patch fails independently in two
species' fits, this is a property of the spatial adjacency graph in that region
rather than of either species' data. See the CAR-field section for what this
does and does not invalidate.

## Occupancy covariate effects

Each model estimates nine covariate effects on occupancy, on the log-odds
scale, so a negative coefficient means a site is less likely to be occupied as
that covariate increases. Nine rather than ten: the three soil-texture
fractions (clay, silt, sand) are compositional and sum to a constant, making
them mutually redundant, so `soil_sand` is dropped fleet-wide and acts as the
implicit reference level.

Counting effects whose 95% credible interval excludes zero, bobcat has 5 of 9,
moose 4 of 9, and white-tailed deer 3 of 9.

**Directions, ordered from most negative to most positive:**

- **Bobcat** — human population −0.38, soil clay −0.13, terrain ruggedness
  +0.07, vegetation greenness (NDVI) +0.12, deciduous forest +0.22.
- **Moose** — human population −3.77, vegetation greenness −0.49, mixed forest
  +0.51, evergreen forest +1.02.
- **White-tailed deer** — human population −0.50, terrain ruggedness −0.04,
  evergreen forest +0.14.

Two patterns are worth naming. First, **human population is negative and
statistically clear for all three species**, and is the only covariate that is
significant in every fit — but its magnitude for moose (−3.77) is roughly
eight-fold larger than for bobcat or deer, which is why the three panels of the
covariate figure use different x-axis scales. Moose avoidance of settled
landscapes is by far the strongest single covariate effect anywhere in this
analysis. Second, **the forest-type coefficients separate the species in an
ecologically sensible way**: moose respond strongly and positively to evergreen
and mixed forest, bobcat to deciduous forest, and deer weakly to evergreen —
consistent with each species' known habitat associations, which is a useful
sanity check that the covariate model is picking up real structure rather than
noise.

**A robustness result worth stating:** across all 36 species-by-covariate
combinations, the national_scalar and ecoregion parameterizations agree on the
sign of every single coefficient (36/36), with a maximum disagreement in
magnitude of 0.15. Covariate inference is therefore insensitive to which
spatial trend structure is assumed.

## iNaturalist-versus-camera congruence, measured the way the source method measures it

The model links the two data streams through
`log(mu) = theta0 + theta1 * log(sum(lambda))`, where `mu` is expected
iNaturalist count intensity and `lambda` is camera-derived abundance. The
exponent `theta1` is therefore the scaling relationship between the two data
types, and it is the quantity the original method (Goldstein et al., bioRxiv
2025.01.17.633640) uses to characterise how far the datasets correspond. We
report it in preference to an information-criterion measure such as WAIC
because it answers the question directly and is comparable to the source
method's own reported diagnostic.

**Every fit is sub-linear.** `theta1` ranges from 0.41 (white-tailed deer) to
0.61 (bobcat), and in all eight fits the 95% credible interval excludes 1.0.
Proportional scaling — where doubling true abundance doubles iNaturalist
counts — is therefore rejected everywhere. iNaturalist counts rise with
abundance, but with diminishing returns, which is what one would expect if
observer attention saturates: the tenth deer in a neighbourhood is much less
likely to be photographed and uploaded than the first.

**An unexplained result we are reporting as unexplained.** `theta1` and the
trend robustness indicator rank the four species-tracks in exactly inverse
order (Spearman ρ = −1.00; exact two-sided permutation p = 0.083 at n = 4).
Bobcat has the highest cross-sectional congruence and the lowest temporal
corroboration; white-tailed deer the reverse.

The two metrics measure different things, and that alone means the pattern is
not a contradiction: `theta1` asks whether iNaturalist counts track camera
abundance across space at a point in time, while the trend indicator asks
whether the two streams agree about change over time. A species can score high
on one and low on the other without inconsistency. But mere distinctness would
predict a scatter, not a perfect monotone reversal across every fit, so this
explanation accounts for the absence of a contradiction without accounting for
the pattern actually present. We tested the obvious alternative — that `theta1`
is attenuated toward zero where abundance has little spatial contrast — and it
is refuted, and refuted backwards: white-tailed deer has the *largest* spatial
contrast in log counts (sd 1.79) and the *lowest* `theta1`, while bobcat has the
smallest contrast (1.32) and the highest. `theta1` is also not monotone in
record count or in occupied-cell count.

With four fits across three species, ρ = −1.00 is not strong evidence of a
mechanism, and we decline to assert one. We flag it as an open question worth
revisiting when a fourth and fifth species enter the fleet, at which point the
pattern will either strengthen into something requiring explanation or dissolve.

## Simulation study: which estimator to trust as abundance rises

Occupancy models are known to lose information once a species is common enough
that nearly every site is occupied — the detection history saturates, and
presence/absence stops discriminating. Royle-Nichols models instead use
detection *frequency*, which continues to carry information past that point.
Because this project's species span a wide abundance range, we tested how much
this matters using simulated data with a known true trend.

**Design.** Synthetic datasets were generated at three abundance levels
anchored to real measured detection rates — bobcat-like (low), intermediate,
and deer-like (high) — and fit with four estimators: occupancy and
Royle-Nichols, each at the individual-camera and the camera-array level. Thirty
replicate datasets per abundance level per estimator. The true trend was
−0.1795 throughout.

**Result: the difference between these estimators is accuracy, not detection.**
Both occupancy estimators become steadily more biased as abundance rises. At
deer-like abundance, camera occupancy overstates the true decline by roughly
3.6-fold and array occupancy by roughly 5.6-fold, while both Royle-Nichols
estimators remain within a few percent of the true value at every level.
Royle-Nichols credible intervals are also 3- to 9-fold narrower. Occupancy's
interval coverage degrades from near-nominal at low abundance to 0.80 at
deer-like abundance, while Royle-Nichols coverage stays in the 0.83–0.97 band
without a monotone decline.

**The trap this creates.** At deer-like abundance camera occupancy flagged the
declining trend more often than camera Royle-Nichols did (66.7% versus 53.3%),
despite its point estimate being 3.6-fold too negative. A biased estimator can
appear to detect a trend more reliably than an unbiased one, because its own
overconfident intervals make rejecting the null easier. Any estimator
comparison that ranks methods by how often they detect a trend will therefore
prefer the more biased method at high abundance. Detection rate is the wrong
criterion; accuracy is.

**One diagnostic never fired.** The camera-corroboration indicator was zero in
every arm, at every abundance level, in all 720 simulated datasets. This is
consistent with the low values seen in the real fits, and suggests the
indicator as currently defined is conservative to the point of rarely
triggering rather than that camera data never corroborates.

> **PENDING — false-positive rate.** The sweep above correctly measures bias,
> precision, and coverage, but a separate defect meant its "null" (no-trend)
> scenario never actually set the trend to zero — it removed only the spatial
> deviation field, leaving the true national trend present in both scenarios.
> No arm's false-positive rate has therefore been measured, and an earlier
> "false-positive exceeds power" finding derived from it has been retracted. A
> corrected re-run with the trend genuinely zeroed is in progress. Do not read
> the absence of a false-positive number here as evidence that it is low; it
> has simply not been measured correctly yet.

> **PENDING — 5-year and 10-year trend windows.** All trend figures in this
> report are 18-year (2008–2025). Camera-trap data begins in 2019, so the
> 18-year window is carried substantially by iNaturalist alone in its early
> years and is not the intended headline. Windowed refits over 2021–2025 and
> 2016–2025 are the intended primary results; they require input bundles
> rebuilt on the current mask and covariate design and have not yet been run.
## Data provenance caveat — moose v2b numbers are from a superseded extraction

Not all eight fits are of equal evidential standing, and the report must say so.

Bobcat, white-tailed deer, and moose_v1fix9 posteriors were all extracted on
2026-08-28 from **single continuous MCMC chains**, which structurally cannot be
affected by the chunked-checkpoint burn-in defect (Issue 4) that contaminated
the earlier fleet — resumed chunks skipped re-adaptation burn-in, inflating
credible intervals and distorting R-hat at every resume boundary.

**`moose_v2b`'s posteriors were never re-extracted from single-shot chains.**
The only available `moose_v2b` extraction dates to 2026-08-20, in the chunked
era. Its numbers are therefore subject to the Issue 4 defect, and two
extractions from that date disagree slightly with each other
(`trend_robust_indicator` 0.2448 versus 0.2551; `total_var_beta` −0.1162 versus
−0.1192), consistent with different chain subsets being pooled.

Practical consequence: the moose_v2b-versus-moose_v1fix9 comparison — the whole
point of running two moose tracks, to separate the range-mask fix from the
covariate fix — is currently **not a clean comparison**, because the two tracks
come from different MCMC regimes as well as different designs. A single-shot
re-extraction of moose_v2b is needed before that comparison can carry weight.
Until then, moose_v1fix9 is the more trustworthy of the two moose tracks on
provenance grounds alone.
