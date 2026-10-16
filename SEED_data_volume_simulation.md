# Seed: how much camera and iNaturalist data does a trend need?

Design notes for a future simulation study. Nothing here is run yet. Written
2026-10-13 from failure modes found in the real fleet fits, each of which the
existing estimator sweep cannot produce because it holds the relevant quantity
fixed.

The question: **at what volume of camera and iNaturalist data can a population
change of a given size be detected, and where does the integration stop
helping?**

---

## What the existing sweep varies, and what it holds fixed

`SEED_array_level_project.md` varies estimator arm (4), abundance level (3) and
spatial scenario (2), at 30-60 replicates. It holds fixed: data volume, range
size, detection probability at given abundance, the temporal window of each
stream, the size and growth of the opportunistic record pool, and camera site
composition over time.

Every failure mode below lives in something it holds fixed. That is why none of
them could have been anticipated from the sweep, and why they are worth a
dedicated design rather than more replicates of the old one.

---

## Tier 1 — the core factorial

### A. Camera volume x iNaturalist volume

The question as asked. A 2-D grid rather than a single axis, because the point
of integration is that the two streams substitute for each other, and the
substitution rate is the unknown. Power as the outcome, at 2-3 effect sizes
spanning what the fleet actually shows (roughly 2, 5 and 15 %/yr).

Anchor the grid on real numbers so cells are interpretable as "a project like
X": the fleet spans 1,806 to 26,517 camera deployments per species and 2.5 % to
80 % detection rates, with iNaturalist records from a few thousand to over a
million.

### B. Range size (number of spatial units)

**Strongest empirical result we have: Spearman -0.61, p=0.005 between range size
and how far the joint fit moves from the camera-only fit.** Species under 100
cells of 100 km shift a median 17.1 pp/yr; those above shift 1.7. The mechanism
is that the CAR fields are declared on the cell100 grid, so a 30-80 cell species
gives the spatial field almost no structure to identify, the iNaturalist side
sets the surface, and the camera trend absorbs the mismatch.

Simulate 25 / 50 / 100 / 300 / 900 occupied cells. Three real species sit below
the threshold, so the threshold's location is barely constrained by data and
this is the cheapest thing a simulation can pin down that observation cannot.

### C. Camera window as a fraction of the iNaturalist window

Not varied anywhere in the current design, and structural in the real data:
**every species has cameras 2019-2025 and iNaturalist 2008-2025.** The trend is
camera-anchored but the cameras see only the last 7 of 18 years.

Simulate camera coverage at 100 / 60 / 40 / 20 % of the opportunistic window,
including the realistic case where camera years sit at the END of the record.
Expect an interaction with A: more cameras over a short window may not substitute
for fewer over a long one.

---

## Tier 2 — mechanisms diagnosed in real fits, never simulated

### D. Detection probability separated from abundance

The current sweep varies abundance, which moves detection rate with it, so the
two are confounded by construction. In the fleet they are not: badger detects at
2.47 % and groundhog at 2.56 %, yet badger's camera interval is 138 pp wide and
groundhog's is 6.7. Body size, trail affinity, nocturnality and home-range size
relative to the detection zone all move detection rate at fixed abundance.

Vary per-visit detection probability at fixed abundance. This separates "few
animals" from "animals rarely seen", which the real data cannot.

### E. An opportunistic stream with its own temporal process

The most important gap, because it is the gap under the paper's headline claim.
**The sweep generates both streams from one truth, so stream disagreement can
only arise from noise.** In the fleet, 6 of 18 species have camera and
iNaturalist trends of opposite sign. The earlier note that the corroboration
indicator "cannot fire in any arm because the generating truth has the two
components in opposite signs" is the same gap seen from the other side.

Generate data where the opportunistic stream carries a temporal signal that is
NOT population change (observer growth, shifting taxon attention, platform
effects), and ask the question that matters: does the model route it into the
sink parameter, or does it contaminate the camera-anchored trend? Vary the size
of that spurious signal relative to the true trend.

---

## Tier 3 — realism refinements, one-at-a-time

### F. Record-pool growth and a mis-specified effort offset

Real growth was 2,798 to 369,355 records, roughly 130-fold, over the same period
the trend is estimated. The offset is meant to absorb that. Test the growth rate
at which it stops absorbing, and separately test an offset that is a noisy or
biased proxy for true effort rather than the exact one.

### G. Camera site turnover

Composition drift was diagnosed in the real data as able to inject a trend, and
the sweep has fixed sites. Vary the fraction of sites replaced per year, and
whether replacement is random or habitat-biased.

---

## Outcome measures

Report all of these per cell, not detection rate alone. **Detection rate
conflates bias with precision** and was already shown to favour the more biased
estimator at high abundance.

1. Bias in the camera-anchored trend parameter.
2. Interval width. Our real problem is 138-pp intervals, not wrong point estimates.
3. Coverage of the nominal interval.
4. **Shift between the joint fit and a camera-only fit of the same data.** This is
   the diagnostic now computed on every real fit, so simulating it makes the two
   directly comparable and calibrates what shift counts as large.
5. Whether the posterior of the congruence parameter is unimodal. Porcupine and
   badger each had two solutions, one with the streams decoupled, and the
   disconnected solution produced +41 %/yr for porcupine. Worth knowing what data
   conditions create that.

---

## Feasibility

Seven axes is not one factorial. Proposed shape: a full grid on Tier 1 (A x B x C)
at 2-3 effect sizes, then Tier 2 and 3 as one-at-a-time perturbations off a
single baseline cell chosen to match a mid-fleet species. Reuse the existing
driver and metrics code; the arms do not need to change, since this study varies
data rather than estimator.

Two lessons from the last sweep to carry over: confirm each axis is live at
runtime with a smoke test before the full launch, since the abundance ladder was
inert for a whole run; and define the null by zeroing the trend parameters, not
the spatial deviation.
