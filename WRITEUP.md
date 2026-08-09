# Prototyping a spatially-varying occupancy trend in the bobcat iSDM

## Motivation

The integrated SDM currently represents the annual occupancy trend as a single
national scalar: one number for how bobcat occupancy is changing per year,
shared identically across the entire country. That's almost certainly hiding
real regional structure — populations plausibly increasing in some parts of the
range and declining in others, not moving uniformly. We prototyped an extension
that lets the trend vary spatially, and validated it via simulation before
touching the real, already-converged bobcat fit.

## What changed

The model already has spatially-varying coefficients (SVCs) for two temperature
covariates (`MWMT_effect`, `MCMT_effect`), each a conditional-autoregressive
(CAR) field over a ~900-cell, 100km grid. We added a third CAR field,
`year_effect`, using the identical template, representing a zero-mean spatial
*deviation* from the national trend. The national trend (`year_beta`) stays as
the model's global mean; `year_effect[cell]` lets individual regions deviate
from it. The change is small and surgical — two additions to the existing model
code, reusing the CAR machinery already in place — and confined to the iNat
half of the model (the camera-trap occupancy submodel's trend term is
unchanged in this first pass; extending it there is a natural but separate
next step).

## Why simulate before refitting

Before fitting this to the real, expensive bobcat data (a multi-day process:
the existing national-trend model needed 16 rounds of chunked MCMC restarts to
converge), we validated the extension against data simulated from a known
truth — the standard way to check whether an estimator can actually recover
what it claims to, before trusting it on real data.

A full-scale simulation wasn't tractable: profiling the real model shows
roughly 7.5 seconds per MCMC iteration, driven by the sheer size of the
iNaturalist prediction surface (~300,000 grid subcells × 18 years). Fifty
replicate fits at that scale would need over 1,000 CPU-hours. Instead, we kept
the real spatial structure exactly (the actual ~900-cell adjacency graph
`year_effect` is estimated on — this is cheap regardless of data volume) and
subsampled the *data volume* — fewer camera sites and iNaturalist grid cells,
still carrying their real covariate and effort values, only the response data
simulated fresh from the model's own likelihood for a known, controlled trend
pattern (positive in the Northeast, negative in the West, near-zero elsewhere).

## Results so far

**Stage 1 (50 replicates, ~20% of the real spatial cells directly informed):**
coverage was already good — the model's 95% credible intervals contained the
true value 99.5% of the time, meaning its uncertainty is honest, not
overconfident. But point estimates showed real bias toward zero in both
informative regions, and only 80% of replicates recovered the correct
Northeast-vs-West sign.

**Stage 2 (20 replicates, ~48% of cells informed, a confirmatory follow-up):**
more than doubling the informed-cell count cut overall bias by roughly
two-thirds and pushed sign-recovery to a clean 100%, with coverage unchanged
(already good). The improvement wasn't perfectly uniform across regions.

| | Stage 1 (20% informed) | Stage 2 (48% informed) |
|---|---|---|
| bias, all cells | −0.040 | −0.015 |
| RMSE, all cells | 0.289 | 0.167 |
| coverage, all cells | 99.5% | 99.4% |
| sign recovered | 40/50 (80%) | 20/20 (100%) |

**Interpretation:** this is a consistent, clear trend in the direction you'd
expect if the original bias were an artifact of deliberately thinning the data
for computational tractability, rather than a flaw in the model formulation
itself — the real deployment will have on the order of 2x the data density we
tested even in the confirmatory run. That's encouraging but not proof: we
cannot test at literal full scale without paying the same computational cost
we built this whole reduced-simulation approach to avoid.

## Proposed follow-up: what grain of heterogeneity can we actually detect?

The validation above tested recovery of one specific spatial pattern — a broad,
roughly continental-scale East/West contrast. It says nothing about whether
finer heterogeneity (state-sized clusters, or smaller) would be detectable
given how the real data are actually distributed in space. That's a distinct
and important question: the model's CAR field has the *capacity* to vary
cell-by-cell, but the CAR prior's smoothing will blur out any true pattern
finer than what the data density can actually resolve — and we haven't yet
mapped out where that resolving power breaks down.

We found that an existing simulation framework already built for this project
addresses a closely related pair of questions — (1) the relative contribution
of camera-trap vs. iNaturalist data to the fitted model, and (2) how much
additional data of each type would be needed for more statistical power — using
an idealized synthetic landscape (not the real geography) and a Gaussian
random field with a controllable spatial-autocorrelation "range" parameter for
generating spatially patterned covariates. Rather than build a new framework
from scratch, we're extending that one: reusing its established data-source
axis and its GRF/range mechanism for controlling spatial grain, applied to the
real bobcat geography (real cell locations, real adjacency graph) instead of a
synthetic grid, and to the new trend field specifically, since the existing
framework's validated scenarios never included a spatially-varying trend at
all.

The proposed design has three phases, run sequentially rather than as one huge
factorial crossing (which would be computationally infeasible — see below):

**Phase 1 — spatial grain.** Simulate true `year_effect` fields at a range of
spatial autocorrelation scales, from the current continental-scale contrast
down through state-cluster-sized patches to fine, local patches, holding data
density fixed at a realistic level. This traces out a "minimum detectable
patch size" curve given current data density.

**Phase 2 — camera vs. iNat contribution.** At the grain identified as
just-detectable in phase 1, re-run the existing framework's five-level
camera/iNat data-source-mix axis (already used in the group's prior validation
work) to see how much each source individually contributes to detecting the
*spatial* trend specifically, not just the existing national scalar.

**Phase 3 — how much more data would help.** Vary camera and iNaturalist data
volume separately (holding the other fixed), at five levels each drawn from
the *real, empirically observed* density distribution across the current
dataset — not arbitrary multipliers. Camera-trap coverage turns out to be far
patchier than iNaturalist: roughly two-thirds of the real spatial grid cells
have zero camera coverage at all, versus only about 5% for iNaturalist, and
both sources are heavily right-skewed (a few well-studied areas, most cells
modest). The five levels for each source are the 10th/30th/50th/70th/90th
percentile of real per-cell density among cells that do have coverage,
translated to a per-cell-per-year rate. This directly answers "how much more
data of each type would be needed for additional power," anchored to
realistically achievable densities (e.g., "extend camera coverage in a sparse
region to what our best-covered regions already have") rather than an abstract
sensitivity curve.

This staged approach keeps the total cost manageable: roughly 4 + 5 + (5+5) =
19 scenarios rather than a full cross of every axis (which would run into the
hundreds), each scenario run at a modest replicate count (10-20) sufficient to
see a clear directional trend, with the option to follow up at higher
precision on whichever scenario turns out to be most consequential — the same
"confirm the direction cheaply, then invest in the interesting case" approach
used for the validation above.

## Where things stand

The spatially-varying-trend extension has passed simulation validation with a
consistent, improving trend, but not yet at full real-data scale. The grain/
data-source/data-volume follow-up study above is designed but not yet built.
The real, multi-day bobcat refit (fitting the extended model to the actual
data) has **not** been started — that's a real compute-time commitment worth
discussing with the group before committing to it, and is a natural thing to
decide once (or if) the follow-up study above is judged worthwhile to run
first.
