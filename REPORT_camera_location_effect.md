# Should the temporal model carry a camera-location grouping term?

A question for Arielle and Krishna, with the supporting counts.

## Take-home points

1. **The non-temporal model groups camera deployments that share a physical
   location; the temporal model does not.** Goldstein et al. include a random
   intercept shared by all deployments occurring at the same unique location in
   any year. Our model loops over deployments with an independent occupancy
   draw for each one and no term linking repeat visits to the same place.
2. **The location identifier exists in our prep and is dropped at the model
   boundary.** `prep_data_for_spoccupancy.R` builds `loc_ID` from
   `distinct(longitude, latitude)` and passes it as spOccupancy's `grid.index`.
   It is not among the constants our NIMBLE model receives.
3. **97.8% of the variance in deployment year is between locations, 2.2% within
   them.** The trend coefficient multiplies a per-deployment year covariate, so
   almost everything identifying the trend is a comparison between different
   places sampled in different years.
4. **Only 1,862 of 21,206 camera locations were sampled in two or more years.**
   Those locations hold 6,234 of 26,776 deployments, 23.3%. They are the part of
   the data that can speak to change at a fixed place.
5. **Adding the term would widen the credible intervals and could unresolve some
   trends.** Four deployments at one trail currently enter as four independent
   statements about occupancy there. Grouping them makes the intervals reflect
   the number of places observed rather than the number of deployments.
6. **The repeat-location data are concentrated in recent years and in under half
   the projects.** 5,752 of the 6,234 repeat-location deployments are from 2019
   onward against 482 before, and they come from 27 of 60 projects.

## What the two models do

**Goldstein et al. (2025), non-temporal.** The camera detection submodel
contains a grouping random effect shared by detections associated with all
deployments occurring at the same unique location in any year. It sits in the
detection linear predictor, which is separate from the intensity predictor, so
it absorbs site idiosyncrasy without competing with abundance.

**Our model, with the trend added.** The camera submodel loops over deployments.
For each one it builds a detection probability from day of year, its square,
canopy height and distance to road; builds an occupancy probability from the
cell-level spatial intercept, the two climate effects, the occupancy covariates,
and the trend term `year_beta * year_occ[i]`; and evaluates the detection
history with `dOcc_v`. The index runs over deployments: 20,531 for bobcat,
21,559 for white-tailed deer, 1,654 for moose. A location deployed in 2019 and
again in 2022 contributes two rows with two independent occupancy values,
differing in their year covariate, and nothing ties them together.

## Definitions

**Location.** A distinct pair of coordinates. Matched here at five decimal
places, about 1 m, which treats a camera re-hung at the same tree as the same
location and two cameras 50 m apart as different ones.

**Deployment.** One camera placement with a start and end date, identified by
project and deployment identifier. There are 26,776 of them.

**Location grouping random effect.** A normally distributed intercept, one per
location, shared by every deployment at that location. It states that
deployments at one place are more like each other than like deployments
elsewhere, for reasons the covariates do not capture.

**Within-location and between-location year contrast.** Split the year of each
deployment into its location mean year and the deviation from that mean. The
within-location part is what a revisit provides. The between-location part is
the difference between places sampled at different times.

**Pseudo-replication.** Treating correlated observations as independent, which
narrows credible intervals below what the data support.

## Evidence

![Repeat-location structure in the camera data]({{artifact:art_PLACEHOLDER}})

### How much repeat sampling exists

| Years sampled | Locations | % of locations | Deployments | % of deployments |
|---|---|---|---|---|
| 1 | 19,344 | 91.22 | 20,542 | 76.72 |
| 2 | 1,033 | 4.87 | 2,500 | 9.34 |
| 3 | 371 | 1.75 | 1,359 | 5.08 |
| 4 | 224 | 1.06 | 1,025 | 3.83 |
| 5 | 134 | 0.63 | 694 | 2.59 |
| 6 | 76 | 0.36 | 481 | 1.80 |
| 7 | 24 | 0.11 | 175 | 0.65 |

Across 26,776 deployments there are 21,206 distinct locations, 1.26 deployments
per location. The 1,862 locations sampled in two or more years span a median of
2 years and at most 11.

### How the trend is identified

Decomposing the variance in deployment year gives 2.2% within locations and
97.8% between them. Restricting to the 6,234 deployments at repeat locations
raises the within-location share to 22.0%.

The practical reading: if the places sampled in 2024 differ from the places
sampled in 2012 in any way the covariates do not capture, that difference enters
the trend coefficient as though it were change over time. The term that would
prevent this is the location intercept, which absorbs between-place differences
and leaves the trend to be carried by revisits.

### Where the repeat data sit

The repeat-location deployments come from 27 of the 60 projects. 5,752 of them
are from 2019 onward and 482 from 2018 or earlier, so the within-location signal
is concentrated in the recent part of the record.

## The question for Arielle and Krishna

Two readings are available and they lead to different actions.

**If the omission was a considered choice**, the argument would be that with 91%
of locations sampled once, a location intercept is estimable for only 8.8% of
them and the rest contribute a prior draw each, adding 21,206 parameters for
little gain. There may also have been a convergence or runtime reason. That
belongs in the methods, because a reviewer who knows the Goldstein model will
ask why the term is absent.

**If it was incidental**, lost when the temporal structure was added, then the
trend is currently identified from a contrast we would not choose, and the fix
is one index vector the prep already computes.

A middle option exists: group at the array or project level rather than the
location level. That estimates far fewer intercepts, absorbs the
project-specific protocol differences that are plausibly the largest source of
between-place variation, and is identifiable for every deployment.

## Proposed test

Bobcat and moose at the 18-year window, each fitted with and without the
location grouping term.

Bobcat has 20,531 deployments and dense repeat coverage, so it shows the largest
change if pseudo-replication is inflating precision. Moose has 1,654 deployments
in a sparse northern footprint, so it shows whether the term is estimable when
repeat locations are few. Both are fitted and converged, so each comparison runs
against a known baseline.

Three quantities decide the answer:

- **The change in the trend estimate.** A shift means between-location
  differences were loading onto the trend.
- **The change in its interval width.** A widening is the pseudo-replication
  being corrected.
- **The standard deviation of the location effect.** If it is indistinguishable
  from zero, locations are not idiosyncratic once covariates are accounted for
  and the current model stands. If it is large, the current intervals are too
  narrow.

Two chains per fit at current settings, four fits, roughly a week of queue time
alongside the windowed refits.

## Limitations

The variance decomposition uses deployment year, which is what the trend
covariate encodes. It describes how the trend is identified and does not measure
how much the estimate would move, which requires the refit.

Locations are matched on exact coordinates at about 1 m. A camera moved 50 m
between years counts as two locations, so the repeat-sampling figures are a
lower bound. Matching within a radius would raise both the count of repeat
locations and the within-location share.

The detection submodel was not examined for the same issue. Whether detection
also needs a location term is a separate question.

Nothing here tests the opposite risk, that a location intercept absorbs real
local abundance variation the spatial fields should carry. With the intensity at
50 km and locations at metres these operate at different scales, but only the
refit would show whether the spatial fields change.
