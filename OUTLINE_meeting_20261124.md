# Slide outline, 24 November meeting

Krishna, Arielle, Ben. Each slide below gives the one point it makes, what goes
on it, and whether the material exists yet.

## Take-home points

1. **Nineteen species should be fitted by Tuesday morning, and only a minority
   resolve a direction.** On the eleven finished now, three resolve, all
   declines. The eight species the abstract names land overnight on Monday, and
   they are the ones expected to carry the increases.
2. **The honest agency comparison is much thinner than the record count
   suggests.** Of 829 agency records, 332 state a direction, and only 111 of
   those are unflagged survey estimates. Those 111 cover 12 species. Twenty-two
   of our 34 have no survey-grade agency trend at all.
3. **Window choice moves the headline.** White-tailed deer gives +4.00%/yr over
   the full record and +0.46%/yr over ten years, on the same data and model.
4. **The trend we report is the best-converged parameter in the model**,
   clean in 10 of 11 species, while the spatial fields are clean in 1 of 11.
   That ordering decides which results carry weight and which are pattern.
5. **The stream-agreement metric is informative where it converges and silent
   where it does not.** Three of eleven five-year fits give a usable theta1.
   Three others fail for a structural reason that more computing cannot fix.

---

## 1. Title and the question

One sentence: population trends are unknown for most North American mammals,
and no existing framework combines camera traps with iNaturalist to estimate
them.

## 2. Why this is hard, and the shape of the answer

The two streams have complementary weaknesses. Camera traps measure effort but
cover little ground; iNaturalist covers the continent with unknown effort.
Fitting both to one shared surface lets each cover for the other.

Figure: `fig_integration.png`, the methods schematic.
Status: ready.

## 3. The model, in one slide

The shared surface, the two observation models, and the one sentence that
matters for reading every later slide: the trend we report comes from the
camera submodel, and the iNaturalist side has its own increment that absorbs
reporting change.

Status: ready, needs a simple diagram drawn from the methods text.

## 4. What changed from Arielle's original, and why

Short slide, since they have the document and the code. Three columns: the
change, the reason, the evidence. The point to make out loud is that the model
equations are untouched. Our fork of the production `nimbleCode` block is
byte-identical to the original; every change sits upstream in data preparation
or downstream in which parameter we report.

Figure: a six-row table lifted from `MODEL_SPEC_AND_CHANGES.docx`.
Status: ready. Keep to one slide and let them pull on it if they want to.

---

## 5. Where the fleet stands

Count of species fitted, converged, resolving a direction, and still running.
Say plainly that the species finished first are the small-range fast ones, so
the current set is not representative of the fleet.

Table: `table_fleet_win5.csv`, extended as the overnight species land.
Status: eleven species now; nineteen expected Tuesday morning.

## 6. The trends

The forest plot, five-year window, every fitted species with its interval,
data-limited species marked, non-converged fits drawn separately and excluded
from every count.

Figure: `fig_fleet_win5.png`, regenerated Tuesday morning.
Status: regenerates with one command when the overnight species land.

## 7. Where the animals are

Relative abundance at the start and end of the window plus the trend map, three
panels per species. The two year panels share a colour scale, so the difference
between them is the trend and can be read directly.

Figures: `fig_abundance_three_p1.png` through `p4.png`.
Status: ready for eleven species. Pick three or four species for the slide.

## 8. Where the cameras say it changed

Change in detection probability per camera array, zoomed and scaled per
species, species whose trend spans zero drawn in grey.

Figures: `fig_array_change_p1.png` through `p3.png`.
Status: ready. One slide, two or three species.

Worth saying out loud: the model applies one national trend, so every array for
a species moves the same direction. What varies between arrays is how far, and
that comes from where each array starts on the detection curve.

---

## 9. Does the model work? Convergence by parameter family

The slide a reviewer will ask for, and the answer is better than expected in
one specific way: **the trend we report is the best-behaved parameter in the
model.** It converges in 10 of 11 species. Nothing else comes close.

| Parameter family | Species passing both thresholds |
|---|---|
| **Trend (reported)** | **10 of 11** |
| Overdispersion | 8 of 11 |
| Detection | 7 of 11 |
| iNaturalist trend terms | 7 of 11 |
| Habitat covariates | 3 of 7 monitored |
| Stream coupling | 3 of 11 |
| **Spatial fields** | **1 of 11** |

That ordering is not an accident. `year_beta` sits directly in the camera
likelihood where the detection histories identify it; the spatial fields are
hundreds of correlated cell-level parameters informed by sparse data.

Figure: `fig_model_performance.png`, a species-by-family R-hat grid.
Status: ready for eleven species, regenerates with the rest.

The consequence to state plainly: the trend table is reportable, the spatial
fields are pattern rather than established local effects, and we say so on the
figures themselves rather than in a footnote.

## 10. Habitat associations

Nine covariates across every species, with intervals, and the effects from
non-converged chains drawn separately and excluded from the counts.

Figure: `fig_covariate_effects.png`.
Status: ready for eleven species.

What it shows on the current set: **61 of 99 effects come from converged
chains, and 27 of those 61 have an interval excluding zero.** All counts below
are on converged chains only; the other 38 effects are drawn in the figure but
excluded from every number, since an interval from chains that did not mix is
not an interval.

The covariates that resolve most often are mixed forest (5 of 7 converged
fits), then agriculture, human population and soil clay (4 each).

Three are one-directional wherever they resolve, which is worth saying out loud
because the consistency is the finding: **agriculture, human population and
soil clay are each negative in all four species that resolve them.** The three that split are mixed
forest at 3 negative to 2 positive, and evergreen forest and soil silt at 1
negative to 2 positive each. The remaining three covariates resolve in fewer
than three species, so their direction carries no weight.

## 11. The spatial fields

Baseline occurrence and the climate responses, per species, clipped to the
cells each species occupies.

Figures: `fig_spatial_baseline_win5.png`, `fig_climate_mwmt_win5.png`.
Status: ready for eleven species.

Pair this directly with slide 9. These fields converged in one of eleven fits,
so the honest framing is that they carry spatial pattern worth investigating
and not established local effects. The proportion of cells whose interval
excludes zero ranges from 4% for collared peccary to 100% for kit fox, which is
itself informative about where the data are thick.

---

## 12. Do the two data streams agree?

theta1, the exponent linking iNaturalist records to the camera-side surface,
with both thresholds marked. Zero means the records carry no spatial
information about abundance; one means they track it proportionally, which no
observer behaviour produces.

Figure: `fig_stream_agreement_win5.png`.
Status: ready for eleven species.

The three bands on the current eleven: five in the healthy range, three whose
interval includes zero, three whose interval reaches one.

## 13. When the coupling cannot be estimated, and why

Eight of eleven fits miss the mixing threshold for theta1, and they are two
different problems. Five are slow and would clear with two to five times the
iterations. Three are structural: the likelihood constrains only the product of
theta1 and the iNaturalist trend, so opposite-sign chains fit equally well and
more computing changes nothing.

The three structural cases are exactly the three whose theta1 interval includes
zero, and they have the lowest camera spatial contrast in the set. A species
spread evenly across its range gives the coupling nothing to latch onto.

Figure: `table_theta1_mixing_diagnosis.csv` as a table, plus the contrast
correlation.
Status: ready.

This is a methods contribution in its own right, so it is worth a slide rather
than a footnote: it says when integration can work.

## 14. Window sensitivity

Three windows for the two species fitted on more than one. White-tailed deer
goes from a resolved increase to no detectable change.

Figures: `fig_window_trends.png` and `fig_win10_windows.png`.
Status: ready. Deer's five-year row fills in Tuesday morning; bobcat's
five-year fit is among the held chains and will not exist by Tuesday.

Decision to put to the group: how many species get a ten-year fit.

---

## 15. What state agencies publish

The map of how many of our 34 species each state publishes a trend for.
Minnesota leads at 13; the median state publishes 4; no state reaches half.

Figure: `fig_agency_infographic.png`.
Status: ready.

## 16. What counts as a data-supported agency trend

The critical slide, and the one to spend time on. Starting from 829 records:

| | n | share of the 332 |
|---|---|---|
| Records that state a direction at all | 332 of 829 | |
| **A. Population or survey estimate, unflagged** | **111** | **33%** |
| B. Survey estimate, flagged as weakened | 28 | 8% |
| C. Harvest index, not abundance | 80 | 24% |
| C. Other index, not a direct estimate | 81 | 24% |
| D. Audit found an invalidating problem | 24 | 7% |
| D. Qualitative or expert opinion | 8 | 2% |

The sentence to land: **the 111 survey-grade records cover 12 species.** They
are almost entirely big game and large carnivores: black bear, brown bear,
bobcat, elk, fisher, grey wolf, moose, mule deer, peccary, pronghorn, puma,
white-tailed deer. **Twenty-two of our 34 species have no survey-grade agency
trend anywhere in North America.**

Also worth naming: 338 flags record a harvest or effort index being used as
abundance, 142 an observation index, and 99 cases where the agency itself
discloses the metric is not abundance.

Table: `table_agency_evidence_tiers.csv`.
Status: ready.

## 17. Model against agency, where a comparison is possible

Only where an agency publishes a survey-grade direction for a species we have
fitted. Give the count first, then the agreements, then the disagreements with
their reasons.

Status: **not ready.** The existing concordance was scored on the full-record
fits and a different window. It needs rescoring once the overnight species
land, and it should be restricted to tier A and B records.

Say what it will be rather than showing the old numbers.

## 18. A student is verifying the agency records

One slide on the dashboard: 323 records across 207 source documents, priority
ordered, six questions each, exports a revised spreadsheet.

Figure: a screenshot of `agency_verification_dashboard.html`.
Status: ready, not yet started by the student.

---

## 19. What is still coming

| Result | When |
|---|---|
| The eight abstract species, five-year | Monday night into Tuesday morning |
| The remaining 15 species | Held; days once released |
| Camera-only control fits | Not run on the five-year window |
| Urban-rural screen | Not run on the five-year window |
| Ten-year fits beyond two species | Not started |
| Student verification of agency records | Not started |

Be explicit that the camera-only control and the urban-rural screen are both
required by our own workflow and neither exists for the reported window. The
California ground squirrel is the live case: it resolves a decline but its
theta1 interval sits above one, which is exactly the state that requires a
control fit before the trend is reported.

## 20. Decisions for the group

1. How many species get a ten-year fit, and does the ten-year become a second
   reported window or stay a sensitivity check.
2. Whether to spend queue time on camera-only controls for the species that
   need them, before or after the remaining 15 species.
3. Whether the structural non-identification result is a paper section or a
   supplement.
4. Whether to report the three species whose theta1 cannot be estimated, with
   the camera trend alone and no congruence statement.
5. Ben's proposal on freezing numbered releases for review.
