# Population trends from integrated camera-trap and iNaturalist data: the five-year window

## Take-home points

1. **18 of the 34 species in the fleet have been fitted so far, and 17 of those converged.** The remaining 16 are still running. Every count below is out of the 17 converged fits, not out of 34, so none of them is a statement about the fleet.

2. **8 species resolve a direction**, meaning the 95% credible interval for the trend excludes zero. They are brown bear, California ground squirrel, grey wolf, black-tailed jackrabbit, grey fox, Northern raccoon, Eastern gray squirrel and coyote. The remaining 9 fits have intervals spanning zero, so the data leave their directions open.

3. **2 of the resolved species sit below the range size where these estimates hold up.** Brown bear and California ground squirrel occupy 30 and 74 grid cells, against a fleet median of 425. At -59% per year, brown bear would fall to 1% of its starting abundance over the five years. For brown bear, over five years, 8 of 14 records give a direction (6 stable and 2 increasing); over ten years, 14 of 14 records give a direction (12 increasing and 2 stable). None of the 4 agency records for California ground squirrel states a direction, so no external assessment is available to compare against.

4. **Setting those aside, 6 species with ranges above 100 cells resolve a direction**: grey wolf at -25.4% per year [-45.7, -0.6]; black-tailed jackrabbit at -19.2% per year [-25.9, -11.9]; grey fox at -5.5% per year [-10.8, -0.0]; Northern raccoon at +4.8% per year [+2.0, +7.7]; Eastern gray squirrel at +5.7% per year [+2.7, +8.8]; coyote at +9.1% per year [+5.6, +12.4].

5. **Stream agreement stays unconfirmed for 9 of the 17 converged fits, on the joint fit alone.** 1 needs a camera-only control fit because the exponent linking iNaturalist records to abundance has an interval reaching 1 or including 0, and 8 have not yet settled: their 12 chains still disagree on the exponent, so its interval cannot be trusted whether it is narrow or wide.

## What is being measured

**Trend.** The change in a species' relative abundance per year, expressed as a percentage. A value of −5% per year means the modelled abundance in an average grid cell falls by five percent of its value each year, so after five years about 77% of the starting abundance remains. The model estimates how abundance changes from year to year, on a relative scale. Converting that into a count of animals would need a separate calibration we have not done.

**Credible interval.** The range containing the true value with 95% probability, given the model and the data. When it excludes zero we say the direction is resolved. When it spans zero, the data leave the direction open: a wide interval spanning zero is a statement about what we can see, and the population may still be rising or falling.

**The five-year window.** The fits here use detections from the most recent five years only. Camera-trap effort grew steeply over the study period, and few cameras were running in the earliest years, so those years add more noise than signal.

**Stream agreement (the exponent θ₁).** The model ties the two data streams together with one number: the expected count of iNaturalist records in a grid cell is proportional to modelled abundance raised to the power θ₁. At θ₁ = 1 records track abundance proportionally; at θ₁ = 0.5 doubling the animals brings about 41% more records; at θ₁ = 0 the records say nothing about where the animals are. When the fitted interval reaches 1 or includes 0, the iNaturalist half of the fit may be tracking where people looked instead of where the animals are, and the species needs a camera-only control fit before we report its joint trend.

**Data-limited.** A species occupying fewer than 100 of the 50 km grid cells. Across the species with both a joint and a camera-only fit, those below this size sat a median of 17 percentage points per year away from their camera-only estimate, against 1.7 points above it. The threshold comes from that measurement.

**Convergence.** Each fit runs several independent chains. R-hat compares the spread between chains with the spread within them, and should be at or below 1.01; effective sample size is how many independent draws the chains are worth, and should be at least 400. We report only the fits that pass both.

## Methods

**The model.** Camera-trap detections and iNaturalist records are fitted to one shared surface of relative abundance on a 50 km grid, each stream through its own observation model. We know how long each camera ran, so a camera that recorded nothing tells us as much as one that recorded an animal, and we anchor the trend to the camera stream. We do not know how much time iNaturalist observers spent looking. We control for that using the total number of mammal records in the same cell and year as a measure of how much observing happened there.

**How the trend is computed.** The model fits a slope against a standardised calendar year, so the raw coefficient gives change per standard deviation of year. Dividing by that standard deviation, exponentiating and subtracting one converts it to proportional change per year. The same scaling factor is read from each fit's own record, because the standard deviation of year differs between species with different numbers of sampled years.

**Which fits are reported.** A fit enters the results only if both R-hat and effective sample size pass for the trend parameter. Fits that fail appear on the figure as open markers on a dotted interval, and every count in the title excludes them.

## Results

![Five-year trend for every fitted species]({{artifact:art_FIGID}})

*The 17 converged fits on the five-year window. Each bar is the 95% credible interval and each dot the posterior mean. Species in amber occupy fewer than 100 grid cells, with their cell count beside the name. Open markers on dotted intervals did not converge and are excluded from the counts in the title.*

| Species                    |   Range (cells) |   Trend (%/yr) |   Lower |   Upper | Direction          |   θ₁ |   θ₁ lower |   θ₁ upper | Stream agreement                 |
|:---------------------------|----------------:|---------------:|--------:|--------:|:-------------------|-----:|-----------:|-----------:|:---------------------------------|
| Brown bear                 |              30 |         -59.13 |  -79.14 |  -37.41 | decline            | 0.32 |       0.08 |       0.84 | theta1 below reporting threshold |
| California ground squirrel |              74 |         -37.73 |  -46.57 |  -27.57 | decline            | 1.54 |       1.28 |       1.83 | needs camera-only control        |
| Grey wolf                  |             179 |         -25.35 |  -45.68 |   -0.64 | decline            | 0.44 |       0.23 |       0.7  | theta1 below reporting threshold |
| Black-tailed jackrabbit    |             424 |         -19.2  |  -25.94 |  -11.93 | decline            | 0.72 |       0.62 |       0.84 | clear                            |
| Moose                      |             258 |         -10.47 |  -20.64 |    0.83 | no change detected | 0.52 |       0.43 |       0.62 | clear                            |
| Groundhog                  |             406 |          -8.19 |  -17.72 |    1.56 | no change detected | 0.21 |      -0.4  |       0.41 | theta1 below reporting threshold |
| Fisher                     |             176 |          -7.7  |  -19.42 |    3.2  | no change detected | 0.96 |       0.57 |       1.31 | theta1 below reporting threshold |
| Collared peccary           |             168 |          -7.48 |  -20.84 |    5.42 | not reportable     | 0.42 |      -0.68 |       0.74 | theta1 below reporting threshold |
| Kit fox                    |             167 |          -6.32 |  -24.4  |   12.95 | no change detected | 0.68 |       0.41 |       1    | theta1 below reporting threshold |
| Grey fox                   |             813 |          -5.5  |  -10.85 |   -0.02 | decline            | 0.5  |       0.45 |       0.56 | clear                            |
| Nine-banded armadillo      |             336 |          -2.88 |   -7.65 |    1.7  | no change detected | 0.42 |       0.35 |       0.49 | clear                            |
| Mule deer                  |             469 |          -2.37 |   -7.25 |    2.3  | no change detected | 0.48 |       0.43 |       0.54 | clear                            |
| White-tailed deer          |             787 |          -0.73 |   -3.49 |    1.9  | no change detected | 0.5  |       0.47 |       0.54 | clear                            |
| Eastern chipmunk           |             434 |          -0.59 |   -5.14 |    3.94 | no change detected | 0.23 |      -0.2  |       0.41 | theta1 below reporting threshold |
| Elk                        |             359 |           2.71 |   -4.71 |   12.18 | no change detected | 0.36 |       0.29 |       0.45 | theta1 below reporting threshold |
| Northern raccoon           |            1002 |           4.84 |    2.01 |    7.7  | increase           | 0.31 |       0.28 |       0.34 | clear                            |
| Eastern gray squirrel      |             489 |           5.66 |    2.71 |    8.77 | increase           | 0.51 |       0.47 |       0.55 | clear                            |
| Coyote                     |            1032 |           9.1  |    5.62 |   12.44 | increase           | 0.53 |       0.46 |       0.59 | theta1 below reporting threshold |

### What the steepest declines rest on

Brown bear and California ground squirrel return the two steepest declines in the table, and are also the smallest ranges in the fleet. The model sees a species in a few dozen grid cells across five years of detections, and fits the steepest slope those detections allow. Both species also failed the anchoring check, which compares the joint estimate against a camera-only fit of the same data.

## What happens next

**16 species are still being fitted.** The species finished so far include several of the fleet's smallest ranges, so a fleet-level statement waits for the rest.

**9 of the converged fits need a camera-only control fit** before their joint trend can be reported, because stream agreement could not be established from the joint fit.

**1 fit did not converge** (Collared peccary) and is excluded from all results.

## Limitations

**These species finished first because they are small and quick to fit**, which is the same property that makes their estimates least reliable. The fleet average will look different.

**The trends are proportional change.** Reading them as population sizes or animal counts would need a calibration step we have not done.

**Stream agreement is unconfirmed for most of these fits.** Where θ₁ came back unsettled, the iNaturalist half of the fit may be tracking where people looked instead of where the animals are. A camera-only control fit settles it.

**The agency comparison is still to come.** It needs the state and provincial records matched to the same five-year window.

## Glossary

**Chain.** One independent run of the fitting algorithm. Several are run per species so that agreement between them can be checked.

**Credible interval.** See *What is being measured*.

**Effective sample size.** The number of independent draws the chains are worth once their internal correlation is accounted for.

**Grid cell.** A 50 km square. Range size is reported as the number of cells a species occupies.

**R-hat.** A convergence statistic comparing variation between chains with variation within them. One means perfect agreement.

**Relative abundance.** How many animals are in a place compared with another place or another year, without estimating the absolute number.

**θ₁ (theta-one).** The exponent linking iNaturalist record counts to modelled abundance. See *What is being measured*.
