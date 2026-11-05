# Organised reporting campaigns and what they do to our iNaturalist effort measure

## Take-home points

1. **Five campaign-and-species combinations shifted reporting enough to matter.**
   Two Twin Cities coyote projects, Backyard Bobcats in Southern California,
   the Grit City Carnivore Project in Tacoma, and the California Ground Squirrel
   Census in the East Bay. Each raised its target species' records faster than
   every other mammal's in the same places over the same years, by 18% to 88%,
   with 95% intervals excluding no effect.
2. **Four of our species are affected, in 58 of about 3,300 grid cells.** Coyote,
   bobcat, northern raccoon and California ground squirrel. The windows run from
   2014, 2018, 2019 and 2020 respectively, all still open.
3. **Coyote is the species that matters most.** Two overlapping Twin Cities
   projects give the two largest effects in the set, 1.88 and 1.66, and coyote is
   one of the eight species whose increase our camera-only refits corroborate.
4. **Campaigns did not bring new people to mammal reporting.** In campaign
   regions the number of people reporting any mammal grew slightly more slowly
   than in the rest of the country. What changed was which species existing
   observers photographed.
5. **The correction belongs in the opportunistic reporting submodel, applied to
   the 58 cells in their campaign years, with the effect size supplied as a prior
   rather than estimated freely.** Three of the five footprints cover only three
   to five cells, too few to estimate a coefficient from the fitted data alone.
6. **Of twelve campaign-and-species combinations tested, five show an increase,
   one shows a significant decrease, four cannot be distinguished from no effect,
   and two are range expansions we deliberately retain.** The decrease is red fox
   under a campaign that also names coyote, where reporting moved toward coyote
   and away from red fox.

## What a campaign is, and why it threatens the effort measure

**Organised reporting campaign.** A project that asks the public to report a
named species or group in a named region, run by a park district, museum,
university or agency. We counted only campaigns that collect records through
iNaturalist, because an agency sighting form creates no iNaturalist record and
therefore cannot touch our data.

**Effort measure.** For each 50 km cell and year, the model scales the expected
count of one species by the total number of research-grade mammal records in
that cell and year. The assumption is that a cell and year with many mammal
records of any kind had many people looking, so one species' count there should
be judged against that background.

**Why a campaign breaks it.** A campaign raises records of its target without
raising records of other mammals proportionally. The denominator barely moves
while the numerator rises, so the model reads the excess as abundance. This is
the one bias mechanism we have found that the effort measure cannot absorb.

**Selectivity.** How much faster the target species' records grew than every
other mammal's in the same cells over the same years, divided by that same ratio
in all cells outside the campaign region. A value of 1 means the target grew
exactly as other mammals did and the campaign changed nothing about the mix. A
value of 1.5 means the target gained half again as much ground as its
neighbours. Below 1 means the target grew more slowly than other mammals, which
some campaigns did.

**Exposure.** The fraction of a 50 km cell's area that falls inside a campaign's
region. Marin County is 1,407 square kilometres against a cell's 2,500, so a
Marin campaign spreads across three cells at about a fifth of each. Exposure
records that partial coverage instead of treating a cell as wholly in or out.

**Range-expansion campaign.** A campaign whose target is new to the region, such
as the eastern fox squirrel in Utah. People report a species they have not seen
before whether or not a project asks them to, so the campaign and the novelty
cannot be separated.

## How we found and measured them

### Finding the campaigns

A systematic search of iNaturalist projects, state and provincial agency pages,
museum and park-district programs, university extension pages and conservation
nonprofits across the United States and Canada returned 61 campaigns targeting
North American mammals, 56 confirmed from the organisation's own page and 5
described only by a secondary source. Of these, 36 collect records through
iNaturalist.

Each campaign was recorded with its target species, whether it targets one
species or a group, its region as specifically as the source states, its start
year, and a URL.

### Bounding each campaign in space and time

Campaign regions were translated into county boundaries, then intersected with
the 50 km model grid to give each cell an exposure fraction. Cells at least a
quarter covered were treated as inside the campaign. Time windows run from the
campaign's stated start year to its stated end, or to 2025 where it continues.

### Measuring the effect

For each campaign paired with each species it targets, we compared the four
years before the campaign's start year with the four years after, counting
records of the target species and records of all other mammals, inside the
campaign's cells and in every cell outside them. The start year itself was
excluded. The resulting ratio of ratios is the selectivity, and its 95% interval
comes from the eight counts that form it.

A pair entered the screen only if the campaign's cells held at least 25 records
of the target species in each window.

### Classifying the result

A pair counts as a strong effect when its 95% interval excludes 1 and its target
is not expanding into the region. Range-expansion pairs were identified from the
campaign's own description of its purpose.

## What we found

### The campaigns that shifted reporting

![Campaign selectivity with 95% intervals]({{artifact:art_SEL}})

| Campaign | Species | Years | Cells | Selectivity | 95% interval |
|---|---|---|---|---|---|
| Coyote Sightings in the Twin Cities | Coyote | 2018–2025 | 4 | 1.88 | 1.40 – 2.53 |
| Twin Cities Coyote and Fox Project | Coyote | 2019–2025 | 5 | 1.66 | 1.33 – 2.09 |
| Backyard Bobcats | Bobcat | 2014–2025 | 47 | 1.47 | 1.05 – 2.06 |
| Grit City Carnivore Project | Northern raccoon | 2018–2025 | 3 | 1.35 | 1.04 – 1.74 |
| California Ground Squirrel Census | California ground squirrel | 2020–2025 | 3 | 1.18 | 1.07 – 1.29 |

Read the largest as follows: in the four Twin Cities cells, coyote records grew
88% faster between the two windows than other mammals did in the same cells,
after removing the growth that coyote and other mammals showed everywhere else
over the same years.

The two Twin Cities projects cover nearly the same counties and name the same
species. They are one effect measured twice, not two independent effects.

### Where the effects fall

![The 58 affected grid cells]({{artifact:art_MAP}})

Southern California carries 47 of the 58 cells, through a bobcat campaign
covering six counties. The other three footprints are three to five cells each:
the East Bay, the Minneapolis-St. Paul metropolitan area, and Pierce County,
Washington.

### Campaigns changed what people reported, not how many people reported

Across twelve measurable campaigns, the number of people reporting any mammal in
the campaign's states grew more slowly than in the rest of the country in eight
of twelve, with a median ratio of 0.89. The number reporting the target species
specifically was indistinguishable from the rest of the country, median 1.05.

So the campaigns did not expand the pool of mammal observers in a way our effort
measure would benefit from. They redirected attention within the existing pool,
which is the pattern that biases us.

### The combinations with no increase

Four of the twelve pairs have intervals bracketing 1 and we exclude them from any
correction: red fox under the Twin Cities MNRRA project, California ground
squirrel and western gray squirrel under the Southern California Squirrel Survey,
and coyote under the Marin project.

One pair moved significantly the other way. Red fox under the Twin Cities Coyote
and Fox Project has a selectivity of 0.54 with an interval of 0.44 to 0.66,
entirely below 1, so red fox reporting grew only about half as fast as other
mammals in those cells. That campaign names coyote, red fox and grey fox
together, and coyote is one of the two largest increases in this screen, so
within a single campaign attention appears to have moved toward one target at
another's expense. A campaign can therefore bias a species' records downward as
well as upward, and the same correction applies with the opposite sign.

### Range expansions, and why the species stay in the analysis

Two pairs were set aside: eastern fox squirrel under the Utah campaign, where
selectivity reaches 11.78, and the same species under the Southern California
Squirrel Survey at 2.79. Both are the largest values in the whole screen.

We keep these species and their records in the model and apply no campaign
correction to them, for two reasons. A newly arrived species draws reports
because it is unfamiliar, so the surge would occur without any campaign, and a
correction keyed to the campaign would claim credit for an effect the campaign
did not cause. More importantly, the underlying increase may be a real range
expansion, which is a population change of exactly the kind this project exists
to detect. A correction here risks deleting the signal rather than a bias.

The cost of that decision is that trends for a species expanding into a region
carry a reporting-driven component we cannot separate, and that belongs in the
limitations of any result for such a species.

## Proposal: where this goes in the model

### The term

The campaign effect is a property of human reporting, so it belongs in the
opportunistic submodel, alongside the two parameters that already describe
reporting, and nowhere near the camera likelihood or the abundance surface:

```
log(mu[g,t]) = theta0 + theta1 * log(intensity[g,t]) + gamma * exposure[g] * active[t]
```

`exposure[g]` is the cell's area fraction inside the campaign region, `active[t]`
is 1 in the campaign's years and 0 otherwise, and both are supplied as data.
`gamma` is the campaign's effect on log reporting intensity.

Because each species is fitted in its own model, `exposure` and `active` are a
single vector of cell-year values per species, zero everywhere outside that
species' campaign footprint.

### Why this term is identifiable when a free regional effect is not

A regional coefficient estimated freely has to find both the regions and their
values in the same spatial pattern that abundance also explains, and in cells
without cameras the two are indistinguishable. This term is different: its
support in space, time and species is fixed by an external document. For
abundance to mimic it, a population would have to step up in the exact year a
museum launched a project, in the exact counties it covers, for the exact species
it names, and nowhere else, while the camera stream in those same cells showed
no matching change.

### Supply the effect size as a prior, do not estimate it freely

Three of the five footprints cover three to five cells. A coefficient estimated
from three cells will be weakly identified, and a vague prior would let it drift
to absorb whatever those cells do.

We propose instead a normal prior on `gamma` centred on the estimate this
analysis produced, with a standard deviation from its 95% interval. For the
Twin Cities coyote effect that is a prior centred on log(1.88) with a standard
deviation near 0.15. The external estimate then constrains the term, the fitted
data can move it if they disagree strongly, and the correction cannot run away.

### Sequence

Fit each affected species twice, once with the term and once without, and report
both. The comparison is the same shape as the camera-only refits already in the
supplement: it shows whether the trend depends on the correction rather than
asserting that the corrected number is right. Coyote is the first case, because
it has the largest effect and because its increase is one of our headline
results.

### What to do if the term is rejected

If collaborators prefer not to add a parameter, the fallback is a sensitivity
refit with the 58 affected cell-years dropped from the opportunistic stream for
the species concerned. That needs no new parameters and answers the same
question, at the cost of discarding data rather than modelling it. We would
report the result either way.

## Limitations

Twelve campaign-and-species pairs is a small screen. The selectivity estimates
carry wide intervals, and the finding that campaigns did not recruit new
observers rests on medians across twelve cases rather than on a significance
test.

Campaign regions were translated into counties by hand from the geography each
source states. A campaign described as covering a metropolitan area was given
that area's core counties, which may be wider or narrower than where the campaign
actually reached.

Three of the five campaign regions are smaller than a single model cell or close
to it, so their exposure is spread thinly across three to five cells. The
correction cannot be bounded more tightly than the grid the model runs on.

We cannot attribute individual records to individual campaigns. iNaturalist
project membership is not in our download, so the link between each campaign and
the shift we measure rests on the campaign's documented region, species and start
year coinciding with a change in the data, not on knowing which records came
through the project.

The four-year windows either side of each start year mean a campaign whose effect
grew slowly, or faded, is measured as an average over that window rather than as
a trajectory.

Campaigns active before 2012 or starting after 2021 were not measurable, because
one of the two windows falls outside our record or has too few records to
compare.

## Glossary

**Active window.** The years a campaign ran, from its stated start to its stated
end or to 2025 where it continues.

**Cell.** A 50 km grid square, the spatial unit of the model.

**Effort measure.** The count of all research-grade mammal records in a cell and
year, used to scale the expected count of any one species there.

**Exposure.** The fraction of a cell's area inside a campaign region.

**Range-expansion campaign.** A campaign targeting a species new to its region,
where reporting interest and genuine population change cannot be separated.

**Research grade.** An iNaturalist record with a photograph, date, coordinates
and community agreement on the identification.

**Selectivity.** Target-species growth divided by other-mammal growth in the
campaign's cells, divided by the same ratio outside them. 1 means no change in
the mix.
