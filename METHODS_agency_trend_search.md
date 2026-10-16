## Comparison with state agency population trends

We compared our modelled trends against population trends published independently
by state and provincial wildlife agencies. Agency trends rest on different data
and different methods from ours, principally aerial and ground surveys, mark
recapture, sightability models and mandatory harvest reporting, so agreement
between the two is concordance between independent estimates rather than
validation of either.

### Defining the jurisdictions where each species occurs

We first established, for each species, the set of jurisdictions in which it
occurs, because a state that reports no trend for a species it holds is a gap in
monitoring rather than a missing record. We intersected each species' fitted
range mask, the union of its IUCN extant-and-native range polygon and the 100 km
grid cells holding research-grade iNaturalist records (see Range masks above),
with the boundaries of the 48 conterminous states and the District of Columbia
(US Census Bureau 2022 cartographic boundaries, 1:20,000,000). We discarded
intersections covering less than 1% of a species' total range, which arise where
a grid cell clips a state corner. This yielded 745 species-by-jurisdiction occurrence
records across 34 species, with a median of 21 occupied jurisdictions per species
(range 4 to 39), reduced from 1,017 before the 1% rule.

### Search and extraction

For each species and occupied jurisdiction we searched agency publications for a
stated population trend, including annual big game and furbearer program reports,
job completion reports, species management plans, harvest summaries, commission
regulatory filings and agency species accounts. For each record we extracted the
reported direction of change over approximately ten years and over approximately
five years, the most recent reporting year, the metric the agency used, the
quoted basis for the trend, and a citation with its URL.

We recorded a direction only where a source states one, and did not derive
direction ourselves from a reported population series. This keeps the comparison
independent: the agencies' own classifications, not our reading of their numbers,
are the quantity we compare against. Where an agency reported by herd unit, data
analysis unit or management zone rather than statewide, we recorded the trend at
that level and retained the reporting level as a field. Where a species occupies a
jurisdiction and we found no reported trend, we recorded that jurisdiction with
direction unknown rather than omitting it.

### The comparison window

We use the approximately five-year agency window for all comparisons. Our camera
trap record spans 2019 to 2025, whereas the ten-year agency window covers
approximately 2014 to 2024, and the two windows can disagree in direction. For
pronghorn (*Antilocapra americana*), agency reports give 5 of 17 jurisdictions
(states and provinces) increasing and 1 decreasing over ten years, against 2 increasing and 5 decreasing
over five. Montana's statewide estimates illustrate the same reversal directly:
159,923 animals in 2019, 168,821 in 2020, then 132,359 in 2021. Several agencies
attribute declines in this period to harsh winters, poor fawn survival and
disease.

### Metric type

We classified every record as either a population or survey estimate or a harvest
index, and we do not pool the two. Harvest indices confound change in animal
abundance with change in hunter and trapper effort and in pelt prices, a
limitation that several agencies state in their own reports. Agency definitions
also vary between jurisdictions in reporting window, survey method, and what
counts as stable, so we compare direction categories and do not compare
magnitudes. Where percentage change is reported it is too sparse to use: among
state level records for moose (*Alces alces*) and white-tailed deer (*Odocoileus
virginianus*), 2 of 20 and 7 of 46 rows respectively carry a five-year
percentage.

### Automated compilation and its verification

We assembled these records using automated large language model agents rather
than human literature review. Each agent received a written instruction
specifying the species, the jurisdictions to search, the fields to populate, the
schema to return, and the rule that a direction be recorded only where a source
states it. We supply the full instructions as Supplementary Material S_n, which
serves the same role as a published search string and allows the compilation to
be repeated.

Automated extraction carries a failure mode that human review does not: a
language model can return a plausible citation that does not exist, or attribute
a real statement to the wrong source. We encountered both while assembling the
reference list for this paper, in one case receiving the correct title with the
wrong authors, year and digital object identifier. We therefore required every
record to carry a resolvable URL and treated the compiled table as unverified
until checked against its sources by an author. [AUTHOR] verified N of M records,
correcting P; verification status is recorded per record in the `url_check` field
of the supplied data file, and records failing verification are [retained with a
flag / excluded]. We report the provenance because it bears on the weight the
comparison can carry. The agency panel is an independent reference compiled by
automated search and checked by hand, and is not a systematic review.

### Coverage

Coverage varies by species and we report it per species rather than pooled
(Table_n). Monitoring is close to complete for the most intensively managed
species and sparse for the rest: 35 of 36 occupied jurisdictions carry a
five-year direction for white-tailed deer, against 14 of 18 for moose, 11 of 16
for pronghorn and 8 of 11 for grey wolf (*Canis lupus*), leaving 3.5%, 11.0%,
22.0% and 19.1% of each species' range respectively without a published trend.
For N of the 34 species in this study no agency compiles a population trend
anywhere in the range, because they are not managed as game or furbearers, and no
agency comparison is possible for them.
