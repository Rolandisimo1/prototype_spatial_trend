# iNaturalist reporting behaviour for the California ground squirrel, and what it means for our effort measure

## Take-home points

1. **Observer specialisation does not explain the California ground squirrel's
   misalignment.** Across 16,644 people who have posted a research-grade
   California ground squirrel record, only 3 have ten or more mammal records of
   which this species is more than 80%, and those three hold 42 records between
   them, 0.09% of the species total. Removing every specialist observer in
   the entire dataset changes the species' share of records in 1 of 481 western
   grid cells by more than 2 percentage points. A filter on observer name would
   not move this species' result.
2. **The specialisation pattern is real, but it sits on other species.** Of the
   39,042 people with ten or more research-grade mammal records, 873 (2.2%)
   have one species accounting for more than 80% of their records. They hold
   31,264 records, 1.39% of the eligible total. The species they concentrate on
   are the eastern gray squirrel (379 observers) and the eastern fox squirrel
   (150), not the California ground squirrel (3).
3. **Two targeted reporting programs do operate in exactly the regions where
   records concentrate.** The California Ground Squirrel Census, run by the East
   Bay Regional Park District, and the Southern California Squirrel Survey, run
   by the Natural History Museum of Los Angeles County since 2011, both recruit
   the public to report squirrels. The Bay Area holds 23% of California's
   records for this species and the southern coast holds 46%. The programs
   recruit observers who also report other species, so their effect falls on
   where effort goes, not on who reports.
4. **Domestic animals enter the effort denominator unfiltered.** The pipeline
   drops records flagged captive, but no record of a domestic species in this
   download carries that flag: all 32,015 domestic cat records, 7,399 horse
   records and 7,460 Norway rat records are research grade and count as mammal
   observation effort. Together, domestic and commensal species contribute
   59,666 records, 1.98% of the download. Two of the largest specialist
   observers post domestic cats exclusively, 801 and 502 records each.
5. **We apply no spatial thinning and no per-observer deduplication.** Effort is
   the raw count of every record in a cell and year. 6.6% of this species'
   records fall at the 150 point locations holding ten or more records, the
   busiest a single observer with 178 records at one coordinate.
6. **Thinning would make this species' misalignment worse, not better.**
   Capping each observer at one record per cell per year retains 56.2% of
   California ground squirrel records but only 40.4% of the effort records in
   the same cells, raising the species-to-effort ratio by 39%. The species'
   records are less clustered per observer than mammal records in general, so
   any cap removes more denominator than numerator.

## Background and definitions

**Effort measure.** For each 50 km grid cell and each year, the model offsets
the opportunistic record count by the total number of research-grade mammal
records in that same cell and year. The reasoning is that a place and year with
many mammal records of any kind had many people looking, so a species' record
count there should be judged against that background. Everything in this report
tests whether that background is a fair measure of looking.

**Spatial coupling exponent (theta_space).** The exponent linking expected
iNaturalist records to camera-derived abundance in a cell, on a log-log scale.
A value of 1 means records rise in proportion to abundance. Below 1, records
under-respond. Above 1, records over-respond: a place with twice the animals
gets more than twice the records. The California ground squirrel sits at 2.24,
with a 95% interval of 1.60 to 3.16, the highest of the 25 fitted species. At
2.24, a tenfold difference in abundance produces a 174-fold difference in
records.

**Specialist observer.** For this report, a person with at least ten
research-grade mammal records whose single most-reported species accounts for
more than 80% of them. The threshold follows the question as posed. Both parts
matter: the ten-record floor excludes people with one or two records, for whom a
100% share carries no information.

**Spatial thinning.** Discarding records that are close together in space, or
close in space and time, so that a dense cluster contributes less than its raw
count. It is standard in species distribution modelling, where the aim is to
stop heavily visited places from dominating a fitted habitat relationship.

**Research grade.** An iNaturalist record with a photograph, a date, coordinates
and community agreement on the identification. All records in this analysis are
research grade.

## Data and methods

### Records

The analysis reads the project's iNaturalist download of North American mammals,
3,014,977 records. Filtering to research grade, class Mammalia, coordinates
present, and observation years 2008 to 2025 leaves 2,967,475 records from
371,805 observers covering 1,323 species. Applying the pipeline's
positional-accuracy filter, which keeps records with a stated accuracy better
than 1,000 m, leaves 2,351,075.

For the California ground squirrel, *Otospermophilus beecheyi*, the download
holds 51,461 research-grade records; 45,351 survive the accuracy filter. 49,958
of the unfiltered records are in California, 927 in Baja California and 557 in
Nevada. 2,489 have obscured coordinates.

### Assigning records to cells

Records within the western extent, longitude -125 to -112 and latitude 30 to 46,
were assigned to their nearest 50 km cell centroid in an equal-area projection,
accepting matches within 40 km. 488,674 of 492,435 western records matched a
cell. Cells with fewer than 20 mammal records were excluded from cell-level
statistics, leaving 481 cells, of which 115 hold at least one record of this
species.

### Observer specialisation screen

For every observer, the analysis counted research-grade mammal records by
species, took the most-reported species, and divided its count by the observer's
total. Observers with at least ten records and a dominant share above 80% were
flagged. The screen ran across all 1,323 species and all 371,805 observers, not
only on the California ground squirrel.

To test whether these observers drive the result, the analysis recomputed each
western cell's species share with all flagged observers' records deleted from
both the numerator and the effort denominator, and compared it to the share
computed from the full data.

### Thinning scenarios

Thinning was evaluated as a cap on records per observer, per cell, per year. For
caps of 1, 3, 5 and 10, the analysis computed how many records the species
retains and how many the effort denominator retains, and formed the ratio of the
two retention rates. A ratio above 1 means the cap removes proportionally more
effort than species records, which raises the apparent record rate.

### Programs

Reporting programs were identified by web search. No membership information
exists in the download, so the analysis cannot attribute individual records to
individual projects.

## Results

### Where the records are

![California ground squirrel records and density]({{artifact:art_PLACEHOLDER}})

Three regions hold 91% of California's records: the southern coast around Los
Angeles and San Diego with 46%, the Bay Area with 23%, and the Central Valley
floor with 22%. The northeast of the state holds 9 records across Modoc, Lassen
and Siskiyou counties, which is effectively unsampled. The Mojave and the desert
southeast hold 3.4%.

The density grid shows the same concentration at a coarser grain. Of 775
occupied hexagons of roughly 30 km across, 445 hold fewer than ten records,
while 88 hold at least a hundred and carry 78% of all records. The densest
hexagon holds 2,146.

This spread is what the coupling exponent of 2.24 measures. The model sees cells
where the species is recorded hundreds of times beside cells where it is
recorded once or twice, and the camera data do not show a comparable spread in
abundance between those same cells.

### Who reports them

The records come from 16,644 observers. The ten most prolific account for 11.2%
of records and the fifty most prolific for 23.8%. The single most prolific
observer has 1,771 records of this species out of 3,681 mammal records, a 48%
share.

| Observer | Records of this species | All mammal records | Share of own records |
|---|---|---|---|
| mooseandsquirrel | 1,771 | 3,681 | 48% |
| dalemeister | 678 | 2,173 | 31% |
| bewickwren | 548 | 6,736 | 8% |
| dgrimmphd | 496 | 1,839 | 27% |
| bethpratt1 | 408 | 10,467 | 4% |

These are prolific generalists. The most prolific reporter of this species posts
roughly half their mammal records as other species, and the fifth most prolific
posts 96% other species. Three observers meet the specialist definition.
They hold 42 records of this species between them, 0.09% of the 45,351 that pass
the accuracy filter, from totals of 11, 23 and 15 mammal records each. A further 45
records of this species come from observers who specialise on something else, so all
873 specialists together account for 87 records, 0.19%.

Concentration does appear at cell level. In 6 of the 89 cells holding at least
ten records of this species, a single observer contributed more than half of
them. The most extreme is a cell where one observer contributed 88.8% of 233
records. These dominant observers are generalists by the 80% test, so a
specialist filter does not reach them.

### What a specialist filter would do

Deleting all 873 specialist observers from both numerator and denominator
removes 81 of 44,732 records of this species, 0.18%, and 5,526 of 488,674
effort records, 1.13%. The median cell's species share does not move. One cell
of 481 moves by more than 2 percentage points, and the largest movement in any
cell is 2.3 percentage points.

Across the fleet, the specialist share of a species' records runs from 0.01% for
the black-tailed jackrabbit to 3.39% for the eastern gray squirrel:

| Species | Records | From specialists | Share |
|---|---|---|---|
| Eastern gray squirrel | 251,566 | 8,523 | 3.39% |
| Eastern fox squirrel | 115,730 | 3,102 | 2.68% |
| American black bear | 38,730 | 1,022 | 2.64% |
| Coyote | 80,507 | 1,536 | 1.91% |
| White-tailed deer | 252,077 | 4,203 | 1.67% |
| California ground squirrel | 45,351 | 87 | 0.19% |

The two eastern tree squirrels are where a specialist filter would have any
measurable effect, and even there it would move about 3% of records.

### The reporting programs

The California Ground Squirrel Census asks the public to report where ground
squirrel colonies are in East Bay Regional Park District parks, on the stated
grounds that the species is a keystone of grassland habitats and that its
distribution is poorly documented. The Southern California Squirrel Survey,
started in 2011 at the Natural History Museum of Los Angeles County, catalogues
the occurrence of all squirrel species across the greater Southern California
region, and compares those records to the museum's specimen collection.

Both programs sit in the regions that hold the most records: the East Bay inside
the Bay Area's 23% and the Los Angeles program inside the southern coast's 46%.
The Southern California survey also covers the western gray squirrel and the
eastern fox squirrel, both in our fleet.

The programs recruit observers who report other species as well, which is why
the specialist screen finds so few. Their effect on our effort measure works
through geography: they raise reporting of this species in two regions without
a matching rise in reporting of other mammals there. A cell-level observer
filter cannot separate that from genuine local abundance.

### Thinning

The pipeline applies no spatial thinning and no per-observer deduplication. The
effort offset for a cell and year is the count of every record in it.

Duplication is present but moderate at the record level. The 44,732 records in
the western extent sit at 35,048 distinct coordinates rounded to four decimal
places, about 11 m, which is 1.28 records per location. Concentration is
stronger at the observer level: the busiest observer-cell-year combination holds
392 records of this species, against a median of 1 and a 95th percentile of 4.
The equivalent figure for all mammals is a maximum of 1,918.

| Cap per observer, cell and year | Records of this species kept | Effort records kept | Effect on the record-to-effort ratio |
|---|---|---|---|
| 1 | 56.2% | 40.4% | ×1.39 |
| 3 | 73.6% | 60.8% | ×1.21 |
| 5 | 79.3% | 68.4% | ×1.16 |
| 10 | 85.5% | 76.9% | ×1.11 |

Every cap raises the ratio. The effort denominator is built from all mammal
records, which include the heavily repeated reporting of common urban species,
so capping removes more denominator than numerator. For this species a cap of 1
would raise the apparent record rate by 39%, pushing the coupling exponent
further above 1.

Thinning across years raises a second problem. Camera effort grew from 1
deployment in 2008 to 3,985 in 2025, and the iNaturalist record pool grew
roughly 154-fold over the same period. A cap applied within each year removes
nothing from a sparse early year and a great deal from a dense recent year, so
it imposes a trend on the effort series. Applying a single cap across the whole
period instead discards the temporal structure the trend is estimated from.

## What to do

**Exclude domestic and commensal species from the effort denominator.** Domestic
cat, dog, horse, cattle, sheep, goat, domestic pig, European rabbit, house mouse
and the two commensal rats contribute 59,666 records that count as mammal
observation effort. None are flagged captive in this download. A person
photographing feral cats in a city is not measuring effort toward wild mammals.
This is a one-line change to the species list and affects 1.98% of records
nationally, more in cities.

**Do not add an observer-name filter for this species.** The test above puts its
effect at 0.18% of records and no change to the cell-level pattern.

**Do not thin.** Both the within-year and across-year forms make the quantity we
report worse.

**Treat the two programs as a regional caveat.** The honest statement is that
this species' records are inflated relative to other mammals in the Bay Area and
Southern California by organised reporting campaigns, that we cannot measure the
size of that inflation from the download, and that its direction is the same as
the direction of the species' coupling anomaly.

**If a correction is wanted, work on geography rather than observers.** A term
allowing the record-to-abundance relationship to differ between regions would
absorb a program footprint. That is a model change, and the simulation study
would need an arm that generates region-specific reporting intensity to show
whether it recovers the trend.

## Limitations

The analysis cannot attribute records to projects. iNaturalist project
membership is not in the download, so the link between the two programs and the
record concentration is geographic coincidence plus the programs' stated aims,
not a measured attribution. Querying the iNaturalist API for project membership
would close this.

The 80% and ten-record thresholds were taken from the question as posed and were
not varied. A lower threshold would flag more observers, and the screen does not
show how the conclusion changes across thresholds.

The cell-level test deletes specialist records from both numerator and
denominator, which is what a filter in the pipeline would do. It does not refit
the model, so it shows the change in the input quantity and not the change in
the estimated trend.

Camera-derived abundance was not recomputed here. The statement that the camera
data show no comparable spread between high-record and low-record cells rests on
the fitted coupling exponent of 2.24, not on a fresh comparison.

The domestic-species list covers eleven species chosen as the clear cases. Other
introduced or human-associated mammals in the download were not assessed.

## Glossary

**Cell.** A 50 km grid square, the spatial unit of the model.

**Effort denominator.** The count of all research-grade mammal records in a cell
and year, used to scale the expected record count for one species.

**Research grade.** An iNaturalist record with photo, date, coordinates and
community agreement on identification.

**Specialist observer.** Here, someone with at least ten research-grade mammal
records whose top species is more than 80% of them.

**Spatial thinning.** Discarding records clustered in space, or in space and
time, so dense places contribute less than their raw count.

**theta_space.** The exponent linking expected records to abundance, on a
log-log scale. 1 is proportional, above 1 over-responsive, 0 no spatial
correspondence.
