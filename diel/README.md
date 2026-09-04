# Continental variation in mammal diel activity

Analysis code for a study of how the daily activity schedules of 13 North American
mammal species vary across the conterminous United States, using compiled
camera-trap records from 49 archived projects.

## Status

**UNVALIDATED as a standalone package.** These files were written and run inside a
working analysis session against data paths and helper modules that were live at the
time. They are committed here as the provenance record of what was actually run, not
as a turnkey pipeline. Nothing in this repository has been tested from a clean
checkout. Expect hard-coded paths and missing glue.

The scientific results these files produced were checked — reproduction gates,
cross-validation against nulls, and simulated-null tests are described in the
analysis reports — but that is a different claim from "this code runs."

## Data

Source records are **not** in this repository. They live outside it (roughly 550 MB
of deployment and sequence tables) and are excluded by `.gitignore` along with all
intermediate CSV, parquet and raster outputs.

## Layout

    lib/          Reusable modules: solar time, harmonic curve fitting, spatial
                  statistics, cross-validation, multivariate inference,
                  counting-noise simulation
    pipeline/     Analysis stages in order: binning, curve fitting, variance
                  partition, cross-validation, prediction surfaces
    covariates/   Earth Engine extraction, national prediction grid, seam repair
    figures/      Figure builders
    reports/      HTML and memo document generators
    viewer/       Interactive map page builder and its verifiers
    methods/      Methods text as submitted, including covariate sources

## Method summary

Detections are filtered to independent events, placed on sun-anchored time (sunrise
and sunset at fixed anchors), and binned into 48 half-hour bins. A negative binomial
model with two harmonic terms is fitted per species and site, carrying an exact
per-bin effort offset -- effort is *not* uniform in sun-anchored time, and treating
it as uniform biases curve shape.

Sites are camera clusters at 25 km, not individual cameras. Among cameras that ran
simultaneously, the median distance to the nearest concurrent neighbour in the same
site is 324 m (IQR 204-765). Ignoring dates gives 18 m, but that is dominated by the
same locations re-occupied in different years, which are not simultaneous samples.

Timestamp quality control is a substantial part of the work. Clock errors are found
by screening with fixed-schedule reference species against locally-correct sunrise
and sunset, confirmed by requiring diurnal and nocturnal species to move in opposite
directions at the same camera, and verified against photographs. Confirmed failures
are **excluded** in this pipeline rather than corrected by an offset, because fitting
an offset per camera improves apparent fit even on cameras known to be sound. Any
timestamp corrections applied upstream by the data owners are a separate matter.

## Model specification

The primary inference is one multivariate model per species and measure with eight predictors, each
mapping to a single mechanism:

    human disturbance    population density (1 km, log)
    agriculture          cropland share (5 km)
    vegetation           tree canopy cover (1 km)
    terrain              ruggedness (5 km)
    summer heat          max temperature of hottest month (1 km)
    predator richness    carnivore species detected on camera at the site
    predator abundance   summed carnivore detection rate (log)
    relative abundance   focal species detections per camera-hour (log)

Max VIF 1.50, no predictor pair above r = 0.5. 61 of 65 species-measure combinations are fitted;
four fail a reliability threshold of 0.5 and are reported as not assessed.

Regression weights are inverse squared standard error of each site's measure, CAPPED at the 90th
percentile within each model. Uncapped, 10 of 61 models had an effective sample size below 10 sites
because the standard error is smallest where a bounded measure carries least information.

Nighttime lights and developed-land share were dropped from the study: they correlate with population
density at 0.93 and 0.88 and are one gradient measured three ways, so disturbance by human presence
cannot be separated from artificial light with these covariates.

Mapping uses only the five environmental covariates, refitted on grid-extracted values at matching
spatial supports. Predator and abundance covariates are not on the grid, and relative abundance is a
property of the survey rather than of the landscape.

## Three traps documented here because they cost time

**Solar time versus clock time.** Computing solar sunrise and comparing it against a clock timestamp
is wrong by the site's longitude offset within its time zone plus daylight saving, up to an hour.
`lib/solarclock.py` returns sunrise and sunset in local clock hours using true IANA zones and agrees
with published times to within 14 minutes across four zones. An earlier version of this helper had
the bug and its results were retracted.

**Circular quantities.** Time of peak activity is circular. A linear difference across midnight is
wrong by up to 12 hours and inflates spatial statistics. This produced two separate bugs here, one in
variogram pairwise differences and one in a cross-validation error metric.

**Inverse-variance weight concentration.** For a measure that approaches a bound at many sites, the
standard error is smallest where the value carries least information about variation among sites, so
raw weights collapse onto a few sites. Capping at the 90th percentile fixes it and leaves
coefficients essentially unchanged (r = 0.97, no sign reversals). Applied from the eight-predictor
run onward.

## Related

Covariate extraction methods and full data-source citations are in the accompanying
methods documents rather than in this repository.
