# Materials and Methods

## Camera-trap records

We compiled camera-trap records from 49 archived projects, principally the Snapshot USA coordinated
national survey from 2019 through 2025 (Cove et al. 2021; Kays et al. 2022; Shamon et al. 2024;
Rooney et al. 2025) and North Carolina's Candid Critters, accessed through Wildlife Insights and the
eMammal archive. Camera trapping conventions and their known biases follow Burton et al. (2015). We
restricted the analysis to 1 August through 31 October, the only window with continental coverage in
every year, and to the conterminous United States. The final dataset comprises 1,210,419 independent
detections of 13 species across 18,552 deployments at 489 sites. We treated repeated photographs of
the same species at the same camera as one detection when separated by less than 30 minutes.

Species were retained when they had enough sites to fit site-level activity curves. We screened for
records falling outside a species' published longitude range, which removed 5 detections in total:
3 mule deer (*Odocoileus hemionus*) records east of -96 degrees, from urban Tampa, Florida, the South
Carolina Piedmont, and coastal Rhode Island, and 2 eastern chipmunk (*Tamias striatus*) records west
of -104 degrees. No white-tailed deer (*Odocoileus virginianus*) record was affected, and the two
deer species were not otherwise treated differently. Given the scale, this screen has no bearing on
any result and we report it for completeness.

## Timestamp quality control

Camera clocks fail, and a shifted clock produces an activity curve that no downstream model can
identify as wrong. We screened every camera using species with strongly stereotyped daily schedules
as internal references, comparing their detections against locally correct sunrise and sunset.
Sunrise and sunset were computed in local clock hours using the IANA time zone database, because
comparing solar sunrise against a clock timestamp is wrong by the site's longitude offset within its
time zone plus daylight saving, which reaches an hour at zone edges. Our implementation agrees with
published times to within 14 minutes across four zones.

A clock error displaces every species at a camera together, so we required a candidate camera to
show diurnal and nocturnal species shifted in opposite directions before flagging it, then verified
flagged cameras against the photographs themselves. We excluded 72 confirmed deployments, 0.51
percent of detections. Nine cameras retained after photographic review still show 20 to 67 percent
of strictly-daytime detections at night; they contribute 1,021 detections, 0.084 percent of the
analysis window.

Removing the 72 deployments moved 105 of 1,929 curves by more than half a percentage point of night
activity. Of these, 48 of 49 diurnal-species curves became less nocturnal and all 36
nocturnal-species curves became more nocturnal (sign test P = 1.8e-13), confirming that the excluded
clocks were displaced rather than merely noisy.

## Spatial unit of analysis

Cameras in this dataset are deployed in clusters, so a camera is not an independent spatial
replicate. We defined sites by complete-linkage clustering of deployments at 25 kilometres, giving
489 sites with a maximum extent of 24.9 kilometres and a minimum separation of 8.5 kilometres. Where
a survey re-occupied the same locations in more than one year, we pooled those deployments into the
one site.

Restricting attention to cameras that operated simultaneously, the median distance from a camera to
its nearest concurrent neighbour within the same site is 324 metres (interquartile range 204 to 765
metres). Ignoring dates gives a median of 18 metres, but that figure is dominated by the same
physical locations re-occupied in different years, which are not simultaneous samples of the
landscape. Either way the cameras within a site are close enough that they sample overlapping
animal populations and shared habitat, so we treat the site as the unit of replication throughout.

## Diel activity curves

We placed each detection on sun-anchored time, rescaling each day so that sunrise and sunset fall on
fixed anchors (Vazquez et al. 2019), and binned detections into 48 half-hour bins. For each species
and site we fitted a negative binomial generalized linear model with two harmonic terms, following
the hierarchical trigonometric approach of Iannarilli et al. (2024), using mgcv 1.9 (Wood 2011,
2017) in R 4.4 (R Core Team 2024).

Every model carries an exact analytic per-bin effort offset. That offset is necessary rather than
cosmetic: per-bin effort is uniform in clock time but not in sun-anchored time, because rescaling to
fixed anchors compresses and stretches the day and night segments unequally, with within-camera
ratios of maximum to minimum per-bin effort reaching 3.3. Omitting the offset biases the fitted
curve shape, and because day length covaries with latitude and season the bias would appear as a
plausible latitudinal gradient.

## Activity measures and their reliability

From each fitted curve we derived five measures: the share of activity at night, within the dawn and
dusk windows, and within two hours of solar noon; the concentration of activity, computed as the
resultant length of the fitted rate curve on the 24-hour circle; and the clock time of the activity
peak. Peak time is a circular quantity, so all differences, means, variances and intervals on it use
circular statistics (Fisher 1993; Jammalamadaka and SenGupta 2001).

Each site's measure is an estimate, not an observation, and carries a delta-method standard error
from the fitted harmonic coefficients. Two sites can therefore differ because their animals differ or
because one site had fewer photographs from which to estimate the curve. We separated these by
partitioning the observed variance among sites. Writing the observed between-site variance of a
measure as var(obs), and the mean squared standard error across sites as var(err), the variance
attributable to real differences among sites is var(true) = var(obs) - var(err), and we define
reliability as var(true) / var(obs). A reliability of 1 means every difference among sites is real;
a reliability of 0 means the apparent spread is entirely estimation error.

We required reliability above 0.5 for a species-measure combination to enter the analysis. Median
reliability across the 13 species was 0.95 for peak time, 0.85 for night share, 0.77 for midday
share, 0.68 for dawn and dusk share, and 0.68 for concentration. Overall activity level, in the sense
of Rowcliffe et al. (2014), reached a median reliability of only 0.31 and passed in 1 of 13 species,
so we excluded it from all spatial analysis. Four further combinations fell below the threshold and
are reported as not assessed rather than fitted: American black bear peak time (reliability 0.34,
71 sites), Virginia opossum midday share (0.32, 152 sites), eastern chipmunk night share (0.48, 82
sites), and wild turkey night share (0.48, 118 sites).

For design guidance, the standard error of every measure falls below half the typical between-site
spread once a site accumulates roughly 100 to 200 detections of the focal species. At fewer than 50
detections the median standard error is 4.0 percentage points for night share against a typical
between-site spread of 7.8 points, so sites below that level contribute little information about
spatial variation. We therefore recommend 100 detections per species per site as a working minimum
for this class of analysis, and we retained sites with at least 25.

## Environmental covariates

We extracted all covariates in Google Earth Engine (Gorelick et al. 2017) at each camera location
and as the mean within 1 and 5 kilometre radius buffers, matching land cover, population and
nighttime lights to each deployment's own epoch and computing temperature against its own active
date window. We aggregated camera-level values to the site by taking the unweighted arithmetic mean
across all cameras assigned to that site, with a median of 13 cameras contributing per site.

Table 1 lists the covariates used. We report only those entering the final models; covariates
extracted during screening and subsequently dropped are not described. For each covariate we used a
single spatial support, chosen as the one closest to the scale at which the mechanism is expected to
act: 1 kilometre for population density and tree canopy, and 5 kilometres for cropland and
ruggedness. We did not enter the same layer at two supports.

**Table 1. Environmental covariates entering the models.**

| Mechanism | Covariate | Support | Product | Citation |
|---|---|---|---|---|
| Human disturbance | Population density | 1 km | GHS-POP R2023A | Schiavina et al. 2023 |
| Agriculture | Cropland share | 5 km | NLCD 2021 | Dewitz 2023 |
| Vegetation | Tree canopy cover | 1 km | NLCD TCC 2021 | Coulston et al. 2012 |
| Terrain | Ruggedness | 5 km | Copernicus DEM GLO-30 | European Space Agency 2021 |
| Summer heat | Max temperature, hottest month | 1 km | Daymet V4 | Thornton et al. 2022 |

Table 2 gives the pairwise correlations among all model predictors.

**Table 2. Pairwise Pearson correlations among the eight model predictors, on the scale at which each
enters the model, over the 1,747 site-species rows complete on all eight.**

| | Pop | Crop | Canopy | Rugged | Heat | PredRich | PredAbund |
|---|---|---|---|---|---|---|---|
| Cropland | -0.11 | | | | | | |
| Tree canopy | -0.32 | -0.26 | | | | | |
| Ruggedness | -0.15 | -0.29 | 0.16 | | | | |
| Summer heat | -0.12 | 0.04 | -0.14 | -0.34 | | | |
| Predator richness | -0.16 | -0.23 | 0.26 | 0.34 | -0.13 | | |
| Predator abundance | 0.36 | -0.06 | -0.28 | 0.08 | 0.11 | 0.08 | |
| Relative abundance | 0.13 | 0.12 | -0.08 | -0.12 | 0.01 | -0.25 | 0.07 |

Pooled across all rows the largest absolute pairwise correlation is 0.36, between population density
and predator abundance, and the maximum variance inflation factor is 1.50. Because each
species-measure model is fitted on its own subset of sites, we also report identification per model:
the median maximum variance inflation factor is 1.78 and the largest is 3.55, and the median largest
pairwise correlation is 0.48 with a maximum of 0.78. All remain below the conventional threshold of 5
for variance inflation (Dormann et al. 2013), but the per-model figures are the relevant ones and are
appreciably higher than the pooled summary, so individual partial coefficients in the smaller models
are less sharply identified than the pooled diagnostics alone would suggest.

The candidate covariate set originally contained three measures of human pressure: population
density, nighttime lights and developed-land share. These correlate at r = 0.88 to 0.93, carry
variance inflation factors of 11.2, 8.3 and 7.6 when fitted together, and their effect vectors
across species correlate at 0.94 to 0.97. They are one gradient measured three ways. We represent
human disturbance by population density alone and dropped the other two; the distinction between
disturbance by human presence and the effect of artificial light at night is therefore not
addressed here.

## Predator and relative-abundance covariates

We derived three count-based covariates from the camera records themselves.

*Predator richness* is the number of the eight large and medium carnivores in the dataset detected
at least once at that site: American black bear, bobcat, coyote, grey fox, grey wolf, grizzly bear,
puma and red fox. It ranges from 1 to 7 with a median of 3.

*Predator abundance* is the sum of the per-site detection rates of those same eight species, in
detections per 100 camera-days. It ranges from 0.13 to 123.6 with a median of 6.3. Richness and
abundance are nearly independent across sites (r = 0.05), so they represent distinct quantities:
how many kinds of predator are present, and how many predators are encountered.

*Relative abundance of the focal species* is the natural logarithm of that species' own detections
per camera-hour at the site, and serves as the density proxy in tests of density dependence.

Because all three are built from detection counts, and because the response measures are also
estimated from detection counts, effects involving them require the specific control described under
Counting-noise null below.

## Statistical model

For each species and measure we fitted a single multivariate model with eight predictors, each
representing one mechanism: human disturbance, agriculture, vegetation, terrain, summer heat,
predator richness, predator abundance, and the focal species' own relative abundance. We used
weighted least squares with weights equal to the inverse squared standard error of each site's
measure, standardized all predictors, and log-transformed population density, predator abundance and
relative abundance. Models were fitted with statsmodels 0.14 (Seabold and Perktold 2010) on NumPy
2.1 (Harris et al. 2020) and SciPy 1.14 (Virtanen et al. 2020).

We report each partial coefficient with its 95 percent confidence interval and its partial
R-squared, so that the relative contribution of mechanisms is comparable within a model.

With 13 species and 5 measures there are 65 possible species-measure combinations. Four fail the
reliability threshold described above, so we fitted 61.

Inverse-variance weighting requires one adjustment. For a measure that approaches a bound at many
sites, such as midday activity in a strongly nocturnal species, the standard error is smallest
exactly where the measure carries least information about variation among sites, and unmodified
weights then concentrate on those sites. Across the 61 fitted models, 7 had an effective sample size below 10 sites under raw
inverse-variance weights, computed as the inverse of the sum of squared normalized weights over the
rows each model actually uses. The most extreme was Virginia opossum night share, where 144 sites
gave an effective sample size of 1.0 because a single site carried 99.8 percent of the weight; the
others were eastern cottontail midday share (1.2), red fox midday share (1.5), American black bear
midday share (2.3), red fox night share (2.6), northern raccoon night share (3.4) and coyote midday
share (4.0). One further combination, white-tailed deer midday share, has an effective sample size of
1.7 over rows complete on the measure alone but 17.8 over the rows the model uses; we report the
latter throughout, since that is the sample the coefficients are estimated from. No sites were
dropped in any of these models; the weights simply became degenerate.

We therefore capped weights at their 90th percentile within each model, which raises the effective
sample size in those 7 models to between 12.6 and 109.4 sites and lifts the median across all 61
models from 38.5 to 55.3. Elsewhere the cap leaves coefficients essentially unchanged: across the 427
effects shared with an uncapped fit the coefficients correlate at 0.97 with 93.7 percent sign
agreement and no sign reversals, while 11 effects gain and 29 lose significance at P = 0.05. The cap
is the larger of the two changes made in this specification, larger than adding the eighth predictor.

## Spatial confounding

Environmental gradients in North America are geographically organized: temperature, vegetation and
human settlement all vary smoothly with latitude and longitude. A covariate can therefore appear to
explain a response simply because both vary along the same continental gradient, without the
covariate acting on the animals at all (Hodges and Reich 2010; Paciorek 2010). This matters here
because we interpret partial coefficients mechanistically.

To bound that risk we refitted every model with a low-rank radial-basis spatial smooth added, using
4 to 15 space-filling centres depending on sample size, with the bandwidth set to the median
nearest-neighbour distance among centres. A covariate effect that persists when position is in the
model is not attributable to smooth geographic position alone. An effect that disappears may still
be real, but it cannot be separated from everything else that varies across the continent, and we
therefore give it no mechanistic interpretation. We report both fits so a reader can see which
effects depend on the choice. The smooth complexity was not tuned, so these labels are conditional
on one specification.

## Counting-noise null

This control applies only to the three count-derived covariates named above, predator richness,
predator abundance and relative abundance. It is not relevant to the environmental covariates, which
are measured independently of the camera records.

The response measures are estimated from binned detection counts, and estimates from sparse counts
are biased as well as noisy: with few detections, random clumping among bins makes a fitted curve
appear more sharply concentrated than the underlying rhythm. Regressing such a measure on a
count-derived covariate can therefore produce a coefficient with no biological content. In
simulation from a single fixed rhythm, concentration read 0.395 at 20 detections against a true
value of 0.338, and a regression of concentration on detection count reached P < 0.05 in 18 percent
of runs when the true effect was zero.

We therefore built a null in which spatial variation is absent by construction. For each species we
held curve shape fixed at its pooled national harmonic shape, drew counts at each site from a
negative binomial scaled to that site's own total detections, own per-bin effort and own fitted
dispersion, refitted the harmonic curve, recomputed the measure, refitted the same multivariate
model and read off the same partial coefficient. Any coefficient arising in this null comes solely
from sites differing in how many animals were photographed. We used 400 replicates per species and
report, for each observed effect, the fraction of it the null reproduces. Effects that the null
reproduces in large part are not interpretable as behaviour.

## Spatial cross-validation and prediction

We tested out-of-sample prediction with leave-one-block-out cross-validation on contiguous 400
kilometre spatial blocks, which prevents a held-out site from being predicted by its immediate
neighbours (Roberts et al. 2017). Models were fitted by ridge regression (Hoerl and Kennard 1970)
with an unpenalized intercept and the penalty selected by nested inner-block cross-validation on
training blocks only. All design-matrix parameters, including covariate centring and scaling and the
positions of basis centres, were estimated from training rows alone, so a held-out block contributes
nothing to the model that predicts it. Skill was measured against the null of a single nationwide
value per species, with block bootstrap confidence intervals from 1,000 resamples of blocks, and a
combination counts as beating the null only when the lower confidence bound exceeds zero.

The set of combinations beating the null is sensitive to the bootstrap draw at its boundary. Point
skills are exactly reproducible, but 48 of 183 competitor-combination rows have a bootstrap lower
bound within 0.02 of zero, and membership of the beats-null set changes with the resampling seed at
that margin. We report the seed and resample count, and we treat the skill estimate and its interval
as the result rather than the count of combinations clearing a threshold.

Prediction surfaces use a reduced model containing only the five environmental covariates, because
predator richness, predator abundance and relative abundance are not available on the prediction
grid. Relative abundance is a property of the survey rather than of the landscape and cannot be
mapped in principle. Surfaces were masked to the species range where a range mask exists and to the
convex covariate envelope of the fitted sites. For each surface we report the mean absolute
difference between adjacent grid cells against the total spread across sites, because a surface
whose cell-to-cell variation approaches its total spread reflects overfitting rather than signal.

Site-aggregated and grid-extracted covariates are different quantities and they disagree: elevation
agrees between them at r = 0.97 and tree canopy at 0.75, whereas population density agrees at
r = 0.55 untransformed and 0.66 on the log scale, because a site value is a mean over camera
locations that are not randomly placed within a 25 kilometre cell whereas a grid value is a cell
average. We refitted on grid-extracted covariates
before predicting rather than applying site-fitted coefficients to grid values, and we report the
disagreement.

## Software

R 4.4 (R Core Team 2024) with mgcv 1.9 (Wood 2011) for harmonic curve fitting; Python 3.11 with
statsmodels 0.14 (Seabold and Perktold 2010), NumPy 2.1 (Harris et al. 2020), SciPy 1.14 (Virtanen
et al. 2020) and pandas 2.2 (McKinney 2010) for the multivariate models, cross-validation, null
simulations and prediction; the Google Earth Engine Python API (Gorelick et al. 2017) for covariate
extraction. Analysis code is archived at [repository].

## References

Burton, A.C., Neilson, E., Moreira, D., Ladle, A., Steenweg, R., Fisher, J.T., Bayne, E., & Boutin,
S. (2015). Wildlife camera trapping: a review and recommendations for linking surveys to ecological
processes. *Journal of Applied Ecology*, 52, 675-685.

Coulston, J.W., Moisen, G.G., Wilson, B.T., Finco, M.V., Cohen, W.B., & Brewer, C.K. (2012). Modeling
percent tree canopy cover: a pilot study. *Photogrammetric Engineering and Remote Sensing*, 78,
715-727.

Cove, M.V., et al. (2021). SNAPSHOT USA 2019: a coordinated national camera trap survey of the United
States. *Ecology*, 102, e03353.

Dewitz, J. (2023). *National Land Cover Database (NLCD) 2021 Products*. U.S. Geological Survey data
release.

Dormann, C.F., et al. (2013). Collinearity: a review of methods to deal with it and a simulation
study evaluating their performance. *Ecography*, 36, 27-46.

European Space Agency (2021). *Copernicus DEM: Global and European Digital Elevation Model, GLO-30*.

Fisher, N.I. (1993). *Statistical Analysis of Circular Data*. Cambridge University Press.

Gorelick, N., Hancher, M., Dixon, M., Ilyushchenko, S., Thau, D., & Moore, R. (2017). Google Earth
Engine: planetary-scale geospatial analysis for everyone. *Remote Sensing of Environment*, 202,
18-27.

Harris, C.R., et al. (2020). Array programming with NumPy. *Nature*, 585, 357-362.

Hodges, J.S., & Reich, B.J. (2010). Adding spatially-correlated errors can mess up the fixed effect
you love. *The American Statistician*, 64, 325-334.

Hoerl, A.E., & Kennard, R.W. (1970). Ridge regression: biased estimation for nonorthogonal problems.
*Technometrics*, 12, 55-67.

Iannarilli, F., Gerber, B.D., Erb, J., & Fieberg, J.R. (2024). A 'how-to' guide for estimating animal
diel activity using hierarchical models. *Journal of Animal Ecology*, 93, 1-13.

Jammalamadaka, S.R., & SenGupta, A. (2001). *Topics in Circular Statistics*. World Scientific.

Kays, R., et al. (2022). SNAPSHOT USA 2020: a second coordinated national camera trap survey of the
United States during the COVID-19 pandemic. *Ecology*, 103, e3775.

McKinney, W. (2010). Data structures for statistical computing in Python. *Proceedings of the 9th
Python in Science Conference*, 56-61.

Paciorek, C.J. (2010). The importance of scale for spatial-confounding bias and precision of spatial
regression estimators. *Statistical Science*, 25, 107-125.

R Core Team (2024). *R: A Language and Environment for Statistical Computing*. R Foundation for
Statistical Computing, Vienna.

Roberts, D.R., et al. (2017). Cross-validation strategies for data with temporal, spatial,
hierarchical, or phylogenetic structure. *Ecography*, 40, 913-929.

Rooney, B., et al. (2025). SNAPSHOT USA 2019-2023: the first five years of data from a coordinated
camera trap survey of the United States. *Global Ecology and Biogeography*, 34, e13941.

Rowcliffe, J.M., Kays, R., Kranstauber, B., Carbone, C., & Jansen, P.A. (2014). Quantifying levels of
animal activity using camera trap data. *Methods in Ecology and Evolution*, 5, 1170-1179.

Schiavina, M., Freire, S., Carioli, A., & MacManus, K. (2023). *GHS-POP R2023A: GHS population grid
multitemporal (1975-2030)*. European Commission, Joint Research Centre.

Seabold, S., & Perktold, J. (2010). statsmodels: econometric and statistical modeling with Python.
*Proceedings of the 9th Python in Science Conference*, 92-96.

Shamon, H., et al. (2024). SNAPSHOT USA 2021: a third coordinated national camera trap survey of the
United States. *Ecology*, 105, e4318.

Thornton, M.M., Shrestha, R., Wei, Y., Thornton, P.E., Kao, S.-C., & Wilson, B.E. (2022). *Daymet:
daily surface weather data on a 1-km grid for North America, version 4 R1*. ORNL DAAC.

Vazquez, C., Rowcliffe, J.M., Spoelstra, K., & Jansen, P.A. (2019). Comparing diel activity patterns
of wildlife across latitudes and seasons: time transformations using day length. *Methods in Ecology
and Evolution*, 10, 2057-2066.

Virtanen, P., et al. (2020). SciPy 1.0: fundamental algorithms for scientific computing in Python.
*Nature Methods*, 17, 261-272.

Wood, S.N. (2011). Fast stable restricted maximum likelihood and marginal likelihood estimation of
semiparametric generalized linear models. *Journal of the Royal Statistical Society B*, 73, 3-36.

Wood, S.N. (2017). *Generalized Additive Models: An Introduction with R*, 2nd edition. CRC Press.

### Citation verification note

Only the GHS-POP citation has been checked against a publisher source; it matches the form given in
the Joint Research Centre data catalogue. The Rooney et al. (2025) citation was supplied directly by
the senior author, and the Snapshot USA project-to-survey-year mapping was corrected by him against
the deployment dates. Every other reference here, including the NLCD, Daymet, Copernicus DEM and
tree-canopy data releases and all software citations, is the standard form from memory and has NOT
been verified against its provider's current recommended text. All of them require checking before
submission, and the data products especially, since their recommended citations change between
versioned releases.
