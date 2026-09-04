# Environmental covariates: methods

## Extraction

We extracted all environmental covariates in Google Earth Engine (Gorelick et al. 2017) at each camera
deployment's coordinates, verifying every asset identifier and band name against the live catalogue before
extraction rather than relying on remembered identifiers. The delivered table holds 26,678 deployments and
60 covariate columns.

We extracted each covariate at three spatial supports: at the camera location itself, and as the mean
within 1 km and 5 km radius buffers. Reductions used `reduceRegions` in 40-feature batches; the 5 km
buffers were reduced at 100 m rather than 30 m resolution to remain inside Earth Engine compute limits,
which affects the precision of small-area percentages but not their expectation. Time-invariant layers
were extracted once per unique coordinate and joined back, since the 26,678 deployments occupy 21,595
distinct coordinates.

We matched each layer to the deployment period rather than using a single epoch. The National Land Cover
Database is published in nine epochs, and we assigned each coordinate to its nearest epoch by median
deployment year; deployments span 2008 to 2025, so all 2022 to 2025 deployments (43% of the sample) draw
on the 2021 release and their land cover is up to four years stale. Nighttime lights and population were
matched to the nearest available annual product, and tree canopy is capped at 2021, the last release
available. The epoch used is recorded per row so that it can be tested as a nuisance term.

Daily maximum temperature was computed server-side against each deployment's own active date window,
giving both the mean daily maximum over the deployment period and the mean daily maximum of the hottest
calendar month within it. Climate normals were extracted separately and kept in separate columns, since
normals cannot substitute for date-matched values in a thermal analysis.

## Scope and coverage

We restricted covariates to the conterminous United States, assigning region from administrative
boundaries rather than a latitude-longitude box. This mattered: a bounding box alone misclassified 62
Mexican deployments as domestic. Of 26,798 deployments, 120 fell outside the lower 48 (51 in Jalisco, 43
in Alaska, 11 in Baja California, 10 in Hawaii, 5 in Puerto Rico) and were dropped.

Land cover coverage is complete for all 26,678 retained deployments. Two residual gaps were filled rather
than left missing: 127 rows at 88 coastal coordinates fell on ocean pixels in the 1 km climate normals
grid and were filled with a 5 km buffer mean, and 12 rows on barrier islands fell outside the daily
climate land grid and were likewise filled and flagged. Nineteen rows lack a date-matched temperature
value because their deployment dates are missing or invalid in the source file.

## Buffer overlap and the effective unit of a covariate

Cameras in this dataset are deployed in tight clusters, and this determines what a buffer covariate can
mean. The median nearest-neighbour distance between deployments is 51 m; 65.5% have a neighbour within
200 m and 93.7% within 1 km. A buffer of radius R overlaps a neighbour's buffer whenever another
deployment lies within 2R, so the typical 1 km buffer shares its footprint with 21 other cameras' buffers
and the typical 5 km buffer with 63.

Buffer covariates are therefore array-level rather than camera-level quantities, and the effective sample
size for a covariate effect is the number of sites, not the number of cameras. Treating cameras as
independent would inflate the degrees of freedom roughly thirtyfold. We consequently aggregated covariates
to the site before fitting, and we do not interpret camera-level covariate residual variation. Point-scale
covariates are the only ones carrying genuine within-array contrast, since a camera in forest and one in a
clearing 60 m away differ at the point scale but share both buffers.

## Keeping agriculture and urbanization separate

We deliberately held agriculture and urbanization as separate axes throughout and constructed no composite
human-footprint index at any stage. The two mechanisms act differently on animals, and the data support
treating them separately: across all 9 agriculture by 12 urbanization column pairs, the maximum absolute
Pearson correlation is 0.17 (mean 0.07). The strongest association is between pasture and population
(Spearman 0.44), which is interpretable as exurban landscapes genuinely mixing pasture with dispersed
housing, and is far below the level at which the two would be inseparable. A single index would have
averaged two nearly independent predictors into one term measuring neither.

Two qualifications apply. Both axes are strongly zero-inflated, which depresses Pearson correlations, so
we report Spearman alongside. And weak correlation is not balanced coverage: cells combining high cropland
with high impervious surface are sparse, so interactions between the two axes are weakly identified even
though their main effects are separable.

## Collinearity within the urbanization axis

Nighttime lights and population density are collinear at r = 0.895 on the log scale (Spearman 0.80),
which is close to the impervious-surface-versus-lights pair (r = 0.92) that we had already excluded for
that reason. Fitted jointly at site level, lights reaches a variance inflation factor of 11.2
and population 8.3, and developed-land share is a third member of the same cluster at 7.6; retaining any
one of the three alone brings it below 1.4.

This has a direct consequence for interpretation. A linear diagnostic at site level showed the lights
coefficient for coyote reversing sign between the joint and single-predictor fits, and the deer
coefficient inflating roughly twenty-sevenfold when fitted jointly. Opposite signs on a strongly collinear
pair are the signature of collinearity rather than of two distinct mechanisms. We therefore fitted three
variants per species, each carrying exactly one of the three human-pressure measures, and report all
three; we do not fit two of them together, and we do not describe a mechanism for a human-pressure effect whose sign does not survive being fitted alone. For prediction the
collinearity is far less damaging than for inference, since a collinear pair can predict jointly while
being unable to attribute the effect between its members; surfaces built from the joint variant are
labelled predictive rather than mechanistic.

This collinearity is confined to the urbanization axis. Agriculture terms carry variance inflation factors
of 1.1, and cropland and pasture are themselves near-independent, so the separation of agriculture from
urbanization is unaffected.

## Covariates entering the models

After removing within-layer redundancy (point and 1 km elevation correlate at 0.999; monthly and bioclim
temperature extremes at 0.93 to 0.99) we assembled a candidate set of ten site-level covariates: human
population, nighttime lights, developed-land share, cropland, agriculture as crop plus pasture, tree
canopy cover, forest share, elevation, ruggedness, and the date-matched maximum temperature of the hottest
month. Site values are means over the cameras at each site, with a median of 13 cameras contributing per
site.

That full candidate set is not estimable as a single model. Computed over the 484 sites with complete
covariates, three terms exceed a variance inflation factor of 5: nighttime lights at 11.2, human
population at 8.3, and developed-land share at 7.6. All three measure the same gradient, and their
pairwise correlations on the log scale run from 0.90 to 0.92. Two further terms are internally redundant:
the agriculture sum correlates with its own cropland component at 0.79, since it is crop plus pasture.

We therefore fitted a seven-covariate model carrying exactly one human-pressure term and one agriculture
term, in three variants distinguished by which human-pressure measure enters. Every variant clears the
conventional threshold comfortably:

| Human-pressure term retained | Maximum VIF | Terms above 5 |
|---|---|---|
| Human population | 4.24 | none |
| Nighttime lights | 4.27 | none |
| Developed-land share | 4.33 | none |

In all three the largest inflation is on the vegetation pair (forest share 4.2 to 4.3, tree canopy 4.0),
which is expected since both measure woody cover, and the retained human-pressure term itself carries a
factor near 1.4. We report all three variants rather than selecting one, because the choice of
human-pressure measure is not resolvable from these data and the three are not interchangeable in
interpretation even where they are near-interchangeable in fit.

Point-scale covariates are not redundant with their buffers (impervious surface at the point versus 1 km,
r = 0.59; forest, r = 0.66), which supports retaining one point-scale term alongside one buffer term for
within-array questions.

## Prediction grid, and a caveat that bounds the maps

For mapping we extracted the same layers onto a 25 km grid over the conterminous United States (14,069
cells). An earlier version of this grid was missing an entire column of cells at the projection's central
meridian, an indexing fault rather than absent source data; the corrected grid adds 97 cells.

Site-aggregated and grid-extracted covariates are different quantities, and they disagree more than we
expected. Elevation agrees well between the two (r = 0.96) and tree cover moderately (r = 0.75), but
nighttime lights agree at only r = 0.64 and human population at r = 0.14. A site covariate is a mean over
camera locations, which are not placed at random within a 25 km cell, whereas a grid covariate is a cell
average; for a strongly skewed layer such as population these are simply not the same measurement.
Neither can be declared correct without an external reference.

We therefore refitted every prediction surface on grid-extracted covariates before mapping, rather than
applying site-fitted coefficients to grid values. This reduced the mappable set from six species-by-measure
combinations to four (coyote night share, midday share and concentration; white-tailed deer crepuscular
share), and it is the honest constraint on what these maps can claim. A developed-land layer exists in the
site covariates but not on the grid, so urbanization results cannot be projected onto the map at all.

## Data sources and citations

| Covariate | Product | Version / release | Native resolution | Citation |
|---|---|---|---|---|
| Land cover, impervious surface | National Land Cover Database | 2021 and 2019 releases | 30 m | Dewitz (2023); Dewitz & USGS (2021) |
| Tree canopy cover | NLCD Tree Canopy Cover | 2021, v2021-4 | 30 m | Coulston et al. (2012) |
| Human population | Global Human Settlement Layer, GHS-POP | R2023A (EE asset P2023A) | 100 m | Schiavina et al. (2023) |
| Nighttime lights | VIIRS Day/Night Band annual composites | V2.1 (2013-2021), V2.2 (2022-2025) | 500 m | Elvidge et al. (2021) |
| Daily maximum temperature | Daymet | V4 | 1 km | Thornton et al. (2022) |
| Temperature normals | WorldClim monthly and bioclimatic | V1 | 1 km | Hijmans et al. (2005) |
| Elevation (primary) | Copernicus DEM GLO-30 | 2021 release | 30 m | European Space Agency (2021) |
| Elevation (cross-check) | USGS 3D Elevation Program | 10 m | 10 m | U.S. Geological Survey (2019) |
| Administrative boundaries | Global Administrative Unit Layers | 2015, simplified 500 m | vector | FAO (2015) |
| Extraction platform | Google Earth Engine | - | - | Gorelick et al. (2017) |

### Full references

Coulston, J.W., Moisen, G.G., Wilson, B.T., Finco, M.V., Cohen, W.B., & Brewer, C.K. (2012). Modeling
percent tree canopy cover: a pilot study. *Photogrammetric Engineering and Remote Sensing*, 78(7),
715-727.

Dewitz, J. (2023). *National Land Cover Database (NLCD) 2021 Products*. U.S. Geological Survey data
release. https://doi.org/10.5066/P9JZ7AO3

Dewitz, J., & U.S. Geological Survey (2021). *National Land Cover Database (NLCD) 2019 Products* (ver.
2.0, June 2021). U.S. Geological Survey data release. https://doi.org/10.5066/P9KZCM54

Elvidge, C.D., Zhizhin, M., Ghosh, T., Hsu, F.-C., & Taneja, J. (2021). Annual time series of global VIIRS
nighttime lights derived from monthly averages: 2012 to 2019. *Remote Sensing*, 13(5), 922.

European Space Agency (2021). *Copernicus DEM: Global and European Digital Elevation Model, GLO-30*.
https://doi.org/10.5270/ESA-c5d3d65

Food and Agriculture Organization of the United Nations (2015). *Global Administrative Unit Layers
(GAUL)*, 2015 edition.

Gorelick, N., Hancher, M., Dixon, M., Ilyushchenko, S., Thau, D., & Moore, R. (2017). Google Earth Engine:
planetary-scale geospatial analysis for everyone. *Remote Sensing of Environment*, 202, 18-27.

Hijmans, R.J., Cameron, S.E., Parra, J.L., Jones, P.G., & Jarvis, A. (2005). Very high resolution
interpolated climate surfaces for global land areas. *International Journal of Climatology*, 25(15),
1965-1978.

Schiavina, M., Freire, S., Carioli, A., & MacManus, K. (2023). *GHS-POP R2023A: GHS population grid
multitemporal (1975-2030)*. European Commission, Joint Research Centre [Dataset].
https://doi.org/10.2905/2FF68A52-5B5B-4A22-8F40-C41DA8332CFE

Thornton, M.M., Shrestha, R., Wei, Y., Thornton, P.E., Kao, S.-C., & Wilson, B.E. (2022). *Daymet: daily
surface weather data on a 1-km grid for North America, version 4 R1*. ORNL DAAC, Oak Ridge, Tennessee.
https://doi.org/10.3334/ORNLDAAC/2129

U.S. Geological Survey (2019). *3D Elevation Program 10-Meter Resolution Digital Elevation Model*.

### Citation verification status

The Global Human Settlement population citation above matches the publisher's official form, verified
against the Joint Research Centre data catalogue and the Earth Engine dataset page. The remaining product
citations are the standard forms for these datasets and should be checked against each provider's current
recommended citation before submission, since several are versioned data releases whose citation text is
revised between releases. The Earth Engine asset identifiers and band names in the table were each
verified live against the catalogue at extraction time and are authoritative.
