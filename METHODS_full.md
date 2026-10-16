# Methods

## Spatial scales

We reconciled three data resolutions by defining two spatial scales, referred to
throughout. **Scale 1** is a 50 km grid on which the intensity process operates
and over which iNaturalist observations are aggregated, indexed by $s$.
**Scale 2** is a 100 km grid used to define the conditional autoregressive (CAR)
spatial fields, indexed by $g$. Each Scale 2 cell contains up to four Scale 1
cells, and the map from a Scale 1 cell to its parent is written $g(s)$; for a
camera deployment $i$ the parent cell is $g(i)$.

Both grids are rasters in a common equal-area projection covering North America.
Our Scale 1 therefore merges the separate intensity and aggregation grids used by
Goldstein et al. (in review), which placed the intensity process at 5 km and
aggregated iNaturalist counts at 50 km.

## Notation

| Symbol | Meaning |
|---|---|
| $s = 1,\dots,S$ | Scale 1 cells (50 km) inside the species' range mask |
| $g = 1,\dots,G$ | Scale 2 cells (100 km), the resolution of the spatial fields |
| $t = 1,\dots,T$ | years, 2008 to 2025 ($T = 18$) |
| $i = 1,\dots,N$ | camera deployments |
| $j = 1,\dots,J_i$ | 10-day detection windows at deployment $i$ |
| $r = 1,\dots,R$ | EPA Level I ecoregions ($R = 8$), regional parameterisation only |
| $\lambda_{st}$ | latent intensity |
| $\psi_i$, $p_{ij}$ | occupancy and detection probability |
| $N_{st}$, $E_{st}$ | iNaturalist count and effort offset |
| $\alpha_g$ | spatial intercept field |
| $\gamma^{\text{W}}_g$, $\gamma^{\text{C}}_g$ | spatially varying climate coefficients |
| $\boldsymbol{\beta}$ | nine covariate effects, shared by both submodels |
| $\theta_0$, $\theta_1$ | correspondence between intensity and iNaturalist counts |
| $\beta^{\text{yr}}$ | camera-identified trend, the reported quantity |
| $\beta^{\text{var}}$ | iNaturalist-identified trend increment |
| $\beta^{\text{tot}}$ | $\beta^{\text{yr}} + \beta^{\text{var}}$, the intensity trend |
| $\delta^{\text{reg}}_r$ | regional trend deviation |
| $\phi$ | negative binomial overdispersion |

## Species and spatial extent

[Fleet definition and selection criteria.] We resolved species names against the
Mammal Diversity Database and maintain a per-source crosswalk, because the
accepted binomial, the name an agency publishes under, and the name the camera
imagery platform uses do not always agree. Three crosswalk entries are
non-identity, and each records its evidence and authority.

We restricted the extent per species using a range mask formed as the union of
the species' IUCN extant-and-native range polygon and the Scale 2 cells
containing research-grade iNaturalist records. We modelled only those camera
deployments and Scale 1 cells falling inside that mask. The union is used rather
than either source alone because a records-only mask discards cells where a
species is absent but was surveyed, while a polygon-only mask excludes
populations documented outside the published range.

For the separate question of whether a state agency publishes a trend where a
species occurs (see *Comparison with agency trends*), the union mask is
inappropriate because it is deliberately inclusive and admits cells arising from
vagrant or captive records. That comparison uses a stricter occurrence rule,
counting a jurisdiction as occupied if it falls inside the IUCN polygon or holds
more than a threshold number of records. [This rule is specified but not yet
applied; coverage figures are currently computed on the union mask.]

## Camera trap data

[Project, deployment and survey-night counts; inclusion rules.] Camera detection
histories were prepared from Wildlife Insights sequence and deployment records by
a single upstream pipeline, which returns one detection-history matrix per
species with site and observation covariates. Detections were discretised into
consecutive 10-day detection windows, following the convention used for the same
camera dataset by Goldstein et al. (in review). Every camera quantity in this
study derives from that pipeline, and we audited it by comparing its count of
deployments with at least one detection against the raw sequence records
(Results).

## iNaturalist data

We used research-grade iNaturalist records of mammals in North America. In each
Scale 1 cell $s$ and year $t$ we recorded the number of records of the target
species as the count $N_{st}$, and the total number of mammal records of any
species as the cell-year effort $E_{st}$. The model treats $E_{st}$ as a known
offset, so the estimated intensity is a rate per unit observer effort.

## Covariates

We used nine covariates in the shared linear predictor, two climate covariates
with spatially varying coefficients, and three covariates for detection
probability only. Each is listed in Table 1 with its scale, temporal resolution,
transformation and source. All covariates were aggregated to the Scale 1 grid and
standardised before fitting.

Covariate selection proceeded in two documented reductions from an initial set of
seventeen. The first removed seven variables (climatic water deficit, NDVI
standard deviation, impervious surface, Palmer drought severity index,
precipitation, mean temperature and elevation). The second removed one of the
three soil texture fractions, which are compositional and sum to a constant; we
retained clay and silt and dropped sand, and we apply this single reduction
uniformly across species so that the covariate design is identical for every
species and estimates remain comparable.

## Integrated species distribution model

We fitted single-species integrated models in which camera trap detections and
iNaturalist counts inform one shared latent intensity field, extending the
spatial formulation of Goldstein et al. (in review) with a temporal trend.

### Intensity process

The log intensity in Scale 1 cell $s$, year $t$, is

$$\log \lambda_{st} = \alpha_{g(s)} + \mathbf{x}_{st}^{\top}\boldsymbol{\beta} + w^{\text{W}}_{st}\gamma^{\text{W}}_{g(s)} + w^{\text{C}}_{st}\gamma^{\text{C}}_{g(s)} + \beta^{\text{tot}} \, y_t$$

where $\alpha_g$ is a spatially structured intercept, $\mathbf{x}_{st}$ is the
vector of nine covariates with spatially constant effects $\boldsymbol{\beta}$,
$w^{\text{W}}_{st}$ and $w^{\text{C}}_{st}$ are the two climate covariates with
spatially varying coefficients $\gamma^{\text{W}}_g$ and $\gamma^{\text{C}}_g$,
and $y_t$ is the standardised year with coefficient $\beta^{\text{tot}}$.

All three spatial fields are intrinsic conditional autoregressive on the Scale 2
adjacency graph,

$$\alpha_{1:G} \sim \text{CAR}(0, \tau_{\alpha}^{-1}), \qquad \gamma^{\text{W}}_{1:G} \sim \text{CAR}(0, \tau_{\text{W}}^{-1}), \qquad \gamma^{\text{C}}_{1:G} \sim \text{CAR}(0, \tau_{\text{C}}^{-1})$$

with $\tau_{\alpha}, \tau_{\text{W}}, \tau_{\text{C}} \sim \text{Gamma}(0.1, 0.1)$.
In plain language, each Scale 2 cell has its own intercept and its own response to
summer and winter temperature, smoothed toward the values of its neighbours.
Because the climate fields carry no separate global mean, each is a relative
surface and only contrasts within a field are interpretable.

### Camera submodel

Each camera deployment contributes a detection history. The $j$th observation at
deployment $i$ is

$$y_{ij} \sim \text{Bernoulli}(z_i \, p_{ij}), \qquad z_i \sim \text{Bernoulli}(\psi_i)$$

where $z_i$ is the latent state indicating whether the area surveyed by the
camera was used by the species during the deployment, and $\psi_i$ is the
probability of that use. We adopt the asymptotic interpretation of occupancy,
treating a deployment as occupied if its surveyed area was used at least once
during the deployment. The likelihood is evaluated in its marginalised form,
summing over $z_i$.

Occupancy is linked to the intensity field by the complementary log-log relation

$$\text{cloglog}(\psi_i) = \alpha_{g(i)} + \mathbf{x}_i^{\top}\boldsymbol{\beta} + w^{\text{W}}_i\gamma^{\text{W}}_{g(i)} + w^{\text{C}}_i\gamma^{\text{C}}_{g(i)} + \beta^{\text{yr}} \, \tilde{y}_i$$

Since $\text{cloglog}(\psi) = \eta$ implies $\psi = 1 - \exp(-e^{\eta})$, the
quantity $e^{\eta}$ is an intensity and the linear predictor is on the
log-abundance scale. The camera submodel therefore shares the intercept field,
the nine covariate effects and both climate fields with the intensity process,
and differs only in its trend coefficient, $\beta^{\text{yr}}$.

Detection probability on occasion $j$ is

$$\text{cloglog}(p_{ij}) = \eta_p + \rho_1 d_{ij} + \rho_2 d_{ij}^{2} + \rho_3 h_i + \rho_4 \log r_i$$

with $d_{ij}$ the standardised day of year, $h_i$ canopy height, $r_i$ distance
to the nearest road, $\eta_p = \text{cloglog}(\pi_0)$ and $\pi_0 \sim
\text{Unif}(0,1)$, and $\rho_{1:4} \sim \mathcal{N}(0, 5)$.

### iNaturalist submodel

Counts are negative binomial with the effort offset entering the mean,

$$N_{st} \sim \text{NegBinomial}(\mu_{st} E_{st},\ \phi), \qquad \log \mu_{st} = \theta_0 + \theta_1 \log\!\Big(\sum_{\ell \in s} \lambda_{\ell t}\Big)$$

where $\phi \sim \text{Gamma}(\text{shape}=1,\text{scale}=5)$ is the
overdispersion parameter and the sum runs over the covariate rows $\ell$ falling
within cell $s$. The parameters $\theta_0$ and $\theta_1$ control the
correspondence between the intensity field and the iNaturalist data:
$\theta_0 \sim \mathcal{N}(0,10)$ is an intercept for the species' reporting
rate, constant in space, and $\theta_1 \sim \mathcal{N}(1,1)$ governs the
strength of that correspondence, with $\theta_1 = 1$ denoting proportionality.

This formulation treats the camera data as the gold-standard dataset directly
informing the latent intensity, and allows the iNaturalist data to influence it
more weakly where the two diverge, reflecting that the camera data are collected
under a standardised protocol with known effort.

### Temporal trend

The trend enters the two submodels differently. The camera submodel carries
$\beta^{\text{yr}}$; the intensity process, and hence the iNaturalist counts,
carries the sum

$$\beta^{\text{tot}} = \beta^{\text{yr}} + \beta^{\text{var}}$$

$$\beta^{\text{yr}} \sim \mathcal{N}(0, \sigma_{\text{yr}}), \qquad \beta^{\text{var}} \sim \mathcal{N}(0, \sigma_{\text{var}}), \qquad \sigma_{\text{yr}}, \sigma_{\text{var}} \sim \text{Unif}(0, 2)$$

Because $\beta^{\text{var}}$ appears only in the iNaturalist likelihood, it is
identified by that stream alone and absorbs temporal signal specific to it,
leaving $\beta^{\text{yr}}$ identified by the cameras. **We report
$\beta^{\text{yr}}$ as the population trend**, on the grounds that the cameras
are the designed survey with known effort.

The two submodels standardise their year covariate on different distributions:
$y_t$ on a fixed sequence of the 18 study years, and $\tilde{y}_i$ on each
species' own camera site-years. One unit therefore represents a different number
of calendar years in each submodel, and the two coefficients are not directly
comparable. We convert each to percent change per calendar year as
$100\,(e^{\beta/\Delta} - 1)$, where $\Delta$ is the number of years per unit for
that submodel and species, before any comparison or reporting.

### Regional parameterisation

We also fitted a parameterisation in which the iNaturalist trend carries a
partially pooled regional deviation. The intensity trend coefficient becomes
$\beta^{\text{tot}} + \delta_{g}$, with

$$\delta_g = \delta^{\text{reg}}_{r(g)}, \qquad \delta^{\text{reg}}_r \sim \mathcal{N}(0, \sigma_{\text{reg}}), \qquad \sigma_{\text{reg}} \sim \text{Exp}(1)$$

for $r$ indexing the eight EPA Level I ecoregions and $r(g)$ the region
containing Scale 2 cell $g$. The camera submodel is unchanged and retains the
global $\beta^{\text{yr}}$, so regional deviations are identified by the
iNaturalist stream alone and no regional trend in these fits is
camera-corroborated.

### Priors

$\boldsymbol{\beta} \sim \mathcal{N}(0, 5)$ independently, with a Laplace
alternative, $\beta_m \sim \text{Laplace}(0, \lambda_{k(m)})$, available for
grouped shrinkage. All remaining priors are stated with their parameters above.

## Fitting

[Software versions, chain count, iterations, burn-in, thinning.] Models were
fitted by Markov chain Monte Carlo in NIMBLE. Each chain was run as a single
continuous sample with no resume boundary, because a chunked-checkpoint scheme
used earlier omitted re-adaptation burn-in on resumption and so contaminated
draws at every chunk boundary.

We assess convergence on the trend coefficients and the covariate effects, and
exclude derived quantities from convergence summaries, since a ratio with a
near-zero denominator yields diagnostics that reflect the ratio rather than
sampler behaviour. Where the posterior of $\theta_1$ was multimodal we refitted
from dispersed initial values and selected between modes by comparing the
pointwise log predictive density of chains sampled in each mode, judging the
difference against the spread among chains within a mode.

## Camera-only comparison fits

For every species we refitted the same bundle with the iNaturalist likelihood
removed and nothing else altered. The displacement between the joint and
camera-only posterior means of $\beta^{\text{yr}}$ measures how far adding the
iNaturalist stream moves the camera-identified trend, and we use it to decide
whether a joint estimate remains camera-anchored.

## Comparison with agency trends

Described in the accompanying methods text, including the search and extraction
rules, the separation of population estimates from harvest indices, the choice of
comparison window, and the disclosure that the compilation was assembled by
automated agents and requires author verification.
