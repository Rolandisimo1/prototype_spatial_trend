# The model, written out in full, and every change we made to Arielle's original

A companion to the results reports. Part 1 states the model. Part 2 lists what
we changed, why, and what evidence settled each one.

## Take-home points

1. **The model equations are unchanged.** Our fork of the production
   `nimbleCode` block is byte-identical to Arielle's original, confirmed by
   direct comparison: 142 lines, one trailing blank line of difference and
   nothing else. Every change we made sits upstream of the model, in how the
   data reach it, or downstream, in which parameter we report.
2. **One change alters a published number: which trend we report.** The model
   contains two trend parameters. We report the one the cameras identify,
   $\beta^{\text{yr}}$. Earlier work reported their sum. This is a change of
   interpretation, and the model formulation was correct throughout.
3. **Four changes correct defects in the data pipeline**, each caught by a check
   against real data: a year-ordering fault in the effort matrix, a range mask
   that discarded documented populations, three soil covariates that are nearly
   collinear, and a chain-resumption fault that skipped re-adaptation.
4. **Three changes are deliberate design decisions**, each with a stated
   justification: a five-year reporting window, exclusion of cell-years touched
   by organised reporting campaigns, and a centred form of the iNaturalist link
   used for sampling.
5. **One piece of the current production code has not been inspected locally.**
   The centred variant runs on Hazel. We have its specification in the methods
   text and its arithmetic is an exact reparameterisation, but the file itself
   has not been read on this machine.

---

# Part 1. The model

## Notation

| Symbol | Meaning |
|---|---|
| $s = 1,\dots,S$ | Scale 1 cells (50 km) inside the species' range mask |
| $g = 1,\dots,G$ | Scale 2 cells (100 km), the resolution of the spatial fields |
| $t = 1,\dots,T$ | years |
| $i = 1,\dots,N$ | camera deployments |
| $j = 1,\dots,J_i$ | 10-day detection windows at deployment $i$ |
| $\lambda_{st}$ | latent intensity, the modelled relative abundance |
| $\psi_i$, $p_{ij}$ | occupancy and detection probability |
| $N_{st}$, $E_{st}$ | iNaturalist count of the species, and the effort offset |
| $\alpha_g$ | spatial intercept field |
| $\gamma^{\text{W}}_g$, $\gamma^{\text{C}}_g$ | spatially varying climate coefficients |
| $\boldsymbol{\beta}$ | nine covariate effects, shared by both submodels |
| $\theta_0$, $\theta_1$ | link between intensity and iNaturalist counts |
| $\beta^{\text{yr}}$ | camera-identified trend, **the reported quantity** |
| $\beta^{\text{var}}$ | iNaturalist-identified trend increment |
| $\beta^{\text{tot}}$ | $\beta^{\text{yr}} + \beta^{\text{var}}$, the intensity trend |
| $\phi$ | negative binomial overdispersion |

## 1.1 The intensity process

One latent surface of relative abundance underlies both data streams. In Scale 1
cell $s$ and year $t$,

$$\log \lambda_{st} = \alpha_{g(s)} + \mathbf{x}_{st}^{\top}\boldsymbol{\beta}
+ w^{\text{W}}_{st}\gamma^{\text{W}}_{g(s)} + w^{\text{C}}_{st}\gamma^{\text{C}}_{g(s)}
+ \beta^{\text{tot}} \, y_t$$

where $y_t$ is the standardised year. The three spatial fields are intrinsic
conditional autoregressive on the Scale 2 adjacency graph,

$$\alpha_{1:G} \sim \text{CAR}(0, \tau_{\alpha}^{-1}), \qquad
\gamma^{\text{W}}_{1:G} \sim \text{CAR}(0, \tau_{\text{W}}^{-1}), \qquad
\gamma^{\text{C}}_{1:G} \sim \text{CAR}(0, \tau_{\text{C}}^{-1})$$

with $\tau_{\alpha}, \tau_{\text{W}}, \tau_{\text{C}} \sim \text{Gamma}(0.1, 0.1)$.
In plain language, each 100 km cell has its own baseline level and its own
response to summer and winter temperature, each smoothed toward its neighbours.
These fields carry no separate global mean, so only contrasts within a field are
interpretable.

## 1.2 The camera submodel

Each deployment contributes a detection history over 10-day windows,

$$y_{ij} \sim \text{Bernoulli}(z_i \, p_{ij}), \qquad z_i \sim \text{Bernoulli}(\psi_i)$$

with $z_i$ the latent state indicating whether the species used the camera's
surveyed area during the deployment. The likelihood is evaluated in marginalised
form, summing over $z_i$. Occupancy links to the intensity field through the
complementary log-log,

$$\text{cloglog}(\psi_i) = \alpha_{g(i)} + \mathbf{x}_i^{\top}\boldsymbol{\beta}
+ w^{\text{W}}_i\gamma^{\text{W}}_{g(i)} + w^{\text{C}}_i\gamma^{\text{C}}_{g(i)}
+ \beta^{\text{yr}} \, \tilde{y}_i$$

Since $\text{cloglog}(\psi) = \eta$ gives $\psi = 1 - \exp(-e^{\eta})$, the
quantity $e^{\eta}$ is an intensity and the linear predictor sits on the
log-abundance scale. The camera submodel shares the intercept field, the nine
covariate effects and both climate fields with the intensity process. It differs
in one term only: its trend coefficient is $\beta^{\text{yr}}$, while the
intensity process carries $\beta^{\text{tot}}$.

Detection probability on occasion $j$ is

$$\text{cloglog}(p_{ij}) = \eta_p + \rho_1 d_{ij} + \rho_2 d_{ij}^{2}
+ \rho_3 h_i + \rho_4 \log r_i$$

with $d_{ij}$ the standardised day of year, $h_i$ canopy height, $r_i$ distance
to the nearest road, $\eta_p = \text{cloglog}(\pi_0)$, $\pi_0 \sim \text{Unif}(0,1)$
and $\rho_{1:4} \sim \mathcal{N}(0, 5)$.

## 1.3 The iNaturalist submodel

Counts are negative binomial with the effort offset in the mean,

$$N_{st} \sim \text{NegBinomial}(\mu_{st} E_{st},\ \phi), \qquad
\log \mu_{st} = \theta_0 + \theta_1 \log\!\Big(\sum_{\ell \in s} \lambda_{\ell t}\Big)$$

with $\phi \sim \text{Gamma}(\text{shape}=1,\text{scale}=5)$, and the sum running
over the covariate rows inside cell $s$. Here $\theta_0 \sim \mathcal{N}(0,10)$ is
the species' reporting rate, constant in space, and $\theta_1 \sim \mathcal{N}(1,1)$
sets how strongly record counts follow abundance. At $\theta_1 = 1$ records track
abundance proportionally; at $\theta_1 = 0$ the records say nothing about where
the animals are.

For fitting we centre the abundance term,

$$\log \mu_{st} = \theta_0^{c} + \theta_1 \Big(\log \sum_{\ell \in s} \lambda_{\ell t} - C\Big),
\qquad \theta_0 = \theta_0^{c} - \theta_1 C, \qquad \theta_0^{c} \sim \mathcal{N}(\theta_1 C, 10)$$

with $C$ a constant fixed per fit. The likelihood and the joint prior on
$(\theta_0, \theta_1)$ are unchanged, so this leaves the posterior the same and
changes only how efficiently the sampler explores it.

## 1.4 The temporal trend

The two submodels carry the trend differently. The cameras carry
$\beta^{\text{yr}}$; the intensity process, and through it the iNaturalist
counts, carries the sum

$$\beta^{\text{tot}} = \beta^{\text{yr}} + \beta^{\text{var}}$$

$$\beta^{\text{yr}} \sim \mathcal{N}(0, \sigma_{\text{yr}}), \qquad
\beta^{\text{var}} \sim \mathcal{N}(0, \sigma_{\text{var}}), \qquad
\sigma_{\text{yr}}, \sigma_{\text{var}} \sim \text{Unif}(0, 2)$$

$\beta^{\text{var}}$ appears only in the iNaturalist likelihood, so that stream
alone identifies it and it absorbs temporal signal specific to iNaturalist. That
leaves $\beta^{\text{yr}}$ identified by the cameras, and $\beta^{\text{yr}}$ is
what we report as the population trend.

## 1.5 Converting the coefficient to a yearly percentage

The model fits the slope against a standardised year, so the raw coefficient is
change per standard deviation of year. With $\sigma_y$ the standard deviation of
the year covariate in that fit,

$$\text{percent change per year} = 100 \left( e^{\,\beta^{\text{yr}} / \sigma_y} - 1 \right)$$

$\sigma_y$ differs between fits because species differ in how many years they
were sampled, so it is read from each fit's own record.

## 1.6 Derived quantities in the model that we do not report

The original code computes a signal-to-noise ratio
$\text{snr} = \beta^{\text{yr}} / \beta^{\text{var}}$ and a robustness indicator
$\text{step}(\text{snr} - 1)$. Both remain in the model and both are monitored.
We exclude them from conclusions: the ratio has no finite posterior mean when
$\beta^{\text{var}}$ can approach zero, and the convergence statistics for that
family describe the ratio, leaving the trend's own diagnostics elsewhere.

---

# Part 2. What we changed, and why

| # | Change | Where it sits | What prompted it |
|---|---|---|---|
| 1 | Order the iNaturalist effort matrix by calendar year | Data preparation | Columns sorted as text, displacing years |
| 2 | Range mask: union of IUCN polygon and record-bearing cells | Data preparation | Polygon-only mask discarded documented populations |
| 3 | Drop one of three soil texture fractions | Data preparation | Near-collinear for narrow-range species |
| 4 | Resolve species names against the Mammal Diversity Database | Data preparation | Source databases disagree on the accepted name |
| 5 | Exclude campaign-affected cell-years from the iNaturalist stream | Data preparation | Organised campaigns break the effort assumption |
| 6 | Single continuous chains, no chunked resumption | Fitting | Resumed chunks skipped re-adaptation |
| 7 | Centred form of the iNaturalist link | Fitting | $\theta_0$ and $\theta_1$ trade off along a ridge |
| 8 | Report $\beta^{\text{yr}}$ as the population trend | Reporting | $\beta^{\text{tot}}$ includes an iNaturalist-only term |
| 9 | Five-year reporting window | Reporting | Camera effort grew steeply across the record |
| 10 | Camera-only control fits | Validation | Tests whether iNaturalist drives a result |

## 2.1 The iNaturalist effort matrix was ordered as text, not by year

**What was happening.** The effort matrix is built by pivoting a long table into
one column per year. The pivot named the columns and sorted them as text, so a
species whose smallest in-range cell lacked records in some years received effort
values attached to the wrong years.

**Why it matters.** The effort offset tells the model how much observing happened
in a cell and year. Misaligning it against the counts corrupts the fitted
temporal trend, which is the quantity we report.

**Evidence.** Confirmed displaced years for moose (2008, 2010, 2011) and
white-tailed deer (2010). Bobcat ran on the same code and happened to escape,
because its cells had full coverage.

**What we did.** Forked the two affected functions with the sort corrected, and
added an assertion that the count and effort matrices carry identical column
order before either reaches the model. The fork is `integration_helper_fix1.R`.
Arielle's original file is untouched.

## 2.2 The range mask discarded documented populations

**What was happening.** The original masked each species to its IUCN range
polygon by setting counts outside the polygon to missing. A missing cell drops
out of the likelihood entirely, so records of animals genuinely present outside
the published polygon were removed from the fit, along with the cells holding
them.

**Why it matters.** Range expansions are exactly the signal a trend model should
detect, and the mask removed them. The effect was largest for species with small
or patchy ranges.

**What we tried, in order.** Dropping the mask entirely was tested on real moose
data: the grid went from 382 to 3,322 cells, recovering 3,145 moose records in
Colorado alone, but only 422 of the 3,322 cells held any records and the model
build failed at that size. A presence mask keeping only cells with records gave
422 cells. Lowering the record threshold from three to one moved the count from
about 290 to about 422, both far below the size that broke the build, so the
looser threshold cost nothing in tractability and kept genuine edge populations.

**What we settled on.** The fitting mask is the union of the IUCN polygon and the
cells holding research-grade records. A records-only mask discards cells where a
species is absent but was surveyed, which is real information. A polygon-only
mask discards documented populations.

**A second mask, for a different question.** Asking whether an agency publishes a
trend where a species occurs needs a stricter rule, because the union mask is
deliberately generous. The moose union mask places moose in Nebraska, New Mexico,
Nevada, Arizona and Iowa, where the cells come from vagrant or captive records.
The occurrence rule counts a jurisdiction as occupied if it falls inside the IUCN
polygon or holds more than a threshold number of records.

## 2.3 Three soil texture fractions are nearly collinear

**What was happening.** Clay, silt and sand fractions all sat in the covariate
set. Where they sum to a constant, the three columns are linearly dependent and
their individual effects are not separately identified.

**Evidence, and two corrections along the way.** The first diagnosis reasoned
from composition alone and named the wrong fraction to drop. A direct check on
real data corrected it. The dependence is not a universal identity: regressing
sand on clay and silt gives an $R^2$ of 0.9995 for moose, whose camera sites sit
in a narrow northern footprint, and much less for species spanning more soil
types. Recovering the raw fractions showed the sum is bimodal, with a spike at
exactly one and a tail below it, so the problem is concentrated in species whose
range overlaps few soil types.

**What we did.** Retain clay and silt, drop sand, and apply the same reduction to
every species. A fleet-wide check across eight species confirmed one uniform
design works for all of them, which keeps estimates comparable between species.

## 2.4 Species names disagree between sources

**What was happening.** The accepted binomial, the name an agency publishes
under, and the name the camera imagery platform files detections under do not
always agree. One species returned an empty grid because the pipeline keyed it on
a superseded binomial.

**What we did.** Resolve names against the Mammal Diversity Database and keep a
per-source crosswalk. Three entries are non-identity, and each records its
evidence and authority. Camera detection histories are keyed on the name the
platform uses, because that is the only name under which those detections exist.

## 2.5 Organised reporting campaigns break the effort assumption

**What was happening.** The effort offset assumes that total mammal records in a
cell and year measure how much observing happened there. A campaign that asks
people to report one species raises that species' records without raising the
others, so the offset understates the attention that species received.

**Evidence.** A difference-in-differences across campaign footprints, measured
against cells outside them, showed most campaigns shifted reporting toward their
target species. Campaigns did not recruit new observers and barely raised overall
reporting, which places the effect on composition instead of volume.

**What we did.** Exclude iNaturalist data from campaign-affected cell-years for
bobcat (47 Scale 1 cells), coyote (5), red fox (5), California ground squirrel
(3) and northern raccoon (3). Camera data in those cells are retained. We exclude
all years, including those before a campaign began, as a conservative choice, and
the methods text says so.

**What we decided against.** Using the measured campaign effect as a prior on a
model term. Ben Goldstein pointed out that the measured effect contains both the
reporting shift and any real local abundance change over the same period, so
using it as a prior would transfer local ecology into a reporting parameter. We
withdrew that proposal.

## 2.6 Resumed chains skipped re-adaptation

**What was happening.** Long fits ran as a sequence of chunks, saving and
restoring sampler state between them. The restore path omitted a call present in
the fresh-start path, so each resumed chunk began without re-adapting its
samplers and the first portion of its draws carried a restart transient.

**Why it matters.** It affected every species and every model at every chunk
boundary, contaminating roughly 7% of each round's draws.

**What we did.** A first fix was drafted and tested, and the test showed it had
the library's reset semantics wrong: it would have run without visible failure
while doing nothing. We abandoned the chunked framework and run single
uninterrupted chains. Several convergence diagnoses made under the old framework
failed to reproduce under the new runs and were set aside.

## 2.7 The iNaturalist link is centred for sampling

**What was happening.** In the uncentred form, $\theta_0$ and $\theta_1$ trade off
along a narrow ridge, and $\theta_1$ explores it slowly.

**What we did.** Subtract a constant $C$ from the abundance term inside the link
and adjust the intercept and its prior to match, as in section 1.3. The
likelihood and the joint prior are unchanged, so the posterior is identical. We
set $C$ to the slope of the ridge, estimated from the chain means of $\theta_0$
and $\theta_1$ in an earlier uncentred fit of the same species and window, or from
that species' full-record draws where the chain means were too weakly correlated
to define a ridge.

**Evidence that it works, and where it does not.** Two species refitted on the
same full-record window with centring reproduced their published estimates to
within 0.5 and 0.1 percentage points per year. For one species the sampler
improved by two orders of magnitude. For species whose chains sit in two separate
modes, centring addresses the wrong problem and those fits remain unusable.

## 2.8 Which parameter we report

**What was happening.** The model contains two trend parameters. Earlier work
reported $\beta^{\text{tot}}$, the sum, as the population trend.

**Why it matters.** $\beta^{\text{var}}$ appears only in the iNaturalist
likelihood, so that stream alone identifies it. Reporting the sum mixes a
camera-identified trend with an increment identified by opportunistic records of
unknown effort.

**What we did.** Report $\beta^{\text{yr}}$, the camera-identified trend. This is
a change in which quantity we read out of the fit. The model formulation and the
design that produced it were correct, and $\beta^{\text{var}}$ serves a purpose
in that form: it gives iNaturalist-specific temporal signal somewhere to go,
which protects $\beta^{\text{yr}}$ from absorbing it.

## 2.9 The five-year reporting window

**What was happening.** Fits used the full record from 2008.

**Why we changed it.** Camera effort grew steeply across the record, and few
cameras were running in the earliest years, so those years add more noise than
signal.

**What this costs.** The window materially moves estimates. Two species refitted
on the same window moved 0.1 and 0.5 percentage points per year, while four
species refitted on a shorter window moved 1.1, 3.5, 12.2 and 35.1 points. The
size of the shift tracks range size, with the largest shifts in the smallest
ranges. We report the five-year window as the headline, the ten-year window as
secondary, and keep the full record as a check.

## 2.10 Camera-only control fits

**What this adds.** Fitting the same bundle with the iNaturalist likelihood
removed gives a direct measurement of how much the opportunistic stream moves a
result. Where the joint and camera-only estimates agree, the joint trend stands.
Where they disagree, the iNaturalist stream is driving the result.

**What it showed.** Across species with both fits, the gap between them is
predicted by range size (Spearman $-0.61$, $p = 0.005$). Species below 100 grid
cells sat a median of 17.1 percentage points per year from their camera-only
estimate; species above sat 1.7 points away. That measurement is the basis for
flagging species below 100 cells as data-limited.

---

## What has not been inspected locally

The centred model code runs on Hazel. Section 1.3 states its specification and
the arithmetic showing it leaves the posterior unchanged, and the refit evidence
in section 2.7 is consistent with that. The file itself has not been read on this
machine, so the specification here rests on the methods text and on the behaviour
of the fits, not on a line-by-line reading of the production code.

## Glossary

**Byte-identical.** Two files whose contents match exactly, character for
character. Used here for the comparison between our model fork and Arielle's
original.

**CAR, conditional autoregressive.** A way of giving each spatial cell its own
value while pulling it toward the values of its neighbours, so the surface varies
smoothly.

**Centring, reparameterisation.** Rewriting a model in different parameters that
describe the same thing, so that the answer is unchanged and the sampler works
more easily.

**Collinear.** Two or more covariates carrying nearly the same information, so
the model cannot separate their individual effects.

**Complementary log-log, cloglog.** A link function connecting a probability to a
quantity on the log scale. Here it connects occupancy probability to log
abundance.

**Difference-in-differences.** Comparing the change inside a treated area against
the change outside it over the same period, so that anything affecting both
cancels.

**Effort offset.** A known quantity entered into the model to account for how
much observing happened, so that the estimate is a rate per unit effort rather
than a raw count.

**Marginalised likelihood.** Evaluating the likelihood with a latent variable
summed out, instead of sampling that variable.

**Negative binomial.** A count distribution allowing more variability than a
Poisson, controlled by the overdispersion parameter $\phi$.

**Scale 1 and Scale 2.** The 50 km grid on which intensity and counts are
modelled, and the 100 km grid on which the spatial fields vary.
