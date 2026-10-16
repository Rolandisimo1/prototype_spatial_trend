# Methods and Results outline

One spine, followed by both sections, so each Results subsection maps onto the
Methods subsection that produced it.

**A** Data assembled -> **B** Trends -> **C** Agreement between the two streams ->
**D** When integration is reliable -> **E** Habitat and climate -> **F** Comparison
with agency trends

---

## Methods

### 1. Data
1.1 Study species and taxonomy: the fleet, the name crosswalk, and the three
    species where accepted, agency and platform names disagree.
1.2 Camera-trap data: sequences, deployments, detection histories, occasion
    length, years covered, per-species detection totals.
1.3 iNaturalist data: research-grade records aggregated to cell-year counts, the
    effort offset, and the share-trend check against circularity.
1.4 Range masks: the union mask used for fitting, and the stricter corroborated
    mask used to define occurrence.
1.5 Covariates: the nine retained, the soil-fraction reduction and the
    collinearity that forced it, standardisation, and the absence of a year
    dimension.

### 2. Model
2.1 Spatial structure: intensity at cell50, the three smoothed fields at
    cell100, and the adjacency graph.
2.2 Integrated likelihood: the camera submodel, the opportunistic count
    submodel, and the link between them.
2.3 Trend structure and the reported quantity: the decomposition into a
    camera-anchored component and an opportunistic-only component, why we report
    the former, and the per-submodel year standardisation that makes raw
    coefficients non-comparable.
2.4 Parameterisations: national slope, and national slope with partially pooled
    ecoregion deviations.

### 3. Fitting and assessment
3.1 Implementation and sampling: software and versions, chains, iterations,
    burn-in, single continuous runs with no resume boundary.
3.2 Convergence criteria, and which parameter families are excluded from them.
3.3 Camera-only comparison fits as the anchoring diagnostic.

### 4. External comparison
4.1 Agency trend compilation, including the automated-agent disclosure and its
    verification (see METHODS_agency_trend_search.md).
4.2 Occurrence denominator and how monitoring coverage is measured.

---

## Results

### 1. Data assembled  <- Methods 1
Deployments, survey nights, records, per-species detection rates, and the effort
imbalance between the two streams in both space and time.

### 2. Population trends  <- Methods 2, 3
2.1 Which fits are reportable and which do not resolve a direction.
2.2 National trends per species, as percent change per calendar year.
2.3 Regional trends for the species fitted under the ecoregion parameterisation.

### 3. Agreement between data streams  <- Methods 2.3
The sign disagreements between camera-anchored and opportunistic trends, and the
fact that reporting the opportunistic quantity would reverse them.

### 4. When integration is reliable  <- Methods 3.3
Range size predicts loss of the camera anchor; camera effort does not.

### 5. Habitat and climate associations  <- Methods 1.5, 2.1
Covariate effects and the climate response surfaces, and their consistency
across parameterisations.

### 6. Comparison with agency trends  <- Methods 4
6.1 Concordance for the intensively managed species.
6.2 The monitoring gap, in three tiers.

---

## Ordering decisions

- **Stream agreement (3) precedes reliability (4).** The reader needs to see that
  the streams disagree before being told when the integration can be trusted;
  reversed, section 4 has no motivation.
- **The agency comparison goes last**, not beside the trends. It is the external
  check and it carries the monitoring-gap result, both of which land harder once
  the trends and their caveats are established. It also closes the paper on what
  the method provides that nothing else does.
- **Habitat and climate sit at 5, as secondary.** This assumes the paper is about
  trends. If it is to read as a distribution-and-trend paper, that section moves
  to 2 and the spine shifts.
