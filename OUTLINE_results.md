# Results outline

Follows the Methods spine. Items marked **[from Methods]** were removed from the
Methods draft under the no-results-in-methods rule; items marked **[drop]** were
removed and judged not worth a line here either.

## 1. Data assembled

- Camera deployments, survey nights, projects, years; iNaturalist records and the
  growth of the record pool over the study period.
- Per-species detection rate and detecting-deployment counts. These set up the
  species for which a camera-identified trend can be resolved at all.
- **[from Methods]** The upstream camera-prep audit: the pipeline retains a
  median of 0.901 of raw detecting deployments across 25 species (IQR
  0.866-0.919), with one species at 0.511, attributable to a genus rename rather
  than data loss. One sentence, as a data-integrity check.
- **[from Methods]** Share of the total mammal record pool for the largest
  contributors, descriptively and without raising circularity: white-tailed deer
  9.7% and eastern gray squirrel 8.9% of all mammal records. One sentence.
- Camera effort is concentrated in recent years and the detection record for most
  species begins in 2019, against an iNaturalist record from 2008. This is the
  reason the camera-identified and iNaturalist-identified trends cover different
  periods, and it should be stated here with the figures rather than asserted.

## 2. Population trends

- 2.1 Which fits are reportable, which do not resolve a direction, and which are
  excluded. **[from Methods]** Includes the species whose $\theta_1$ posterior
  was multimodal: the near-proportional mode was favoured by 40 to 75 times the
  within-mode chain spread for two species, while for a third the better-fitting
  mode implies a negative $\theta_1$, which has no interpretation, so that
  species is not reported.
- 2.2 National trends per species in percent change per calendar year, both
  streams.
- 2.3 Regional trends for the species fitted under the regional
  parameterisation.

## 3. Agreement between the two data streams

- Sign disagreements between the camera-identified and iNaturalist-identified
  trends, and the consequence that reporting the iNaturalist quantity would
  reverse them.
- The distribution of $\theta_1$ across the fleet, and that it is below one for
  most species, indicating sub-linear correspondence.

## 4. When integration is reliable

- Range size predicts displacement between the joint and camera-only fits;
  camera effort does not. The three species below the threshold and the mechanism
  in the Scale 2 field resolution.

## 5. Habitat and climate associations

- Covariate effects per species, and their consistency between
  parameterisations.
- The two spatially varying climate surfaces and the proportion of cells at which
  each resolves.

## 6. Comparison with agency trends

- 6.1 Concordance for the species agencies monitor intensively.
- 6.2 The monitoring gap in three tiers, including the species for which no
  agency publishes a trend anywhere in the range.
- **[from Methods]** The states the union range mask admits that hold no
  population, which is the evidence for using a stricter occurrence rule in this
  comparison. Belongs here, with the rule stated in Methods.

## Dropped entirely

- **[drop]** That the chunked-checkpoint defect left point estimates largely
  unaffected while inflating credible intervals. Inconsequential: the affected
  fits were discarded and none are reported, so the comparison has no bearing on
  any result presented.
- **[drop]** That the soil texture collinearity stalled convergence for the
  narrowest-range species. The reduction is a design decision; the convergence
  behaviour of a superseded fit is not a finding worth a line.
- **[drop]** The robustness indicator and the numbers behind not using it. We
  simply do not report it, so it needs no defence unless a reviewer raises it.
