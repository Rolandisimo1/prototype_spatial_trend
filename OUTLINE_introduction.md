# Introduction outline — continental mammal trends from cameras and iNaturalist

*Revised from the working outline. Additions and changes are flagged; the
original spine (P1 importance, P2 birds, P3 mammal gap, P4 our study) is kept
because it funnels correctly.*

---

## P1 — Why population trends matter, and why they are hard at scale

- Knowing which populations are growing and which are declining is a first-order
  conservation priority; it is what lets effort move from documentation to
  diagnosis and action.
- Example: a local, small-scale monitoring programme that detected a decline and
  triggered a management response. **Use a mammal example if one is available**,
  so the taxon is established before P2 goes to birds and has to come back.
- **ADD — what managers actually have for North American mammals.** State agency
  trend products are overwhelmingly harvest-based: catch-per-unit-effort from
  furbearer harvest surveys, bowhunter observation surveys, spotlight counts.
  Agencies say plainly in their own reports that harvest indices confound
  population change with effort change (fewer trappers, shifting pelt prices),
  and that coverage is biased toward game and furbearer species. This is the
  practical gap the paper fills, and naming it here makes P4's contribution
  concrete rather than methodological.
- The scaling problem: doing this annually, across a continent, for many species
  at once is the hard version of the problem.

## P2 — Birds show it can be done

- Citizen science has produced unprecedented biodiversity data, most of all for
  birds.
- Statistical advances have turned those data into continental trend estimates
  (eBird Status and Trends).
- Example of on-the-ground conservation action stemming from eBird results.
- **ADD — the hinge sentence, one line, that P3 depends on.** eBird works for
  trends because checklists are *semi-structured*: observers report duration,
  distance and protocol, and declare whether the list is complete. That makes a
  detection-nondetection likelihood available and the observation process
  modellable. Without this sentence, P3's problem has no setup.
- Optional: Rosenberg et al. 2019 for the scale of what trend estimation can
  demonstrate; Johnston et al. 2025 if we want to foreshadow testing whether
  declines concentrate in abundant species.

## P3 — Why mammals are harder, and what is already in place

- Mammals are more secretive and harder to sample than birds.
- iNaturalist mammal observations are growing rapidly (Jensen et al. 2026),
  **but carry no effort field at all** — not merely "less clear" than eBird.
  State this as the specific structural difference, since it is what forces the
  count-with-offset formulation later.
- **ADD — the published caution.** Soberon and Christen (2026) is the one paper
  directly on using iNaturalist to estimate trends. Richness per unit effort
  trends artifactually negative as effort grows; abundance per unit effort is
  the safer metric; and record counts approximate distinct individuals sighted
  when populations are large. Citing this sharpens the gap and pre-empts the
  obvious reviewer objection, while supporting our choice of response variable.
- Camera-trap networks now provide systematic mammal data (Snapshot USA, and
  internationally Snapshot Global including national monitoring of 35 Chilean
  protected areas).
- However, camera arrays remain low-density samples relative to the continent.
- **ADD — the closest precedent, which must be cited here.** Jachowski et al.
  (2021) assessed three North American weasel species using four streams:
  historical harvest, museum records, iNaturalist, and a US-wide trail camera
  survey. Same data, same taxa, same continent. **Omitting this invites "this
  has already been done."** The difference to state plainly: they compare
  streams descriptively and triangulate; we fit them jointly to a shared
  surface and estimate one trend.
- Goldstein et al. showed the potential of integrating camera traps and
  iNaturalist, with cameras as the reference stream, and showed the two streams
  agreed in ecological signal and spatial pattern.
- However, no temporal trend was included.
- **ADD — the logical step that makes P4 a finding rather than an increment.**
  Spatial agreement between two streams does not imply temporal agreement. The
  processes that generate iNaturalist records change over time in ways that
  have nothing to do with the animals: more observers, better phone cameras,
  shifting attention. Whether the streams still agree once a time axis is added
  is an open and testable question, not a safe assumption.

## P4 — Hypothesis and predictions

**ADD this as its own short paragraph (or the first half of P5).** The current
outline has no hypothesis, and the competing-predictions structure is both the
house style and genuinely available here.

- If iNaturalist records carry population signal over time, the two streams
  should produce concordant trends, and the integrated estimate should simply be
  more precise than either alone.
- If iNaturalist's temporal signal is dominated by observer dynamics, the
  streams should diverge in time even while agreeing in space, and the choice of
  which model quantity to report becomes the substantive decision.
- These make different predictions about the *sign agreement* between streams,
  which we can test species by species.

## P5 — What we did and what we found

- We add a temporal dimension to the Goldstein et al. framework and fit it
  across a fleet of North American mammals.
- **Numbers to use (verified):** 26,517 camera deployments across 60 projects,
  966,225 survey nights, 2008-2025; 1,053,349 research-grade iNaturalist records
  across the fleet over the same period, within species-specific range masks.
  The iNaturalist record pool grew from 2,798 records in 2008 to 369,355 in
  2025, which is both the opportunity and the problem.
- **Species count — decide before drafting.** The fleet is 34 species; 26 have
  jointly converged; 18 are currently reportable. "~30" is defensible only if
  Tier B and C land. Write the number the results section will actually carry.
- Spatially varying coefficients for climate response, and the shared spatial
  field. **Keep this to one sentence unless an SVC result is a headline**;
  otherwise it is Methods.
- Range masks defined as the union of published range polygons and iNaturalist
  presence, which keeps documented populations outside static range maps and
  retains true-absence cells. **Worth one clause** — it is a real methodological
  choice with consequences, not bookkeeping.
- **Power analysis / simulation — resolve its status.** If the estimator-choice
  simulation is part of this paper, it needs a sentence of motivation in P3
  (which camera-side observation model, and whether array-level data suffice).
  If it is the companion paper, cite it here and drop it from the outline.
- **"One sentence about what we discover" — make it the real finding.** The
  camera-anchored and iNaturalist-side trends disagree in sign for 6 of 18
  species, including raccoon, black bear, bobcat, red fox and coyote, where
  cameras show increases and iNaturalist shows declines. Reporting the
  iNaturalist-side quantity would have reversed those conclusions. Across the
  fleet, 14 of 18 species are increasing and 4 declining on the camera-anchored
  estimate.

---

## Three framing notes

1. **Lead the contribution with the disagreement, not with "we added time."**
   Adding a time axis is the method; discovering that the streams diverge in
   time while agreeing in space is the result, and it is what makes the paper
   matter to anyone using opportunistic records for trends.
2. **State the honest counterweight once, early (end of P3 or in Discussion).**
   The bird literature solved a version of this with better data, because eBird
   carries effort and iNaturalist does not. Our answer is the camera stream, and
   the paper should say that iNaturalist contributes spatial coverage rather
   than temporal signal.
3. **A fleet is a claim, not just a sample size.** Fitting ~30 species on a
   common footprint lets us ask which species *can* be monitored this way, which
   a single-species paper cannot. Worth claiming explicitly in P5.
