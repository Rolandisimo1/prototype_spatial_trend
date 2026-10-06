# STATUS — camera/iNaturalist abundance fleet

Live state of the primary deliverable. Update whenever a species moves from
fitted to reportable, or a side-quest opens or closes. Read this when resuming
after a gap rather than reconstructing project state from conversation.

**Last updated:** 2026-10-06

---

## Primary deliverable

Continental population trends for the mammal fleet, camera-anchored, as a
report for Krishna and Arielle.

| | Count | |
|---|---|---|
| **Reportable** | **4** | bobcat, coyote, moose, white-tailed deer |
| **Fitted, not extracted** | **14** | black bear, jackrabbit, chipmunk, fox squirrel, gray squirrel, elk, grey fox, mule deer, armadillo, raccoon, red fox, red squirrel, skunk, opossum |
| Fleet once extracted | 18 | all of Tier A plus moose |
| Not yet fitted | 16 | Tier B (9) and Tier C (8), minus moose |

**Current blocker:** extraction, not compute. All 14 have finished, converged
`v2u_national_scalar` chains sitting on Hazel. Nothing is queued behind them.
Requested 2026-10-06: one table matching `trend_decomp_union_vs_presence.csv`
for all 18, with R-hat on the trend parameters in the same table.

**Behind that:** union ecoregion fits (moose, bobcat, deer) and shared-ecoregion
fits (deer, bobcat) are finished and also unextracted. Those carry the regional
trends and the shared-pathway comparison.

---

## Settled, and not to be reopened

- **Report `year_beta`, not `total_var_beta`.** The camera-anchored parameter is
  the population trend. For bobcat the two differ in sign.
- **Union mask** (IUCN polygon OR iNaturalist records) is the design. The
  presence-only and polygon-only variants are superseded and their comparison
  figures deleted.
- **Nine-covariate design**, dropping `soil_sand` fleet-wide.
- **Single-shot chains**, no chunked resume.
- Convert to percent per calendar year before comparing streams; the two
  submodels standardise their year covariate differently.

---

## Parked side-quests

Each is real, each has changed something, none is the deliverable. Listed so
they can be resumed deliberately rather than by whichever report lands next.

| Item | State |
|---|---|
| Coyote `year_pop` / `year_pop_var` tier-2 test | chains finished 10/4, unextracted; leverage-concentration check requested with it |
| Project fixed effect, 4 camera-gradient species | done for chipmunk, fox squirrel, mule deer — unreported; gray squirrel segfaulted twice (exit 139), needs rerun with more memory |
| Urban-rural stream screen | ran on 16 species; camera arm resolves the gradient for only 3, so most verdicts are "not assessed" |
| Fisher range mask | built; 4 Tennessee cells removed per domain judgement, recorded in `manual_range_overrides.csv` |
| Tier B/C bundles | not built; fisher's override lives outside the pipeline and must be applied when Tier B is bundled |

---

## Known open questions

- **Bobcat's stream divergence is unexplained.** Cameras +1.0%/yr, iNaturalist
  -3.4%/yr. Not effort mis-correction, not an urban gradient — the urban
  explanation was tested and withdrawn.
- **The camera arm is the binding constraint on every stream comparison**,
  3-23x wider intervals than iNaturalist. Most comparisons cannot be made.
- `theta1` and the trend-robustness indicator rank species in inverse order,
  n=4, unexplained.
