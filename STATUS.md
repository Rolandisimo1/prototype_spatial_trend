# STATUS — camera/iNaturalist abundance fleet

Live state of the primary deliverable. Update whenever a species moves from
fitted to reportable, or a side-quest opens or closes. Read this when resuming
after a gap rather than reconstructing project state from conversation.

**Last updated:** 2026-10-08

---

## Primary deliverable

Continental population trends for the mammal fleet, camera-anchored, as a
report for Krishna and Arielle.

| | Count | |
|---|---|---|
| **Reportable** | **18** | all Tier A plus moose — `trend_table/trend_table_v2u_national_scalar_18spp.csv` |
| — increasing | 9 | gray squirrel, red squirrel, elk, black bear, opossum, raccoon, deer, chipmunk, coyote |
| — declining | 3 | jackrabbit, mule deer, grey fox |
| — no change detected (tight interval spanning zero) | 5 | armadillo, fox squirrel, red fox, bobcat, skunk |
| — not resolvable (interval too wide to state anything) | 1 | moose |
| **Fitted, awaiting R-hat + extraction** | **7** | badger, puma, porcupine, snowshoe hare, fisher*, western gray squirrel, California ground squirrel, collared peccary |
| Fleet once extracted | **26** | |
| Held | 8 | Tier C — none projects a usable trend interval |

All 18 converged; worst trend R-hat 1.010 (gray squirrel).

**Current blocker:** extraction of the 7 finished Tier B species, requested
2026-10-08 along with their R-hat. *Fisher launching separately (bundle built
late, flag accepted 10/8) — 3 chains, ~12 h each.

**Headline method finding, only visible at fleet scale.** Regressing the
iNaturalist trend on the camera trend across the 18 gives a slope of
**0.28 [0.04, 0.62]** — iNaturalist recovers about a quarter of the
camera-measured trend and compresses six species past zero into the wrong sign
(black bear, raccoon, coyote, red fox, bobcat, mule deer). This is not
regression dilution: an errors-in-variables (Deming) fit using the posterior
intervals moves OLS 0.271 to only 0.285. It dissolves the bobcat "divergence"
puzzle — bobcat is one point on this line, not an anomaly. Mechanism test
requested: `theta1` attenuation vs `sigma_year_var` shrinkage vs already
present in raw records per unit effort. The raw-counts version is decisive.

**Regional trends:** use the shared-ecoregion fits (deer, bobcat converged at
R-hat <= 1.003). **Do not report union ecoregion bobcat** — R-hat 2.48 on the
trend, 5.57 on `year_region`, regional values to +27%/yr; the centred
`year_region` prior reproduces the deer funnel. Coyote shared ecoregion still
awaiting extraction.

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
| Coyote `year_pop` / `year_pop_var` tier-2 test | extraction + leverage job 1111765 ran 10/5, results not yet read |
| Project fixed effect, 4 camera-gradient species | CLOSED: chipmunk, fox squirrel, mule deer done; gray squirrel UNCOVERED (3 lme4 segfaults at 2.5 GB peak, not memory) — no further action |
| Convergence sweep, whole-model R-hat | array 1111766 (21 runs) ran 10/5, results in `~/isdm/fleet_v2u/rhat_sweep/`, not yet read |
| Urban-rural stream screen | ran on 16 species; camera arm resolves the gradient for only 3, so most verdicts are "not assessed" |
| Fisher range mask | built; 4 Tennessee cells removed per domain judgement (failed introduction), recorded in `manual_range_overrides.csv`; override now applied IN-pipeline by `build_union_grid_fleet.R` via `RANGE_OVERRIDES`. **Verified on disk 2026-10-08** — 176 cells, zero in TN/NC/KY, the 5 VA/WV reintroduction cells retained. The dropped cells held 246 camera site-years (177 after the covariate filter) with 0 fisher detections, but **208 of the 246 are from 2013 alone**, the other 38 from 2022–2025, median 2 occasions per site. So the camera record is *consistent with* the introduction having failed and does **not** independently confirm it — the override rests on domain judgement, not on the cameras. |
| Camera taxonomy splits | **FOUND 2026-10-08.** One common name carrying two scientific names, in the camera data: peccary (*Pecari*→*Dicotyles tajacu*, clean break 2023/24, 47.9% of the record lost if keyed to the old name — **this is the −54%/yr**), long-tailed weasel (same pattern, 50.3%), fisher, domestic pig — all one MDD id, union the names. **Elk is two accepted species** (*elaphus* 1006313 / *canadensis* 1006312); 855 NA records reassigned to *canadensis* per domain judgement, recorded in `fleet_union_prep/manual_taxonomy_decisions.csv`. Elk's +6.59%/yr needs rechecking after reassignment. Crosswalk cause: one name slot per source cannot hold a mid-series rename, and `mdd_binomial` was never validated against the MDD accepted list (peccary's was a synonym). |
| Tier C bundles | not built; held deliberately — none of the 8 projects a usable trend interval, see the deliverable table above. Fisher's override is now applied IN-pipeline by `build_union_grid_fleet.R` via `RANGE_OVERRIDES` (default `fleet34_isdm/manual_range_overrides.csv`), so it no longer needs applying by hand. |

---

## Known open questions

- **Bobcat's stream divergence is unexplained.** Cameras +1.0%/yr, iNaturalist
  -3.4%/yr. Not effort mis-correction, not an urban gradient — the urban
  explanation was tested and withdrawn.
- **The camera arm is the binding constraint on every stream comparison**,
  3-23x wider intervals than iNaturalist. Most comparisons cannot be made.
- `theta1` and the trend-robustness indicator rank species in inverse order,
  n=4, unexplained.
