# STATUS — camera/iNaturalist abundance fleet

Live state of the primary deliverable. Update whenever a species moves from
fitted to reportable, or a side-quest opens or closes. Read this when resuming
after a gap rather than reconstructing project state from conversation.

**Last updated:** 2026-10-05 (late)

---

## Primary deliverable

Continental population trends for the mammal fleet, camera-anchored, as a
report for Krishna and Arielle.

| | Count | |
|---|---|---|
| **Extracted (camera trend in table)** | **18** | all Tier A plus moose — `trend_table/trend_table_v2u_national_scalar_18spp.csv` (2026-10-05) |
| — camera direction resolved | 12 | increase 9 (deer, coyote, black bear, chipmunk, gray squirrel, elk, raccoon, red squirrel, opossum); decline 3 (jackrabbit, grey fox, mule deer) |
| — NOT RESOLVABLE (95% CrI spans 0, wider than ±2 %/yr) | 6 | moose, bobcat, fox squirrel, armadillo, red fox, skunk (skunk +0.09 [−2.01, +2.25] is effectively a precise near-zero, at the rule's edge) |
| Not yet fitted | 16 | Tier B (9) and Tier C (8), minus moose |

All 18 converged on the trend parameters (max R-hat 1.0098, gray squirrel). The 4 previously
reported species reproduce `trend_decomp_union_vs_presence.csv` exactly.
Streams opposite in direction with both resolved: black bear, raccoon, mule deer → report camera,
flag iNaturalist as questioned.

**Regional (ecoregion) fits extracted** — `trend_table/trend_table_ecoregion_{national,regional}.csv`:
- shared ecoregion, deer + bobcat: converged (R-hat ≤ 1.003) — these are the usable regional trends
- union ecoregion, moose: converged (1.011), almost every region not resolvable
- union ecoregion, deer: marginal (1.065)
- **union ecoregion, bobcat: NOT CONVERGED** (trend R-hat 2.48, year_region 5.57) — do not report
- Coyote shared ecoregion: chains c1/c2 finish 10/5 night; R-hat job 1111767 queued on them; then extract.

**Current blocker:** Tier B/C bundles (16 species) not built. Next deliverable step after
the write-up of the 18.

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
| Coyote `year_pop` / `year_pop_var` tier-2 test | extraction + leverage job 1111765 ran 10/5, results not yet read |
| Project fixed effect, 4 camera-gradient species | done for chipmunk, fox squirrel, mule deer — unreported; gray squirrel segfaulted a 3rd time at 128G (job 1111759, peak 2.5 GB) — NOT memory, lme4 crash in first project-FE fit; species uncovered |
| Convergence sweep, whole-model R-hat | array 1111766 (21 runs) ran 10/5, results in `~/isdm/fleet_v2u/rhat_sweep/`, not yet read |
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
