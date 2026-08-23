# Raw camera data — verification and revised estimator plan

Files: `/Users/rwkays/claude_code/data_integration_arielle/raw_cam_data/`
- `combined_deployments_all.csv` — 26,798 rows x 11 cols (4.0 MB)
- `combined_sequences_all.csv` — 2,517,928 rows x 17 cols (561 MB)

## 1. Does it match what we've been using? Yes — it is a superset

| quantity | raw data | `umflist.RDS` siteCovs (in use) | read |
|---|---|---|---|
| deployment rows | 26,798 | 24,869 | raw is superset; prep filters ~7% |
| `project_id` | 60 | 58 | consistent |
| `subproject_name` | 962 | 943 | consistent |
| year range | 2008–2025 (17 distinct) | 2008–2025 | matches |
| deployments spanning >1 yr | 240 of 26,517 | 209 of 246 dup rows | same phenomenon, same scale |

**Linkage is perfect:** every one of the 26,517 unique `deployment_id`s in the
sequences file appears in the deployments file and vice versa — zero orphans in
either direction.

**Effort cross-check passes.** `survey_nights` median 31 days -> 3.1 ten-day
windows, against the prep pipeline's measured mean `J` of 3.22. Independent
confirmation that the window discretization we've been using is faithful to the
raw deployment durations.

**Conclusion: safe to use.** Same projects, same subprojects, same years, same
deployments, consistent effort. The prep filtering (~7% of rows) is expected
and should be characterized before a real fit, but nothing here contradicts
what the models have been fitting.

## 2. Does this enable N-mixture? Not usefully — and this is the key finding

`group_size` exists, which was the missing ingredient. But it is close to
degenerate for our focal species:

| species | sequences | `group_size` NA | mean group | % groups >1 | max |
|---|---|---|---|---|---|
| Moose | 1,510 | 5.3% | **1.08** | 6.9% | 4 |
| Bobcat | 10,462 | 41.4% | **1.04** | 2.7% | 4 |
| White-tailed deer | 581,083 | 20.8% | 1.22 | 16.7% | 56 |

Moose and bobcat are **solitary** — 93% and 97% of detections are a single
animal. A classic N-mixture model on these counts would be estimating
abundance from a variable that is almost always 1. There is essentially no
group-size signal to exploit for the two species this project has walked back
to.

**Worse, the missingness is structured, not random.** `group_size` is 100% NA
for exactly 2 of 27 project groups (`NS` = 551,297 sequences, `CA` = 131,537)
and 0% NA for the other 25 — zero projects are partially missing. So the gap is
whole-project metadata absence (those projects never recorded group size),
which means dropping NA rows would drop two entire projects and their
geography. `camera_trap_array` shares almost exactly the same missing rows
(682,834 of 682,863 identical), confirming a common provenance rather than
per-record data loss.

Bobcat's 41.4% NA rate is the practical killer: it is driven by which projects
happen to hold bobcat detections, and it is not missing-at-random with respect
to space.

**However — the raw data DOES unlock something better.** Sequence *count* per
deployment is richly graded, independent of `group_size`:

| species | detecting deployments | sequences/deployment (median / mean / p95 / max) | % with >1 |
|---|---|---|---|
| Moose | 338 | 2 / 4.5 / 14 / 146 | **53.8%** |
| Bobcat | 3,548 | 1 / 2.9 / 9 / 178 | **47.5%** |
| White-tailed deer | 19,059 | 13 / 30.5 / 113 / 1,358 | **92.9%** |

This is detection *frequency*, which is exactly what Royle–Nichols exploits
("more animals produce more detections"), and it is available at full coverage
with no structured missingness. **It is a far better abundance signal for
solitary species than group size**, because for moose and bobcat abundance
expresses itself as more visits, not larger groups.

## 3. A third array field exists, and it is better on one axis

`camera_trap_array` in the sequences file is distinct from `subproject_name`:

| field | groups | populated | x year: units | cams/unit med / mean / max | singletons |
|---|---|---|---|---|---|
| `subproject_name` | 962 | **99.9%** | 1,512 | 13 / 17.7 / **168** | 3.6% |
| `camera_trap_array` | 682 | **63.0%** | 1,026 | 14 / 16.4 / **158** | **1.2%** |

`camera_trap_array` is the *ecologically intended* grouping (fewer singletons,
tighter) but only 63% populated, and its missingness is the same structured
whole-project gap. **Neither field solves the size problem on its own** — both
still produce a max unit of ~160 cameras, so the spatial-splitting step and its
sign-off remain necessary regardless of which field is chosen.

Recommendation: use `camera_trap_array` where present, fall back to
`subproject_name` elsewhere, record which source each unit came from, and
report the array-construction diagnostics separately by source so a
field-driven artifact would be visible.

## 4. Revised plan

**Arms: still three, but the third is now better justified.**

1. `camera_occ` — baseline (unchanged)
2. `array_occ` — array-level occupancy (unchanged)
3. `array_rn` — Royle–Nichols on **sequence counts**, not group sizes

**Classic N-mixture stays dropped**, now for a stronger reason than before. In
`sim_estimator_feasibility.md` I dropped it because counts were unavailable.
Counts *are* available — and the honest finding is that they are uninformative
for solitary species (mean group 1.04–1.08) and structurally missing for 27% of
records. That is a better-evidenced rejection, and it should replace the
earlier one.

**One genuine addition the raw data makes possible:** RN can now be fit to
*graded* detection frequency rather than only detection/non-detection. Roughly
half of detecting deployments for moose and bobcat have >1 sequence, so there
is real signal beyond presence/absence. This strengthens arm 3 materially and
is the main reason the raw data is worth incorporating.

**What does not change:** the design invariant (trend block byte-identical,
md5 `69213000bc266620f68b21c71316358b`), the two metric families (power and
discrimination, per the Brazil power/AUC divergence), the abundance crossing,
the per-replicate reseed, no chunk-resume, and the requirement for PI sign-off
on array construction before launch.

**New prerequisite:** a sequence-count matrix per deployment-year must be built
from `combined_sequences_all.csv` and carried into `sim_inputs.RDS` alongside
`site_array`. This is a straightforward aggregation but it is new work, and it
should be characterized against the existing binary `y` (every deployment with
count >= 1 must have `y = 1` in the current data; any mismatch means the
filtering differs and must be understood before use).

## 5. Caveats to carry forward

- **Prep filtering is uncharacterized.** ~7% of raw deployment rows do not
  reach the model. Benign pending explanation, but it should be explained
  rather than assumed before any real fit uses raw-derived quantities.
- **87 deployments have `survey_nights = 0`** and 19 have NA year. Trivial in
  volume; must be handled explicitly rather than silently coerced.
- **`survey_nights` max is 2,688 days** (~7.4 years) — almost certainly a
  data-entry artifact or a deployment that was never closed. Needs a cap or
  exclusion rule, especially since effort now enters the array arms as a
  covariate.
- **Sequence counts are not independent detections.** Sequences are
  algorithmically grouped triggers; a single animal lingering can produce
  several. RN assumes detections accumulate with abundance, so temporal
  clustering within a deployment inflates counts in a way RN reads as more
  animals. A minimum-separation rule (e.g. 30 min between sequences of the
  same species) is the standard mitigation and should be applied and reported.
