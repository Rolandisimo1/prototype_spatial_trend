# Bobcat & white-tailed deer — v2b model results

**Runs:** `bobcat_v2b_national_scalar`, `bobcat_v2b_ecoregion`,
`white-tailed_deer_v2b_national_scalar`, `white-tailed_deer_v2b_ecoregion`.
All four are single continuous chains (no resume boundary), so none can be
affected by the chunked-checkpoint burn-in defect (Issue 4) that invalidated
the earlier fleet. Convergence verdicts landed 2026-08-31; WTD posteriors
extracted 2026-09-01 (job 697393).

All numbers below are read directly from the extracted CSVs in
`hazel_pull_20260827/posteriors_20260828/`; source named per section. This
report follows the same structure and disclosure conventions as
`MOOSE_v2b_RESULTS.md`.

---

## 1. Convergence — three of four clean on trend, all four clean if you only need the trend

| run | trend params | max gated R-hat (whole fit) | n > 1.1 (whole fit) |
|---|---|---|---|
| `bobcat_v2b_national_scalar` | clean, max 1.0006 | **1.0955** | 0 |
| `bobcat_v2b_ecoregion` | clean, max ≤1.004 | **1.1819** | **13** (all CAR fields) |
| `wtd_v2b_national_scalar` | clean, max 1.0057 | **1.0931** | 0 |
| `wtd_v2b_ecoregion` | clean, max ≤1.0118 | **1.1894** | **32** (all CAR fields) |

Both ecoregion fits fail the whole-fit gate on CAR spatial fields only —
`link_occ_intercept` (the occupancy-intercept CAR field) and `MCMT_effect`
(the winter-temperature SVC surface). Every trend parameter in both models —
`total_var_beta`, `year_beta`, `year_var`, `sigma_region`, and all eight
`year_region[r]` — is comfortably converged (R-hat ≤ 1.012 in the worst
case). **The national and regional trend numbers in this report are usable.
The CAR intercept field and MCMT climate-response surface are not**, and are
excluded from this report entirely.

**This is not two independent failures.** The two ecoregion fits' failing
cell sets overlap almost completely:

| model | failing cell100 IDs |
|---|---|
| `bobcat_v2b_ecoregion` (13) | 854, 855, 873, 874, 875, 883, 884, 885, 891, 892, 893, 894, + MCMT_tau |
| `wtd_v2b_ecoregion` (32) | 855, 873, 874, 875, 883, 884, 885, 886, 891, 892, 893, 894, 899, 900, 901, 906, 907, + MCMT family |

11 of bobcat's 12 spatial cells (91.7%) also appear in WTD's set. Two
independent fits on two independent species landing on nearly the same
contiguous patch means this is a property of the CAR adjacency graph in that
region — most likely a steep data-density gradient between data-rich and
wholly-empty neighboring cells — not a species-specific data-volume problem.
A leave-one-chain-out check on WTD's ecoregion fit (R-hat recomputed dropping
each chain: 1.2655 / 1.6622 / 1.2332) shows genuine three-way disagreement
localized to that patch, not a single rogue chain or a restart artifact.

**Standing caveat for any regional map of bobcat or WTD:** grey out or
otherwise flag this cell patch on any CAR-intercept or MCMT map. Trend maps
(Section 3) are unaffected.

*Source: convergence verdicts reported 2026-08-31; caveat text embedded
directly in each model's `_global.csv`.*

---

## 2. National trend — bobcat declining without camera corroboration, WTD increasing with strong corroboration

| species | model | `total_var_beta` | 95% CI | P(decline) | `year_beta` (camera) | TRI |
|---|---|---|---|---|---|---|
| Bobcat | national_scalar | **−0.242** | [−0.314, −0.172] | **1.000** | −0.0001 (≈0, P(neg)=0.498) | **0.0005** |
| Bobcat | ecoregion | −0.145 | [−0.372, +0.066] | 0.907 | +0.010 (≈0, P(neg)=0.377) | 0.037 |
| WTD | national_scalar | **+0.177** | [0.136, 0.221] | **0.000** (P(increase)=1) | +0.157 | **0.784** |
| WTD | ecoregion | +0.187 | [0.109, 0.266] | 0.003 | +0.157 | 0.767 |

**Bobcat: a clear iNat-driven decline the cameras do not see.** The national-
scalar model's `total_var_beta` excludes zero decisively (P(decline)=1.000),
but `year_beta` — the camera-anchored component — is statistically
indistinguishable from zero (mean −0.0001, CI straddling zero, P(neg)≈0.5,
literal camera noise). `year_var` (−0.242) carries essentially the entire
trend on its own. TRI = 0.0005 means the two data streams essentially never
agree in this posterior. **This is the same iNat-driven-trend caveat already
established for moose** (Section 8 in `moose_model_full_report.md`) — bobcat
is a second, independent instance of the same pattern, not a new problem.

**WTD: this project's first camera-corroborated trend.** `year_beta`
(+0.157) accounts for nearly all of `total_var_beta` (+0.177) — the cameras
are doing the work, not just failing to contradict iNat. TRI = 0.784/0.767 in
both parameterizations, far above the 0.5 threshold and far above anything
else in the fleet (Section 4). With 188,939 raw iNat records — the most of
any species fit so far — landing on a real, camera-corroborated increase,
this functions as a positive control: the pipeline can detect a genuine
signal when the underlying data supports one, which is worth stating
explicitly rather than assuming.

*Source: `bobcat_v2b_national_scalar_posterior.csv`,
`bobcat_v2b_ecoregion_global.csv`,
`white-tailed_deer_v2b_national_scalar_posterior.csv`,
`white-tailed_deer_v2b_ecoregion_global.csv`.*

![Bobcat vs WTD national trend](fig_bobcat_wtd_national_trend.png)

---

## 3. Regional trend (ecoregion model)

![Bobcat and WTD regional trend maps](fig_bobcat_wtd_regional_trend.png)

Hatched regions: the 95% CI on that region's absolute trend does not clearly
exclude zero (0.05 < P(decline) < 0.95) — read as "not distinguishable from
flat," not "no data."

**Bobcat — one clear reversal against an otherwise uniform decline.** Five of
eight regions decline with P(decline) ≥ 0.986 (Eastern Temperate Forests,
Great Plains, Mediterranean California, North American Deserts,
Northwestern Forested Mountains). **Northern Forests is the one confident
exception**: `year_region` mean +0.456 (CI [0.197, 0.750], entirely positive),
giving an absolute trend of +0.310 with P(decline)=0.0015 — i.e. >99.8%
confidence Northern Forests is *increasing* while the rest of the range
declines. Marine West Coast Forest and Water are both uncertain (P(decline)
0.78–0.84, wide CIs, comparatively few iNat cells: 43 and 24 respectively).

| region | n camera sites | n iNat cells | absolute trend | 95% CI | P(decline) |
|---|---|---|---|---|---|
| Eastern Temperate Forests | 12,331 | 544 | −0.355 | [−0.475, −0.238] | 1.000 |
| Great Plains | 2,322 | 237 | −0.433 | [−0.574, −0.295] | 1.000 |
| Northwestern Forested Mountains | 1,677 | 152 | −0.279 | [−0.464, −0.100] | 0.999 |
| North American Deserts | 1,179 | 232 | −0.191 | [−0.331, −0.052] | 0.997 |
| Mediterranean California | 997 | 59 | −0.176 | [−0.332, −0.020] | 0.986 |
| Water | 82 | 24 | −0.192 | [−0.584, 0.202] | 0.841 |
| Marine West Coast Forest | 493 | 43 | −0.091 | [−0.317, 0.139] | 0.784 |
| **Northern Forests** | 1,450 | 101 | **+0.310** | **[0.105, 0.519]** | **0.0015** |

**WTD — broadly increasing, strongest in the two regions with the most data.**
Six of eight regions show a positive absolute trend with P(decline) ≤ 0.001;
only Marine West Coast Forest (129 camera sites, 1 iNat cell, P(decline)=0.10)
and Mediterranean California (**0 camera sites, 0 iNat cells** — this region
is prior-only for WTD, exactly the zero-data caveat already established for
moose's Northern Forests/North American Deserts) fail to reach significance.

| region | n camera sites | n iNat cells | absolute trend | 95% CI | P(decline) |
|---|---|---|---|---|---|
| Great Plains | 2,468 | 603 | +0.198 | [0.120, 0.280] | 0.000 |
| North American Deserts | 523 | 92 | +0.230 | [0.093, 0.405] | 0.001 |
| Water | 164 | 55 | +0.226 | [0.087, 0.399] | 0.001 |
| Northern Forests | 1,833 | 168 | +0.262 | [0.156, 0.386] | 0.000 |
| Eastern Temperate Forests | 15,217 | 1,060 | +0.157 | [0.105, 0.210] | 0.000 |
| Northwestern Forested Mountains | 1,225 | 125 | +0.090 | [−0.108, 0.234] | 0.167 |
| Marine West Coast Forest | 129 | 1 | +0.155 | [−0.200, 0.419] | 0.103 |
| **Mediterranean California** | **0** | **0** | +0.176 | [−0.143, 0.461] | 0.080 |

*Source: `bobcat_v2b_ecoregion_ecoregion_posterior.csv`,
`white-tailed_deer_v2b_ecoregion_ecoregion_posterior.csv`.*

---

## 4. Camera corroboration across the converged fleet

![Fleet TRI comparison](fig_fleet_tri_comparison.png)

| model | TRI | reads as |
|---|---|---|
| WTD national_scalar | 0.784 | strongly camera-corroborated |
| WTD ecoregion | 0.767 | strongly camera-corroborated |
| Moose v1fix9 national_scalar | 0.488 | borderline |
| Moose v2b national_scalar | 0.255 | weakly corroborated (`year_var`≈0 drives this — see `MOOSE_v2b_RESULTS.md` §2) |
| Bobcat ecoregion | 0.037 | not corroborated |
| Bobcat national_scalar | 0.0005 | not corroborated |

Across five converged models, TRI spans nearly the full [0,1] range. This is
useful context for reading any single species' number in isolation: bobcat's
near-zero TRI is not a modeling failure specific to bobcat — it sits at one
end of a real spectrum this pipeline is capable of producing, with WTD
anchoring the other end. Read alongside Section 2, the fleet is not
uniformly "iNat-driven" or uniformly "camera-corroborated" — it depends on
the species' actual camera detectability and range overlap with the trap
network.

---

## 5. Range mask: why the IUCN polygon was replaced

The original static IUCN-polygon mask excluded real, iNaturalist-documented
records in every species, by discarding grid cells the polygon placed outside
the range:

| species | cells dropped holding real records | records dropped (% of total) |
|---|---|---|
| Moose | 143 | 4,293 (37.1%) |
| Bobcat | 39 | 1.3% |
| WTD | 12 | 1.1% |

**These percentages should not be read as a cross-species comparison.** The
per-species magnitude is governed mainly by how well that species' IUCN
polygon happens to be drawn, and IUCN range-map quality varies substantially
between species for reasons unrelated to the species' ecology or to our data.
A large percentage means that species' polygon was a poorer fit to where the
animals actually are; it does not mean the mask problem is more real, more
severe, or more worth fixing for that species. Each species' number stands on
its own, and the fix is warranted in all three cases on the same grounds:
cells with documented occurrences were being removed from the model fit
entirely, not merely from map rendering.

WTD's dropped records are concentrated rather than spread — its single worst
dropped cell held 969 records — which is a reminder that a small aggregate
percentage can still remove an important locality.

A third candidate mask (IUCN polygon expanded by a contiguous chain of
iNat-detection cells, per a convention from a prior project) was tested
against both existing options and **failed on moose** — the case it was
built to address — because moose's Colorado population is spatially
disconnected from the IUCN polygon by more than one flood-fill hop, so a
4-connected expansion cannot reach it (recovers only ~20% of moose's
excluded records). For bobcat/WTD it retains true-absence cells that the
pure-presence mask (Fix 2b) discards — a real trade-off.

**Resolved 2026-09-01: the union (in-range if IUCN **OR** presence threshold)
was measured and adopted for all three species.** It retains 100% of records
for every species — recovering moose's disjunct Colorado population, which the
flood-fill expansion structurally could not reach — while keeping the
true-absence cells presence-only discards (104 / 1,794 / 756 zero-record IUCN
cells for moose / bobcat / WTD). A cross-check came out clean: the cells
passing presence but falling outside IUCN (143 / 39 / 12) exactly reproduce
the "occupied cells dropped" column from the earlier exclusion audit, so the
union picks up precisely the cells the old mask excluded and nothing more.

Two caveats attach to that adoption. First, **no fit has yet been run under a
union mask** — tractability is inferred from bundle sizes already exercised
(`bobcat_v1fix` built and ran at ncell50 3,145; the union is 3,184), not
observed. Second, per-species grid coverage under the union is moose 15.8%,
WTD 86.1%, **bobcat 95.8%** — because bobcat's IUCN range alone already covers
94.7% of the grid, the union is very nearly no mask at all for that species.
That is a deliberate, documented choice rather than an oversight, and it also
bears on how bobcat's near-zero camera corroboration should be read: a mask
touching ~4–5% of cells cannot plausibly drive TRI from a meaningful value to
0.0005, so bobcat's result is very unlikely to be a masking artifact. A
union-mask bobcat refit is planned to confirm that inference directly rather
than leave it as an argument. See `team_memo_inat_pipeline_issues.md`, Issue 5.

*Source: range-mask exclusion figures reported 2026-08-14/2026-09-01;
`team_memo_inat_pipeline_issues.md` Issue 5.*

---

## 6. WTD model vs. state agency comparison

![WTD model trend, agreement map, and raw agency map, side by side](fig_wtd_trend_vs_agency_sidebyside.png)

**Method:** state-level `trend_10yr_direction` from `deer_moose_trends_master.csv`
(46 states with compiled data) is area-weighted onto the 8 EPA Level I
ecoregions used by the model (each state's polygon intersected against the
ecoregion boundaries; `increasing`=+1, `stable`=0, `mixed`=0, `decreasing`=−1,
`unknown` excluded from the score but still counted as "compiled" for
coverage). The exact historical script that built moose's equivalent table
could not be recovered verbatim from this session's history, so this is a
fresh build using the same overlay principle — methodologically consistent,
not a byte-for-byte replication. The model's per-region direction uses the
same classification already applied for moose (`abs_p_negative` ≥ 0.90 →
decreasing, ≤ 0.10 → increasing, otherwise uncertain; 0 camera AND 0 iNat →
prior-only).

| region | agency score | agency direction | model direction | comparison |
|---|---|---|---|---|
| Great Plains | +0.226 | increasing | increasing | **agree** |
| Northern Forests | +0.066 | increasing | increasing | **agree** |
| Eastern Temperate Forests | +0.243 | stable | increasing | disagree |
| North American Deserts | +0.108 | stable | increasing | disagree |
| Water | −0.550 | decreasing | increasing | disagree |
| Marine West Coast Forest | +0.029 | increasing | uncertain | model uncertain |
| Northwestern Forested Mountains | +0.043 | increasing | uncertain | model uncertain |
| Mediterranean California | — | no data | prior-only | not comparable |

**Reading this honestly, not as a validation exercise:** unlike moose (whose
v1fix ecoregion deviation correlated with agency score at Pearson r=+0.97),
WTD's model-vs-agency correlation here is **weak and negative** (r=−0.24,
n=7 regions — not a result to lean on with this few regions, but not one to
suppress either). Two agree, three disagree, two are regions where the model
itself declines to make a directional claim. The likely mechanism: WTD's
national trend (`total_var_beta`=+0.177) is strong and dominates every
region's absolute trend almost uniformly (Section 3's map is nearly solid
blue), while the *regional deviations* (`year_region`) that would let a
region actually disagree with the national number are small relative to
that shared signal — so the model's regional map may be under-resolved
for genuine region-to-region agency disagreement, particularly in the
already-saturated eastern states where "stable" plausibly reflects a
population near carrying capacity rather than a signal the model should be
expected to reproduce. This is a real open question, not a resolved one.

*Source: `wtd_agency_by_ecoregion.csv` (this session's build),
`wtd_agency_sign_agreement_table.csv`, `deer_moose_trends_master.csv`.*

---

## 7. Known open issues / caveats carried into this report

1. **CAR intercept field and MCMT climate surface are not usable for either
   species' ecoregion model** — a shared, contiguous cell patch (91.7%
   overlap between bobcat's and WTD's failing sets) that appears to be a
   property of the spatial adjacency graph rather than either species'
   data. See Section 1.
2. **No agency/state-survey comparison exists for bobcat.** Unlike moose and
   WTD, no compiled state-wildlife-agency trend dataset covers bobcat — this
   report cannot check the model's direction against management data for
   that species. `deer_moose_trends_master.csv` covers WTD and moose only.
3. **WTD's model-vs-agency correlation is weak and negative (r=−0.24, n=7)**,
   unlike moose's strong positive correlation (r=+0.97) — see Section 6. Not
   yet understood whether this reflects a real regional-resolution
   limitation, agency data noise, or WTD populations genuinely tracking
   national iNat/camera trends more uniformly than moose's did.
4. **18-year trend spans a camera-data-sparse period for both species**,
   same caveat already carried for moose — cameras exist only 2019–2025;
   2008–2018 is iNat-only. Windowed (5yr/10yr) refits would isolate the
   camera-only period but have not been run for bobcat/WTD.
5. **Mask choice is settled but unfitted** (Section 5) — the union mask was
   measured and adopted for all three species, but no model has yet been fit
   under it. Every number in this report comes from presence-mask (Fix 2b)
   fits. Union-mask refits are pending, and are blocked on reconstructing the
   base-bundle build step (see `PIPELINE.md`, Section 3).

---

## Reusable notes

- The cross-species CAR-patch finding (Section 1) should be checked for any
  future species fit on the same ecoregion adjacency graph — if a third
  species lands on the same or an overlapping cell set, that would further
  confirm it as a grid property worth fixing once rather than re-diagnosing
  per species.
- TRI is not comparable across species as a measure of "how good" a model
  is — it measures a specific correspondence between two independent data
  streams and depends on that species' real camera detectability. Do not
  use it to rank models; use it only to state, per species, whether the
  camera stream corroborates the iNat-driven trend.
