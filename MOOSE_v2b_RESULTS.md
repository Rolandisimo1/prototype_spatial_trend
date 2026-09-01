# Moose — v2b model results and agency comparison

**Run:** `moose_v2b_national_scalar` and `moose_v2b_ecoregion`, converged
2026-08-26. **Single continuous chains — no resume boundary**, so these are the
first moose fits that structurally avoid the chunked-MCMC checkpoint defect that
invalidated the earlier fleet.

All numbers read from result files in `prototype_spatial_trend/`; source named
per section.

---

## 1. Convergence — clean

| run | parameters | max R-hat | n > 1.1 |
|---|---|---|---|
| `moose_v2b_national_scalar` | 2,750 | **1.0312** | 0 |
| `moose_v2b_ecoregion` | 2,759 | **1.0376** | 0 |

For context from the same 8/26 batch: `moose_v1fix9` also converged in both
parameterizations (1.093 / 1.076, 0 bad), and `bobcat_v2b_national_scalar` is
fully clean (max **gated** R-hat 1.0955, **0** offenders).
`bobcat_v2b_ecoregion` is the lone holdout at 1.1819 with 13 gated offenders.

CORRECTION (2026-08-27): an earlier version of this file described
`bobcat_v2b_national_scalar` as "clean (1.1145, 1 offender)". That was wrong.
The 1.1145 value belongs to `trend_robust_indicator`, which is **ungated** and
correctly so: it is `step(snr - 1)`, a binary derived quantity, and
Gelman-Rubin is not a valid convergence statistic for binary variables — the
same reason the RJMCMC indicator work monitored the data log-probability
instead. Counting it as an offender applies a diagnostic to a parameter class
it does not apply to. All 2,749 gated parameters converged.

*Source: `hazel_pull_20260827/convergence/gelman_single_*_all.csv`.*

## 2. National trend — a decline the cameras do not corroborate

`moose_v2b_national_scalar`:

| parameter | mean | median | sd | 95% CI | P(<0) |
|---|---|---|---|---|---|
| `total_var_beta` | **−0.1192** | −0.1179 | 0.0505 | [−0.2223, −0.0240] | **0.994** |
| `year_beta` (camera-anchored) | −0.1170 | −0.1131 | 0.0926 | [−0.3143, +0.0463] | 0.910 |
| `year_var` (iNat-specific) | −0.0022 | −0.0026 | 0.0996 | [−0.1967, +0.2068] | 0.519 |
| `trend_robust_indicator` | **0.255** | 0 | 0.436 | — | — |

The ecoregion model agrees closely on the global terms (`total_var_beta`
−0.1177, `year_beta` −0.1168, `year_var` −0.0008), so the trend estimate is not
sensitive to whether the regional term is included.

**Read this carefully.** The combined trend is negative with 99.4% posterior
probability — but `trend_robust_indicator` = P(`year_beta`/`year_var` > 1) is
**0.255 in national_scalar and 0.241 in ecoregion, both far below the 0.5
threshold.** That flag exists to ask whether the camera stream corroborates an
iNat-driven trend, and **moose fails it in both parameterizations.**

The mechanism is visible in the table: `year_var` is essentially zero with a CI
straddling zero, so the ratio `year_beta/year_var` is unstable by construction.
So the honest statement is *"a national decline of about −0.12 per year on the
link scale, driven by the camera-anchored component, with the robustness flag
not satisfied"* — not *"moose are declining nationally."*

*Source: `moose_v2b_national_scalar_posterior_clean.csv`,
`moose_v2b_ecoregion_global_clean.csv`.*

## 3. Regional structure — present but weak

`sigma_region` = 0.191 (median 0.166, 95% CI [0.029, 0.504]). The lower bound is
above zero, so the regional term is doing *something* — but the interval is wide
and the point estimate is small relative to the regional deviations themselves.

Per-region deviations from the national trend, ordered by camera support:

| region | cameras | iNat cells | `year_region` mean | 95% CI | absolute trend mean | P(abs < 0) |
|---|---|---|---|---|---|---|
| Northern Forests | **925** | 95 | −0.127 | [−0.394, +0.078] | **−0.245** | **0.998** |
| NW Forested Mountains | **602** | 161 | +0.049 | [−0.160, +0.274] | −0.069 | 0.868 |
| Eastern Temperate Forests | 78 | 36 | −0.085 | [−0.438, +0.189] | −0.202 | 0.932 |
| North American Deserts | 43 | 62 | +0.075 | [−0.187, +0.412] | −0.043 | 0.656 |
| Water | 6 | 13 | −0.036 | [−0.374, +0.269] | −0.153 | 0.861 |
| Great Plains | **0** | 55 | +0.156 | [−0.083, +0.522] | +0.038 | 0.435 |
| Marine West Coast Forest | **0** | **0** | −0.0003 | [−0.486, +0.479] | −0.118 | 0.753 |
| Mediterranean California | **0** | **0** | −0.0002 | [−0.481, +0.481] | −0.118 | 0.751 |

**Only one region has a credible absolute trend: Northern Forests** (−0.245,
P(<0) = 0.998), which is also the region carrying 925 of the ~1,650 camera
sites. Every other region's `year_region` CI includes zero.

**Two regions have literally zero data** — Marine West Coast Forest and
Mediterranean California, with no cameras and no iNat cells. Their posteriors are
the prior, and their near-identical values (−0.0003 / −0.0002, CIs [−0.49,
+0.48]) confirm it. **These must be greyed out in any map**, not reported as
"no change." This is exactly the data-support gate the simulation work argued
for.

*Source: `moose_v2b_ecoregion_posterior_clean.csv`.*

## 4. Agency comparison — the central result

State-agency directional assessments compiled from published harvest/survey
reports (categorical directions only — an earlier continuous "agency score" was
retired as an undisclosed translation of categorical source data).

**Important:** the on-disk `moose_agency_sign_agreement_table.csv` was built
from **v1fix** posterior values, not v2b. The table below is recomputed against
the v2b posterior and supersedes it.

| region | cameras | model direction | agency | agreement |
|---|---|---|---|---|
| Northern Forests | 925 | **decreasing** (P=0.997) | stable | **disagree** |
| NW Forested Mountains | 602 | uncertain/stable (P=0.860) | stable | agree |
| Eastern Temperate Forests | 78 | uncertain/stable (P=0.930) | stable | agree |
| North American Deserts | 43 | uncertain/stable (P=0.666) | increasing | disagree |
| Great Plains | 0 | uncertain/stable (P=0.426) | stable | agree |
| Water | 6 | uncertain/stable | mixed | not comparable |
| Marine West Coast Forest | 0 | no data (prior only) | stable | no model data |
| Mediterranean California | 0 | no data (prior only) | — | no agency data |

**Tally: 3 agree, 2 disagree, 3 not assessable.**

The comparison hinges on one region. **Northern Forests is the only region where
the model makes a confident claim, and it is the one that disagrees with the
agencies** — the model says decreasing at P=0.997 where the dominant agency
direction across 9 states (87% area coverage) is stable.

That is the finding to take to the team. It is not obviously a model failure:
agency "stable" designations aggregate across states with genuine internal
variation, the model's estimand is occupancy/intensity rather than harvestable
population, and moose declines in the southern edge of the Northern Forests are
documented. But it is a real discrepancy in the one place the model speaks
clearly, and it should not be smoothed over.

The remaining three "agree" cases are weak agreements — the model says
*uncertain* and the agency says *stable*, which is concordance at the level of
"nothing detected," not mutual confirmation.

*Source: `moose_v2b_agency_comparison.csv`, recomputed to
`moose_v2b_sign_agreement_recomputed.csv`.*

## 5. The `snr_derived` discrepancy — resolved, benign

The 8/26 re-extraction reports `snr_derived` with **n_na = 48** where the 8/20
run reported n_na = 0. Checked: this is **not a silently dropped subset.**

`snr = year_beta / year_var`, and `year_var` has mean −0.0022 with a 95% CI of
[−0.197, +0.207] — it straddles zero. So a minority of posterior draws have a
near-zero denominator and the ratio is undefined or explosive. The 48 NAs are a
guard on those draws. The two runs' `year_var` means differ by 0.0046, far
inside the 0.0996 posterior sd, so nothing about the fit changed.

All other parameters reproduce closely between runs (`total_var_beta` −0.1192
vs −0.1162; `year_beta` −0.1170 vs −0.1186). **The re-extraction validates the
8/20 numbers.**

Worth noting: the fact that `snr` is undefined for some draws is the same
underlying condition that makes `trend_robust_indicator` fail — `year_var` is
indistinguishable from zero. Reporting `trend_robust_indicator` (a probability,
always defined) rather than `snr` (a ratio that can blow up) is the right choice.

*Source: `moose_v2b_ns_validate/moose_v2b_national_scalar_posterior_validate.csv`.*

---

## Bottom line

1. **Both moose v2b fits converged cleanly on single unresumed chains** — the
   first moose results not exposed to the checkpoint defect.
2. **National trend ≈ −0.12/yr, P(<0) = 0.994 — but the robustness flag fails**
   (0.255 vs a 0.5 threshold). The camera stream does not independently
   corroborate the decline.
3. **Regional structure exists but only one region is resolvable:** Northern
   Forests, −0.245, P(<0) = 0.998, on 925 camera sites.
4. **That one confident region disagrees with the agencies** (model decreasing,
   agencies stable across 9 states). Three other regions "agree" only in the
   weak sense of both saying nothing.
5. **Two regions have zero data and must be greyed out**, not shown as stable.
6. **The `snr_derived` NA count is benign** and the re-extraction validates the
   earlier numbers.

**Recommended next step:** the disagreement in Northern Forests is worth a
targeted look before it goes to the team — specifically whether the model's
decline is concentrated in the states where agencies *do* report declines, which
would reconcile the two at sub-region scale. The agency data is state-level and
the model is cell-level, so that comparison is makeable.
