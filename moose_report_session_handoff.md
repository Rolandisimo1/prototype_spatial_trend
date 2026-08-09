# Handoff: moose model results portfolio (continue in new conversation)

## Project context
iSDM temporal-trend project, PI Roland Kays / Pacifici lab, NCSU. NIMBLE
spatially-varying-coefficient occupancy/intensity models fusing camera-trap +
iNaturalist data for North American mammals, on NCSU Hazel HPC. Original
analyst Arielle (ahwaldst) — her files/envs are NEVER modified, always
fork/copy. Claude Science (this agent) drafts specs/figures/analysis; Claude
Code (separate tool, running on the user's laptop with Hazel SSH access)
executes on Hazel — relayed through the user pasting reports back and forth.
Prototype folder (rw host grant): `/Users/rwkays/claude_code/data_integration_arielle/prototype_spatial_trend`.
Analysis env: `geo`. Project `proj_230540b0c90b`.

## Current task
Building a **full results portfolio report for the moose model** — the first
of a template meant to be reused for bobcat, white-tailed deer (WTD), and
eventually the full 35-species fleet. Moose is both a validation species
(agencies track it closely, so results can be checked against real
population data) and the first species with **both** the national-scalar and
ecoregion models fully converged.

## Model background (need to understand before continuing)
- Single latent intensity field λ (`mu`) in each grid cell/year. Cameras see it
  via occupancy (`cloglog(psi) = log(lambda)`, saturates at high abundance);
  iNaturalist sees it via counts (`log(mu) = theta0 + theta1*log(sum(lambda))`,
  doesn't saturate). `theta1` is the camera/iNat coupling strength.
- Trend: `total_var_beta = year_beta + year_var`. `year_beta` is
  camera-anchored (appears in both submodels), `year_var` is iNat-specific
  (only in the iNat pathway). `trend_robust_indicator = step(snr-1)`,
  `snr = year_beta/year_var` — 1 only when year_beta and year_var share sign
  AND |year_beta| >= |year_var|. Indicator < 0.5 means the trend is NOT
  camera-corroborated (true for moose: 0.478 national-scalar, 0.432 ecoregion).
- Ecoregion extension: `year_region[r] ~ dnorm(0, sigma_region)` random effect
  by EPA Level I ecoregion, added ONLY to the iNat pathway. Nests the
  national-scalar model as sigma_region→0.
- Model provenance is fully audited: `model_code_national_scalar.R` is
  byte-identical (token-level, md5-verified) to Arielle's original. Extensions
  add only the documented ecoregion term. See `model_audit_vs_arielle.md`.
- `cell100` (908 cells) carries the CAR spatial fields (link_occ_intercept,
  MWMT_effect, MCMT_effect). `cell50` (382 for moose) carries per-cell-year
  intensity `mu`. Confirmed: `cell100` index = the `grid100.tif` raster value
  directly, no re-indexing needed. Coordinates trustworthy in both
  `moose_v1fix_*_spatial_car_fields.csv` (x/y EPSG:5070 + lon/lat WGS84) and
  `cell50_coordinates.csv`.

## Data-density finding (IMPORTANT — governs what's trustworthy)
Claude Code found: **moose has zero camera trap data 2008-2018** — camera
coverage starts only in 2019. iNat data alone carries the full early period.
Bobcat/WTD are less extreme but still real: ~88-90% of camera effort sits in
2017-2025, near-empty before 2012. This means:
- The "18-year trend" (`total_var_beta`) is NOT uniformly supported across its
  own window — it's disproportionately shaped by the camera-sparse early
  years, especially for moose.
- **User's explicit direction: report 5-year and 10-year trends instead of
  the 18-year figure.** Windows: 5yr = 2021-2025, 10yr = 2016-2025 (same
  calendar years across species; moose's 10yr keeps 3 camera-sparse years
  2016-2018, accepted with a caveat rather than special-cased).
- **NEVER build a change map spanning back to 2008** for moose — this
  directly contradicts the data-density caveat. (I made this mistake once
  this session — built a 2008→2025 % change map before catching that it
  wasn't the map the user actually wanted. Don't repeat it.)

## Windowed refit status (separate parallel thread on Hazel, not yet done)
12 new fits planned: bobcat/moose/WTD x national_scalar/ecoregion x 5yr/10yr.
Same architecture as the current v1fix pipeline (Fix 1 sorting fix + 10-
covariate KEEP list + current range mask) — a data-subsetting change, not a
new methodology. As of this handoff:
- The `monitors2`/`thin2` patch (to directly monitor `mu`/`psi` in the
  windowed fits, avoiding the reconstruction workaround needed for the
  already-converged 18-year moose models) was implemented and **smoke-tested
  successfully on moose** (mechanism proven: checkpoint save/rbind across
  chunk boundaries works correctly).
- **Next step (in progress, not yet reported back):** build the bobcat-5yr
  run directory with the patch applied and retest at bobcat's actual scale
  (27,572 extra columns vs moose's 8,530) before fanning out to build all 12
  run directories and launching. User approved this order.
- The other 11 preps and all 12 run directories/sbatch scripts do not exist
  yet — only the windowed input bundles exist so far.

## The `mu`/`psi` reconstruction (why maps needed extra work)
Neither model directly monitored the deterministic per-cell intensity `mu` or
per-site occupancy `psi` — only top-level stochastic parameters. Claude Code
reconstructed both **exactly** (not approximated) from monitored posterior
draws, thinned to every 25th draw (6,000 of 150,000 pooled draws) for
tractability. `mu` reconstruction uses `theta0, theta1, total_var_beta` +
covariate data (NOT `overdisp_inat` — that's NB dispersion in the observation
layer, doesn't enter mu). `psi` reconstruction uses `link_occ_intercept,
occ_beta, MWMT_effect, MCMT_effect, year_beta`. Confirmed exact year↔index
mapping: t=1→2008, t=6→2013, t=12→2019, t=18→2025 (18-year model; changes for
windowed bundles). Camera sites only exist 2019/2025 for moose — psi is NOT
available for 2008/2013 (zero sites, not approximated as anything).

## Range mask — DO NOT use the IUCN polygon (established bug, multiple times)
`moose_range.geojson` (circumboreal *Alces alces* IUCN range, includes
Europe/Russia — must filter to North America bounds even if used) is
UNRELIABLE as a mask: it excludes real moose activity (documented Colorado
range expansion, see `fig15_moose_range_mask_review.png`,
`moose_cell_counts_unmasked.csv`). **Do not clip any map to this polygon.**
The model's own `mu`/`psi` output already only exists for cells with real
camera/iNat data (382 cells for moose, all conterminous US, no Canada/Alaska)
— use THAT as the natural extent, zoomed to the data's own bounding box, no
external range product needed. A team decision to replace the range mask
entirely (drop IUCN polygon, add captive-record + spatial-isolation filters
instead) is designed but NOT YET APPLIED to any converged fit — see
`team_memo_inat_pipeline_issues.md` from earlier work.

## Agency comparison — needs REBUILDING (flagged, not yet fixed)
An earlier moose-vs-agency comparison figure (`fig21_moose_model_vs_agency.png`,
saved before this session) has THREE known problems, all explicitly flagged
by the user and NOT YET FIXED:
1. No range mask clipping (ironically the opposite problem — needs SOME
   sensible extent, just not the bad IUCN polygon; use the model's own data
   footprint as with the intensity maps).
2. Lower rendering resolution than the earlier agency maps (`fig9`, `fig12`,
   `fig13`) — should reuse that established rendering approach, not rebuild
   from scratch.
3. **The `agency_score` continuous number is an undisclosed data translation**
   — built by mapping categorical agency labels (increasing/stable/
   decreasing/mixed) to numbers (+1/0/-1/0) and area-weighting across states
   per ecoregion. User's explicit decision: **drop this, compare categorically
   only** (region-by-region sign-agreement table against the plain agency
   categories, not a correlation coefficient). This rebuild has NOT been done
   yet — it's real remaining work.

## Report structure established this session (the template)
See `moose_model_full_report.md` (artifact `5602182a-3f02-4fe4-8231-af62495a19b2`)
for the current draft with all sections, some marked `[PENDING HAZEL]`
(now mostly resolved — see completed artifacts below) and the agency section
marked as needing the rebuild above. Structure: (1) Convergence diagnostics,
(2) National trend, (3) Regional trend, (4) Model vs agency comparison
[NEEDS REBUILD], (5) Spatial intensity over time, (6) Climate SVC maps
[NOT YET BUILT], (7) Camera/iNat congruence — temporal [DONE] + spatial
[NOT YET BUILT], (8) Caveats. Ends with reusable-template notes for
bobcat/WTD (e.g.: add `mu`/`psi` to monitor list before launching future
fits, so the reconstruction workaround isn't needed again — already being
done for the windowed refits via `monitors2`).

## Manuscript reference for methodology
Goldstein et al., bioRxiv 2025.01.17.633640 ("Mammal niches are not conserved
over continental scales") — same lab/lineage, same NIMBLE camera+iNat
integration approach. Two borrowed methods:
1. **Model comparison via RJMCMC indicator**, not WAIC (used earlier in the
   simulation-validation phase of this project, not directly relevant to this
   report).
2. **Camera/iNat congruence via `theta1` posterior** — their paper's own
   quote: "we inspected the posterior estimates of theta_1 to characterize
   the frequency at which the two datasets corresponded." Our model's
   `theta1` plays the identical role. `theta1`'s R-hat is already extracted
   (in the `_rhat_by_family.csv` files, iNat family) but its full posterior
   summary (mean/median/CI) for the congruence writeup is not yet pulled into
   a report table — small, cheap addition, same samples matrix as everything
   else already extracted.

## Completed artifacts this session (all under project proj_230540b0c90b)
- `fig_moose_trend_results.png` → version `2732181d-bfe9-4730-a586-404dc7807278`
  — national vs regional trend forest plot, all 3 caveats stated in-figure.
- `fig_moose_temporal_congruence.png` → version `46d23e11-1ec0-4cac-9e48-6f71a92c1fe7`
  — year_beta vs year_var CI comparison + trend_robust_indicator bars, correct
  `step(snr-1)` definition stated in caption.
- `fig_moose_rhat_by_family.png` → version `40784d2d-eba1-4bc0-a82a-db73de01c112`
  — R-hat by parameter family, both models; shows ecoregion's single miss
  (`link_occ_intercept[7]`) is isolated to CAR fields, strengthening the
  accept-as-converged call.
- `fig_moose_intensity_snapshots.png` → version `44810f73-d747-448e-8843-148731102612`
  — 4-panel mu snapshot maps (2008/2013/2019/2025), zoomed to the model's real
  382-cell footprint, NO range mask applied (correct approach).
- `moose_national_scalar_trend_table.csv` → version `8537fae4-ff47-4b60-be06-c6cf0b7d5e88`
- `moose_ecoregion_regional_trend_table.csv` → version `115d23f0-1699-4d8a-b82a-4fe03a8046ac`
- `moose_convergence_summary_partial.csv` → version `7fa60c09-d85e-4319-abf5-b9eaf7077f93`
  (superseded by the fuller R-hat-by-family data, now in
  `moose_v1fix_{national_scalar,ecoregion}_rhat_by_family.csv` on disk)
- `moose_model_full_report.md` → version `5602182a-3f02-4fe4-8231-af62495a19b2`
  (DRAFT, needs updating — several `[PENDING HAZEL]` sections are now
  resolved by data already on disk, and the agency section needs the rebuild
  described above)
- `pending_hazel_extraction_notes.md` → version `2202ba9f-9e08-4a4b-8d2b-d5cd8785d5e6`
  (mostly resolved — theta1 by-family R-hat landed with item 4; only theta1's
  full posterior mean/CI table is still a small open item)

## Raw data files on disk (prototype folder, not yet all turned into figures)
- `moose_v1fix_national_scalar_spatial_car_fields.csv` (908 rows) — CAR fields,
  NOT YET mapped (climate SVC maps, step 6 of report, not started).
- `moose_v1fix_ecoregion_spatial_car_fields.csv` (908 rows) — same, ecoregion
  model, not yet mapped.
- `moose_v1fix_national_scalar_mu_snapshots.csv` (1528 rows = 382 cells x 4 yrs)
  — used for intensity snapshot maps (done) and change map (DROPPED, see
  above) — NOT yet used for spatial congruence.
- `moose_v1fix_national_scalar_psi_snapshots.csv` (410 rows = 179 sites 2019 +
  231 sites 2025) — NOT yet used for anything; needed for spatial congruence
  map (compare against mu in same cells/years).
- `moose_v1fix_national_scalar_rhat_by_family.csv` (5 rows),
  `moose_v1fix_ecoregion_rhat_by_family.csv` (6 rows) — used for the R-hat
  figure; theta1's row is in the iNat family but its posterior mean/CI (not
  just R-hat) hasn't been pulled into a table yet.
- `moose_agency_by_ecoregion.csv`, `deer_moose_trends_master.csv` — needed for
  the agency comparison rebuild.
- `fig9_moose_agency_trend_map.png`, `fig12_deer_combined_trend_map.png`,
  `fig13_moose_combined_trend_map.png` — the higher-resolution rendering
  approach to reuse for the agency comparison rebuild.
- `cell50_coordinates.csv` — reusable cell50→coordinate crosswalk.

## Immediate next steps for the new conversation
1. **Rebuild the agency comparison section** (step 4 of the report): drop
   `agency_score` (categorical-only comparison), apply the model's own data
   footprint as extent (not the IUCN polygon), reuse the `fig9`/`fig12`/`fig13`
   rendering style.
2. **Build climate SVC maps** (step 6): map `MWMT_effect_mean`/
   `MCMT_effect_mean` from the spatial CAR fields CSVs, cell100-level,
   zoomed to the data footprint (join cell100→cell50 via the site/cell
   mapping already used in `psi_snapshots.csv` if needed, or use cell100's
   own lon/lat directly).
3. **Build spatial congruence map** (step 7b): compare `mu` (iNat-side,
   cell50) against `psi` (camera-side, per-site, only 2019/2025) — this needs
   a join strategy since they're on different spatial units (grid cell vs.
   point site); consider aggregating psi to cell50 by mean, or plotting both
   as separate but adjacent panels rather than forcing a single bivariate map.
4. **Pull theta1's full posterior summary** (small addition) for the
   congruence section's primary metric, per the Goldstein et al. precedent.
5. **Update `moose_model_full_report.md`** to reflect all of the above once
   built, removing stale `[PENDING HAZEL]` markers.
6. Continue monitoring the parallel windowed-refit thread on Hazel (bobcat
   5yr checkpoint test at scale, then fan-out to all 12) — report back to the
   user once Claude Code's next update arrives; this is independent of the
   report-building work above and will supply the real 5yr/10yr trend maps
   that should eventually replace/supplement the 18-year figures in this
   report.
7. Once bobcat and WTD converge, repeat this same report structure for them
   using `moose_model_full_report.md` as the template.
