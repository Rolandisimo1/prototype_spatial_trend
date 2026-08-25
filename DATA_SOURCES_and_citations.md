# Data sources and citations

**Status: INCOMPLETE — this is a working inventory, not a finished reference
list.** Read §5 before using it in a manuscript.

Compiled by inspecting the data files themselves, plus the one method citation
recorded in the project. **No citation list existed in this project before
this file**, and several sources below could not be attributed from the data on
hand. Everything stated here is either read directly from a file (and the file
is named) or flagged as unverified.

---

## 1. Camera trap data

**File:** `raw_cam_data/combined_deployments_all.csv` (26,798 deployment rows),
`raw_cam_data/combined_sequences_all.csv` (2,517,928 sequence-class rows).
**Extent:** 60 `project_id` values, 962 subprojects, 2008–2025.

This is an aggregation of many independent camera trap projects. **40 of the 60
project IDs are human-readable; 20 are opaque numeric database keys covering
15,789 deployments — 59% of all deployments — whose contributing institution
cannot be recovered from these files.** Full inventory of the named projects,
with deployments/years/trap-nights each, is in
`camera_projects_inventory.csv`.

Largest named contributors:

| project | deployments | years | trap-nights |
|---|---|---|---|
| North Carolina's Candid Critters | 2,610 (+192) | 2016–2019 | 52,269 |
| Recreation Effects on Wildlife Project | 1,894 | 2012–2013 | 45,438 |
| Snapshot USA 2020 | 1,454 | 2020 | 50,634 |
| Snapshot USA 2019 | 958 | 2019 | 37,004 |
| Triangle Camera Trap Survey Project | 324 | 2014–2016 | 7,174 |

ID prefixes suggest the aggregator of origin — `NS_` (29 ids, 8,560
deployments), `CA_` (6 ids, 1,368), `Snapshot USA` (2 ids, 2,412), and bare
numeric (23 ids, 14,458). **The expansion of these prefixes is an inference
from their form, not something the files state.** Confirm with Roland before
publishing.

**Citations required but NOT yet resolved (see §5):**

- **Snapshot USA** — the annual national survey; 2019 and 2020 seasons are
  present. Each season has its own data paper in *Ecology*. The specific
  citations for the seasons used must be added.
- **North Carolina's Candid Critters** — a statewide citizen-science camera
  survey with its own publication.
- The **remaining ~35 named projects and all 20 numeric-ID projects** — an
  aggregated multi-project dataset normally requires either a citation to the
  aggregate data paper or an acknowledgement of contributing projects. Roland
  will know which applies.

---

## 2. iNaturalist data

**File:** `raw_cam_data/inat_combo_nam_mams.csv` — 3,014,977 records,
North American mammals, observation dates spanning 1908-11-01 to 2026-01-31
(the model window is 2008–2025).

**Verified from the file:**

- **Quality grade: 3,014,976 of 3,014,977 records are `research` grade**
  (one `needs_id`). So this is effectively a research-grade-only export.
- **Licensing is mixed, and 21% is not openly licensed:**

| license | records | share |
|---|---|---|
| CC-BY-NC | 1,968,231 | 65.3% |
| **none / all rights reserved** | **633,565** | **21.0%** |
| CC-BY | 227,137 | 7.5% |
| CC0 | 84,290 | 2.8% |
| CC-BY-NC-ND | 48,682 | 1.6% |
| CC-BY-NC-SA | 23,079 | 0.8% |
| CC-BY-SA | 18,830 | 0.6% |
| CC-BY-ND | 11,163 | 0.4% |

**Two things to check before publication.** First, the standard iNaturalist
citation depends on how the export was obtained — a direct iNaturalist export
and a GBIF-mediated download are cited differently, and a GBIF download has its
own DOI that must be quoted. **Which route produced this file is not recorded
anywhere in the project** and needs confirming. Second, the 21% all-rights-
reserved fraction is fine for analysis but constrains redistribution — do not
republish those records as a supplementary data file without checking.

---

## 3. Environmental covariates

The model uses **10 ordinary occupancy covariates**, plus temperature entering
separately as a spatially-varying coefficient (its own CAR field, not in the
covariate list):

`Human_pop`, `NDVI_mean`, `Ag`, `Deciduous`, `Evergreen`, `Mixed`,
`terrain_ruggedness`, `soil_clay`, `soil_silt`, `soil_sand`
— plus **MWMT** (mean warmest month temperature) and **MCMT** (mean coldest
month temperature) as SVCs.

This 10-covariate set is the reviewed final decision, recorded in
`make_reduced_input.R`'s header: keep the SVC on MWMT+MCMT, drop the ordinary
temperature covariate and the collinear set (CWD, NDVI_sd, Impervious, PDSI,
PPT, elevation, Temp).

**Source datasets are NOT documented anywhere in this project.** The covariate
rasters were prepared upstream by the original analyst (Arielle / ahwaldst) and
only the extracted values reach the modelling code. `MWMT`/`MCMT` are ClimateNA
variable names by convention, but that is an inference from nomenclature, not a
verified provenance. **Every entry in this section needs its source and version
supplied by Arielle before it can be cited.** This matters beyond bookkeeping:
a raster source drift once silently changed the covariate count from 10 to 17
and broke a model fit, so the exact source and version are load-bearing.

---

## 4. Ancillary spatial data — verified

**EPA / CEC North America Level I Ecoregions.** Used for the regional trend
term. Layer `NA_CEC_Eco_Level1`, field `NA_L1`. Downloaded in-session from:
`https://dmap-prod-oms-edc.s3.us-east-1.amazonaws.com/ORD/Ecoregions/cec_na/na_cec_eco_l1.zip`
(the older `ORD/NA_CEC_Eco_Level1.zip` path now 404s). Produced by the
Commission for Environmental Cooperation; cite CEC, not EPA, as originator.
Cached at `/Users/rwkays/claude_code/geospatialdata`.

**US state boundaries.** US Census Bureau TIGER/Line cartographic boundary
file `cb_2022_us_state_20m` — used for mapping only, not modelling.

---

## 5. Method and software citations

**Verified — recorded in the project:**

- **Goldstein et al.**, continental mammal SVC integrated-SDM preprint,
  bioRxiv **2025.01.17.633640**. Source of the reversible-jump Bernoulli
  inclusion-indicator approach that replaced WAIC as the structure-detection
  gate, and of the FP/FN validation design.
- **eBird double-machine-learning trend paper**, doi **10.1111/2041-210X.14186**
  (code: Zenodo 8092408). **Not adopted** — reviewed and deliberately declined;
  four design lessons were borrowed. See `ebird_dml_lessons.md`.

**Software:** NIMBLE (version **1.4.3** in the local `isdm-sim` environment)
and nimbleEcology. NIMBLE has a standard citation (de Valpine et al., *JCGS*)
plus a version-specific package citation — use `citation("nimble")` on the
machine that produced the results, so the version matches.

**Also cite, if used:** the Royle–Nichols model formulation (Royle & Nichols
2003, *Ecology*) for the RN arm, and the occupancy-model formulation
(MacKenzie et al. 2002) for the occupancy arms.

---

## 6. What is missing — resolve before any manuscript

This list is the point of the document; treat it as a to-do rather than a
formality.

1. **Snapshot USA season citations** for 2019 and 2020 (each season has a
   distinct *Ecology* data paper).
2. **North Carolina's Candid Critters** publication.
3. **Attribution for the 20 numeric-ID projects** — 59% of deployments, no
   recoverable provider from these files. Needs the source database export or
   Roland's knowledge.
4. **How the iNaturalist export was obtained** — direct export vs GBIF
   download; if GBIF, the download DOI.
5. **Covariate raster sources and versions**, from Arielle (§3).
6. **Confirmation of the `NS_` / `CA_` prefix meanings**, which are currently
   an inference from their form.

**Do not paper over any of these with a plausible-looking citation.** Several
sources here are genuinely unrecoverable from the files on hand, and a wrong
citation is worse than an acknowledged gap.
