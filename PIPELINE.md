# Pipeline of record

Written 2026-09-01, after a search for "which script built the v2b bundles?"
came up empty and cost a round trip. The point of this file is that no future
session should have to reconstruct that answer from commit archaeology.

**Read this first if you are about to build a new input bundle.**

---

## 1. Real-data model pipeline (bobcat / white-tailed deer / moose)

Run in this order. Steps marked **[GAP]** are not currently captured as a
committed script — see Section 3.

| # | Step | Script | Produces |
|---|---|---|---|
| 1 | iNat grid prep, with range mask | `prep_inat_data_grid_v2.R` | masked per-species cell100/cell50 counts |
| 2 | Effort/count matrix build (Fix 1) | `integration_helper_fix1.R` | `make_inat_cell_year_matrix()`, `make_inat_effort_matrix()` with `names_sort=TRUE` |
| 3 | Prep driver — wires 1 and 2 in order | `run_data_prep_v2.R` | `inat_y`, `inat_effort`, `inat_cell_summary` **in memory only** |
| 4 | **[GAP]** Base bundle write | *unknown / not committed* | `HPC/<species>/input_data_<species>.RDS` |
| 5 | Fork into two parameterizations | `00b2_join_ecoregion_real.R` | `input_data_<tok>_ecoregion.RDS`, `input_data_<tok>_national_scalar.RDS` |
| 6 | Covariate reduction (soil fix) | `make_reduced_input_v2.R` | 9-covariate bundle (drops `soil_sand`) |
| 7 | Fit, single-shot MCMC | `HPC/<run>/…_chain<N>.R` | chain output on Hazel |
| 8 | Extract posteriors | `extract_posteriors_batch.R` | per-model posterior CSVs |

**Critical ordering constraint in step 3:** `run_data_prep_v2.R` sources
Arielle's unmodified `integration_helper.R` *first*, then
`integration_helper_fix1.R` (which overwrites two matrix builders), then
`prep_inat_data_grid_v2.R` (which defines the masked prep under a new name so
the original stays callable). Changing that order silently reverts Fix 1.

**What step 5 does and does not do:** `00b2_join_ecoregion_real.R` *forks an
existing base bundle* — it reads `input_data_<species>.RDS`, attaches
`ecoregion_of_cell100`/`nregion`, and writes the two variants. It never touches
`inat_y`/`inat_effort`. It is not the consumer of step 3's in-memory output.
It does carry useful hard assertions: `adj`/`num` must be byte-identical to
`prepped_sim_inputs.RDS` and `ncell100` must match, so a wrong CAR graph fails
loudly rather than silently.

## 2. Range masks — four variants, one adopted

| variant | rule | built by |
|---|---|---|
| IUCN-only | IUCN extant+native polygon | original pipeline |
| presence-only ("Fix 2b") | cell has ≥ `PRESENCE_MIN_RECORDS` iNat records, after isolation/captive filters; IUCN discarded | `prep_inat_data_grid_v2.R` |
| IUCN + expansion | IUCN base, grown by a 4-connected chain of cells clearing a >2-detection threshold | `build_iucn_expanded_mask.R` |
| **union (adopted)** | in-range if IUCN **OR** presence threshold | `build_iucn_expanded_mask.R` |

Union was adopted for all three species (2026-09-01). It retains 100% of
records for every species — recovering moose's disjunct Colorado population,
which the flood-fill expansion structurally could not reach — while keeping the
true-absence cells that presence-only discards.

Per-species grid coverage under the union: moose 15.8%, WTD 86.1%, bobcat 95.8%.
**Note for bobcat:** its IUCN range already covers 94.7% of the grid, so the
union is very nearly no mask at all there. That is a deliberate, documented
choice, not an oversight.

**No fit has yet been run under a union mask** — tractability is inferred from
bundle sizes already exercised (`bobcat_v1fix` built and ran at ncell50 3,145;
union is 3,184), not observed.

## 3. The gap, stated plainly

Step 4 above — the script that takes `run_data_prep_v2.R`'s in-memory
`inat_y`/`inat_effort` and writes the base `input_data_<species>.RDS` — is not
in this repo and has not been located. `run_data_prep_v2.R` says so in its own
closing message: it instructs the user to "wire these into the existing
bundle-forking step … to produce a real `input_data_<species>_v2.RDS`."

Consequences:

- The Aug 18 v2b bundles that every currently-converged fit was built from
  cannot be regenerated from committed code alone.
- Union-mask bundles cannot be built until this step is reconstructed.
- `run_data_prep_v2.R`'s header states "neither fix has been run against real
  data as of this commit" (Aug 10) while the v2b bundles date to ~Aug 18, so
  roughly a week of interactive patching happened that this repo did not
  capture. Two runtime fixes from that period are recorded in the header
  (missing `library()` calls; presence-mask NA rows not being dropped, which
  silently reproduced the full 3,322-cell grid).

**Reconstruction protocol, when it happens:** rebuild a *presence-mask* (v2b)
bundle first and diff it field-for-field against the real Aug 18 `moose_v2b`
bundle. Validate the reconstruction against a bundle whose correct answer is
already known, before using it for anything new.

## 4. Simulation pipeline

Numbered `00*` (prep) → `01*` (run) → `02*`/`01j` (collect). Estimator sweep
arms: `array_occ`, `array_rn`, `camera_occ`, `camera_rn` (rows 1-180, 181-360,
361-540, 541-720 respectively under the 4-arm design).

`isdm_sim_codebase/` is an **August snapshot bundle**, deliberately left out of
version control (see `.gitignore`). Do not treat it as authoritative: it is
stale for at least one file (`01i_camera_rn_sbatch.sh` — the bundle copy is a
broken earlier draft; the repo-root copy is correct). Two files existed *only*
there and have been recovered into the repo with provenance headers:
`01j_collect_estimator_sweep.R` and `HPC_reference/bobcat/integration_helper.R`.

**Known defect, not yet confirmed fixed locally:** the sweep's "null" scenario
zeroed only the spatial deviation amplitude, never `year_beta`/`year_var`, so
the true national trend was present in both scenarios. No false-positive rate
has ever been measured; the earlier "false-positive > power" finding was
retracted. See the header of `01j_collect_estimator_sweep.R`.

## 5. Conventions worth not rediscovering

- **Version control flow.** Changes are drafted and tested locally, committed
  and pushed by the user from GitHub Desktop, and pulled by Claude Code from
  `origin/main`. Files are never hand-delivered to the cluster outside git —
  that is what produced the untracked hand-edited cluster copies now sitting in
  `HPC_reference/`.
- **Hazel sbatch preamble.** `source ~/.bashrc; conda activate <env>` does not
  work: `compileNimble()` shells out to `R CMD SHLIB`, so bare `R` must be on
  `PATH`. Use the `export PATH="$DST/bin:$PATH"` form from a known-working
  wrapper. `--qos=short` caps at 02:00:00; longer jobs need
  `--partition=compute --qos=long`.
- **`*.RDS` is gitignored** — bundles are regenerable pipeline intermediates,
  not source. That is also precisely why step 4 going missing was costly.
