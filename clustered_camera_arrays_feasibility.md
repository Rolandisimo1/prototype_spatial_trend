# Phase 2 feasibility: Snapshot-USA-style clustered camera arrays

**Status: scoped only, per instructions. Not built.**

## What "clustered" would change

The current model treats every camera site as an independent occupancy draw:
`y[i, 1:J[i]] ~ dOcc_v(probOcc = psi[i], probDetect = p[i, 1:J[i]], len = J[i])`,
with `psi[i]` a function of that site's own covariates (plus the shared
cell100-level CAR intercept/SVCs). Snapshot USA-style deployment instead
places ~10-20 cameras within a small radius ("array"), replicated at many
array locations. Cameras within an array are not independent replicates of
occupancy — they are repeat detectors of essentially the SAME true
occupancy state, and the true number of independent spatial samples is the
**array count**, not the camera count.

This matters specifically for this project because the entire point of the
effort/abundance sweeps is resolving **regional** trend structure, and
regional resolvability is governed by the number of independent spatial
samples within a region, not raw camera count. A region with 3 arrays × 15
cameras (45 cameras, 3 independent samples) is much worse-resolved than 45
independently-scattered sites (45 independent samples) — even though "45
cameras" looks identical in a naive count. If clustered data were fed
through the CURRENT model unchanged (independent-site likelihood, just
relabeling grouped cameras as if they were independent sites), the model
would silently overstate resolvability — a wrong answer, not a conservative
one.

## What real data is available

Checked `input_data_bobcat.RDS` directly: `constants_list` has no
array/cluster/deployment field (the only similarly-named field,
`interaction_group`, is unrelated — it's the Laplace-prior grouping for
`occ_beta`, not a spatial camera grouping). `real_data` has `coords` (site
lat/lon) and the existing `grid.index`/`grid100` cell indexing, nothing
array-specific. **There is no real Snapshot-USA-style array metadata in this
dataset to build from.** Any clustered design would need to be synthetic
(e.g., k-means or fixed-radius clustering of the existing 20,531 real site
coordinates into invented "arrays"), which is a real modeling assumption to
defend, not a faithful reconstruction of an actual deployment.

## What building it would actually take

1. **Model fork** (moderate effort, ~1-2 hrs given established patterns in
   this project): add an array-level `psi_array[a] ~ ...` node and a
   deterministic pass-through `psi[i] <- psi_array[array_of_site[i]]` —
   architecturally identical to the `year_effect[c] <- year_region[ecoregion_of_cell100[c]]`
   pattern already built and validated twice this session (CAR model →
   ecoregion model → switch model). Not a new capability, a known pattern.
2. **A clustering scheme** to assign the real 700-site (or scaled) subsample
   into synthetic arrays — needs a defensible, documented choice of cluster
   radius / cameras-per-array, and should be sensitivity-tested since the
   choice is itself an assumption with no real-data anchor.
3. **A third design axis**: cameras-per-array (e.g., 5/10/15/20) at fixed
   total camera count is the actual severity parameter for the independence
   violation — this is what would need to be swept, not just camera count
   again. Crossed with the existing camera × iNat grid, this multiplies the
   design further.
4. **Compute**: given the main effort sweep is already 320-480 replicates,
   a full clustered-array third axis would plausibly 3-4x the total compute
   again on top of that.

## Recommendation

Feasible in principle, not advisable to attempt now. The blocking issue
isn't the model code (that part is straightforward given this project's
existing patterns) — it's that there's no real array metadata to ground the
clustering assumption in, and the honest version of this study needs a new
design axis (cameras-per-array) that meaningfully adds to an already large
compute commitment. Recommend treating this as a genuinely separate
follow-up, gated on: (a) whether the main effort sweep's results actually
create a reason to care about clustering specifically, and (b) checking with
Arielle/the PI whether real deployment/array metadata exists somewhere
outside this prototype's `input_data_bobcat.RDS` bundle (e.g., in the raw
camera-trap records before they were aggregated into cell100/cell50) before
committing to a synthetic clustering scheme.
