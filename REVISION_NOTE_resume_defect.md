# Revision note: chunked MCMC resume does not restore chain state

**Date:** 2026-08-23
**Status:** defect confirmed; no fix validated or applied
**Affects:** every multi-chunk model fit produced by this pipeline to date

## Summary

The chunked "pause and resume" mechanism in
`HPC_run_model_chunks_chain{1,2,3}.R` does not restore MCMC chain state. Every
resumed chunk silently restarts from the model's hardcoded initial values.
Fits recorded as long single chains are in fact several short chains, each
started from the same point, concatenated.

All model results produced by this pipeline should be treated as
provisional pending a re-fit.

## The mechanism

The original save/restore pair:

```r
save_chain_state <- function(file, Cmcmc, samples_matrix, iter_total) {
  chain_list <- list(mcmc_state = Cmcmc$mvSaved, ...)
  saveRDS(chain_list, file = file)
}

restore_chain_state <- function(file, Cmcmc) {
  obj <- readRDS(file)
  Cmcmc$mvSaved <- obj$mcmc_state
  return(obj)
}
```

`Cmcmc$mvSaved` is a compiled `CmodelValues` object backed by an external
pointer. `saveRDS` writes the R-level shell; the pointer does not survive
deserialization (touching the restored object raises `Sextptr is not a valid
external pointer`). The assignment therefore transfers no values — and does
not error — so sampling in the new chunk begins from the model's inits.

This is a NIMBLE-specific consequence of compiled objects not being
serializable. The NIMBLE developers have confirmed on the nimble-users list
that model building and compilation must be redone in each R session, and
distinguish saving model variables from saving internal sampler parameters —
consistent with this diagnosis.

## Evidence

All 164 `chain_*.RDS` files on Hazel were read (206 GB, 58 fit directories).
Boundaries occur at each multiple of `chunk_iter = 10000`, verified against
every chunk script.

| | fits | chains | boundaries |
|---|---|---|---|
| Restart at all boundaries | 52 | 154 | 1,234 |
| Clean resume at any boundary | 0 | 0 | 0 |
| Never resumed (single chunk, N/A) | 6 | 12 | — |
| Unreadable / indeterminate | 0 | 0 | — |

Of 1,234 boundaries, 1,132 were auto-confirmed against the fit's own inits
bundle; the remaining 102 were confirmed by inspection against the
hardcoded init constants (their directories carry no inits bundle, so the
automated check returned blank rather than false).

The signature is the post-resume draw landing exactly (< 1e-9) on an init:

```
Converged/coyote             det_intercept  pre=0.37963007  post=0.5  init=0.5
                             intercept_tau  pre=0.35114758  post=0.3  init=0.3
                             overdisp_inat  pre=0.50513761  post=0.1  init=0.1
Converged/white-tailed_deer  det_intercept  pre=0.82166334  post=0.5  init=0.5
moose_v1fix_national_scalar  det_intercept  pre=0.28489817  post=0.5  init=0.5
bobcat_v1fix_ecoregion       overdisp_inat  pre=0.72564112  post=0.1  init=0.1
```

Supporting symptoms: runs of up to 634 identical consecutive draws
immediately after a boundary (samplers rejecting while the chain re-climbs);
pre-to-post jumps with a median of 3.2 within-round SD, p90 11.9, max 282.8;
`intercept_tau` collapsing to fresh-start values (e.g. 0.37 -> 4411.97).

Affected groups: the 19-species `Converged/` fleet, `Poor_mixing/`
(grey_fox, grey_wolf), older single-model fits (american_black_bear, bobcat,
brown_bear, mule_deer, north_american_porcupine), the v0 trend prototype
fits, all six `v1fix` fits, and the `v2b` fits.

The 12 windowed (5yr/10yr) bundles were built but never launched, so they
carry no affected output — they would be affected on launch.

## Why this was not caught by convergence checks

Convergence was gated on Gelman-Rubin R-hat across chains. R-hat asks whether
independently-started chains agree. It was being computed on chains that had
each been chopped up and restarted internally, which is a different question.
A restart transient shared across all three chains can *deflate* R-hat,
making convergence appear better than it is.

Consequently no fit in this project currently carries a trustworthy
convergence assessment.

## Consequence for run length

Because `set.seed()` and the once-only burn-in both sit in the fresh-start
branch, each resumed chunk is an *independently seeded* 10,000-iteration run
from identical initial values. R-hat computed across chunks within a chain is
therefore available at no compute cost, and answers a question the project
could not otherwise answer: how long do these models actually need?

Across 12 models / 36 chains:

| max cross-chunk R-hat | chains |
|---|---|
| > 1.05 | 35 / 36 |
| > 1.1 | 29 / 36 |
| > 1.2 | 19 / 36 |
| > 1.5 | 9 / 36 |

On the headline trend parameter `total_var_beta`: 14/36 chains above 1.1,
median 1.077, max 7.718.

This test is one-directional. These runs share initial values, which makes
agreement *easier* than for genuinely independent chains — so disagreement
proves 10,000 iterations is insufficient, while agreement would have proven
nothing. The true requirement is at least this large. Two extreme values
(12.5, 9.8) come from already-known runaway chains and are inflated by
divergence rather than chunk length; excluding them the bulk still sits at
1.1-1.3.

## What has and has not been done

**Done:** the audit above; all 37 queued and running jobs cancelled; a fork
point committed (orphan branch `arielle-original`, tag `arielle-fork-point`)
recording the inherited code exactly as received.

**Drafted, NOT validated, NOT applied:** `checkpoint_numeric_values.R`
replaces the save/restore pair with one that persists plain numeric node
values, writes them back into the rebuilt model, and *verifies* continuity
instead of assuming it. Its header states the validation plan. No production
fit uses it.

**Explicitly retracted:** an earlier diagnosis in
`run_chain_chunk_fix_resume_burnin.R` attributed the transient solely to lost
adaptive-sampler tuning and asserted that the checkpoint/restore mechanism
itself was sound. That assertion was wrong. It rested on a probe that read
its verdict off a conjugate-sampled node, which draws from its full
conditional in one step from any starting value and therefore cannot
distinguish a restored state from a fresh start. Both that file's header and
Issue 4 of `team_memo_inat_pipeline_issues.md` have been corrected in place
rather than deleted, with this explanation.

Lost sampler adaptation is a real and separate effect, and it remains
unaddressed: NIMBLE provides no supported way to serialize adaptive sampler
state. It is much smaller than restarting from inits.

## Open question

Model building (`nimbleModel()` graph construction) costs ~8 h before
sampling begins, and is re-paid on every chunk. Two paths:

1. If build cost can be reduced substantially (e.g. by vectorizing the
   `mu[g,t]` / `y_inat` loops into nimbleFunctions, collapsing ~56,610
   scalar deterministic nodes to one per year), run each chain as a single
   unresumed job and retire chunking entirely — the defect then cannot recur.
2. If it cannot, keep chunking but with a validated checkpoint, and accept
   long runs on the long queue.

A scaling measurement is queued to decide this: `nimbleModel()` timed alone
at nyear = 5, 10, 18 with `nsite` held fixed, so build time can be attributed
to node count without the site-count confound present in existing timings.
Note that published NIMBLE efficiency gains of 25x-100x refer to MCMC
*sampling* efficiency, not graph-construction time; no literature figure
applies to this bottleneck.

A caveat for path 1: vectorizing `mu` coarsens the dependency graph, so
updating one `link_occ_intercept[k]` would dirty an entire year-column of
`mu` rather than only the cells mapping to that `cell100`. That reduces build
time and may increase per-iteration sampling cost.
