# Email draft — chunked-resume defect (to Arielle + Krishna)

Krishna, Arielle — following up on the burn-in question. Checking it turned up
something bigger, and I need your read before we go further.

**On the burn-in point itself: you're both right.** Burn-in should happen once.
A checkpointed chain that restores its state is a genuine continuation, and
re-running burn-in on it would be wrong. No disagreement there.

**The problem is that the restore isn't happening.** We checked the actual
chain output rather than reasoning about the code. At every resume boundary,
the first post-resume draw sits *exactly* on the model's hardcoded initial
values — `det_intercept` -> 0.5, `overdisp_inat` -> 0.1, `intercept_tau` ->
0.3 — against pre-resume values like 0.38, 0.51, 0.37. Not close to the
inits; identical, to within 1e-9.

We read all 164 chain files on Hazel (206 GB, 58 fit directories):
**1,234 of 1,234 resume boundaries restart. Zero clean resumes.**

**The mechanism**, in `HPC_run_model_chunks_chain1.R`: `save_chain_state()`
does `saveRDS(Cmcmc$mvSaved, ...)`; `restore_chain_state()` does
`Cmcmc$mvSaved <- obj$mcmc_state`. Our hypothesis — not yet proven — is that
`Cmcmc$mvSaved` is a compiled object backed by an external pointer. `saveRDS`
writes the R wrapper, the pointer is dead on read-back, and the assignment
silently no-ops instead of erroring. The chain then starts the chunk from the
model's inits.

**What it means practically:** a fit recorded as "50,000 iterations,
converged" is really five ~10,000-iteration runs from the same fixed starting
point, concatenated. Effective sample size is overstated and R-hat is
distorted — and R-hat was our convergence gate. Supporting symptoms: runs of
up to 634 identical consecutive draws right after a boundary (samplers
rejecting while they re-climb), and `intercept_tau` jumping 0.37 -> 4,412 at
one boundary.

**What we're proposing — tell us if you'd go a different way.** Rather than
write a new checkpoint mechanism to repair this, we'd drop chunk-resume
entirely: run each chain as a single job, and get additional samples by adding
chains rather than extending them. That removes the broken mechanism instead
of patching it, and needs no new serialization code to validate.

**The question that most affects whether that works, Arielle:** we assume
chunk-resume exists because you hit Hazel's walltime ceiling. Do you know
roughly how many iterations these models actually need to converge? If it's
well beyond what one job can hold, "one job per chain" may not be viable and
we'd need a real checkpoint scheme after all — you'd know that far better than
we would.

**One heads-up rather than a request:** everything in `Converged/` (the 19
species) is affected. It needs re-running regardless — those fits also predate
the iNat effort-matrix sorting fix — so we're not asking anyone to reconstruct
what came out of them. But if any of those numbers went to anyone outside the
project, they should know they'll change.

We've paused all fitting and haven't modified any of your code. Happy to share
the full per-boundary audit table (10,811 rows) if it's useful.
