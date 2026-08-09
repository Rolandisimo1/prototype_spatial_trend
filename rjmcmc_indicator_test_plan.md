# Replacing the WAIC gate: RJMCMC indicator test for the ecoregion trend

## Why (from the corrected sweep)
The reseeded sweep showed WAIC cannot distinguish the null scenario (no
regional trend) from the varying scenario (real regional trend) — the two
WAIC distributions overlap at every abundance level (p = 0.80 / 0.91 / 0.61).
So WAIC is unusable as the "does this species need the ecoregion term" gate.

## The method we're adopting (Goldstein et al., bioRxiv 2025.01.17.633640)
Same lab lineage, same camera+iNat iSDM, same NIMBLE stack. They did NOT use
an information criterion. Instead they put model comparison INSIDE one model
as Bayesian variable selection:
- each spatially-varying component carries a Bernoulli "switch" indicator gamma
- gamma = 1 turns the component on, gamma = 0 turns it off
- fit with NIMBLE's reversible-jump MCMC (RJMCMC) sampler on the indicators
- "support" = posterior mean of gamma; thresholds 0.5 (moderate) and 0.9 (strong)

Validation was a simulation: simulate data with a spatial effect present or
absent, fit, and read the FALSE-POSITIVE and FALSE-NEGATIVE rates of the
indicator at the 0.5 and 0.9 thresholds. Published benchmark (single effect):
FP 9.1% / FN 14.9% at 0.5; FP 0.7% / FN 38.1% at 0.9.

## Our test (single indicator, minimal — chosen scope)
Add ONE Bernoulli indicator on the ecoregion trend deviation:
- gamma = 1  -> year_region[r] ~ dnorm(0, sigma_region)  (spatially-varying trend)
- gamma = 0  -> year_region[r] = 0 for all r  (national scalar trend only)
Fit this single "switch" model to the SAME 180 simulated datasets already
generated (null + varying x 3 abundance x 30 reps). Then:
- null scenario:   P(gamma=1) over 30 reps = false-positive rate
- varying scenario: 1 - P(gamma=1) over 30 reps = false-negative rate
Report FP/FN by abundance level at thresholds 0.5 and 0.9 — a direct,
interpretable replacement for the WAIC column, comparable to the 9.1%/14.9%
benchmark.

## Two gotchas carried from their Methods
1. Convergence CANNOT use Gelman-Rubin on gamma (Bernoulli variables sit at
   0 or 1 for long stretches even when mixing well). Monitor the log-prob of
   the observed data (y and y_inat) and compute Gelman-Rubin + ESS on THAT.
2. RJMCMC needs the indicator sampler configured explicitly in NIMBLE with a
   prior inclusion probability pi (start pi = 0.5). This is a real but
   well-trodden setup.

## Deliverable
`indicator_test_summary.csv` (FP/FN by abundance x threshold) + a figure
showing FP and FN vs abundance at both thresholds, with bobcat's rung marked
and the paper's benchmark drawn as a reference — the honest "can we even tell
if a species needs this term" answer WAIC couldn't give.
