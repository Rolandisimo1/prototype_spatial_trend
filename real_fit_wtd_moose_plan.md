# Real-data ecoregion fit — white-tailed deer + moose

FIRST real-data refit of the project. Crosses the staging gate held through
all simulation work. Chosen because both are closely tracked by state/provincial
agencies -> external truth to validate recovered regional trends against.

## What already exists (prototype)
- model_code_national_scalar.R = Arielle's model, token-identical (audited).
- model_code_ecoregion_trend.R = + year_region random effect (iNat pathway).
- 00b_prep_ecoregion.R = downloads EPA/CEC NA_CEC_Eco_Level1, assigns
  ecoregion_of_cell100 (login node; EPA S3 source, gaftp 404s).
- Arielle real-fit template: HPC_run_model_chunks_chain{1,2,3}.R reads
  input_data_<species>.RDS (constants_list, real_data$y, inat_y, inat_effort,
  inits_list) -> nimbleModel -> configureMCMC(enableWAIC) -> 3 chains.

## The GAP (must be built, name it to Claude Code)
There is NO real-data fit runner in the prototype -- only simulation drivers.
A real ecoregion fit needs:
1. input_data_white-tailed_deer.RDS CONFIRMED on Hazel (prep-only exists per
   prior notes). input_data_moose.RDS -- STATUS UNKNOWN, may need prep first.
2. ecoregion_of_cell100 (FULL length-ncell100) + nregion appended to the real
   bundle's constants_list, joined on the bundle's OWN cell100 identity (NOT
   the sim's reduced 210-cell set). This is the real integration point of
   00b into a real bundle.
3. A fit runner = fork of Arielle's HPC_run_model_chunks_chain*.R that swaps
   model_code -> model_code_ecoregion_trend, adds year_region/sigma_region to
   monitors, keeps everything else (WAIC, 3 chains, checkpointed chunks to
   R-hat<1.1). Plus the national_scalar baseline fit for comparison.

## Slurm reality (Arielle's real fits)
70h walltime/chain, ~110GB RAM, 3 chains, 5000 burnin + 10000-iter
checkpointed chunks until Gelman-Rubin R-hat < 1.1. This is NOT the ~40-min
reduced sim fit -- it is the full-scale model. Budget accordingly.

## Validation (the point)
Compare recovered per-ecoregion trend (year_region posterior) against
agency-known WTD and moose trends. Gate display on the data-support /
posterior-width diagnostic (gray out under-supported regions). Report
national scalar trend + per-region deviations + which regions are
data-supported.

## Guardrails
- Never modify Arielle's originals/envs; fork/copy (new tokens, e.g.
  test_wtd_ecoregion). Use nimble_env.
- Confirm bundles exist BEFORE building the runner; if moose bundle absent,
  stop and report -- prep is a separate step.
- Fit BOTH national_scalar and ecoregion_trend per species for comparison.
- Stop after fits + convergence + validation table; no fleet launch.
