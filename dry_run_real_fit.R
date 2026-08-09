#!/usr/bin/env Rscript
# dry_run_real_fit.R -- smoke test only. Builds (does NOT compile or run) the
# uncompiled nimbleModel DAG for both real-fit forks against their real
# bundles, to catch constants/data/inits mismatches before committing to a
# 70h/chain sbatch submission. Not a deliverable itself.
suppressPackageStartupMessages({ library(nimble); library(nimbleEcology) })

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"

# source the chain scripts' model_code blocks by extracting just the
# nimbleCode(...) definition, not the whole run_chain_chunk() driver --
# simplest safe way is to re-source the actual chain1.R files up to the
# model_code assignment, but they immediately proceed to run_chain_chunk()
# at the bottom. Instead, pull the nimbleCode block via the standalone
# model_code_*.R files in this directory (identical bodies, audited).
check_one_standalone <- function(token, local_model_file, model_obj_name) {
  dir <- paste0(PROJ, "/HPC/", token)
  setwd(dir)
  source("integration_helper.R")
  source(paste0(PROJ, "/prototype_spatial_trend/", local_model_file))
  model_code <- get(model_obj_name)

  input_data <- readRDS(paste0("input_data_", token, ".RDS"))

  cat("===", token, "===\n")
  t0 <- Sys.time()
  model <- nimbleModel(model_code,
                       constants = input_data$constants_list,
                       data = list(y = input_data$real_data$y,
                                   y_inat = input_data$inat_y,
                                   inat_effort = input_data$inat_effort),
                       inits = input_data$inits_list,
                       calculate = FALSE)
  cat("nimbleModel() built OK in", round(as.numeric(difftime(Sys.time(), t0, units="secs")),1), "sec\n")
  cat("model nodes:", length(model$getNodeNames()), "\n")
  ll <- tryCatch(model$calculate(), error = function(e) paste("ERROR:", conditionMessage(e)))
  cat("initial logProb (calculate()):", ll, "\n\n")
  invisible(model)
}

check_one_standalone("wtd_national_scalar", "model_code_national_scalar.R", "model_code_national_scalar")
check_one_standalone("wtd_ecoregion",       "model_code_ecoregion_trend.R", "model_code_ecoregion")

cat("DRY RUN DONE\n")
