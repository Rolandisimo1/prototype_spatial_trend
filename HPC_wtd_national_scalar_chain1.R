library(nimbleEcology)
library(tidyverse)
library(MCMCvis)
library(parallel)
library(coda)

species <- "wtd_national_scalar"

setwd(paste0("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/HPC/", species))

project_dir <- getwd()

source('pipeline_helper.R')
source('integration_helper.R')

# Forked from Arielle's HPC_run_model_chunks_chain1.R, token-identical model
# body (see model_audit_vs_arielle.md), object renamed model_code_national_scalar.
model_code_national_scalar <- nimbleCode({

  for (g in 1:ncell50) {
    for (t in 1:nyear) {

      #Using original intensity formula from Ben's code, grid cell index only
      if (hasSVC) {
        mu[g,t] <- calcIntensity_SVC(
          intensity_intercept = link_occ_intercept[inat_cell100[g]],
          theta0 = theta0,
          theta1 = theta1,
          MWMT_effect = MWMT_effect[inat_cell100[g]],
          MCMT_effect = MCMT_effect[inat_cell100[g]],
          beta = occ_beta[1:numOccCovars],
          xdat = xdat_inat[inat_cell50_start[g]:inat_cell50_end[g], 1:numOccCovars, t],
          MWMT_dat = MWMT_inat[inat_cell50_start[g]:inat_cell50_end[g], t],
          MCMT_dat = MCMT_inat[inat_cell50_start[g]:inat_cell50_end[g], t],
          year_dat = year_vals[t],
          total_var_beta = total_var_beta
        )

      } else {

        mu[g,t] <- calcIntensity_noSVC(
          intensity_intercept = link_occ_intercept[inat_cell100[g]],
          theta0 = theta0,
          theta1 = theta1,
          beta = occ_beta[1:numOccCovars],
          xdat = xdat_inat[inat_cell50_start[g]:inat_cell50_end[g], 1:numOccCovars, t],
          year_dat = year_vals[t],
          total_var_beta = total_var_beta
        )

      }
    }
  }

  for (t in 1:nyear) {
    for(k in 1:n_cells_year[t]) {

      #iNat input data also indexed by year
      y_inat[inat_cells_by_year[k,t],t] ~ dnbinom(size = 1 / overdisp_inat,
                                                  prob = 1 / (1 + overdisp_inat * inat_effort[inat_cells_by_year[k,t],t] *
                                                                mu[inat_cells_by_year[k,t],t]))

    }
  }

  #Occ detection model does not change
  for (i in 1:nsite) {
    for (j in 1:J[i]) {
      cloglog(p[i, j]) <- link_det_intercept +
        p_beta[1] * yday[i, j] +
        p_beta[2] * yday[i, j]^2 +
        p_beta[3] * canopy_height[i] +
        p_beta[4] * log_roaddist[i]
    }

    #Year trend added to the occupancy submodel
    if (has_SVC) {
      cloglog(psi[i]) <- link_occ_intercept[cell[i]] +
        MWMT[i] * MWMT_effect[cell[i]] +
        MCMT[i] * MCMT_effect[cell[i]] +
        inprod(occ_beta[1:numOccCovars], occ_covars[i, 1:numOccCovars]) +
        year_beta * year_occ[i]

    } else {
      cloglog(psi[i]) <- link_occ_intercept[cell[i]] +
        inprod(occ_beta[1:numOccCovars], occ_covars[i, 1:numOccCovars]) +
        year_beta * year_occ[i]
    }

    y[i, 1:J[i]] ~ dOcc_v(probOcc = psi[i], probDetect = p[i, 1:J[i]], len = J[i])
  }

  for (i in 1:numOccCovars) {

    if (prior_type == "Normal") {
      occ_beta[i] ~ dnorm(0, sd = 5)
    } else if (prior_type == "Laplace") {
      occ_beta[i] ~ ddexp(0, scale = lambda[interaction_group[i]])
    }

  }

  for (i in 1:4) {
    p_beta[i] ~ dnorm(0, sd = 5)
  }

  det_intercept ~ dunif(0, 1)
  link_det_intercept <- cloglog(det_intercept)

  link_occ_intercept[1:ncell100] ~ dcar_normal(
    adj = adj[1:nadj],
    num = num[1:nnum],
    tau = intercept_tau
  )

  intercept_tau ~ dgamma(0.1, 0.1)

  if (has_SVC) {

    MWMT_effect[1:ncell100] ~ dcar_normal(
      adj = adj[1:nadj],
      num = num[1:nnum],
      tau = MWMT_tau
    )

    MCMT_effect[1:ncell100] ~ dcar_normal(
      adj = adj[1:nadj],
      num = num[1:nnum],
      tau = MCMT_tau
    )

    MWMT_tau ~ dgamma(0.1, 0.1)
    MCMT_tau ~ dgamma(0.1, 0.1)

  }

  theta0 ~ dnorm(0, sd = 10)
  theta1 ~ dnorm(1, sd = 1)

  overdisp_inat ~ dgamma(shape = 1, scale = 5)

  # Temporal trend priors
  year_beta ~ dnorm(0, sd = sigma_year_beta)
  year_var ~ dnorm(0, sd = sigma_year_var)
  total_var_beta <- year_beta + year_var

  sigma_year_beta ~ dunif(0, 2)
  sigma_year_var ~ dunif(0, 2)

  ###############
  # Robustness indicator
  ##############

  # Signal - to - noise ratio
  snr <- year_beta / year_var

  # Probabilistic "robustness" indicator based on snr
  # Posterior prob that annual occ trend dominates annual noise
  trend_robust_indicator <- step(snr - 1) # 1 if >1, 0 if <=1

})

input_data <- readRDS(paste0("input_data_", species, ".RDS"))

#Set the spatial covariate names
annual_covars <- c("CWD","Human_pop","NDVI_mean","NDVI_sd","Ag","Deciduous",
                   "Evergreen","Impervious","Mixed","PDSI","PPT","Temp")

static_covars <- c("terrain_ruggedness",
                   "soil_clay",
                   "soil_silt",
                   "soil_sand",
                   "elevation")

#Set the time-varying spatial covariates
temp_covars <- c("MWMT", "MCMT")

############################ Functions #########################################
save_chain_state <- function(file, Cmcmc, samples_matrix, iter_total) {

  chain_list <- list(mcmc_state = Cmcmc$mvSaved,
                     samples    = samples_matrix,
                     iter_total = iter_total)

  saveRDS(chain_list, file = file)
}

restore_chain_state <- function(file, Cmcmc) {

  obj <- readRDS(file)

  Cmcmc$mvSaved <- obj$mcmc_state

  return(obj)
}

run_chain_chunk <- function(chain_id,
                            input_data,
                            chunk_iter = 10000,
                            burnin_once = 5000) {

  chain_file <- file.path(project_dir,
                          paste0("chain_", species, "_", chain_id, ".RDS"))

  cat("Building model for chain", chain_id, "\n")

  # Must rebuild model each time

  model <- nimbleModel(model_code_national_scalar,
                       constants = input_data$constants_list,
                       data = list(y = input_data$real_data$y,
                                   y_inat = input_data$inat_y,
                                   inat_effort = input_data$inat_effort),
                       inits = input_data$inits_list,
                       calculate = FALSE)

  Cmodel <- compileNimble(model)

  conf <- configureMCMC(model, enableWAIC = TRUE)

  conf$addMonitors("occ_beta",
                   "link_occ_intercept",
                   "year_beta",
                   "year_var",
                   "total_var_beta",
                   "MWMT_effect",
                   "MCMT_effect",
                   "trend_robust_indicator")

  mcmc  <- buildMCMC(conf)
  Cmcmc <- compileNimble(mcmc, project = model)

  # If there is a chain file - restart, otherwise run for first time
  if (file.exists(chain_file)) {

    cat("Restoring previous state...\n")

    obj <- restore_chain_state(chain_file, Cmcmc)

    samples_so_far <- obj$samples
    iter_total     <- obj$iter_total

  } else {

    cat("Starting new chain...\n")

    set.seed(1000 + chain_id)

    # Burn-in once only
    Cmcmc$run(burnin_once)

    samples_so_far <- NULL
    iter_total     <- 0
  }

  # Run next chunk

  cat("Running", chunk_iter, "iterations...\n")

  Cmcmc$run(chunk_iter)

  new_samples <- as.matrix(Cmcmc$mvSamples)

  if (is.null(samples_so_far)) {
    combined_samples <- new_samples
  } else {
    combined_samples <- rbind(samples_so_far, new_samples)
  }

  iter_total <- iter_total + chunk_iter

  # Save updated state
  save_chain_state(chain_file,
                   Cmcmc,
                   combined_samples,
                   iter_total)

  cat("Chain", chain_id, "now at", iter_total, "iterations\n")
}

############################## Run the chain ###################################

run_chain_chunk(chain_id = 1,
                input_data = input_data,
                chunk_iter = 10000,
                burnin_once = 5000)
