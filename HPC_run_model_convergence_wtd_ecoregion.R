library(coda)

species <- "wtd_ecoregion"

setwd(paste0("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/HPC/", species))

project_dir <- getwd()

############################# Assess convergence ###############################
chain_list <- list()

for (chain_id in 1:3) {

  obj <- readRDS(file.path(project_dir,
                           paste0("chain_", species, "_", chain_id, ".RDS")))

  chain_list[[chain_id]] <- mcmc(obj$samples)
}

mcmc_list <- mcmc.list(chain_list)

gd <- gelman.diag(mcmc_list, autoburnin = FALSE, multivariate = FALSE)

# extract PSRF table
rhat <- gd$psrf[,1]

# flag anything problematic
bad <- rhat[rhat > 1.1]

# summary report
cat("\n===== GELMAN SUMMARY =====\n")
cat("Total parameters:", length(rhat), "\n")
cat("Max R-hat:", max(rhat), "\n")
cat("Mean R-hat:", mean(rhat), "\n")
cat("N > 1.1:", length(bad), "\n\n")

# Save to a file
if (length(bad) > 0) {
  cat("Parameters > 1.1:\n")
  print(sort(bad, decreasing = TRUE))
} else {
  cat("No parameters exceed 1.1\n")
}

if (length(bad) > 0) {
  cat("Parameters > 1.1:\n")
  print(sort(bad, decreasing = TRUE))

  bad_df <- data.frame(
    parameter = names(bad),
    rhat = as.numeric(bad)
  )

  write.csv(
    bad_df,
    file = file.path(project_dir, paste0("gelman_bad_rhat_", species, ".csv")),
    row.names = FALSE
  )

} else {
  cat("No parameters exceed 1.1\n")
}
