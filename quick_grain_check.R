#!/usr/bin/env Rscript
# quick_grain_check.R -- cheap, nimble-free sanity check of the new grain
# truth-field generator and design_df row-selection logic, before spending a
# real MCMC smoke test. Not a deliverable.
PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
setwd(PROTO_DIR)
source("sim_helpers.R")

prepped <- readRDS("prepped_sim_inputs.RDS")
cl <- prepped$constants_list

cat("=== truth field checks across grain levels ===\n")
ref <- make_true_year_effect(prepped$cell100_geo, cl$adj, cl$num, amplitude = 0.3)
cat("reference (continental, amplitude=0.3) sd:", round(sd(ref), 4), "\n\n")

for (g in c(0, 1, 3, 8)) {
  ye <- make_true_year_effect_grain(prepped$cell100_geo, cl$adj, cl$num, n_smooth_iter = g)
  cat("grain=", g, ": length=", length(ye), " mean=", round(mean(ye), 6),
      " sd=", round(sd(ye), 4), " range=", paste(round(range(ye), 3), collapse=" to "),
      " any NA/Inf:", any(!is.finite(ye)), "\n")
}

cat("\n=== design_df row-selection check ===\n")
GRAIN_LEVELS <- c(0, 1, 8)
N_REPS <- 15
grain_design_df <- expand.grid(grain = GRAIN_LEVELS, rep = seq_len(N_REPS))
grain_design_df <- grain_design_df[order(grain_design_df$grain, grain_design_df$rep), ]
grain_design_df$row <- seq_len(nrow(grain_design_df))
cat("nrow:", nrow(grain_design_df), "(expect 45)\n")
cat("row 1:", grain_design_df$grain[1], grain_design_df$rep[1], "(expect grain=0 rep=1)\n")
cat("row 15:", grain_design_df$grain[15], grain_design_df$rep[15], "(expect grain=0 rep=15)\n")
cat("row 16:", grain_design_df$grain[16], grain_design_df$rep[16], "(expect grain=1 rep=1)\n")
cat("row 45:", grain_design_df$grain[45], grain_design_df$rep[45], "(expect grain=8 rep=15)\n")

cat("\n=== compute_grain_metrics dry check (fake samples) ===\n")
ncell <- cl$ncell100
fake_samples <- matrix(rnorm(100 * ncell, 0, 0.1), nrow = 100)
colnames(fake_samples) <- paste0("year_effect[", 1:ncell, "]")
ye_true <- make_true_year_effect_grain(prepped$cell100_geo, cl$adj, cl$num, n_smooth_iter = 3)
informed <- sample(1:ncell, 178)
m <- compute_grain_metrics(fake_samples, ye_true, prepped$cell100_geo, informed)
print(m)
cat("\nQUICK CHECK DONE\n")
