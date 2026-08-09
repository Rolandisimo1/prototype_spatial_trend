#!/usr/bin/env Rscript
# debug_collector.R -- isolate why grain=3 (original run) isn't appearing in
# the grain-sweep aggregate table. Not a deliverable.
PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
PROTO_DIR <- paste0(PROJ, "/prototype_spatial_trend")
setwd(PROTO_DIR)

read_results <- function(dir, grain_value = NA) {
  files <- list.files(dir, pattern = "^rep_.*\\.RDS$", full.names = TRUE)
  cat("n files found in", dir, ":", length(files), "\n")
  if (length(files) == 0) return(NULL)
  rows <- lapply(files, function(f) {
    x <- readRDS(f)
    data.frame(
      rep = x$rep, status = x$status,
      grain = if (!is.null(x$grain)) x$grain else grain_value,
      bias_all = if (!is.null(x$bias_all)) x$bias_all else NA,
      rmse_all = if (!is.null(x$rmse_all)) x$rmse_all else NA,
      coverage_all = if (!is.null(x$coverage_all)) x$coverage_all else NA,
      n_informed = if (!is.null(x$n_informed)) x$n_informed else NA,
      spatial_cor = if (!is.null(x$spatial_cor)) x$spatial_cor else NA_real_
    )
  })
  do.call(rbind, rows)
}

cat("=== original ===\n")
orig <- read_results(paste0(PROTO_DIR, "/sim_results"), grain_value = 3)
cat("class(orig):", class(orig), " nrow:", nrow(orig), "\n")
print(head(orig, 3))
cat("table(status):\n"); print(table(orig$status))
cat("table(grain):\n"); print(table(orig$grain))

cat("\n=== new_grain ===\n")
new_grain <- do.call(rbind, lapply(c(0, 1, 8), function(g) {
  read_results(paste0(PROTO_DIR, "/sim_results_grain/grain_", g), grain_value = g)
}))
cat("nrow new_grain:", nrow(new_grain), "\n")

cat("\n=== combined ===\n")
combined <- rbind(new_grain, orig)
cat("nrow combined (before status filter):", nrow(combined), "\n")
combined_ok <- combined[combined$status == "ok", ]
cat("nrow combined (after status=='ok' filter):", nrow(combined_ok), "\n")
cat("table(grain) in combined_ok:\n"); print(table(combined_ok$grain))
