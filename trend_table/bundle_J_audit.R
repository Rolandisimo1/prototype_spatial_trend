#!/usr/bin/env Rscript
# bundle_J_audit.R -- fisher build failure (dOcc_v dimension on y[i, 1:J[i]]) survived trimming y/yday to max(J).
# Compare fisher's bundle with every fitted bundle on the shape facts NIMBLE checks: J type/range/NA/zero,
# y and yday class/dims/storage mode, all-NA rows, rows whose non-NA count != J. Also answers the zero-night
# question fleet-wide: any J == 0 or all-NA y rows in a fitted bundle?
H <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/HPC"
dirs <- list.dirs(H, recursive = FALSE, full.names = FALSE)
dirs <- dirs[grepl("_v2un?_national_scalar$", dirs)]
out <- list()
for (d in dirs) {
  f <- file.path(H, d, paste0("input_data_", d, ".RDS")); if (!file.exists(f)) next
  b <- readRDS(f); cl <- b$constants_list; y <- b$real_data$y; J <- cl$J
  yd <- if (!is.null(b$real_data$yday)) b$real_data$yday else cl$yday
  nn <- if (is.null(dim(y))) NA else rowSums(!is.na(y))
  out[[d]] <- data.frame(run_dir = d, fitted = file.exists(file.path(H, d, sprintf("chain_%s_1_single.RDS", d))),
    y_class = paste(class(y), collapse = "/"), y_mode = storage.mode(y), y_dim = paste(dim(y), collapse = "x"),
    yday_where = if (!is.null(b$real_data$yday)) "data" else if (!is.null(cl$yday)) "constants" else "none",
    yday_dim = paste(dim(yd), collapse = "x"), J_class = class(J), J_len = length(J), n_sites = if (is.null(dim(y))) NA else nrow(y),
    J_min = min(J, na.rm = TRUE), J_max = max(J, na.rm = TRUE), J_NA = sum(is.na(J)), J_zero = sum(J == 0, na.rm = TRUE),
    J_one = sum(J == 1, na.rm = TRUE), J_nonint = sum(J != round(J), na.rm = TRUE),
    allNA_rows = sum(nn == 0), nonNA_ne_J = sum(nn != J, na.rm = TRUE),
    J_in = paste(intersect(c("constants", "data"), c(if (!is.null(cl$J)) "constants", if (!is.null(b$real_data$J)) "data")), collapse = "+"))
  rm(b); invisible(gc())
}
r <- do.call(rbind, out); options(width = 260); print(r, row.names = FALSE)
write.csv(r, "/home/rwkays/isdm/fleet_v2u/bundle_J_audit.csv", row.names = FALSE)
cat("DONE\n")
