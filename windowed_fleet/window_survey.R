#!/usr/bin/env Rscript
# window_survey.R (2026-10-17) -- per species x window, from the FITTED 18-yr bundle (read-only):
# camera deployments (site-years), detecting deployments, DETECTING YEARS (calendar years with >= 1
# detecting deployment), first/last detecting year, iNat records and cells with records, ncell50.
# Exclusion rule, both versions reported so the cost of each is visible before the list is fixed:
#   strict  (rwkays 10/17 first): >= 5 camera detecting years in the window (win5 then needs all 5 of 5)
#   relaxed (rwkays 10/17, ADOPTED): win5 >= 4 of 5, win10 >= 5 of 10. Not estimable -> not fitted.
H <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/HPC"
W <- list(full = 2008:2025, win10 = 2016:2025, win5 = 2021:2025)
dirs <- sort(list.files(H, "_v2u_national_scalar$"))
sp_of <- sub("_v2u_national_scalar$", "", dirs)
src <- ifelse(file.exists(file.path(H, paste0(sp_of, "_v2un_national_scalar"))), paste0(sp_of, "_v2un_national_scalar"), dirs)
out <- do.call(rbind, lapply(seq_along(src), function(i) {
  b <- readRDS(file.path(H, src[i], paste0("input_data_", src[i], ".RDS")))
  yr <- b$real_data$occ.covs$year; det <- rowSums(b$real_data$y, na.rm = TRUE) > 0
  iy <- as.integer(colnames(b$inat_y))
  do.call(rbind, lapply(names(W), function(w) {
    k <- yr %in% W[[w]]; dy <- sort(unique(yr[k & det])); ic <- iy %in% W[[w]]
    data.frame(species = sp_of[i], source = src[i], window = w, ncell50 = b$constants_list$ncell50,
      cam_deps = sum(k), cam_det_deps = sum(k & det), det_years = length(dy),
      first_det = if (length(dy)) min(dy) else NA, last_det = if (length(dy)) max(dy) else NA,
      cam_years = length(unique(yr[k])), J_last_window = if (any(k)) b$constants_list$J[max(which(k))] else NA,
      inat_records = sum(b$inat_y[, ic], na.rm = TRUE), inat_cells_w_rec = sum(rowSums(b$inat_y[, ic, drop = FALSE], na.rm = TRUE) > 0),
      inat_years_w_rec = sum(colSums(b$inat_y[, ic, drop = FALSE], na.rm = TRUE) > 0))
  }))
}))
out$estimable_strict  <- out$det_years >= 5
out$estimable_relaxed <- out$det_years >= ifelse(out$window == "win5", 4, 5)
out$estimable <- out$estimable_relaxed
write.csv(out, "/home/rwkays/isdm/fleet_v2u/windowed/window_survey.csv", row.names = FALSE)
options(width = 250); print(out, row.names = FALSE)
for (w in c("win10", "win5")) for (r in c("strict", "relaxed")) {
  e <- out[[paste0("estimable_", r)]][out$window == w]
  cat(sprintf("\n%s %-7s: %d of %d estimable; NOT ESTIMABLE: %s\n", w, r, sum(e), length(e),
              paste(out$species[out$window == w][!e], collapse = ", ")))
}
cat("\nfits (relaxed, 2 arms x 12 chains): ", sum(out$estimable_relaxed[out$window != "full"]) * 2, "fits,",
    sum(out$estimable_relaxed[out$window != "full"]) * 24, "chain jobs\n")
