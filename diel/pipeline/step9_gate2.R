## harmonic-parameterisation fidelity: metrics from each array's OWN fitted harmonics
## vs metrics computed directly from the raw binned rates. Separates "is the 5-coefficient
## curve a faithful summary" (this) from "do covariates predict those coefficients" (sanity_gate.csv).
source("covlib.R")
H <- fread("array_harmonics_raw.csv")
rows <- list()
for (SP in unique(H$species)) {
  h <- H[species == SP]
  for (mm in c("noct","crep","conc","peak","act")) {
    o <- h[[paste0("raw_", mm)]]; f <- h[[paste0("obs_", mm)]]
    if (is.null(o)) next
    ok <- is.finite(o) & is.finite(f)
    err <- if (mm == "peak") abs(circ_diff_h(f[ok], o[ok])) else abs(f[ok] - o[ok])
    rows[[length(rows)+1]] <- data.table(species = SP, metric = mm, n = sum(ok),
      r_harmonic_vs_raw = cor(o[ok], f[ok]), mae_harmonic_vs_raw = mean(err))
  }
}
G <- rbindlist(rows); fwrite(G, "harmonic_fidelity.csv")
print(dcast(G, species ~ metric, value.var = "r_harmonic_vs_raw"))
print(dcast(G, species ~ metric, value.var = "mae_harmonic_vs_raw"))
