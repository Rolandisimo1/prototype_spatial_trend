## step8_predict.R -- prediction surfaces + extrapolation envelope on the 25 km grid.
## Support: grid covariates are 1 km BUFFER MEANS at the cell centroid, matching the
## camera-side 1 km buffer means (grid_covariates_25km.csv).
## Intervals come from POSTERIOR SIMULATION of the coefficient models (not delta method).
## The recommended variant is the DEFAULT map; other variants are written for comparison.
source("covlib.R")
NSIM <- 400L
REC_VARIANT <- "setA_poponly"     # recommendation justified in the report (sign-stable + best CV)

H  <- fread("array_harmonics_raw.csv"); AC <- fread("array_covariates.csv")
CB <- fread("camera_setb_array.csv")
D  <- merge(H, merge(AC, CB, by = "array_id", all.x = TRUE), by = "array_id")
G  <- fread("grid_covariates_25km.csv")
RMK <- fread(Sys.getenv("RMP"))
SCI <- c("White-tailed Deer"="Odocoileus virginianus","Northern Raccoon"="Procyon lotor",
         "Eastern Gray Squirrel"="Sciurus carolinensis","Coyote"="Canis latrans",
         "American Black Bear"="Ursus americanus")
MODELS <- readRDS("final_models.rds")
CVS <- fread("covariate_model_cv.csv")

## MESS + Mahalanobis + convex-hull-style flag, training arrays vs grid cells
env_calc <- function(train, newd, vars) {
  Mm <- sapply(vars, function(v) {
    tr <- train[[v]]; tr <- tr[is.finite(tr)]; nv <- newd[[v]]
    rng <- range(tr); n <- length(tr)
    f <- vapply(nv, function(z) 100*sum(tr <= z)/n, numeric(1))
    ifelse(nv < rng[1], 100*(nv - rng[1])/max(diff(rng), 1e-9),
    ifelse(nv > rng[2], 100*(rng[2] - nv)/max(diff(rng), 1e-9),
           pmin(f, 100 - f)*2))
  })
  Mm <- matrix(Mm, nrow = nrow(newd), dimnames = list(NULL, vars))
  mu <- colMeans(train[, ..vars], na.rm = TRUE)
  S  <- cov(as.matrix(train[, ..vars]), use = "pairwise.complete.obs")
  Si <- tryCatch(solve(S + diag(1e-8, ncol(S))), error = function(e) diag(ncol(S)))
  md <- function(X) sqrt(pmax(rowSums(((X <- sweep(as.matrix(X), 2, mu, "-")) %*% Si) * X), 0))
  list(mess = apply(Mm, 1, min), maha = md(newd[, ..vars]),
       maha_cut = as.numeric(quantile(md(train[, ..vars]), 0.95)),
       per_cov = Mm)
}

pred_rows <- list(); env_rows <- list(); curve_rows <- list()
for (SP in names(SCI)) {
  keep_cells <- RMK[species == SCI[[SP]], cell_id]
  for (vn in unique(CVS[species == SP, variant])) {
    key <- paste(SP, vn); MM <- MODELS[[key]]
    if (is.null(MM)) next
    vars <- MM$vars
    ## a model containing a camera-derived covariate CANNOT be mapped: heterospecific
    ## detection rate does not exist at an unsampled cell. Skipped for prediction,
    ## reported in CV only.
    if (any(grepl("^rate_", vars))) next
    if (!all(vars %in% names(G))) next
    gg <- G[cell_id %in% keep_cells]
    gg <- gg[complete.cases(gg[, ..vars])]
    if (!nrow(gg)) next
    gv <- apply_trans(gg, vars); tv <- MM$dat

    ev <- env_calc(tv, gv, vars)
    inside <- ev$mess >= 0 & ev$maha <= ev$maha_cut
    env_rows[[length(env_rows)+1]] <- data.table(
      species = SP, variant = vn, cell25 = gg$cell25, lon = gg$lon, lat = gg$lat,
      mess = ev$mess, maha = ev$maha, maha_cut = ev$maha_cut, in_envelope = inside,
      n_cov_outside = rowSums(ev$per_cov < 0))

    if (vn != REC_VARIANT) next            # surfaces only for the recommended variant
    pc <- predict_curves(MM$models, gv, nsim = NSIM, seed = 5)
    SM <- sim_metrics(pc$sims)
    for (mm in MET) {
      q <- apply(SM[[mm]], 1, quantile, c(0.025, 0.5, 0.975), na.rm = TRUE)
      pred_rows[[length(pred_rows)+1]] <- data.table(
        species = SP, variant = vn, metric = mm, cell25 = gg$cell25,
        lon = gg$lon, lat = gg$lat, cell_id = gg$cell_id,
        value = pc$M[[mm]], ci_lo = q[1,], ci_med = q[2,], ci_hi = q[3,],
        in_envelope = inside)
    }
    ## 48-bin predicted curves per cell (quantised to integers per mille of the daily total)
    cur <- pc$curves                        # 48 x ncell
    prop <- sweep(cur, 2, colSums(cur), "/")
    q1000 <- round(prop*1000)
    curve_rows[[length(curve_rows)+1]] <- data.table(
      species = SP, cell25 = rep(gg$cell25, each = NBIN),
      bin = rep(seq_len(NBIN), ncol(q1000)),
      rate_permille = as.integer(as.vector(q1000)),
      in_envelope = rep(inside, each = NBIN))
    cat(sprintf("  [%s] %s: %d cells, inside %.1f%%\n", SP, vn, nrow(gg), 100*mean(inside)))
  }
}
PR <- rbindlist(pred_rows); EN <- rbindlist(env_rows); CU <- rbindlist(curve_rows)
fwrite(PR, "covariate_predictions.csv")
fwrite(EN, "extrapolation_envelope.csv")
fwrite(CU, "covariate_cellcurves.csv")
cat(sprintf("PRED %d  ENV %d  CURVES %d\n", nrow(PR), nrow(EN), nrow(CU)))
print(EN[, .(cells = .N, pct_inside = round(100*mean(in_envelope),1)), by = .(species, variant)][order(species, variant)])
