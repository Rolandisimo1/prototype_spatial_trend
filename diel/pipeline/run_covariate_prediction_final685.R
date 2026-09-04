## run_covariate_prediction_final685.R
##
## Covariate-driven predictive-mapping analysis of camera-trap mammal diel activity,
## keyed on the manuscript's own site definition (`final_array`, 685 sites, 13 species).
##
## SUPERSEDES the prior array_id-keyed covariate-prediction track (covlib.R /
## step1_harmonics.R ... step10_resid_moran.R), which used a different, incompatible
## site definition than the manuscript's main results. Only that track's MODELING
## APPROACH is replicated here (unit of analysis = array, no array random effect,
## GAM harmonic coefficients with select=TRUE shrinkage, inverse-variance weighting,
## posterior-simulation intervals, 400 km spatial-block CV against a nationwide-curve
## null AND a position-only competitor, residual Moran's I diagnostic).
##
## Inputs (see README.md for provenance):
##   array_harmonics_final685.parquet   -- per-array harmonic coefficients (already fit)
##   array_covariates_complete_final.csv -- 5 manuscript covariates, complete for 685 sites
##   array_audit_v4_conus_corrected.csv  -- site locations / deployment counts
##   national_grid_covariates_4cov_25km.csv -- 25 km CONUS prediction grid (4 mappable covariates)
##
## Outputs: table_cv_results_full685.csv, table_mapping_decision.csv,
##          table_residual_moran_full685.csv / _mapped5.csv, national_predictions_final685.csv,
##          fig4_cv_skill_final685.png, fig5_residual_moran_final685.png,
##          fig_covariate_maps_array685.png, fig_national_predictions_mappable5.png
##
## Run: Rscript run_covariate_prediction_final685.R
suppressPackageStartupMessages({library(data.table); library(arrow); library(mgcv)})
source("covlib_final685.R")

## ============================================================================
## 0. Paths (edit if running outside the packaged data bundle)
## ============================================================================
HARM_PARQUET <- "array_harmonics_final685.parquet"
COV_CSV      <- "array_covariates_complete_final.csv"
AUDIT_CSV    <- "array_audit_v4_conus_corrected.csv"
GRID_CSV     <- "national_grid_covariates_4cov_25km.csv"   # built by build_national_grid.py

## ============================================================================
## 1. Load and merge
## ============================================================================
H   <- as.data.table(read_parquet(HARM_PARQUET))
AC  <- fread(COV_CSV)
AUD <- fread(AUDIT_CSV)

D <- merge(H, AC, by = "final_array")
D <- merge(D, AUD[, .(final_array, n_deployments, n_projects, extent_km)], by = "final_array")
for (h in HARM) D[[paste0("v_", h)]] <- D[[paste0(h, "_se")]]^2
D[, obs_noct := pct_noct]; D[, obs_crep := pct_crep]; D[, obs_conc := conc]
D[, obs_peak := peak_h];   D[, obs_act  := act]

cat(sprintf("Merged: %d rows, %d arrays, %d species\n",
            nrow(D), uniqueN(D$final_array), uniqueN(D$species)))
stopifnot(uniqueN(D$final_array) == 489)   # arrays WITH a fitted harmonic curve (>= 30 events);
                                            # AC/AUD cover all 685 final_array sites, D is the modeling subset

SE_COL <- c(noct = "pct_noct_se", crep = "pct_crep_se", conc = "conc_se",
            peak = "peak_h_se",   act = "act_se")

## ============================================================================
## 2. 400 km spatial-block cross-validation: covariate model vs null vs position
##    (5 manuscript covariates: pop_1km, ag_5km, tcc_1km, tmax_hottest_month, rug_5km)
## ============================================================================
mk_folds <- function(dd, block_km = 400, nfold = 5, seed = 7) {
  a <- unique(dd[, .(final_array, x_km, y_km)])
  a[, block := paste0("B", floor(x_km/block_km), "_", floor(y_km/block_km))]
  blk <- a[, .(x_km = mean(x_km), y_km = mean(y_km)), by = block]
  set.seed(seed)
  km <- kmeans(as.matrix(blk[, .(x_km, y_km)]), centers = min(nfold, nrow(blk)),
               nstart = 25, iter.max = 100)
  blk[, fold := km$cluster]
  merge(a[, .(final_array, block)], blk[, .(block, fold)], by = "block")
}

run_species_cv <- function(SP, vars, nfold = 5, block_km = 400, include_position = TRUE) {
  dd <- D[species == SP]
  if (nrow(dd) < 20) return(NULL)
  fo <- mk_folds(dd, block_km, nfold); dd <- merge(dd, fo, by = "final_array")
  dv <- apply_trans(dd, vars)
  res <- list()
  for (k in seq_len(uniqueN(dv$fold))) {
    tr <- dv[fold != k]; te <- dv[fold == k]
    if (nrow(te) < 3 || nrow(tr) < 30) next
    mods <- lapply(setNames(HARM, HARM), function(h) try(fit_coef(tr, h, vars, k = 5), silent = TRUE))
    if (any(sapply(mods, inherits, "try-error"))) next
    posm <- NULL
    if (include_position) {
      posm <- lapply(setNames(HARM, HARM), function(h) {
        w <- 1/pmax(tr[[paste0("v_", h)]], 1e-8); w <- w/mean(w)
        t2 <- copy(tr); t2$.w <- w
        try(gam(as.formula(paste(h, "~ s(x_km, y_km, k=30)")), data = t2,
                weights = .w, method = "REML", select = TRUE), silent = TRUE)
      })
      if (any(sapply(posm, inherits, "try-error"))) posm <- NULL
    }
    nullB <- sapply(HARM, function(h) {
      w <- 1/pmax(tr[[paste0("v_", h)]], 1e-8); sum(w*tr[[h]])/sum(w) })
    pc <- predict_curves(mods, te); pm <- pc$M
    nm <- metrics_mat(curves_from_B(matrix(rep(nullB, each = nrow(te)), nrow(te),
                                           dimnames = list(NULL, HARM))))
    sm <- if (!is.null(posm)) predict_curves(posm, te)$M else NULL
    for (mm in MET) {
      ob <- te[[paste0("obs_", mm)]]
      r  <- skill_of(ob, pm[[mm]], nm[[mm]], mm)
      rs <- if (!is.null(sm)) skill_of(ob, sm[[mm]], nm[[mm]], mm) else list(skill = NA, mae = NA)
      res[[length(res)+1]] <- data.table(species = SP, metric = mm, fold = k,
        n_te = nrow(te), n_tr = nrow(tr), mae_cov = r$mae, mae_null = r$mae_null,
        skill_cov = r$skill, mae_pos = rs$mae, skill_pos = rs$skill)
    }
  }
  rbindlist(res)
}

summarise_cv <- function(cv) {
  Ssum <- cv[, {
    ci  <- boot_ci(.SD, "n_te", "mae_cov", "mae_null")
    cip <- boot_ci(.SD, "n_te", "mae_pos", "mae_null")
    .(cov_skill = uw(n_te, mae_cov, mae_null), ci_lo = ci[1], ci_hi = ci[2],
      pos_skill = uw(n_te, mae_pos, mae_null), pos_ci_lo = cip[1], pos_ci_hi = cip[2],
      n_folds = .N, n_te_tot = sum(n_te))
  }, by = .(species, metric)]
  Ssum[, beats_null     := is.finite(ci_lo) & ci_lo > 0]
  Ssum[, beats_position := is.finite(cov_skill) & is.finite(pos_skill) & cov_skill > pos_skill]
  Ssum[, pos_beats_null := is.finite(pos_ci_lo) & pos_ci_lo > 0]
  Ssum
}

set.seed(20260825)
CV_ALL <- rbindlist(lapply(unique(D$species), run_species_cv, vars = COVARS))
S <- summarise_cv(CV_ALL)
setorder(S, species, metric)
fwrite(S, "table_cv_results_full685.csv")
cat(sprintf("5-covariate CV: %d species x metric combos, %d beat the null (95%% CI > 0)\n",
            nrow(S), sum(S$beats_null)))

## ============================================================================
## 3. Mappability test: which winners still beat null+position on the 4 covariates
##    that exist at an unsampled grid cell? (tmax_hottest_month is a per-deployment,
##    survey-window-matched Daymet value -- undefined off-camera; see README.)
## ============================================================================
MAPPABLE_COVARS <- c("pop_1km", "ag_5km", "tcc_1km", "rug_5km")
winners <- S[beats_null == TRUE & beats_position == TRUE, .(species, metric)]

CV_MAP4 <- rbindlist(lapply(unique(winners$species), run_species_cv, vars = MAPPABLE_COVARS))
S4 <- summarise_cv(CV_MAP4)
fwrite(S4, "cv_summary_mappable4_final685.csv")

map_decision <- merge(S4, winners, by = c("species","metric"))
map_decision[, mappable := beats_null & beats_position]
map_decision[, decision := ifelse(mappable, "MAPPED", "not supported for national mapping")]
fwrite(map_decision, "table_mapping_decision.csv")
cand_species <- unique(map_decision[mappable == TRUE, species])
target_metrics <- setNames(map_decision[mappable == TRUE, metric], map_decision[mappable == TRUE, species])
cat(sprintf("Mappable on the 4 grid-available covariates: %d of %d candidate combos\n",
            sum(map_decision$mappable), nrow(map_decision)))
print(map_decision[mappable == TRUE, .(species, metric, cov_skill, ci_lo, pos_skill)])

## ============================================================================
## 4. Full-data final models (all 13 species, 5-covariate) + the mappable-set
##    final models (5 candidates, 4-covariate) used for prediction
## ============================================================================
Dt  <- apply_trans(D, COVARS)
Dt4 <- apply_trans(D, MAPPABLE_COVARS)

fit_final <- function(SP, dat, vars) {
  dd <- dat[species == SP]
  mods <- lapply(setNames(HARM, HARM), function(h) try(fit_coef(dd, h, vars, k = 5), silent = TRUE))
  if (any(sapply(mods, inherits, "try-error"))) return(NULL)
  list(models = mods, dat = dd)
}
FINAL_MODELS <- setNames(lapply(unique(D$species), fit_final, dat = Dt,  vars = COVARS),         unique(D$species))
FINAL4       <- setNames(lapply(cand_species,      fit_final, dat = Dt4, vars = MAPPABLE_COVARS), cand_species)
saveRDS(FINAL_MODELS, "final_models_685.rds")
saveRDS(FINAL4,       "final_models_mappable4.rds")

## ============================================================================
## 5. Residual spatial autocorrelation diagnostic (Moran's I by distance band)
##    on the mappable-set final models, at the metric each species was mapped for
## ============================================================================
BANDS <- list(c(0,50), c(50,200), c(200,800))
moran_rows <- list()
for (SP in cand_species) {
  MM <- FINAL4[[SP]]; dv <- MM$dat; mm <- target_metrics[[SP]]
  pc <- predict_curves(MM$models, dv)
  obscol <- paste0("obs_", mm)
  res <- if (mm == "peak") circ_diff_h(pc$M[[mm]], dv[[obscol]]) else dv[[obscol]] - pc$M[[mm]]
  w <- 1/pmax(dv[[SE_COL[[mm]]]]^2, 1e-6)
  for (bd in BANDS) {
    raw <- moran_band(dv$x_km, dv$y_km, dv[[obscol]], w, bd[1], bd[2])
    rsd <- moran_band(dv$x_km, dv$y_km, res, w, bd[1], bd[2])
    moran_rows[[length(moran_rows)+1]] <- data.table(species = SP, metric = mm,
      band = sprintf("%d-%d km", bd[1], bd[2]), n = raw[3],
      moran_raw = raw[1], p_raw = raw[2], moran_resid = rsd[1], p_resid = rsd[2],
      pct_clustering_removed = 100*(1 - rsd[1]/raw[1]))
  }
}
RM4 <- rbindlist(moran_rows)
fwrite(RM4, "table_residual_moran_mapped5.csv")

## Full 13x5 diagnostic (5-covariate models), for the supplementary record
moran_rows_full <- list()
for (SP in names(FINAL_MODELS)) {
  MM <- FINAL_MODELS[[SP]]; if (is.null(MM)) next
  dv <- MM$dat
  pc <- predict_curves(MM$models, dv)
  for (mm in MET) {
    obscol <- paste0("obs_", mm)
    res <- if (mm == "peak") circ_diff_h(pc$M[[mm]], dv[[obscol]]) else dv[[obscol]] - pc$M[[mm]]
    w <- 1/pmax(dv[[SE_COL[[mm]]]]^2, 1e-6)
    for (bd in BANDS) {
      raw <- moran_band(dv$x_km, dv$y_km, dv[[obscol]], w, bd[1], bd[2])
      rsd <- moran_band(dv$x_km, dv$y_km, res, w, bd[1], bd[2])
      moran_rows_full[[length(moran_rows_full)+1]] <- data.table(species = SP, metric = mm,
        band = sprintf("%d-%d km", bd[1], bd[2]), n = raw[3],
        moran_raw = raw[1], p_raw = raw[2], moran_resid = rsd[1], p_resid = rsd[2],
        pct_clustering_removed = 100*(1 - rsd[1]/raw[1]))
    }
  }
}
fwrite(rbindlist(moran_rows_full), "table_residual_moran_full685.csv")

## ============================================================================
## 6. National prediction (25 km grid), posterior simulation for 95% CIs,
##    restricted to the 5 species x metric combinations that (a) beat the null,
##    (b) beat the position-only competitor, on the 4 grid-available covariates.
##    Extrapolation envelope (MESS + Mahalanobis) flags cells outside the
##    training arrays' covariate range.
## ============================================================================
G <- fread(GRID_CSV)
Gt <- apply_trans(copy(G), MAPPABLE_COVARS)

env_calc <- function(train, newd, vars) {
  Mm <- sapply(vars, function(v) {
    tr <- train[[v]]; tr <- tr[is.finite(tr)]; nv <- newd[[v]]
    rng <- range(tr); n <- length(tr)
    f <- vapply(nv, function(z) 100*sum(tr <= z)/n, numeric(1))
    ifelse(nv < rng[1], 100*(nv - rng[1])/max(diff(rng), 1e-9),
    ifelse(nv > rng[2], 100*(rng[2] - nv)/max(diff(rng), 1e-9), pmin(f, 100 - f)*2))
  })
  Mm <- matrix(Mm, nrow = nrow(newd), dimnames = list(NULL, vars))
  mu <- colMeans(train[, ..vars], na.rm = TRUE)
  Sg <- cov(as.matrix(train[, ..vars]), use = "pairwise.complete.obs")
  Si <- tryCatch(solve(Sg + diag(1e-8, ncol(Sg))), error = function(e) diag(ncol(Sg)))
  md <- function(X) sqrt(pmax(rowSums(((X <- sweep(as.matrix(X), 2, mu, "-")) %*% Si) * X), 0))
  list(mess = apply(Mm, 1, min), maha = md(newd[, ..vars]),
       maha_cut = as.numeric(quantile(md(train[, ..vars]), 0.95)))
}

circfun <- function(v) {                       # circular 95% CI for peak-time draws
  ang <- v/24*2*pi
  mu  <- atan2(mean(sin(ang)), mean(cos(ang)))
  dif <- ((ang - mu + pi) %% (2*pi)) - pi
  lo  <- mu + unname(quantile(dif, 0.025)); hi <- mu + unname(quantile(dif, 0.975))
  c(lo = (lo %% (2*pi))/(2*pi)*24, hi = (hi %% (2*pi))/(2*pi)*24)
}

NSIM <- 400
pred_rows <- list(); env_rows <- list()
for (SP in cand_species) {
  MM <- FINAL4[[SP]]; mm <- target_metrics[[SP]]
  ev <- env_calc(MM$dat, Gt, MAPPABLE_COVARS)
  inside <- ev$mess >= 0 & ev$maha <= ev$maha_cut
  env_rows[[SP]] <- data.table(species = SP, cell25 = G$cell25, lon = G$lon, lat = G$lat,
                                mess = ev$mess, maha = ev$maha, in_envelope = inside)
  pc <- predict_curves(MM$models, Gt, nsim = NSIM, seed = 5)
  sm <- sim_metrics(pc$sims)[[mm]]
  if (mm == "peak") { q <- apply(sm, 1, circfun); ci_lo <- q["lo",]; ci_hi <- q["hi",] }
  else { q <- apply(sm, 1, quantile, c(0.025, 0.975), na.rm = TRUE); ci_lo <- q[1,]; ci_hi <- q[2,] }
  pred_rows[[SP]] <- data.table(species = SP, metric = mm, cell25 = G$cell25,
    lon = G$lon, lat = G$lat, value = pc$M[[mm]], ci_lo = ci_lo, ci_hi = ci_hi)
  cat(sprintf("  predicted %-24s %-6s inside envelope %.1f%%\n", SP, mm, 100*mean(inside)))
}
PRED <- merge(rbindlist(pred_rows), rbindlist(env_rows)[, .(species, cell25, in_envelope, mess, maha)],
              by = c("species","cell25"))
fwrite(PRED, "national_predictions_final685.csv")

cat("\nDONE. Key outputs: table_cv_results_full685.csv, table_mapping_decision.csv,\n",
    "table_residual_moran_mapped5.csv, national_predictions_final685.csv\n")
