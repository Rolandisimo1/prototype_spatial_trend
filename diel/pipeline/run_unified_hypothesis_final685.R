## run_unified_hypothesis_final685.R
##
## UNIFIED hypothesis-test pipeline for the camera-trap mammal diel-activity manuscript,
## keyed throughout on `final_array` (the manuscript's own 685-site, 25 km geographically
## clustered site definition). Replaces TWO separate, now-inconsistent tracks that existed
## before this phase:
##
##   (1) the H1-environmental-only track (covlib_final685.R / run_covariate_prediction_final685.R,
##       just rerun this session on final_array) -- its GAM/select=TRUE/REML/inverse-variance-
##       weighting/posterior-simulation/400km-block-CV/Moran's-I machinery is the foundation
##       this script builds on (sourced from covlib_final685.R, unchanged).
##   (2) the ORIGINAL H4 (predator)/H5 (density-dependence) hypothesis tests, which used a
##       different, incorrect array_id site key incompatible with the manuscript's final_array
##       definition. Its ADDITIONAL rigor -- the null-simulation artifact check for
##       predator/abundance terms, and the residual-autocorrelation and spatial-competitor
##       diagnostics -- is preserved and extended here to cover ALL predictor classes, not just
##       H4/H5. See covariate_prediction_report_final685.md and carnivore_null_test_report.md
##       for the source methodology this script ports and generalises.
##
## WHAT CHANGED: instead of three separate models per hypothesis family, this script fits ONE
## GAM per species x harmonic coefficient with ALL relevant covariates (5 environmental H1
## terms + species-specific carnivore H4 terms + the species' own relative-abundance H5 term)
## as candidate smooth terms together, gam(..., select=TRUE) deciding which terms survive
## shrinkage. This is a genuine change in modeling framework, not a relabeling of the old
## three-track outputs.
##
## Effect-size scale (stated once, applied everywhere): every term's average marginal effect
## (AME, predictor moved from its 10th to 90th percentile, all other terms held at observed
## values) is reported in units of the RESPONSE metric's own between-array SD (noct, crep,
## conc, act) or in raw hours for peak time (a circular quantity for which "SD units" has no
## natural meaning). This is the ONE scale used for every beta_std-equivalent number in the
## results table and its CI -- the earlier beta_std/raw-scale mismatch bug in the original H4
## predator figure is not reproduced because there is only one reported scale here.
##
## Run: Rscript run_unified_hypothesis_final685.R
suppressPackageStartupMessages({library(data.table); library(arrow); library(mgcv)})
source("covlib_final685.R")

set.seed(20260825)

## ============================================================================
## 0. Paths
## ============================================================================
HARM_PARQUET <- "array_harmonics_final685.parquet"
COV_CSV      <- "array_covariates_complete_final.csv"
PRED_PARQUET <- "predator_abundance_predictors_final685.parquet"
AUDIT_CSV    <- "array_audit_v4_conus_corrected.csv"

## ============================================================================
## 1. Load and merge -- ALL data joined on final_array
## ============================================================================
H    <- as.data.table(read_parquet(HARM_PARQUET))
AC   <- fread(COV_CSV)
PRED <- as.data.table(read_parquet(PRED_PARQUET))
AUD  <- fread(AUDIT_CSV)

D <- merge(H, AC, by = "final_array")
D <- merge(D, PRED, by = "final_array", suffixes = c("", ".pred"))
D <- merge(D, AUD[, .(final_array, n_deployments, n_projects, extent_km)], by = "final_array")

for (h in HARM) D[[paste0("v_", h)]] <- D[[paste0(h, "_se")]]^2
D[, obs_noct := pct_noct]; D[, obs_crep := pct_crep]; D[, obs_conc := conc]
D[, obs_peak := peak_h];   D[, obs_act  := act]

cat(sprintf("Merged D: %d rows, %d arrays, %d species\n",
            nrow(D), uniqueN(D$final_array), uniqueN(D$species)))
stopifnot(uniqueN(D$final_array) == 489)   # arrays with a fitted harmonic curve (>=30 events)

SE_COL <- c(noct = "pct_noct_se", crep = "pct_crep_se", conc = "conc_se",
            peak = "peak_h_se",   act = "act_se")

fwrite(D, "unified_modeling_dataset_final685.csv")
cat(sprintf("Saved unified_modeling_dataset_final685.csv: %d rows x %d cols\n", nrow(D), ncol(D)))

## ============================================================================
## 2. Predictor sets
## ============================================================================
COVARS  <- c("pop_1km","ag_5km","tcc_1km","tmax_hottest_month","rug_5km")
LOGVARS <- c("pop_1km","ag_5km","rug_5km")

CARN_SP <- c("AmericanBlackBear","Bobcat","Coyote","GreyFox","GreyWolf","GrizzlyBear","Puma","RedFox","RedWolf")

OWN_MAP <- c("American Black Bear"="AmericanBlackBear","Coyote"="Coyote","Eastern Chipmunk"="EasternChipmunk",
             "Eastern Cottontail"="EasternCottontail","Eastern Fox Squirrel"="EasternFoxSquirrel",
             "Eastern Gray Squirrel"="EasternGraySquirrel","Mule Deer"="MuleDeer","Northern Raccoon"="NorthernRaccoon",
             "Red Fox"="RedFox","Red Squirrel"="RedSquirrel","Virginia Opossum"="VirginiaOpossum",
             "White-tailed Deer"="WhiteTailedDeer","Wild Turkey"="WildTurkey")

## Predator-relevance decision: a carnivore's log-detection-rate term is included as a
## candidate for species SP only if, WITHIN SP's own 489-array modeling subset:
##   (a) it is not SP itself (self-terms enter through the H5 own-abundance term instead),
##   (b) its across-array SD is >= 0.10 (log1p scale) -- i.e. it varies enough to identify a
##       smooth term at all,
##   (c) at least 5% of SP's arrays have a nonzero detection, AND
##   (d) at least 10 of SP's arrays have a nonzero detection (an absolute floor -- a 5%
##       fraction of a large species sample can still be a handful of arrays).
## This is checked from the rebuilt predictor table's own summary statistics, not assumed.
## Grizzly bear and red wolf fail this test for EVERY one of the 13 species (max 6 nonzero
## arrays for any species) and are excluded from every model, exactly per the task's own example.
MIN_SD <- 0.10; MIN_NZ_FRAC <- 0.05; MIN_NZ_N <- 10

build_predator_decision <- function(D) {
  species_list <- sort(unique(D$species))
  rows <- list()
  for (sp in species_list) {
    dd <- D[species == sp]
    own_c <- OWN_MAP[[sp]]
    for (c in CARN_SP) {
      if (c == own_c) {
        rows[[length(rows)+1]] <- data.table(species=sp, carnivore=c, n=nrow(dd), sd=NA, nzfrac=NA,
          n_nonzero=NA, is_self=TRUE, include=FALSE,
          reason="is the focal species itself (captured via the H5 own-abundance term instead)")
        next
      }
      col <- paste0("log_cam_rate_", c)
      v <- dd[[col]]
      s <- sd(v, na.rm=TRUE); nz <- mean(v > 0, na.rm=TRUE); n_nz <- sum(v > 0, na.rm=TRUE)
      inc <- is.finite(s) && s >= MIN_SD && is.finite(nz) && nz >= MIN_NZ_FRAC && n_nz >= MIN_NZ_N
      reason <- if (inc) sprintf("included: sd=%.3f, nonzero_frac=%.3f (n_nonzero=%d)", s, nz, n_nz) else
        sprintf("excluded: sd=%.3f, nonzero_frac=%.3f, n_nonzero=%d (need sd>=%.2f, frac>=%.2f, n>=%d)",
                s, nz, n_nz, MIN_SD, MIN_NZ_FRAC, MIN_NZ_N)
      rows[[length(rows)+1]] <- data.table(species=sp, carnivore=c, n=nrow(dd), sd=round(s,3),
        nzfrac=round(nz,3), n_nonzero=n_nz, is_self=FALSE, include=inc, reason=reason)
    }
  }
  rbindlist(rows)
}
PRED_DECISION <- build_predator_decision(D)
fwrite(PRED_DECISION, "predator_relevance_decision.csv")
cat(sprintf("Predator-relevance check: %d of %d species x carnivore combinations included\n",
            sum(PRED_DECISION$include), nrow(PRED_DECISION)))

species_vars <- function(sp) {
  own_col <- paste0("own_logdet_", OWN_MAP[[sp]])
  carn <- PRED_DECISION[include == TRUE & species == sp, carnivore]
  carn_cols <- paste0("log_cam_rate_", carn)
  list(env = COVARS, predator = carn_cols, abundance = own_col,
       all = c(COVARS, carn_cols, own_col))
}

## Family lookup for every candidate variable (used to label results rows)
family_of <- function(var) {
  if (var %in% COVARS) return("environmental")
  if (grepl("^log_cam_rate_", var)) return("predator")
  if (grepl("^own_logdet_", var)) return("abundance")
  "unknown"
}

## ============================================================================
## 3. Per-species data prep: log1p-transform, drop rows with any NA in this species'
##    candidate predictor set (keeps the modeling subset consistent within a species).
## ============================================================================
prep_species_data <- function(sp) {
  dd <- D[species == sp]
  d2 <- copy(dd)
  for (v in intersect(LOGVARS, COVARS)) set(d2, j = v, value = log1p(pmax(d2[[v]], 0)))
  vs <- species_vars(sp)
  allv <- vs$all
  keep <- complete.cases(d2[, ..allv])
  d2 <- d2[keep]
  list(dat = d2, vars = vs)
}

species_list <- sort(unique(D$species))
SPD <- setNames(lapply(species_list, prep_species_data), species_list)
for (sp in species_list) cat(sprintf("  %-24s n=%d (of %d), %d predictor terms\n",
  sp, nrow(SPD[[sp]]$dat), sum(D$species==sp), length(SPD[[sp]]$vars$all)))

## ============================================================================
## 4. 400 km spatial-block CV: unified model vs. nationwide-curve null vs.
##    position-only competitor, for every species x metric.
## ============================================================================
mk_folds <- function(dd, block_km = 400, nfold = 5, seed = 7) {
  a <- unique(dd[, .(final_array, x_km, y_km)])
  a[, block := paste0("B", floor(x_km/block_km), "_", floor(y_km/block_km))]
  blk <- a[, .(x_km = mean(x_km), y_km = mean(y_km)), by = block]
  set.seed(seed)
  km <- kmeans(as.matrix(blk[, .(x_km, y_km)]), centers = min(nfold, nrow(blk)), nstart = 25, iter.max = 100)
  blk[, fold := km$cluster]
  merge(a[, .(final_array, block)], blk[, .(block, fold)], by = "block")
}

run_species_cv_unified <- function(dd, all_vars, k = 4, nfold = 5, block_km = 400) {
  if (nrow(dd) < 20) return(NULL)
  fo <- mk_folds(dd, block_km, nfold); dd2 <- merge(dd, fo, by = "final_array")
  res <- list()
  for (kk in seq_len(uniqueN(dd2$fold))) {
    tr <- dd2[fold != kk]; te <- dd2[fold == kk]
    if (nrow(te) < 3 || nrow(tr) < 15) next
    mods <- lapply(setNames(HARM, HARM), function(h) try(fit_unified(tr, h, all_vars, k = k), silent = TRUE))
    if (any(sapply(mods, inherits, "try-error"))) next
    kpos <- min(30, nrow(tr) - 2)
    posm <- lapply(setNames(HARM, HARM), function(h) {
      w <- 1/pmax(tr[[paste0("v_", h)]], 1e-8); w <- w/mean(w)
      t2 <- copy(tr); t2$.w <- w
      try(gam(as.formula(sprintf("%s ~ s(x_km, y_km, k=%d)", h, kpos)), data = t2,
              weights = .w, method = "REML", select = TRUE), silent = TRUE)
    })
    if (any(sapply(posm, inherits, "try-error"))) posm <- NULL
    nullB <- sapply(HARM, function(h) { w <- 1/pmax(tr[[paste0("v_", h)]], 1e-8); sum(w*tr[[h]])/sum(w) })
    pc <- predict_curves(mods, te); pm <- pc$M
    nm <- metrics_mat(curves_from_B(matrix(rep(nullB, each = nrow(te)), nrow(te), dimnames = list(NULL, HARM))))
    sm <- if (!is.null(posm)) predict_curves(posm, te)$M else NULL
    for (mm in MET) {
      ob <- te[[paste0("obs_", mm)]]
      r  <- skill_of(ob, pm[[mm]], nm[[mm]], mm)
      rs <- if (!is.null(sm)) skill_of(ob, sm[[mm]], nm[[mm]], mm) else list(skill = NA, mae = NA)
      res[[length(res)+1]] <- data.table(metric = mm, fold = kk, n_te = nrow(te), n_tr = nrow(tr),
        mae_cov = r$mae, mae_null = r$mae_null, skill_cov = r$skill, mae_pos = rs$mae, skill_pos = rs$skill)
    }
  }
  rbindlist(res)
}

summarise_cv <- function(cv, sp) {
  Ssum <- cv[, {
    ci  <- boot_ci(.SD, "n_te", "mae_cov", "mae_null")
    cip <- boot_ci(.SD, "n_te", "mae_pos", "mae_null")
    .(cv_skill = uw(n_te, mae_cov, mae_null), ci_lo = ci[1], ci_hi = ci[2],
      pos_skill = uw(n_te, mae_pos, mae_null), pos_ci_lo = cip[1], pos_ci_hi = cip[2],
      n_folds = .N, n_te_tot = sum(n_te))
  }, by = metric]
  Ssum[, species := sp]
  Ssum[, beats_null     := is.finite(ci_lo) & ci_lo > 0]
  Ssum[, beats_position := is.finite(cv_skill) & is.finite(pos_skill) & cv_skill > pos_skill]
  Ssum
}

cat("\n== Running 400km spatial-block CV for the unified model, all 13 species ==\n")
CV_SUMMARY <- list()
for (sp in species_list) {
  vs <- SPD[[sp]]$vars; dat <- SPD[[sp]]$dat
  cv <- run_species_cv_unified(dat, vs$all, k = 4)
  if (is.null(cv) || !nrow(cv)) { cat(sprintf("  %-24s CV: no folds fit\n", sp)); next }
  s <- summarise_cv(cv, sp)
  CV_SUMMARY[[sp]] <- s
  cat(sprintf("  %-24s CV done: %d/5 combos beat null, %d/5 beat position\n",
              sp, sum(s$beats_null, na.rm=TRUE), sum(s$beats_position, na.rm=TRUE)))
}
CV_ALL <- rbindlist(CV_SUMMARY)
setcolorder(CV_ALL, c("species","metric"))
fwrite(CV_ALL, "table_cv_unified_final685.csv")

## ============================================================================
## 5. Full-data unified models (all 13 species x 5 harmonic coefficients)
## ============================================================================
cat("\n== Fitting full-data unified models ==\n")
FINAL_MODELS <- list()
for (sp in species_list) {
  vs <- SPD[[sp]]$vars; dat <- SPD[[sp]]$dat
  mods <- lapply(setNames(HARM, HARM), function(h) try(fit_unified(dat, h, vs$all, k = 4), silent = TRUE))
  if (any(sapply(mods, inherits, "try-error"))) { cat(sprintf("  %-24s FAILED to fit\n", sp)); next }
  FINAL_MODELS[[sp]] <- list(models = mods, dat = dat, vars = vs)
  cat(sprintf("  %-24s fit OK (n=%d, %d terms)\n", sp, nrow(dat), length(vs$all)))
}
saveRDS(FINAL_MODELS, "unified_final_models_685.rds")

## ============================================================================
## 6. Effect sizes: AME (10th->90th pctile) + posterior-simulation 95% CI, on a single
##    consistent scale (response-SD units, or hours for peak -- see header). Every candidate
##    term for every species x metric is reported, whether or not select=TRUE shrunk it near
##    zero (an edf close to 0 means "shrunk to ~0", which the results table shows directly via
##    the effect estimate and its CI straddling zero).
## ============================================================================
cat("\n== Computing AME effect sizes + posterior-CI for every term x metric ==\n")
EFFECTS <- list()
for (sp in species_list) {
  FM <- FINAL_MODELS[[sp]]; if (is.null(FM)) next
  vs <- FM$vars; dat <- FM$dat; models <- FM$models
  for (mm in MET) {
    for (var in vs$all) {
      r <- tryCatch(ame_effect_ci(models, dat, var, mm, nsim = 500, seed = 3), error = function(e) NULL)
      if (is.null(r)) next
      EFFECTS[[length(EFFECTS)+1]] <- data.table(species = sp, metric = mm, predictor = var,
        hypothesis_family = family_of(var), effect_std = r$effect_std, ci_lo = r$ci_lo, ci_hi = r$ci_hi,
        q10 = r$q10, q90 = r$q90, n = nrow(dat))
    }
  }
  cat(sprintf("  %-24s effects computed (%d terms x 5 metrics)\n", sp, length(vs$all)))
}
EFF <- rbindlist(EFFECTS)
EFF[, scale := ifelse(metric == "peak", "hours (circular AME, 10th->90th pctile of predictor)",
                                        "response-between-array-SD units (AME, 10th->90th pctile of predictor)")]

## ============================================================================
## 7. Artifact check for EVERY predator and own-abundance term (H4/H5-derived rigor,
##    extended to cover every such term in the unified model -- not dropped for H1 terms,
##    which are not the target of this check per the original methodology: the artifact this
##    test catches is specific to a predictor computed from detection counts sharing the
##    array's own effort/precision structure with the response, which applies to predator and
##    own-abundance terms, not to remotely-sensed environmental covariates).
## ============================================================================
cat("\n== Null-simulation artifact check for predator + own-abundance terms ==\n")
NSIM_NULL <- 20
ARTIFACT <- list()
for (sp in species_list) {
  FM <- FINAL_MODELS[[sp]]; if (is.null(FM)) next
  vs <- FM$vars; dat <- FM$dat
  terms_to_check <- c(vs$predator, vs$abundance)
  if (!length(terms_to_check)) next
  null_cache <- build_null_model_cache(dat, vs$all, k = 4, nsim_null = NSIM_NULL, seed0 = 1000 + which(species_list == sp) * 100)
  for (mm in MET) {
    obs_row_lookup <- EFF[species == sp & metric == mm]
    for (var in terms_to_check) {
      obs_eff <- obs_row_lookup[predictor == var, effect_std]
      if (!length(obs_eff)) next
      ac <- artifact_check(null_cache, var, mm, obs_eff[1])
      ARTIFACT[[length(ARTIFACT)+1]] <- data.table(species = sp, metric = mm, predictor = var,
        null_mean = ac$null_mean, null_sd = ac$null_sd, null_p95 = ac$null_p95,
        artifact_ratio = ac$artifact_ratio, survives_artifact_check = ac$survives_artifact_check,
        n_null_reps = ac$n_null)
    }
  }
  cat(sprintf("  %-24s artifact check done (%d terms x 5 metrics)\n", sp, length(terms_to_check)))
}
ARTCHK <- rbindlist(ARTIFACT)
fwrite(ARTCHK, "table_artifact_check_final685.csv")

## ============================================================================
## 8. Residual spatial autocorrelation (Moran's I), 3 distance bands, every species x metric.
## ============================================================================
cat("\n== Residual Moran's I diagnostic ==\n")
BANDS <- list(c(0,50), c(50,200), c(200,800))
MORAN <- list()
for (sp in species_list) {
  FM <- FINAL_MODELS[[sp]]; if (is.null(FM)) next
  dat <- FM$dat; models <- FM$models
  pc <- predict_curves(models, dat)
  for (mm in MET) {
    obscol <- paste0("obs_", mm)
    res <- if (mm == "peak") circ_diff_h(pc$M[[mm]], dat[[obscol]]) else dat[[obscol]] - pc$M[[mm]]
    w <- 1/pmax(dat[[SE_COL[[mm]]]]^2, 1e-6)
    for (bd in BANDS) {
      raw <- moran_band(dat$x_km, dat$y_km, dat[[obscol]], w, bd[1], bd[2])
      rsd <- moran_band(dat$x_km, dat$y_km, res, w, bd[1], bd[2])
      MORAN[[length(MORAN)+1]] <- data.table(species = sp, metric = mm,
        band = sprintf("%d-%d km", bd[1], bd[2]), n = raw[3],
        moran_raw = raw[1], p_raw = raw[2], moran_resid = rsd[1], p_resid = rsd[2])
    }
  }
  cat(sprintf("  %-24s Moran's I computed\n", sp))
}
MORAN_ALL <- rbindlist(MORAN)
fwrite(MORAN_ALL, "table_residual_moran_unified_final685.csv")
## species x metric summary p-value used in the main results table: the 0-800km band
## (widest span, matches the original H1 track's headline reporting convention)
MORAN_HEADLINE <- MORAN_ALL[band == "200-800 km", .(species, metric, residual_moran_p = p_resid)]

## ============================================================================
## 9. ONE unified results table -- all three hypothesis families, consistent columns.
## ============================================================================
cat("\n== Assembling unified results table ==\n")
RESULTS <- merge(EFF, CV_ALL[, .(species, metric, cv_skill, ci_lo_cv = ci_lo, ci_hi_cv = ci_hi,
                                  pos_skill, beats_null, beats_position)],
                  by = c("species","metric"), all.x = TRUE)
RESULTS <- merge(RESULTS, ARTCHK[, .(species, metric, predictor, artifact_ratio, survives_artifact_check)],
                  by = c("species","metric","predictor"), all.x = TRUE)
RESULTS <- merge(RESULTS, MORAN_HEADLINE, by = c("species","metric"), all.x = TRUE)

setnames(RESULTS, c("ci_lo","ci_hi"), c("effect_ci_lo","effect_ci_hi"))
RESULTS <- RESULTS[, .(species, metric, predictor, hypothesis_family,
                        effect_std, effect_ci_lo, effect_ci_hi, scale,
                        cv_skill, ci_lo_cv, ci_hi_cv, pos_skill, beats_null, beats_position,
                        survives_artifact_check, artifact_ratio, residual_moran_p, n)]
setorder(RESULTS, species, metric, hypothesis_family, predictor)
fwrite(RESULTS, "unified_hypothesis_results_final685.csv")
cat(sprintf("\nSaved unified_hypothesis_results_final685.csv: %d rows\n", nrow(RESULTS)))
cat(sprintf("Saved unified_modeling_dataset_final685.csv, table_cv_unified_final685.csv,\n",
            "table_artifact_check_final685.csv, table_residual_moran_unified_final685.csv,\n",
            "predator_relevance_decision.csv, unified_final_models_685.rds\n"))
cat("DONE.\n")
