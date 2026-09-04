## step3_cv.R -- the decisive test: 400 km spatial-block CV, three competitors.
##  (a) covariate model   (b) position model   (c) nationwide-curve null
## Fold construction is copied from cv_run2.R (400 km blocks, kmeans seed 7, 5 folds)
## so results line up with skill_bootstrap_cis.csv.
## The position competitor is fitted in the SAME harmonic framework (coefficients ~
## s(x_km,y_km)) so the ONLY difference from the covariate model is the predictor set;
## the published bam benchmark is reported alongside in the report.
source("covlib.R")
set.seed(20260814)
BLOCK_KM <- 400; NFOLD <- 5; NB <- 4000L

H  <- fread("array_harmonics_raw.csv")
AC <- fread("array_covariates.csv")
CB <- if (file.exists("camera_setb_array.csv")) fread("camera_setb_array.csv") else NULL
if (!is.null(CB)) AC <- merge(AC, CB, by = "array_id", all.x = TRUE)
D  <- merge(H, AC, by = "array_id")
D[, project_id := sub("\\|.*$", "", array_id)]     # cluster unit for the project bootstrap
cat("merged rows", nrow(D), " projects", uniqueN(D$project_id), "\n")

## ---- fold construction, identical scheme to cv_run2.R ----
mk_folds <- function(dd) {
  a <- unique(dd[, .(array_id, x_km, y_km)])
  a[, block := paste0("B", floor(x_km/BLOCK_KM), "_", floor(y_km/BLOCK_KM))]
  blk <- a[, .(x_km = mean(x_km), y_km = mean(y_km)), by = block]
  set.seed(7)
  km <- kmeans(as.matrix(blk[, .(x_km, y_km)]), centers = min(NFOLD, nrow(blk)),
               nstart = 25, iter.max = 100)
  blk[, fold := km$cluster]
  merge(a[, .(array_id, block)], blk[, .(block, fold)], by = "block")
}

## ---- MESS-style extrapolation diagnostic ----
## per-covariate coverage (negative = outside training range) + Mahalanobis distance
mess_calc <- function(train, newd, vars) {
  M <- sapply(vars, function(v) {
    tr <- train[[v]]; tr <- tr[is.finite(tr)]; nv <- newd[[v]]
    rng <- range(tr); n <- length(tr)
    f <- sapply(nv, function(z) 100*sum(tr <= z)/n)
    ifelse(nv < rng[1], 100*(nv - rng[1])/max(diff(rng), 1e-9),
    ifelse(nv > rng[2], 100*(rng[2] - nv)/max(diff(rng), 1e-9),
           pmin(f, 100 - f)*2))
  })
  M <- matrix(M, nrow = nrow(newd), dimnames = list(NULL, vars))
  mu <- colMeans(train[, ..vars], na.rm = TRUE)
  S  <- cov(as.matrix(train[, ..vars]), use = "pairwise.complete.obs")
  Si <- tryCatch(solve(S + diag(1e-8, ncol(S))), error = function(e) diag(ncol(S)))
  Z  <- sweep(as.matrix(newd[, ..vars]), 2, mu, "-")
  md <- sqrt(pmax(rowSums((Z %*% Si) * Z), 0))
  list(mess = apply(M, 1, min), maha = md, per_cov = M)
}

## ---- evaluation ----
skill_of <- function(obs, pred, null, metric) {
  if (metric == "peak") {
    e  <- abs(circ_diff_h(pred, obs)); en <- abs(circ_diff_h(null, obs))
  } else { e <- abs(pred - obs); en <- abs(null - obs) }
  ok <- is.finite(e) & is.finite(en)
  list(mae = mean(e[ok]), mae_null = mean(en[ok]),
       skill = 1 - mean(e[ok])/mean(en[ok]), n = sum(ok))
}

VARIANTS <- list(
  setA_lightsonly = list(vars = function(sp) setA("A"), tag = "A_lights"),
  setA_poponly    = list(vars = function(sp) setA("B"), tag = "A_pop"),
  setA_both       = list(vars = function(sp) setA("C"), tag = "A_both"),
  setB_map        = list(vars = function(sp) c(setA("A"), SETB_MAP), tag = "B_map"),
  setB_full       = list(vars = function(sp) c(setA("A"), SETB_MAP, HETERO[[sp]]), tag = "B_full"))

run_species <- function(SP) {
  dd <- D[species == SP]
  if (nrow(dd) < 20) { cat("skip", SP, "\n"); return(NULL) }
  fo <- mk_folds(dd); dd <- merge(dd, fo, by = "array_id")
  is_bear <- SP == "American Black Bear"
  res <- list(); envrows <- list()

  for (vn in names(VARIANTS)) {
    vars <- VARIANTS[[vn]]$vars(SP)
    if (is_bear) {
      if (vn %in% c("setB_map","setB_full")) next          # 55 arrays cannot afford Set B
      vars <- intersect(BEAR_SET, vars)
      if (!length(vars)) vars <- BEAR_SET[1]
    }
    vars <- intersect(vars, names(dd))
    vars <- vars[sapply(vars, function(v) sum(is.finite(dd[[v]])) > 0.8*nrow(dd) &&
                                          length(unique(dd[[v]])) > 5)]
    if (!length(vars)) next
    dv <- apply_trans(dd, vars)
    dv <- dv[complete.cases(dv[, ..vars])]

    for (k in 1:NFOLD) {
      tr <- dv[fold != k]; te <- dv[fold == k]
      if (nrow(te) < 3 || nrow(tr) < 30) next
      kk <- if (is_bear) 3 else 5
      lin <- is_bear
      mods <- lapply(setNames(HARM, HARM), function(h)
        try(fit_coef(tr, h, vars, k = kk, linear = lin), silent = TRUE))
      if (any(sapply(mods, inherits, "try-error"))) next
      pos <- lapply(setNames(HARM, HARM), function(h)
        try(fit_coef(tr, h, NULL) , silent = TRUE))       # placeholder, replaced below
      ## position competitor: same framework, s(x_km,y_km) instead of covariates
      posm <- lapply(setNames(HARM, HARM), function(h) {
        w <- 1/pmax(tr[[paste0("v_", h)]], 1e-8); w <- w/mean(w)
        t2 <- copy(tr); t2$.w <- w
        try(gam(as.formula(paste(h, "~ s(x_km, y_km, k=30)")), data = t2,
                weights = .w, method = "REML", select = TRUE), silent = TRUE)
      })
      if (any(sapply(posm, inherits, "try-error"))) posm <- NULL
      ## null: intercept-only weighted mean coefficient (one nationwide curve)
      nullB <- sapply(HARM, function(h) {
        w <- 1/pmax(tr[[paste0("v_", h)]], 1e-8); sum(w*tr[[h]])/sum(w) })

      pc <- predict_curves(mods, te)
      pm <- pc$M
      nm <- metrics_mat(curves_from_B(matrix(rep(nullB, each = nrow(te)), nrow(te),
                                             dimnames = list(NULL, HARM))))
      sm <- if (!is.null(posm)) predict_curves(posm, te)$M else NULL

      ## envelope for this fold's held-out arrays, against this fold's training set
      ms <- mess_calc(tr, te, vars)
      inside <- ms$mess >= 0 & ms$maha <= quantile(mess_calc(tr, tr, vars)$maha, 0.95)

      for (mm in MET) {
        ob <- te[[paste0("obs_", mm)]]
        r  <- skill_of(ob, pm[[mm]], nm[[mm]], mm)
        rs <- if (!is.null(sm)) skill_of(ob, sm[[mm]], nm[[mm]], mm) else list(skill = NA, mae = NA)
        ri <- if (sum(inside) >= 3) skill_of(ob[inside], pm[[mm]][inside], nm[[mm]][inside], mm)
              else list(skill = NA_real_, mae = NA_real_, mae_null = NA_real_, n = sum(inside))
        ro <- if (sum(!inside) >= 3) skill_of(ob[!inside], pm[[mm]][!inside], nm[[mm]][!inside], mm)
              else list(skill = NA_real_, mae = NA_real_, mae_null = NA_real_, n = sum(!inside))
        res[[length(res)+1]] <- data.table(
          species = SP, variant = vn, metric = mm, fold = k,
          n_te = nrow(te), n_tr = nrow(tr), n_inside = sum(inside),
          mae_cov = r$mae, mae_null = r$mae_null, skill_cov = r$skill,
          skill_pos = rs$skill, mae_pos = rs$mae,
          skill_cov_inside = ri$skill, n_inside_eval = ri$n,
          mae_cov_inside = ri$mae, mae_null_inside = ri$mae_null,
          skill_cov_outside = ro$skill, n_outside_eval = ro$n,
          mae_cov_outside = ro$mae, mae_null_outside = ro$mae_null,
          n_vars = length(vars))
      }
      envrows[[length(envrows)+1]] <- data.table(species = SP, variant = vn, fold = k,
        array_id = te$array_id, mess = ms$mess, maha = ms$maha, inside = inside)
    }
    cat(sprintf("  [%s] %-16s vars=%d folds done\n", SP, vn, length(vars)))
  }
  list(cv = rbindlist(res), env = rbindlist(envrows))
}

ALL <- lapply(unique(D$species), run_species)
CV  <- rbindlist(lapply(ALL, `[[`, "cv"))
EN  <- rbindlist(lapply(ALL, `[[`, "env"))
fwrite(CV, "cv_folds_raw.csv"); fwrite(EN, "cv_envelope_folds.csv")
cat("CV rows", nrow(CV), "\n")
