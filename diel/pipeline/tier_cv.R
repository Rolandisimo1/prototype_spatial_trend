## tier_cv.R -- per-held-out-array CV predictions + envelope diagnostics.
## Replicates step3_cv.R fold construction EXACTLY (400 km blocks, kmeans seed 7,
## 5 folds) so results line up with covariate_model_cv.csv / cv_folds_raw.csv.
## New output: one row per species x metric x held-out ARRAY carrying the covariate
## prediction, the nationwide-null prediction, the observed value, and that array's
## envelope diagnostics measured against ITS OWN FOLD'S training set. That is what
## the confidence tiers must be validated on -- a tier assigned from training data
## the array was not part of.
suppressPackageStartupMessages({library(data.table); library(mgcv)})

REF_SR <- 6.2202; REF_SS <- 17.9998; NBIN <- 48
HARM <- c("b0","s1","c1","s2","c2"); MET <- c("noct","crep","conc","peak","act")
bcen <- (seq_len(NBIN) - 0.5) * 0.5
trad <- bcen/24*2*pi
XD   <- cbind(1, sin(trad), cos(trad), sin(2*trad), cos(2*trad))
ang48 <- bcen/24*2*pi
NIGHT <- bcen < REF_SR | bcen >= REF_SS
CREP  <- (abs(bcen - REF_SR) <= 1.5) | (abs(bcen - REF_SS) <= 1.5)

metrics_mat <- function(M) {
  M <- pmax(M, 0); s <- colSums(M)
  P <- sweep(M, 2, ifelse(s > 0, s, NA_real_), "/")
  C <- as.numeric(crossprod(cos(ang48), P)); S <- as.numeric(crossprod(sin(ang48), P))
  data.table(noct = 100*as.numeric(crossprod(NIGHT, P)),
             crep = 100*as.numeric(crossprod(CREP,  P)),
             conc = sqrt(C^2 + S^2),
             peak = bcen[max.col(t(M), ties.method = "first")],
             act  = colMeans(M)/apply(M, 2, max))
}
curves_from_B <- function(B) exp(XD %*% t(as.matrix(B)))
circ_diff_h <- function(a, b) { d <- (a - b) %% 24; ifelse(d > 12, d - 24, d) }

TRANS <- list(pop_1km="log1p", ntl_1km="log1p", nlcd_1k_crop="log1p",
              nlcd_1k_pasture="log1p", rug_1km="log1p")
apply_trans <- function(D, vars) {
  D <- copy(D)
  for (v in intersect(vars, names(D))) {
    tt <- TRANS[[v]]; if (is.null(tt)) next
    if (tt == "log1p") set(D, j = v, value = log1p(pmax(D[[v]], 0)))
  }
  D[]
}
## PRIMARY variant only: setA_poponly (population as the human-pressure term).
## ntl_1km is NOT used; nlcd_1k_impervious appears nowhere; crop and pasture are
## separate terms and separate from the human-pressure term; no composite index.
VARS_POP  <- c("pop_1km","nlcd_1k_crop","nlcd_1k_pasture","tcc_1km",
               "t_warmmonth","t_coldmonth","rug_1km")
## bear: 55 arrays cannot afford the 7-covariate specification. The envelope
## reconstruction (verified to <0.002 Mahalanobis units against the published
## extrapolation_envelope.csv) shows the shipped bear envelope used tcc + t_coldmonth.
## step3_cv.R takes intersect(BEAR_SET, variant_vars) where BEAR_SET is
## c("ntl_1km","tcc_1km","t_coldmonth"); for setA_poponly that intersection is
## tcc_1km + t_coldmonth ONLY -- no human-pressure term enters bear's model at all.
BEAR_ENV_VARS <- c("tcc_1km","t_coldmonth")
BEAR_FIT_VARS <- c("tcc_1km","t_coldmonth")

mk_form <- function(y, vars, k = 5, linear = FALSE) {
  if (!length(vars)) return(as.formula(paste(y, "~ 1")))
  tt <- if (linear) vars else sprintf("s(%s, k=%d)", vars, k)
  as.formula(paste(y, "~", paste(tt, collapse = " + ")))
}
fit_coef <- function(dat, y, vars, k = 5, linear = FALSE) {
  w <- 1/pmax(dat[[paste0("v_", y)]], 1e-8); w <- w/mean(w)
  d2 <- copy(dat); d2$.w <- w
  gam(mk_form(y, vars, k, linear), data = d2, weights = .w, method = "REML", select = TRUE)
}
predict_B <- function(models, newd) {
  B <- sapply(HARM, function(h) as.numeric(predict(models[[h]], newd, type = "response")))
  matrix(B, nrow = nrow(newd), dimnames = list(NULL, HARM))
}

## MESS-style per-covariate coverage + Mahalanobis distance, identical to step8_predict.R
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
  S  <- matrix(S, nrow = length(vars))
  Si <- tryCatch(solve(S + diag(1e-8, ncol(S))), error = function(e) diag(ncol(S)))
  md <- function(X) { X <- sweep(as.matrix(X), 2, mu, "-"); sqrt(pmax(rowSums((X %*% Si) * X), 0)) }
  list(mess = apply(Mm, 1, min), maha = md(newd[, ..vars]),
       maha_cut = as.numeric(quantile(md(train[, ..vars]), 0.95)),
       n_out = rowSums(Mm < 0))
}

H  <- fread("array_harmonics_raw.csv")
AC <- fread("array_covariates_aso.csv")
D  <- merge(H, AC, by = "array_id")
cat("merged", nrow(D), "species-array rows\n")

BLOCK_KM <- 400; NFOLD <- 5
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

out <- list()
for (SP in unique(D$species)) {
  dd <- D[species == SP]
  fo <- mk_folds(dd); dd <- merge(dd, fo, by = "array_id")
  is_bear <- SP == "American Black Bear"
  fitv <- if (is_bear) BEAR_FIT_VARS else VARS_POP
  envv <- if (is_bear) BEAR_ENV_VARS else VARS_POP
  allv <- union(fitv, envv)
  dv <- apply_trans(dd, allv)
  dv <- dv[complete.cases(dv[, ..allv])]

  for (k in 1:NFOLD) {
    tr <- dv[fold != k]; te <- dv[fold == k]
    if (nrow(te) < 3 || nrow(tr) < 30) next
    kk <- if (is_bear) 3 else 5
    mods <- lapply(setNames(HARM, HARM), function(h)
      try(fit_coef(tr, h, fitv, k = kk, linear = is_bear), silent = TRUE))
    if (any(sapply(mods, inherits, "try-error"))) { cat("fit fail", SP, k, "\n"); next }
    ## nationwide-curve null: inverse-variance weighted mean coefficient
    nullB <- sapply(HARM, function(h) {
      w <- 1/pmax(tr[[paste0("v_", h)]], 1e-8); sum(w*tr[[h]])/sum(w) })
    Bc <- predict_B(mods, te)
    Mc <- metrics_mat(curves_from_B(Bc))
    Mn <- metrics_mat(curves_from_B(matrix(rep(nullB, each = nrow(te)), nrow(te),
                                           dimnames = list(NULL, HARM))))
    ev <- env_calc(tr, te, envv)
    for (mm in MET) {
      ob <- te[[paste0("obs_", mm)]]
      e  <- if (mm == "peak") abs(circ_diff_h(Mc[[mm]], ob)) else abs(Mc[[mm]] - ob)
      en <- if (mm == "peak") abs(circ_diff_h(Mn[[mm]], ob)) else abs(Mn[[mm]] - ob)
      out[[length(out)+1]] <- data.table(
        species = SP, metric = mm, fold = k, array_id = te$array_id,
        lon = te$lon, lat = te$lat, events = te$events,
        obs = ob, pred_cov = Mc[[mm]], pred_null = Mn[[mm]],
        ae_cov = e, ae_null = en,
        mess = ev$mess, maha = ev$maha, maha_cut = ev$maha_cut,
        maha_ratio = ev$maha / ev$maha_cut, n_cov_outside = ev$n_out)
    }
  }
  cat(sprintf("  %-24s folds done\n", SP))
}
R <- rbindlist(out)
fwrite(R, "cv_array_tiers.csv")
cat("rows", nrow(R), " arrays", uniqueN(R$array_id), "\n")

## sanity: fold-level skill must reproduce cv_folds_raw.csv (setA_poponly)
FS <- R[, .(n_te = .N, skill = 1 - sum(ae_cov)/sum(ae_null)), by = .(species, metric, fold)]
fwrite(FS, "cv_fold_check.csv")
print(FS[species == "Eastern Gray Squirrel"][order(metric, fold)][1:10])
