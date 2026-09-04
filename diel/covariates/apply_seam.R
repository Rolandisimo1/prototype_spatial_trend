## apply_seam.R -- APPLY the existing fitted covariate models to the updated grid.
## The GAM objects were never serialised (final_models.rds is absent from the store), so
## fit_species() is re-instantiated from the UNCHANGED training data (array_harmonics_raw +
## array_covariates) exactly as dielfix.R defines it. This is a deterministic re-instantiation
## of the SAME fit, not a refit against new data: no training row, formula, weight, basis
## dimension or method changes. Reproduction on unchanged cells is verified downstream.
source("dielfix.R")
H <- fread("array_harmonics_raw.csv"); AC <- fread("array_covariates.csv")
D <- merge(H, AC, by = "array_id"); G <- fread("grid_covariates_25km_seamfixed.csv")
RMK <- fread("range_mask_cells.csv"); SUB <- fread("subset_cells.csv")$cell25
NSIM <- 400L
SCI <- c("White-tailed Deer"="Odocoileus virginianus","Northern Raccoon"="Procyon lotor",
         "Eastern Gray Squirrel"="Sciurus carolinensis","Coyote"="Canis latrans",
         "American Black Bear"="Ursus americanus")
BEAR_ENV <- c("tcc_1km","t_coldmonth")

env_calc <- function(train, newd, vars) {
  Mm <- sapply(vars, function(v) {
    tr <- train[[v]]; tr <- tr[is.finite(tr)]; nv <- newd[[v]]
    rng <- range(tr); n <- length(tr)
    f <- vapply(nv, function(z) 100*sum(tr <= z)/n, numeric(1))
    ifelse(nv < rng[1], 100*(nv-rng[1])/max(diff(rng),1e-9),
    ifelse(nv > rng[2], 100*(rng[2]-nv)/max(diff(rng),1e-9), pmin(f,100-f)*2)) })
  Mm <- matrix(Mm, nrow=nrow(newd), dimnames=list(NULL,vars))
  mu <- colMeans(train[, ..vars], na.rm=TRUE)
  S <- matrix(cov(as.matrix(train[, ..vars]), use="pairwise.complete.obs"), nrow=length(vars))
  Si <- tryCatch(solve(S + diag(1e-8,ncol(S))), error=function(e) diag(ncol(S)))
  md <- function(X){X<-sweep(as.matrix(X),2,mu,"-"); sqrt(pmax(rowSums((X%*%Si)*X),0))}
  list(mess=apply(Mm,1,min), maha=md(newd[, ..vars]),
       maha_cut=as.numeric(quantile(md(train[, ..vars]),0.95)), n_out=rowSums(Mm<0))
}

## posterior simulation of the coefficient models -> metric CIs (as covlib.R predict_curves)
sim_ci <- function(models, newd, nsim, seed = 5) {
  set.seed(seed)
  n <- nrow(newd)
  sims <- array(NA_real_, c(n, 5, nsim))
  for (h in seq_along(HARM)) {
    m <- models[[HARM[h]]]
    Xp <- predict(m, newd, type="lpmatrix"); Vb <- vcov(m)
    L <- tryCatch(chol(Vb + diag(1e-10, ncol(Vb))), error=function(e) NULL)
    if (is.null(L)) { sims[,h,] <- as.numeric(predict(m, newd, type="response")); next }
    Bs <- matrix(rnorm(nsim*ncol(Vb)), nsim) %*% L
    Bs <- sweep(Bs, 2, coef(m), "+")
    sims[,h,] <- Xp %*% t(Bs)
  }
  arr <- lapply(setNames(MET, MET), function(mm) matrix(NA_real_, n, nsim))
  for (i in seq_len(nsim)) {
    MM <- metrics_mat(curves_from_B(sims[,,i]))
    for (mm in MET) arr[[mm]][,i] <- MM[[mm]]
  }
  arr
}

coef_rows <- list(); pred_rows <- list(); conf_rows <- list(); curve_rows <- list()
for (SP in names(SCI)) {
  dd <- D[species==SP]; fitv <- vars_for(SP)
  envv <- if (SP=="American Black Bear") BEAR_ENV else fitv
  allv <- union(fitv, envv)
  dv <- apply_trans(dd, allv); dv <- dv[complete.cases(dv[, ..allv])]
  mods <- fit_species(dv, SP); env <- coef_envelope(dv)

  gg <- G[cell_id %in% RMK[species==SCI[[SP]], cell_id] & cell25 %in% SUB]
  gv <- apply_trans(gg, allv); gv <- gv[complete.cases(gv[, ..allv])]
  if (!nrow(gv)) next

  cn <- constrain_B(predict_B(mods, gv), env)
  coef_rows[[SP]] <- data.table(species=SP, cell25=gv$cell25, b0=cn$B[,1], s1=cn$B[,2],
                                c1=cn$B[,3], s2=cn$B[,4], c2=cn$B[,5], lambda=cn$lambda)
  Mf <- metrics_mat(curves_from_B(cn$B))
  SM <- sim_ci(mods, gv, NSIM)
  for (mm in MET) {
    q <- apply(SM[[mm]], 1, quantile, c(0.025,0.5,0.975), na.rm=TRUE)
    pred_rows[[length(pred_rows)+1]] <- data.table(
      species=SP, variant="setA_poponly", metric=mm, cell25=gv$cell25,
      lon=gv$lon, lat=gv$lat, cell_id=gv$cell_id,
      value=Mf[[mm]], ci_lo=q[1,], ci_med=q[2,], ci_hi=q[3,])
  }
  cur <- curves_from_B(cn$B); prop <- sweep(cur, 2, colSums(cur), "/")
  q1000 <- round(prop*1000)
  curve_rows[[SP]] <- data.table(species=SP, cell25=rep(gv$cell25, each=NBIN),
                                 bin=rep(seq_len(NBIN), ncol(q1000)),
                                 rate_permille=as.integer(as.vector(q1000)))
  ev <- env_calc(dv, gv, envv)
  conf_rows[[SP]] <- data.table(species=SP, cell25=gv$cell25, lon=gv$lon, lat=gv$lat,
    mess=ev$mess, maha=ev$maha, maha_cut=ev$maha_cut, maha_ratio=ev$maha/ev$maha_cut,
    n_cov_outside=ev$n_out, in_envelope=(ev$mess>=0 & ev$maha<=ev$maha_cut),
    seam_interpolated=gv$seam_interpolated)
  cat(SP, nrow(gv), "\n", flush=TRUE)
}
CO <- rbindlist(coef_rows); PR <- rbindlist(pred_rows)
CU <- rbindlist(curve_rows); CF <- rbindlist(conf_rows)
CF[, score := pmax(0, pmin(100, 100 - 50*maha_ratio))]
CF[, tier := fifelse(score>=50,"Interpolation", fifelse(score>=37.5,"Near",
             fifelse(score>=12.5,"Moderate extrapolation","Severe extrapolation")))]
fwrite(CO,"new_coefs.csv"); fwrite(PR,"new_preds.csv")
fwrite(CU,"new_curves.csv"); fwrite(CF,"new_conf.csv")
cat("DONE", nrow(CO), nrow(PR), nrow(CU), nrow(CF), "\n")
