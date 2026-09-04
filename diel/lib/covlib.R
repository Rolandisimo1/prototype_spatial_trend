## covlib.R -- shared machinery for the covariate-driven predictive map.
## Unit of analysis = ARRAY. NO array random effect (the covariate coefficients must
## be free to carry the BETWEEN-array relationship, which is the predictive estimand).
suppressPackageStartupMessages({library(data.table); library(mgcv)})
source("cvlib.R")                       # REF_SR/REF_SS, AEA, curve_metrics, NBIN

HARM  <- c("b0","s1","c1","s2","c2")
MET   <- c("noct","crep","conc","peak","act")
bcen  <- (seq_len(NBIN) - 0.5) * 0.5
trad  <- bcen/24*2*pi
XD    <- cbind(1, sin(trad), cos(trad), sin(2*trad), cos(2*trad))
ang48 <- bcen/24*2*pi
NIGHT <- bcen < REF_SR | bcen >= REF_SS
CREP  <- (abs(bcen - REF_SR) <= 1.5) | (abs(bcen - REF_SS) <= 1.5)

## metrics from a 48 x n matrix of rate curves (columns = curves)
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
curves_from_B <- function(B) exp(XD %*% t(as.matrix(B)))     # B: n x 5 coefficient matrix

## covariate transforms: log1p for the heavily right-skewed / zero-inflated layers.
## Applied identically at cameras and on the grid.
TRANS <- list(ntl_1km = "log1p", pop_1km = "log1p", nlcd_1k_crop = "log1p",
              nlcd_1k_pasture = "log1p", rug_1km = "log1p", dist_water_m = "logkm")
apply_trans <- function(D, vars) {
  D <- copy(D)
  for (v in intersect(vars, names(D))) {
    tt <- TRANS[[v]]
    if (is.null(tt)) next
    if (tt == "log1p") set(D, j = v, value = log1p(pmax(D[[v]], 0)))
    if (tt == "logkm") set(D, j = v, value = log1p(pmax(D[[v]], 0)/1000))
  }
  D[]
}

## ---- covariate sets ----
## tmax_window DROPPED (survey-window quantity, not evaluable at an unsampled cell);
## replaced by WorldClim normals t_warmmonth / t_coldmonth.
SETA_BASE <- c("nlcd_1k_crop","nlcd_1k_pasture","tcc_1km","t_warmmonth","t_coldmonth","rug_1km")
setA <- function(variant) switch(variant,
  A = c("ntl_1km", SETA_BASE),                 # lights only  (drop pop)
  B = c("pop_1km", SETA_BASE),                 # population only (drop lights)
  C = c("ntl_1km","pop_1km", SETA_BASE))       # both -- the collinear version
## Set B additions: mechanistically distinct axes.
## NOTE: hetero rate is camera-derived and does NOT exist at an unsampled cell ->
## it enters the CV signal test only, never the continental map. Flagged everywhere.
SETB_MAP  <- c("snowfrac","dist_water_m","daylength_h")
HETERO <- list("White-tailed Deer" = "rate_Coyote",
               "Eastern Gray Squirrel" = "rate_Coyote",
               "Northern Raccoon" = c("rate_Coyote","rate_AmericanBlackBear"),
               "Coyote" = "rate_WhitetailedDeer",
               "American Black Bear" = character(0))

## bear cannot afford the full specification (55 arrays): reduced, linear terms only
BEAR_SET <- c("ntl_1km","tcc_1km","t_coldmonth")

mk_form <- function(y, vars, k = 5, linear = FALSE) {
  if (!length(vars)) return(as.formula(paste(y, "~ 1")))
  tt <- if (linear) vars else sprintf("s(%s, k=%d)", vars, k)
  as.formula(paste(y, "~", paste(tt, collapse = " + ")))
}

## fit one harmonic coefficient ~ covariates, weighted by inverse coefficient variance.
## select=TRUE gives each smooth an extra shrinkage penalty so unsupported terms can be
## shrunk to zero -- this IS the "where they are strong enough" filter.
fit_coef <- function(dat, y, vars, k = 5, linear = FALSE, select = TRUE) {
  w <- 1/pmax(dat[[paste0("v_", y)]], 1e-8)
  w <- w/mean(w)
  d2 <- copy(dat); d2$.w <- w
  gam(mk_form(y, vars, k, linear), data = d2, weights = .w,
      method = "REML", select = select)
}

## predict all five coefficients on newdata -> curves -> metrics, with posterior
## simulation for intervals (NOT the delta method).
predict_curves <- function(models, newd, nsim = 0, seed = 1) {
  B <- sapply(HARM, function(h) as.numeric(predict(models[[h]], newd, type = "response")))
  B <- matrix(B, nrow = nrow(newd), dimnames = list(NULL, HARM))
  out <- list(B = B, M = metrics_mat(curves_from_B(B)), curves = curves_from_B(B))
  if (nsim > 0) {
    set.seed(seed)
    sims <- array(NA_real_, c(nrow(newd), 5, nsim))
    for (h in seq_along(HARM)) {
      m  <- models[[HARM[h]]]
      Xp <- predict(m, newd, type = "lpmatrix")
      Vb <- vcov(m)
      L  <- tryCatch(chol(Vb + diag(1e-10, ncol(Vb))), error = function(e) NULL)
      if (is.null(L)) { sims[, h, ] <- B[, h]; next }
      Bs <- matrix(rnorm(nsim*ncol(Vb)), nsim) %*% L
      Bs <- sweep(Bs, 2, coef(m), "+")
      sims[, h, ] <- Xp %*% t(Bs)
    }
    out$sims <- sims
  }
  out
}

## metric quantiles across posterior draws
sim_metrics <- function(sims) {
  n <- dim(sims)[1]; ns <- dim(sims)[3]
  res <- vector("list", ns)
  for (i in seq_len(ns)) res[[i]] <- metrics_mat(curves_from_B(sims[,,i]))
  arr <- lapply(MET, function(mm) sapply(res, function(z) z[[mm]]))
  names(arr) <- MET
  arr
}

## circular-aware mean of an angle set (hours)
circ_mean_h <- function(h, w = NULL) {
  a <- h/24*2*pi; if (is.null(w)) w <- rep(1, length(a))
  ok <- is.finite(a) & is.finite(w); a <- a[ok]; w <- w[ok]
  ((atan2(sum(w*sin(a)), sum(w*cos(a)))) %% (2*pi))/(2*pi)*24
}
