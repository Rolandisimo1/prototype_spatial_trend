## covlib_final685.R -- shared machinery for the covariate-driven predictive map,
## rerun on the manuscript's own final_array site key (685-site definition).
## Unit of analysis = ARRAY (final_array). NO array random effect -- the covariate
## coefficients must be free to carry the BETWEEN-array relationship, which is the
## predictive estimand. Model approach ported unchanged from covlib.R / cvlib.R
## (prior array_id-keyed track); only the join key and covariate set differ.
suppressPackageStartupMessages({library(data.table); library(mgcv)})

REF_SR <- 6.2202; REF_SS <- 17.9998
NBIN   <- 48
HARM   <- c("b0","s1","c1","s2","c2")
MET    <- c("noct","crep","conc","peak","act")

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

circ_diff_h <- function(a, b) { d <- (a - b) %% 24; ifelse(d > 12, d - 24, d) }

## ---- the manuscript's 5 covariates (complete for all 685 final_array sites) ----
COVARS  <- c("pop_1km","ag_5km","tcc_1km","tmax_hottest_month","rug_5km")
LOGVARS <- c("pop_1km","ag_5km","rug_5km")     # right-skewed / zero-inflated -> log1p
apply_trans <- function(dt, vars = COVARS) {
  dt <- copy(dt)
  for (v in intersect(LOGVARS, vars)) set(dt, j = v, value = log1p(pmax(dt[[v]], 0)))
  dt[]
}

mk_form <- function(y, vars, k = 5) {
  if (!length(vars)) return(as.formula(paste(y, "~ 1")))
  tt <- sprintf("s(%s, k=%d)", vars, k)
  as.formula(paste(y, "~", paste(tt, collapse = " + ")))
}

## fit one harmonic coefficient ~ covariates, weighted by inverse coefficient variance.
## select=TRUE gives each smooth an extra shrinkage penalty so unsupported terms can be
## shrunk to zero -- this IS the "where they are strong enough" filter.
fit_coef <- function(dat, y, vars, k = 5, select = TRUE) {
  w <- 1/pmax(dat[[paste0("v_", y)]], 1e-8)
  w <- w/mean(w)
  d2 <- copy(dat); d2$.w <- w
  gam(mk_form(y, vars, k), data = d2, weights = .w, method = "REML", select = select)
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
  ns <- dim(sims)[3]
  res <- vector("list", ns)
  for (i in seq_len(ns)) res[[i]] <- metrics_mat(curves_from_B(sims[,,i]))
  arr <- lapply(MET, function(mm) sapply(res, function(z) z[[mm]]))
  names(arr) <- MET
  arr
}

## unit-weighted skill (weight each fold by its number of held-out arrays) + cluster
## bootstrap over folds (matches step4_summary.R RULE 5)
uw <- function(w, m, mn) {
  ok <- is.finite(m) & is.finite(mn) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  1 - sum(w[ok]*m[ok])/sum(w[ok]*mn[ok])
}
boot_ci <- function(dt, wcol, mcol, ncol_, nb = 4000L, seed = 11) {
  set.seed(seed); n <- nrow(dt)
  if (n < 2) return(c(NA_real_, NA_real_))
  b <- replicate(nb, { i <- sample.int(n, n, TRUE)
                       uw(dt[[wcol]][i], dt[[mcol]][i], dt[[ncol_]][i]) })
  b <- b[is.finite(b)]
  if (!length(b)) return(c(NA_real_, NA_real_))
  as.numeric(quantile(b, c(0.025, 0.975)))
}

skill_of <- function(obs, pred, null, metric) {
  if (metric == "peak") {
    e  <- abs(circ_diff_h(pred, obs)); en <- abs(circ_diff_h(null, obs))
  } else { e <- abs(pred - obs); en <- abs(null - obs) }
  ok <- is.finite(e) & is.finite(en)
  list(mae = mean(e[ok]), mae_null = mean(en[ok]),
       skill = 1 - mean(e[ok])/mean(en[ok]), n = sum(ok))
}

## inverse-variance weighted Moran's I in a distance band (binary band weights),
## with a 200-shuffle permutation p-value. Ported unchanged from step10_resid_moran.R.
moran_band <- function(x, y, z, w, lo, hi, nperm = 200, seed = 5) {
  ok <- is.finite(z) & is.finite(w) & w > 0
  x <- x[ok]; y <- y[ok]; z <- z[ok]; w <- w[ok]
  n <- length(z); if (n < 20) return(c(NA_real_, NA_real_, n))
  zc <- z - sum(w*z)/sum(w)
  den <- sum(w*zc^2)
  Dmat <- as.matrix(dist(cbind(x, y)))
  band <- Dmat >= lo & Dmat < hi
  diag(band) <- FALSE
  Wm <- outer(w, w) * band
  S0 <- sum(Wm)
  if (S0 <= 0 || den <= 0) return(c(NA_real_, NA_real_, n))
  num <- sum(Wm * outer(zc, zc))
  I <- (n/S0) * (num/den)
  set.seed(seed); pm <- numeric(nperm)
  for (b in seq_len(nperm)) {
    ord <- sample.int(n)
    zp <- zc[ord]
    pm[b] <- (n/S0) * sum(Wm * outer(zp, zp)) / den
  }
  c(I, mean(abs(pm) >= abs(I)), n)
}

## ============================================================================
## UNIFIED MODEL EXTENSIONS (phase 3: H1 env + H4 predator + H5 own-abundance,
## one GAM per species x harmonic coefficient, select=TRUE shrinkage decides
## which of the candidate terms survive).
## ============================================================================

## fit one harmonic coefficient on an ARBITRARY variable list (env + predator + own-abundance
## terms together). Identical machinery to fit_coef() above, just generalised var naming.
fit_unified <- function(dat, y, vars, k = 4, select = TRUE) {
  w <- 1/pmax(dat[[paste0("v_", y)]], 1e-8); w <- w/mean(w)
  d2 <- copy(dat); d2$.w <- w
  gam(mk_form(y, vars, k), data = d2, weights = .w, method = "REML", select = select)
}

## Null-coefficient simulator for the artifact check (H4/H5 methodology, ported to the
## final_array-keyed, coefficient-level data we have -- no raw per-bin counts are available at
## this key, so the null draws each array's simulated harmonic coefficients from
## N(species inverse-variance-weighted mean coefficient, that array's OWN fitted coefficient
## variance). This preserves exactly the effort/precision structure (array-to-array variation in
## measurement noise, which is what a predator or abundance term could spuriously track) while
## erasing any true between-array biological signal -- anything a GAM "finds" here is an
## artifact of that shared noise structure, not biology.
simulate_null_coefs <- function(dat, seed) {
  set.seed(seed)
  d2 <- copy(dat)
  for (h in HARM) {
    w  <- 1/pmax(dat[[paste0("v_", h)]], 1e-8)
    mu <- sum(w * dat[[h]]) / sum(w)
    se <- sqrt(dat[[paste0("v_", h)]])
    d2[[h]] <- rnorm(nrow(dat), mean = mu, sd = se)
  }
  Bsim <- as.matrix(d2[, ..HARM])
  Msim <- metrics_mat(curves_from_B(Bsim))
  d2$obs_noct <- Msim$noct; d2$obs_crep <- Msim$crep; d2$obs_conc <- Msim$conc
  d2$obs_peak <- Msim$peak; d2$obs_act  <- Msim$act
  d2
}

## Average marginal effect of one predictor, holding all others at their observed values,
## moving the predictor from its 10th to 90th percentile. Point estimate only (nsim=0
## equivalent) -- used inside the null-artifact check where only the statistic (not its CI)
## is needed, many times over.
ame_point <- function(models, dat, var, mm) {
  x <- dat[[var]]
  qlo <- as.numeric(quantile(x, 0.10, na.rm = TRUE)); qhi <- as.numeric(quantile(x, 0.90, na.rm = TRUE))
  if (!is.finite(qlo) || !is.finite(qhi) || qhi <= qlo) return(NA_real_)
  nd_lo <- copy(dat); nd_lo[[var]] <- qlo
  nd_hi <- copy(dat); nd_hi[[var]] <- qhi
  pc_lo <- predict_curves(models, nd_lo, nsim = 0)
  pc_hi <- predict_curves(models, nd_hi, nsim = 0)
  resp_sd <- sd(dat[[paste0("obs_", mm)]], na.rm = TRUE)
  if (mm == "peak") { diff <- circ_diff_h(pc_hi$M[[mm]], pc_lo$M[[mm]]); scale_val <- 1 }
  else              { diff <- pc_hi$M[[mm]] - pc_lo$M[[mm]];            scale_val <- resp_sd }
  mean(diff, na.rm = TRUE) / scale_val
}

## Same AME, but with posterior-simulation draws for a bootstrap 95% CI -- on the SAME
## standardized scale as ame_point(): response-SD units for noct/crep/conc/act, hours
## (circular) for peak. This is the single consistent effect-size scale reported everywhere
## in the unified results table (stated explicitly in the report).
ame_effect_ci <- function(models, dat, var, mm, nsim = 500, seed = 1) {
  x <- dat[[var]]
  qlo <- as.numeric(quantile(x, 0.10, na.rm = TRUE)); qhi <- as.numeric(quantile(x, 0.90, na.rm = TRUE))
  if (!is.finite(qlo) || !is.finite(qhi) || qhi <= qlo)
    return(list(effect_std = NA_real_, ci_lo = NA_real_, ci_hi = NA_real_, q10 = qlo, q90 = qhi))
  nd_lo <- copy(dat); nd_lo[[var]] <- qlo
  nd_hi <- copy(dat); nd_hi[[var]] <- qhi
  pc_lo <- predict_curves(models, nd_lo, nsim = nsim, seed = seed)
  pc_hi <- predict_curves(models, nd_hi, nsim = nsim, seed = seed)
  resp_sd <- sd(dat[[paste0("obs_", mm)]], na.rm = TRUE)
  if (mm == "peak") { obs_diff <- circ_diff_h(pc_hi$M[[mm]], pc_lo$M[[mm]]); scale_val <- 1 }
  else              { obs_diff <- pc_hi$M[[mm]] - pc_lo$M[[mm]];            scale_val <- resp_sd }
  ame_obs <- mean(obs_diff, na.rm = TRUE) / scale_val
  sm_lo <- sim_metrics(pc_lo$sims)[[mm]]; sm_hi <- sim_metrics(pc_hi$sims)[[mm]]
  diffs <- if (mm == "peak") circ_diff_h(sm_hi, sm_lo) else (sm_hi - sm_lo)
  eff_draws <- colMeans(diffs, na.rm = TRUE) / scale_val
  ci <- quantile(eff_draws, c(0.025, 0.975), na.rm = TRUE)
  list(effect_std = ame_obs, ci_lo = unname(ci[1]), ci_hi = unname(ci[2]), q10 = qlo, q90 = qhi)
}

## Build a cache of NSIM_NULL null-simulated datasets + refitted unified models for one
## species x metric's variable set. Reused across every term's artifact check for that
## species so the (expensive) refits happen once, not once per term.
build_null_model_cache <- function(dat, all_vars, k = 4, nsim_null = 20, seed0 = 100) {
  lapply(seq_len(nsim_null), function(s) {
    dnull <- simulate_null_coefs(dat, seed = seed0 + s)
    mods_null <- tryCatch(lapply(setNames(HARM, HARM), function(h) fit_unified(dnull, h, all_vars, k = k)),
                           error = function(e) NULL)
    list(models = mods_null, dat = dnull)
  })
}

## Artifact check: does the observed AME exceed the 95th percentile of |AME| computed on
## NSIM_NULL null-simulated datasets (same predictor set, same model, no true signal)?
artifact_check <- function(null_cache, var, mm, obs_effect) {
  null_effects <- sapply(null_cache, function(nc) {
    if (is.null(nc$models)) return(NA_real_)
    tryCatch(ame_point(nc$models, nc$dat, var, mm), error = function(e) NA_real_)
  })
  null_effects <- null_effects[is.finite(null_effects)]
  if (!length(null_effects))
    return(list(null_mean = NA, null_sd = NA, null_p95 = NA, artifact_ratio = NA, survives_artifact_check = NA, n_null = 0))
  null_p95 <- unname(quantile(abs(null_effects), 0.95, na.rm = TRUE))
  ratio <- if (is.finite(obs_effect) && obs_effect != 0) mean(abs(null_effects)) / abs(obs_effect) else NA_real_
  list(null_mean = mean(null_effects), null_sd = sd(null_effects), null_p95 = null_p95,
       artifact_ratio = ratio,
       survives_artifact_check = is.finite(obs_effect) && is.finite(null_p95) && abs(obs_effect) > null_p95,
       n_null = length(null_effects))
}
