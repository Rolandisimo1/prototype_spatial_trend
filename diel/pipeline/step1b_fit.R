## step1b_fit.R -- fit 5 harmonic coefficients per species x array + observed metrics with SEs.
suppressPackageStartupMessages({library(data.table)})
source("cvlib.R")
S <- readRDS("s1_agg.rds"); Xd <- S$Xd; bcen <- S$bcen
THR <- 30L; NSIM <- 1000L
act_level <- function(v) { v <- pmax(v,0); m <- max(v); if (!is.finite(m) || m <= 0) NA_real_ else mean(v)/m }
all_metrics <- function(v) { m <- curve_metrics(v); m$act <- act_level(v); m }

fit_one <- function(cnt, eff) {
  ## quasipoisson log-link with offset; identical linear predictor for all arrays
  fm <- try(suppressWarnings(glm(cnt ~ 0 + Xd, family = quasipoisson(),
                                 offset = log(eff))), silent = TRUE)
  if (inherits(fm,"try-error") || !fm$converged) return(NULL)
  b  <- coef(fm); if (any(!is.finite(b))) return(NULL)
  V  <- vcov(fm)
  ## dispersion clamped at >=1: quasipoisson can report phi<1 for smooth counts, which
  ## would make the SEs anticonservative relative to the Poisson floor.
  phi <- summary(fm)$dispersion
  if (is.finite(phi) && phi < 1) V <- V / phi
  list(b = b, V = V, phi = phi)
}

## vectorised metrics over a 48 x nsim matrix of rate curves (columns = draws)
ang48 <- bcen / 24 * 2*pi
NIGHT <- bcen < REF_SR | bcen >= REF_SS
CREP  <- (abs(bcen - REF_SR) <= 1.5) | (abs(bcen - REF_SS) <= 1.5)
metrics_mat <- function(M) {
  M <- pmax(M, 0); s <- colSums(M)
  P <- sweep(M, 2, ifelse(s > 0, s, NA_real_), "/")
  C <- as.numeric(crossprod(cos(ang48), P)); S <- as.numeric(crossprod(sin(ang48), P))
  R <- sqrt(C^2 + S^2)
  list(noct = 100*as.numeric(crossprod(NIGHT, P)),
       crep = 100*as.numeric(crossprod(CREP,  P)),
       conc = R,
       peak = bcen[max.col(t(M), ties.method = "first")],
       mean = ((atan2(S, C)) %% (2*pi))/(2*pi)*24,
       act  = colMeans(M)/apply(M, 2, max))
}

metric_sim <- function(b, V, nsim = NSIM) {
  ## posterior simulation of the fitted rate curve -> metric point estimate + SE
  L  <- try(chol(V + diag(1e-10, ncol(V))), silent = TRUE)
  mu <- as.numeric(exp(Xd %*% b))
  m0 <- all_metrics(mu)
  if (inherits(L,"try-error")) return(list(est = m0, se = rep(NA_real_, 6)))
  Bs <- matrix(rnorm(nsim*length(b)), nsim) %*% L
  Bs <- sweep(Bs, 2, b, "+")
  mm <- metrics_mat(exp(Xd %*% t(Bs)))            # 48 x nsim
  M  <- cbind(mm$noct, mm$crep, mm$conc, mm$peak, mm$mean, mm$act)
  ## circular SD for the two timing metrics
  csd <- function(h) { a <- h/24*2*pi; R <- sqrt(mean(cos(a))^2 + mean(sin(a))^2)
                       sqrt(-2*log(pmin(pmax(R,1e-12),1)))/(2*pi)*24 }
  se <- c(sd(M[,1],na.rm=TRUE), sd(M[,2],na.rm=TRUE), sd(M[,3],na.rm=TRUE),
          csd(M[,4][is.finite(M[,4])]), csd(M[,5][is.finite(M[,5])]), sd(M[,6],na.rm=TRUE))
  list(est = m0, se = se)
}

set.seed(42)
res <- list()
for (SP in names(S$out)) {
  ag <- S$out[[SP]]$ag; ev <- S$out[[SP]]$ev
  keep <- ev[events >= THR, array_id]
  ag <- ag[array_id %in% keep]
  rows <- list()
  for (aid in keep) {
    z <- ag[array_id == aid]
    setorder(z, bin)
    if (nrow(z) != 48L) next                      # need the full 48-bin support
    f <- fit_one(z$count, z$eff)
    if (is.null(f)) next
    ms <- metric_sim(f$b, f$V)
    se <- sqrt(pmax(diag(f$V), 0))
    rows[[length(rows)+1]] <- data.table(
      species = SP, array_id = as.character(aid),
      events = sum(z$count), eff_h = sum(z$eff), ndep = z$ndep[1], phi = f$phi,
      b0 = f$b[1], s1 = f$b[2], c1 = f$b[3], s2 = f$b[4], c2 = f$b[5],
      se_b0 = se[1], se_s1 = se[2], se_c1 = se[3], se_s2 = se[4], se_c2 = se[5],
      v_b0 = f$V[1,1], v_s1 = f$V[2,2], v_c1 = f$V[3,3], v_s2 = f$V[4,4], v_c2 = f$V[5,5],
      obs_noct = ms$est$pct_noct, obs_crep = ms$est$pct_crep, obs_conc = ms$est$conc,
      obs_peak = ms$est$peak_h,  obs_mean = ms$est$mean_h,   obs_act  = ms$est$act,
      se_noct = ms$se[1], se_crep = ms$se[2], se_conc = ms$se[3],
      se_peak = ms$se[4], se_mean = ms$se[5], se_act = ms$se[6])
  }
  r <- rbindlist(rows)
  ## raw (model-free) metrics straight from the binned rates, as a fit check
  raw <- ag[, { o <- count/eff; m <- all_metrics(o)
                .(raw_noct = m$pct_noct, raw_crep = m$pct_crep, raw_conc = m$conc,
                  raw_peak = m$peak_h, raw_act = m$act) }, by = array_id]
  raw[, array_id := as.character(array_id)]
  r <- merge(r, raw, by = "array_id")
  res[[SP]] <- r
  cat(sprintf("[%s] fitted %d arrays | median phi %.2f | cor(fitted,raw) noct %.3f conc %.3f act %.3f\n",
      SP, nrow(r), median(r$phi, na.rm=TRUE),
      cor(r$obs_noct, r$raw_noct, use="complete.obs"),
      cor(r$obs_conc, r$raw_conc, use="complete.obs"),
      cor(r$obs_act,  r$raw_act,  use="complete.obs")))
}
H <- rbindlist(res)
H <- merge(H, S$arr_ll[, .(array_id = as.character(array_id), lon, lat, x_km, y_km, cell_id)],
           by = "array_id")
fwrite(H, "array_harmonics_raw.csv")
cat(sprintf("TOTAL rows %d\n", nrow(H)))
