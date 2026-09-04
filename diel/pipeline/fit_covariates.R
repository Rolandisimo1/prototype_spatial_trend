#!/usr/bin/env Rscript
# Refit per-species covariate models to recover FITTED covariate effect curves for the
# viewer panels. PI-mandated covariate set exactly: ntl_1km (urbanisation, NOT impervious),
# nlcd_1k_crop and nlcd_1k_pasture kept SEPARATE from each other and from urbanisation,
# plus pop_1km, tcc_1km, tmax_window, elev. No composite human-footprint index anywhere.
suppressPackageStartupMessages({library(mgcv); library(data.table); library(arrow); library(jsonlite)})

ARGS <- commandArgs(trailingOnly = TRUE); SPEC <- ARGS[1]; OUTPFX <- ARGS[2]
NDRAW <- as.integer(Sys.getenv("NDRAW", "400"))
set.seed(20240815)
P <- fromJSON("work/paths.json")
SUNRISE <- 6.2202; SUNSET <- 17.9998; CREP_H <- 1.5
bins <- seq(0.25, 23.75, by = 0.5); nB <- 48L

COVS <- c("ntl_1km", "nlcd_1k_crop", "nlcd_1k_pasture", "pop_1km", "tcc_1km", "tmax_window", "elev")
LOGC <- c("ntl_1km", "pop_1km")            # heavy right skew -> model on log1p, report on raw scale
stopifnot(!("nlcd_1k_impervious" %in% COVS))   # PI: impervious must NOT be fitted

## ---- data: binned counts with the sun-time effort offset (RULE 1) ----
st <- as.data.table(read_parquet(P$suntime,
        col_select = c("dep_key", "species", "bin", "bin_center_h", "count",
                       "effort_h_aso", "log_effort_aso", "array_id")))
st <- st[species == SPEC]
dc <- fread(P$depcov, select = c("deployment_id", "project_id", COVS))
dc[, dep_key := paste0(project_id, "|", deployment_id)]
d <- merge(st, dc[, c("dep_key", COVS), with = FALSE], by = "dep_key")
d <- d[complete.cases(d[, COVS, with = FALSE]) & is.finite(log_effort_aso) & effort_h_aso > 0]
for (cv in LOGC) d[[paste0(cv, "_t")]] <- log1p(d[[cv]])
for (cv in setdiff(COVS, LOGC)) d[[paste0(cv, "_t")]] <- d[[cv]]
d[, array_id := factor(array_id)]
cat(sprintf("[%s] rows=%d deps=%d arrays=%d\n", SPEC, nrow(d), uniqueN(d$dep_key), nlevels(d$array_id)))

## ---- model: diel shape allowed to vary with each covariate SEPARATELY ----
tv <- paste0(COVS, "_t")
f <- as.formula(paste0(
  "count ~ s(bin_center_h, bs='cc', k=10) + ",
  paste(sprintf("s(%s, k=5)", tv), collapse = " + "), " + ",
  paste(sprintf("ti(bin_center_h, %s, bs=c('cc','tp'), d=c(1,1), k=c(10,5))", tv), collapse = " + "),
  " + s(array_id, bs='re')"))
t0 <- Sys.time()
m <- bam(f, family = nb(), data = d, offset = d$log_effort_aso, discrete = TRUE,
         nthreads = 6, knots = list(bin_center_h = c(0, 24)))
cat(sprintf("[%s] fitted in %.1f min, dev.expl=%.3f\n", SPEC,
            as.numeric(difftime(Sys.time(), t0, units = "mins")), summary(m)$dev.expl))

## ---- population block: exclude the array random effect by coefficient index (RULE 2/3) ----
labs <- vapply(m$smooth, function(s) s$label, character(1))
re_idx <- which(vapply(m$smooth, function(s) inherits(s, "random.effect"), logical(1)))
stopifnot(identical(labs[re_idx], "s(array_id)"))
p_last <- min(vapply(m$smooth[re_idx], function(s) s$first.para, numeric(1))) - 1L
KP <- seq_len(p_last)
fixed_idx <- setdiff(seq_along(m$smooth), re_idx)
stopifnot(max(vapply(m$smooth[fixed_idx], function(s) s$last.para, numeric(1))) == p_last)

B <- rmvn(NDRAW, coef(m)[KP], m$Vp[KP, KP, drop = FALSE])

ov <- function(a, b) { lo <- pmax(bins - 0.25, a); hi <- pmin(bins + 0.25, b); pmax(hi - lo, 0) / 0.5 }
w_day  <- ov(SUNRISE, SUNSET)
w_crep <- ov(SUNRISE - CREP_H, SUNRISE + CREP_H) + ov(SUNSET - CREP_H, SUNSET + CREP_H) +
          ov(SUNRISE - CREP_H + 24, 24) + ov(0, SUNSET + CREP_H - 24)
ang <- 2 * pi * bins / 24
metrics_of <- function(R) {
  p <- sweep(R, 2, colSums(R), "/")
  Cb <- colSums(p * cos(ang)); Sb <- colSums(p * sin(ang))
  k <- max.col(t(R), ties.method = "first")
  list(pct_nocturnal = 100 * (1 - colSums(p * w_day)), crepuscular = 100 * colSums(p * w_crep),
       concentration = sqrt(Cb^2 + Sb^2), peak_hour = (bins[k] - SUNRISE) %% 24,
       activity_level = colMeans(R) / apply(R, 2, max))
}
MK <- c("pct_nocturnal", "crepuscular", "concentration", "peak_hour", "activity_level")

NX <- 21L
med <- sapply(tv, function(v) median(d[[v]]))
out <- list()
for (ci in seq_along(COVS)) {
  cv <- COVS[ci]; tvv <- tv[ci]
  qs <- quantile(d[[tvv]], c(0.10, 0.90), names = FALSE)
  if (!(qs[2] > qs[1])) next
  xg <- seq(qs[1], qs[2], length.out = NX)
  xraw <- if (cv %in% LOGC) expm1(xg) else xg
  nd <- as.data.frame(matrix(rep(med, each = NX * nB), NX * nB, length(tv)))
  names(nd) <- tv
  nd[[tvv]] <- rep(xg, each = nB)
  nd$bin_center_h <- rep(bins, NX)
  nd$array_id <- d$array_id[1]                       # dummy; RE columns are dropped below
  Xp <- predict(m, nd, type = "lpmatrix")[, KP, drop = FALSE]
  ETA <- Xp %*% t(B); RATE <- exp(ETA)
  dr <- lapply(MK, function(k) matrix(NA_real_, NX, NDRAW)); names(dr) <- MK
  for (j in seq_len(NDRAW)) {
    mm <- metrics_of(matrix(RATE[, j], nB, NX))
    for (k in MK) dr[[k]][, j] <- mm[[k]]
  }
  for (k in MK) {
    M <- dr[[k]]; q <- t(apply(M, 1, quantile, probs = c(0.05, 0.95), names = FALSE))
    out[[length(out) + 1]] <- data.table(species = SPEC, metric = k, covariate = cv,
      x = xraw, y = rowMeans(M), lo = q[, 1], hi = q[, 2])
  }
  cat(sprintf("  [%s] %s done\n", SPEC, cv))
}
fwrite(rbindlist(out), sprintf("%s_covpanels.csv", OUTPFX))
fwrite(data.table(species = SPEC, dev_expl = summary(m)$dev.expl, n_rows = nrow(d),
                  n_arrays = nlevels(d$array_id)), sprintf("%s_covfit.csv", OUTPFX))
