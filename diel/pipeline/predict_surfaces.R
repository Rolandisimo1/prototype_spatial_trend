#!/usr/bin/env Rscript
# Population-level (all random effects excluded) predictions on the 100 km viewer grid,
# with posterior-simulation credible intervals, for five diel metrics.
suppressPackageStartupMessages({library(mgcv); library(data.table)})

ARGS <- commandArgs(trailingOnly = TRUE)
SPEC <- ARGS[1]; RDS <- ARGS[2]; OUTPFX <- ARGS[3]
NDRAW <- as.integer(Sys.getenv("NDRAW", "500"))
set.seed(20240815)

SUNRISE <- 6.2202; SUNSET <- 17.9998; CREP_H <- 1.5

PATHS <- jsonlite::fromJSON("work/paths.json")
cen  <- fread("work/grid_centroids_5070km.csv")
mask <- fread(PATHS$range_mask)
l48  <- fread(PATHS$l48)$cell_id
SCI <- c("White-tailed Deer"="Odocoileus virginianus","Northern Raccoon"="Procyon lotor",
         "Eastern Gray Squirrel"="Sciurus carolinensis","Coyote"="Canis latrans",
         "American Black Bear"="Ursus americanus")
keep <- intersect(intersect(mask[species == SCI[[SPEC]], unique(cell_id)], l48), cen$cell_id)
cen <- cen[cell_id %in% keep][order(cell_id)]
nC <- nrow(cen); stopifnot(nC > 0)

m <- readRDS(RDS)

## ---- RULE 2: derive random-effect labels from the fitted smooths, never hardcode ----
labs <- vapply(m$smooth, function(s) s$label, character(1))
re_idx <- which(vapply(m$smooth, function(s) inherits(s, "random.effect"), logical(1)))
re_labs <- labs[re_idx]
expected <- c("s(array_id)", "s(sin1,array_id)", "s(cos1,array_id)", "s(sin2,array_id)",
              "s(cos2,array_id)")
stopifnot(setequal(re_labs, expected))                 # every label matches a fitted smooth
fixed_idx <- setdiff(seq_along(m$smooth), re_idx)
stopifnot(identical(labs[fixed_idx],
                    c("s(bin_center_h)", "s(X,Y)", "ti(bin_center_h,X,Y)")))
# coefficient block for intercept + the three population smooths; must be contiguous from 1
p_last <- max(vapply(m$smooth[fixed_idx], function(s) s$last.para, numeric(1)))
p_first <- min(vapply(m$smooth[fixed_idx], function(s) s$first.para, numeric(1)))
stopifnot(p_first == 2L)                               # 1 = intercept
re_first <- min(vapply(m$smooth[re_idx], function(s) s$first.para, numeric(1)))
stopifnot(re_first == p_last + 1L)                     # all RE coefs sit strictly after
KP <- seq_len(p_last)
cat(sprintf("[%s] fixed coefs 1:%d of %d; RE blocks %s\n", SPEC, p_last, length(coef(m)),
            paste(re_labs, collapse = " ")))

## ---- design matrix: 48 sun-anchored bins x cells, all RE columns absent by construction ----
bins <- seq(0.25, 23.75, by = 0.5); nB <- length(bins)
nd <- data.table(cell_id = rep(cen$cell_id, each = nB), X = rep(cen$X, each = nB),
                 Y = rep(cen$Y, each = nB), bin_center_h = rep(bins, nC))
ndf <- as.data.frame(nd)
Xp <- matrix(0, nrow = nrow(nd), ncol = p_last)
Xp[, 1] <- 1
for (i in fixed_idx) {
  s <- m$smooth[[i]]
  Xp[, s$first.para:s$last.para] <- PredictMat(s, ndf)
}
stopifnot(all(is.finite(Xp)))

## ---- posterior draws from Vp on the population block ----
B <- rmvn(NDRAW, coef(m)[KP], m$Vp[KP, KP, drop = FALSE])   # NDRAW x p_last
ETA <- Xp %*% t(B)                                          # (nB*nC) x NDRAW ; offset = 0
RATE <- exp(ETA); rm(ETA); gc()

## ---- weights for the sun-anchored day / crepuscular windows, fractional bin overlap ----
ov <- function(a, b) { lo <- pmax(bins - 0.25, a); hi <- pmin(bins + 0.25, b); pmax(hi - lo, 0) / 0.5 }
w_day  <- ov(SUNRISE, SUNSET)
w_crep <- ov(SUNRISE - CREP_H, SUNRISE + CREP_H) + ov(SUNSET - CREP_H, SUNSET + CREP_H) +
          ov(SUNRISE - CREP_H + 24, 24) + ov(0, SUNSET + CREP_H - 24)
stopifnot(all(w_day >= 0 & w_day <= 1), all(w_crep >= 0 & w_crep <= 1))
ang <- 2 * pi * bins / 24

metrics_of <- function(R) {                    # R: nB x nC matrix of rates
  s <- colSums(R); p <- sweep(R, 2, s, "/")
  pct_noct <- 100 * (1 - colSums(p * w_day))
  crep     <- 100 * colSums(p * w_crep)
  Cb <- colSums(p * cos(ang)); Sb <- colSums(p * sin(ang))
  conc <- sqrt(Cb^2 + Sb^2)
  # primary peak: argmax bin refined by circular parabolic interpolation
  k <- max.col(t(R), ties.method = "first")
  km <- ifelse(k == 1, nB, k - 1); kp <- ifelse(k == nB, 1, k + 1)
  ii <- seq_len(ncol(R))
  y0 <- R[cbind(km, ii)]; y1 <- R[cbind(k, ii)]; y2 <- R[cbind(kp, ii)]
  den <- (y0 - 2 * y1 + y2)
  d <- ifelse(abs(den) < 1e-12, 0, 0.5 * (y0 - y2) / den)
  d <- pmax(pmin(d, 0.5), -0.5)
  tpk <- (bins[k] + d * 0.5) %% 24
  peak <- (tpk - SUNRISE) %% 24
  act <- colMeans(R) / apply(R, 2, max)
  list(pct_nocturnal = pct_noct, crepuscular = crep, concentration = conc,
       peak_hour = peak, activity_level = act)
}

MK <- c("pct_nocturnal", "crepuscular", "concentration", "peak_hour", "activity_level")
draws <- lapply(MK, function(k) matrix(NA_real_, nC, NDRAW)); names(draws) <- MK
for (d in seq_len(NDRAW)) {
  mm <- metrics_of(matrix(RATE[, d], nB, nC))
  for (k in MK) draws[[k]][, d] <- mm[[k]]
}

circ_summ <- function(M) {                      # rows = cells, cols = draws, values in [0,24)
  a <- 2 * pi * M / 24
  mu <- atan2(rowMeans(sin(a)), rowMeans(cos(a)))
  dev <- ((a - mu + pi) %% (2 * pi)) - pi
  q <- t(apply(dev, 1, quantile, probs = c(0.025, 0.975), names = FALSE))
  hr <- function(x) (x * 24 / (2 * pi)) %% 24
  list(value = hr(mu), lo = hr(mu + q[, 1]), hi = hr(mu + q[, 2]),
       width = (q[, 2] - q[, 1]) * 24 / (2 * pi))
}

surf <- rbindlist(lapply(MK, function(k) {
  M <- draws[[k]]
  if (k == "peak_hour") {
    s <- circ_summ(M)
    data.table(species = SPEC, metric = k, cell_id = cen$cell_id, value = s$value,
               ci_lo = s$lo, ci_hi = s$hi, ci_width = s$width, in_mask = TRUE)
  } else {
    q <- t(apply(M, 1, quantile, probs = c(0.025, 0.975), names = FALSE))
    data.table(species = SPEC, metric = k, cell_id = cen$cell_id, value = rowMeans(M),
               ci_lo = q[, 1], ci_hi = q[, 2], ci_width = q[, 2] - q[, 1], in_mask = TRUE)
  }
}))
fwrite(surf, sprintf("%s_surfaces.csv", OUTPFX))

## ---- per-cell curves: posterior-mean rate, normalised to cell mean = 1 ----
Rmean <- matrix(rowMeans(RATE), nB, nC)
Rlo   <- matrix(apply(RATE, 1, quantile, 0.025, names = FALSE), nB, nC)
Rhi   <- matrix(apply(RATE, 1, quantile, 0.975, names = FALSE), nB, nC)
sc <- colMeans(Rmean)
Nm <- sweep(Rmean, 2, sc, "/"); Nlo <- sweep(Rlo, 2, sc, "/"); Nhi <- sweep(Rhi, 2, sc, "/")
ci_rel <- colMeans((Nhi - Nlo) / (2 * pmax(Nm, 1e-9)))
cc <- data.table(species = SPEC, cell_id = rep(cen$cell_id, each = nB),
                 bin = rep(seq_len(nB) - 1L, nC), sun_hour = rep(bins, nC),
                 rate = as.vector(Nm), ci_rel = rep(ci_rel, each = nB))
fwrite(cc, sprintf("%s_cellcurves.csv", OUTPFX))

## ---- featured curves: greedy max-min distance on shape-normalised curves ----
P <- sweep(Nm, 2, colSums(Nm), "/")
sel <- which.max(apply(P, 2, function(v) max(v) - min(v)))
for (j in 2:5) {
  dmin <- apply(P, 2, function(v) min(sqrt(colSums((P[, sel, drop = FALSE] - v)^2))))
  dmin[sel] <- -Inf
  sel <- c(sel, which.max(dmin))
}
cur <- rbindlist(lapply(sel, function(j) {
  lab <- sprintf("%s (%.1f\u00b0N, %.1f\u00b0W)", cen$cell_id[j], cen$lat[j], -cen$lon[j])
  data.table(species = SPEC, location = lab, bin = seq_len(nB) - 1L, sun_hour = bins,
             rate = Nm[, j], lo = Nlo[, j], hi = Nhi[, j])
}))
fwrite(cur, sprintf("%s_curves.csv", OUTPFX))

## ---- sanity gate + activity-level summary ----
gate <- data.table(species = SPEC, n_cells = nC,
                   fitted_mean_pct_noct = mean(surf[metric == "pct_nocturnal", value]),
                   activity_level_mean = mean(surf[metric == "activity_level", value]),
                   activity_level_min = min(surf[metric == "activity_level", value]),
                   activity_level_max = max(surf[metric == "activity_level", value]))
fwrite(gate, sprintf("%s_gate.csv", OUTPFX))
cat(sprintf("[%s] cells=%d fitted_mean_pct_noct=%.2f act_level=%.3f\n",
            SPEC, nC, gate$fitted_mean_pct_noct, gate$activity_level_mean))
