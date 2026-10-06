#!/usr/bin/env Rscript
# Population gradient in the effort-corrected iNat trend vs the camera trend, per species,
# union (v2u) national_scalar bundles. Reads existing data only; no refit. (2026-10-01)
#
# Unit for both streams = cell50 (camera rows are single-year deployments, so per-site trends
# are not estimable; deployments are pooled within their cell50).
#   iNat  : per cell, Poisson GLM  count ~ yr + offset(log(inat_effort)) over years with effort > 0
#           -> slope = trend in log(count/effort) per year.  Needs >= 5 effort-years and >= 5 records.
#   Camera: per cell, binomial GLM (cloglog) of detection-occasions / valid occasions ~ yr over the
#           cell's deployments -> slope = trend in log detection rate per year.
#           Needs deployments in >= 3 distinct years and >= 5 detections.
#   Cell Human_pop = model's standardised covariate (occ_beta is shared, so xdat_inat and occ_covars
#           are on one scale): iNat = mean of xdat_inat[,"Human_pop",] in the cell; camera = mean of
#           occ_covars[,"Human_pop"] over the cell's deployments.
#   Gradient = random-effects meta-regression of per-cell slopes on Human_pop
#           (weights 1/(se^2 + tau^2), tau^2 by method of moments), reported x100 as
#           percentage points of %/yr trend per +1 SD Human_pop. 95% CI by bootstrap over cells (1000).
#   Cells with non-converged fits or slope SE > 1 (separation) are dropped and counted.
suppressMessages({library(data.table); library(terra)})
PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
H <- file.path(PROJ, "HPC"); OUT <- "/home/rwkays/isdm/fleet_v2u"
g50 <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/grid50.tif")
set.seed(20261001)

FIT_ERR <- new.env(); FIT_ERR$n <- 0L; FIT_ERR$first <- NULL
slope_fit <- function(fml, data, fam) {
  f <- tryCatch(suppressWarnings(glm(fml, family = fam, data = data)), error = function(e) {
    FIT_ERR$n <- FIT_ERR$n + 1L; if (is.null(FIT_ERR$first)) FIT_ERR$first <- conditionMessage(e); NULL })
  if (is.null(f) || !f$converged) return(c(NA_real_, NA_real_))
  s <- summary(f)$coefficients
  if (!"yr" %in% rownames(s)) return(c(NA_real_, NA_real_))
  as.numeric(c(s["yr", 1], s["yr", 2]))
}
metareg <- function(b, se, x) {           # returns gradient (slope on x) and tau2
  X <- cbind(1, x); w <- 1 / se^2
  fitw <- function(w) solve(crossprod(X, w * X), crossprod(X, w * b))
  beta <- fitw(w); Q <- sum(w * (b - X %*% beta)^2); k <- length(b); p <- 2
  P <- diag(w) - (w * X) %*% solve(crossprod(X, w * X)) %*% t(w * X)
  tau2 <- max(0, (Q - (k - p)) / sum(diag(P)))
  w2 <- 1 / (se^2 + tau2); beta2 <- fitw(w2)
  c(grad = beta2[2], tau2 = tau2)
}
metareg_fast <- function(b, se, x) {      # same estimator without the k x k matrix (for bootstrap)
  X <- cbind(1, x); w <- 1 / se^2
  A <- crossprod(X, w * X); beta <- solve(A, crossprod(X, w * b))
  Q <- sum(w * (b - X %*% beta)^2); k <- length(b)
  trP <- sum(w) - sum(diag(solve(A, crossprod(X, w^2 * X))))
  tau2 <- max(0, (Q - (k - 2)) / trP); w2 <- 1 / (se^2 + tau2)
  solve(crossprod(X, w2 * X), crossprod(X, w2 * b))[2]
}
grad_ci <- function(d) {
  if (nrow(d) < 10) return(c(NA_real_, NA_real_, NA_real_))
  est <- metareg_fast(d$b, d$se, d$hp)
  bs <- replicate(1000, { i <- sample.int(nrow(d), replace = TRUE); metareg_fast(d$b[i], d$se[i], d$hp[i]) })
  100 * c(est = est, lo = unname(quantile(bs, .025)), hi = unname(quantile(bs, .975)))
}
fmt <- function(v) if (is.na(v[1])) "too few cells (<10)" else sprintf("%+.2f [%+.2f, %+.2f]", v[1], v[2], v[3])

res <- list(); cellrows <- list()
for (sp in c("moose", "bobcat", "white-tailed_deer", "coyote")) {
  tok <- paste0(sp, "_v2u_national_scalar")
  b <- readRDS(file.path(H, tok, paste0("input_data_", tok, ".RDS"))); cl <- b$constants_list
  Y <- b$inat_y; E <- b$inat_effort; cells <- as.integer(rownames(Y)); yrs <- as.integer(colnames(Y))
  stopifnot(identical(dimnames(Y), dimnames(E)))
  G <- cl$ncell50; row_g <- integer(dim(cl$xdat_inat)[1])
  for (g in 1:G) row_g[cl$inat_cell50_start[g]:cl$inat_cell50_end[g]] <- g
  hp_inat <- as.vector(tapply(rowMeans(cl$xdat_inat[, "Human_pop", ]), row_g, mean))

  # ---- iNat per-cell corrected trend ----
  ib <- matrix(NA, G, 2)
  for (g in 1:G) {
    ok <- !is.na(Y[g, ]) & !is.na(E[g, ]) & E[g, ] > 0
    if (sum(ok) < 5 || sum(Y[g, ok]) < 5) next
    ib[g, ] <- slope_fit(cnt ~ yr + offset(log(eff)), data.frame(cnt = Y[g, ok], eff = E[g, ok], yr = yrs[ok] - 2016), poisson)
  }
  inat <- data.table(cell50 = cells, b = ib[, 1], se = ib[, 2], hp = hp_inat)
  n_inat_try <- sum(!is.na(ib[, 1]) | TRUE)
  inat_ok <- inat[!is.na(b) & se < 1]

  # ---- camera per-cell trend ----
  oc <- b$real_data$occ.covs; y <- b$real_data$y
  xy <- project(as.matrix(oc[, c("longitude", "latitude")]), from = "EPSG:4326", to = crs(g50))
  site_c50 <- terra::extract(g50, xy)[, 1]
  ndet <- rowSums(y > 0, na.rm = TRUE); nocc <- rowSums(!is.na(y))
  cam_sites <- data.table(cell50 = site_c50, yr = as.integer(oc$year) - 2016, det = ndet, occ = nocc,
                          hp = cl$occ_covars[, "Human_pop"])[!is.na(cell50) & cell50 %in% cells & occ > 0]
  cb <- cam_sites[, {
    if (uniqueN(yr) < 3 || sum(det) < 5) list(b = NA_real_, se = NA_real_, hp = mean(hp), nsite = .N)
    else { s <- slope_fit(cbind(det, occ - det) ~ yr, .SD, binomial(link = "cloglog")); list(b = s[1], se = s[2], hp = mean(hp), nsite = .N) }
  }, by = cell50]
  cam_ok <- cb[!is.na(b) & se < 1]

  # ---- gradients ----
  gi_all <- grad_ci(inat_ok)
  inat_m <- inat_ok[cell50 %in% cam_ok$cell50]
  gi_match <- grad_ci(inat_m)
  gc <- grad_ci(cam_ok)
  cat(sprintf("\n==== %s\n  iNat: %d cells estimated (of %d), %d dropped for SE>1; matched-to-camera %d\n  camera: %d cells estimated from %d deployments (of %d cells with any deployment; %d deployments in footprint)\n",
              sp, nrow(inat_ok), G, sum(!is.na(inat$b) & inat$se >= 1), nrow(inat_m),
              nrow(cam_ok), sum(cam_ok$nsite), nrow(cb), nrow(cam_sites)))
  cat("  gradient (pp of %/yr per +1 SD Human_pop):  iNat all", fmt(gi_all), "| iNat matched", fmt(gi_match), "| camera", fmt(gc), "\n")
  cat("  glm errors so far:", FIT_ERR$n, if (!is.null(FIT_ERR$first)) paste("first:", FIT_ERR$first) else "", "\n")
  if (nrow(inat_ok) == 0 && nrow(cam_ok) == 0) stop("no per-cell fits succeeded for ", sp, " -- code fault, not data")
  res[[sp]] <- data.table(species = sp,
    inat_gradient_all = fmt(gi_all), inat_cells_all = nrow(inat_ok),
    inat_gradient_matched = fmt(gi_match), inat_cells_matched = nrow(inat_m),
    camera_gradient = fmt(gc), camera_cells = nrow(cam_ok), camera_deployments = sum(cam_ok$nsite),
    inat_mean_trend_pct = round(100 * (exp(weighted.mean(inat_ok$b, 1 / inat_ok$se^2)) - 1), 2),
    camera_mean_trend_pct = if (nrow(cam_ok)) round(100 * (exp(weighted.mean(cam_ok$b, 1 / cam_ok$se^2)) - 1), 2) else NA_real_,
    hp_sd_inat_cells = round(sd(inat_ok$hp), 3), hp_sd_camera_cells = round(sd(cam_ok$hp), 3))
  cellrows[[sp]] <- rbind(inat_ok[, .(species = sp, stream = "inat", cell50, slope = b, se, human_pop = hp, n_deploy = NA_integer_)],
                          cam_ok[, .(species = sp, stream = "camera", cell50, slope = b, se, human_pop = hp, n_deploy = nsite)])
  rm(b, cl, Y, E, y, oc); invisible(gc())
}
res <- rbindlist(res); options(width = 250); print(res)
fwrite(res, file.path(OUT, "pop_gradient_summary.csv"))
fwrite(rbindlist(cellrows), file.path(OUT, "pop_gradient_percell.csv"))
cat("POP GRADIENT DONE\n")
