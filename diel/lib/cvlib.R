## cvlib.R — shared machinery for Phase 1 track 2 (validation)
suppressPackageStartupMessages({
  library(data.table); library(arrow); library(mgcv); library(sf)
})

REF_SR <- 6.2202; REF_SS <- 17.9998        # sun-time anchors (phase 0)
KN     <- list(bin_center_rad = c(0, 2*pi))
AEA    <- "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=37.5 +lon_0=-96 +datum=WGS84 +units=km"
NBIN   <- 48

## ---- model specification (matches the spec described in the project context) ----
## tensor space x time-of-day + array random diel curves (intercept + 1st/2nd harmonic
## random slopes) + deployment random intercept + offset(log_effort_aso), nb()
f_spatial <- count ~ te(bin_center_rad, x_km, y_km, bs = c("cc","ds"), d = c(1,2),
                        k = c(10, 25), xt = list(list(), list(max.knots = 1500, seed = 11))) +
  offset(log_effort_aso) +
  s(array_id, bs = "re") +
  s(array_id, s1, bs = "re") + s(array_id, c1, bs = "re") +
  s(array_id, s2, bs = "re") + s(array_id, c2, bs = "re") +
  s(dep_key, bs = "re")

## null baseline: one continental diel curve, no spatial term at all
f_null <- count ~ s(bin_center_rad, bs = "cc", k = 10) + offset(log_effort_aso) + s(dep_key, bs = "re")

## naive deployment-only variant (for the effective-sample-size / CI comparison)
f_naive <- count ~ te(bin_center_rad, x_km, y_km, bs = c("cc","ds"), d = c(1,2),
                      k = c(10, 25), xt = list(list(), list(max.knots = 1500, seed = 11))) +
  offset(log_effort_aso) + s(dep_key, bs = "re")

SPEC_STR <- paste(
  "count ~ te(bin_center_rad, x_km, y_km, bs=c('cc','ds'), d=c(1,2), k=c(10,25))",
  "+ s(array_id, bs='re') + s(array_id, s1, bs='re') + s(array_id, c1, bs='re')",
  "+ s(array_id, s2, bs='re') + s(array_id, c2, bs='re') + s(dep_key, bs='re')",
  "+ offset(log_effort_aso); family=nb(); bam(discrete=TRUE, method='fREML')")

## ---- data prep ----
load_species <- function(d, sp) {
  x <- d[species == sp & aso_window == TRUE & is.finite(log_effort_aso) & effort_h_aso > 0]
  x[, `:=`(s1 = sin(bin_center_rad), c1 = cos(bin_center_rad),
           s2 = sin(2*bin_center_rad), c2 = cos(2*bin_center_rad))]
  p  <- sf::st_as_sf(unique(x[, .(dep_key, longitude, latitude)]),
                     coords = c("longitude","latitude"), crs = 4326)
  co <- sf::st_coordinates(sf::st_transform(p, AEA))
  x  <- merge(x, data.table(dep_key = p$dep_key, x_km = co[,1], y_km = co[,2]), by = "dep_key")
  x[, `:=`(dep_key = factor(dep_key), array_id = factor(array_id))]
  x[]
}

## thin cameras within arrays: arrays are the replicate unit for spatial inference,
## so capping cameras/array keeps ~all spatial information at a fraction of the
## coefficient count (p is dominated by the deployment random effect).
thin_deps <- function(x, cap = 8, seed = 42) {
  set.seed(seed)
  keys <- unique(x[, .(dep_key, array_id, count_tot = 0)])[, count_tot := NULL]
  tot  <- x[, .(tot = sum(count)), by = dep_key]
  keys <- merge(keys, tot, by = "dep_key")
  keys[, r := runif(.N)]
  ## keep the cameras with detections preferentially, then fill at random
  keys[, ord := order(-(tot > 0), r), by = array_id]
  sel  <- keys[, .SD[ord <= cap], by = array_id]$dep_key
  y <- x[dep_key %in% sel]
  y[, `:=`(dep_key = droplevels(dep_key), array_id = droplevels(array_id))]
  y[]
}

## ---- fitting ----
fit_bam <- function(form, dat, sp = NULL, theta = NULL, nthreads = 7) {
  fam <- if (is.null(theta)) nb() else nb(theta = theta)
  ## HARD GUARD: every fit must carry the exact per-bin sun-time effort offset
  stopifnot("offset(log_effort_aso) missing from model formula" =
              grepl("offset(log_effort_aso)", paste(deparse(form), collapse=""), fixed = TRUE))
  stopifnot(all(is.finite(dat$log_effort_aso)))
  bam(form, family = fam, data = dat, knots = KN,
      sp = sp, discrete = TRUE, method = "fREML", nthreads = nthreads,
      gc.level = 1, select = FALSE)
}

## ---- diel-shape metrics on a 48-bin curve ----
## normalized curve (sums to 1); night = sun-bins outside [REF_SR, REF_SS]
bc      <- (seq_len(NBIN) - 0.5) * 0.5
is_night <- bc < REF_SR | bc >= REF_SS
is_crep  <- (abs(bc - REF_SR) <= 1.5) | (abs(bc - REF_SS) <= 1.5)
ang      <- bc / 24 * 2*pi

curve_metrics <- function(v) {
  v <- pmax(as.numeric(v), 0)
  s <- sum(v)
  if (!is.finite(s) || s <= 0) return(list(pct_noct = NA_real_, pct_crep = NA_real_,
                                           conc = NA_real_, peak_h = NA_real_))
  p <- v / s
  C <- sum(p * cos(ang)); S <- sum(p * sin(ang))
  R <- sqrt(C^2 + S^2)
  mu <- (atan2(S, C) %% (2*pi)) / (2*pi) * 24
  list(pct_noct = 100 * sum(p[is_night]),
       pct_crep = 100 * sum(p[is_crep]),
       conc     = R,                      # 1 - circular variance
       peak_h   = bc[which.max(v)],
       mean_h   = mu)
}

circ_diff_h <- function(a, b) {           # signed-magnitude circular difference in hours
  d <- (a - b) %% 24
  ifelse(d > 12, d - 24, d)
}

## ---- random-effect labels, READ FROM THE FITTED MODEL ----
## mgcv relabels s(array_id, s1, bs="re") as "s(s1,array_id)". Hardcoding the
## as-written label makes predict()'s `exclude` silently no-op, which leaks the
## reference array's harmonic deviations into every prediction and corrupts the
## very curve shape this validation measures. Always derive labels, then assert.
re_labels <- function(m) {
  labs <- vapply(m$smooth, function(s) s$label, character(1))
  hit  <- grepl("array_id|dep_key", labs)
  stopifnot("no random-effect smooths found" = any(hit))
  labs[hit]
}
assert_excluded <- function(m, excl) {
  labs <- vapply(m$smooth, function(s) s$label, character(1))
  miss <- setdiff(excl, labs)
  stopifnot("exclude labels do not match fitted smooths" = length(miss) == 0)
  invisible(TRUE)
}
## surface-only prediction on new data: all grouping terms excluded, verified
pred_surface <- function(m, newd, tr) {
  excl <- re_labels(m); assert_excluded(m, excl)
  nd <- copy(newd)
  if ("array_id" %in% names(nd)) nd$array_id <- factor(levels(tr$array_id)[1], levels = levels(tr$array_id))
  nd$dep_key <- factor(levels(tr$dep_key)[1], levels = levels(tr$dep_key))
  p <- as.numeric(predict(m, nd, type = "response", exclude = excl, discrete = FALSE))
  ## independence check: prediction must not depend on which reference level we plug in
  if (length(levels(tr$array_id)) > 1) {
    nd2 <- copy(nd); nd2$array_id <- factor(levels(tr$array_id)[2], levels = levels(tr$array_id))
    p2 <- as.numeric(predict(m, nd2, type = "response", exclude = excl, discrete = FALSE))
    stopifnot("exclude failed: prediction depends on reference array level" =
                isTRUE(all.equal(p, p2, tolerance = 1e-8)))
  }
  p
}
