#!/usr/bin/env Rscript
# =============================================================================
# screen_stream_gradient.R -- STANDARD per-species screen (required step, 2026-10-01): does either
# data stream's trend change along the human-population gradient, and do the streams agree?
# Reads the species' union (v2u) national_scalar bundle only; no refit. Usage: SPECIES=<token>.
#
# Covariate (same as the tier-2 urbantrend model): hp = the model's standardised Human_pop at the
# BASELINE year 2008, averaged over each cell100's iNat sub-cell rows, CENTRED on the footprint mean.
# One value per cell100, used by BOTH streams (camera deployment -> cell[i]; iNat cell50 -> inat_cell100).
#
# HEADLINE (per stratum: national + each ecoregion):
#   camera: glmer cbind(det, occ - det) ~ yr * hp + (1 + yr | cell50), cloglog
#   iNat  : glmer cnt ~ yr * hp + offset(log(eff)) + (1 + yr | cell50), Poisson
#   gradient = 100 * coef(yr:hp): percentage points of %/yr trend per +1 unit of hp (Wald 95% CI).
#   Random year-slopes per cell mean the gradient is estimated across cells' own trends, so a
#   single high-deployment cell cannot drive it (a random intercept alone does not prevent that).
# Tertile contrast and per-cell precision-weighted slopes are NOT used (camera-inat-stream-agreement
#   skill: both failed on real data -- coyote -0.6 / +3.8, bobcat +1.6 vs pooled +13.8).
# Fitting: yr in decades and hp in SD units internally (the raw scales made glmer fail with
#   'Downdated VtV is not positive definite'); the gradient is converted back to the reporting units.
#   If the full fit still fails, a nAGQ = 0 fit is used and labelled.
# CONCENTRATION (always): the ten cells with most camera deployments in the stratum -- deployments,
#   detecting deployments, years, projects, own trend -- and the camera headline refitted without them.
# AGREEMENT (camera-inat-stream-agreement section 2; thresholds set by rwkays 2026-10-01):
#   diff = iNat gradient - camera gradient, se = sqrt(se_cam^2 + se_inat^2) (streams are separate data)
#   OPPOSE        diff interval excludes 0 and the gradients have opposite signs
#   DIFFER        diff interval excludes 0, same sign (magnitude)
#   NOT ASSESSED  diff interval contains 0 AND reaches -10 or +10 gradient units (~6 pp/yr urban-rural,
#                 about the smallest difference that would matter). ABSOLUTE bar: an earlier relative bar
#                 (scaled to the larger gradient) gave larger estimates a laxer test.
#   AGREE         diff interval contains 0 and lies strictly inside (-10, +10)
# DEGENERATE FIT: a half-width < 2% of an |estimate| > 2 is a BUG SIGNAL (elk camera +18.76 +/- 0.03), not a
#   precise estimate -- printed as a loud warning naming the species and recorded in column fit_alarm.
# CONFOUNDS reported per stratum: projects and projects spanning >= 2 years (only those identify a
#   within-project trend); year coverage = Spearman corr(deployment year, hp) and median year of the
#   most- vs least-urban third of deployments.
# Region names come from prepped$ecoregion_levels (IDs are alphabetical factor levels), checked
# against each region's mean lon/lat.
# =============================================================================
suppressMessages({library(data.table); library(lme4)})
SP <- Sys.getenv("SPECIES"); stopifnot(nzchar(SP))
PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
H <- file.path(PROJ, "HPC"); OUT <- "/home/rwkays/isdm/fleet_v2u/stream_screen"; dir.create(OUT, showWarnings = FALSE)
# National stratum minimums; REGIONAL strata need more (thin regions gave degenerate glmer intervals,
# e.g. a camera gradient of -50.9 [-51.1, -50.6] from 42 detecting deployments, which the agreement
# rule would then have classified). MIN_HP_SD: across-cell spread of hp needed to estimate a gradient.
MIN_CAM_DEP <- 300; MIN_CAM_DET <- 40; MIN_CAM_CELLS <- 10; MIN_CAM_YEARS <- 4; MIN_INAT_CELLS <- 30
REG_MIN_CAM_DET <- 200; REG_MIN_CAM_CELLS <- 30; REG_MIN_INAT_CELLS <- 60; MIN_HP_SD <- 0.05
fmt <- function(e, se) if (is.na(e)) NA_character_ else sprintf("%+.2f [%+.2f, %+.2f]", e, e - 1.96 * se, e + 1.96 * se)
pct <- function(b) 100 * (exp(b) - 1)

tok <- paste0(SP, "_v2u_national_scalar")
b <- readRDS(file.path(H, tok, paste0("input_data_", tok, ".RDS"))); cl <- b$constants_list
prepped <- readRDS(file.path(PROJ, "prototype_spatial_trend", "prepped_sim_inputs.RDS"))
stopifnot(identical(cl$adj, prepped$constants_list$adj), identical(cl$num, prepped$constants_list$num))
E_OF <- prepped$ecoregion_of_cell100; LEV <- prepped$ecoregion_levels

# ---- baseline, centred cell100 Human_pop (identical for both streams) ----
yrs <- as.integer(colnames(b$inat_y)); stopifnot(yrs[1] == 2008)
G <- cl$ncell50; row_g <- integer(dim(cl$xdat_inat)[1])
for (g in 1:G) row_g[cl$inat_cell50_start[g]:cl$inat_cell50_end[g]] <- g
hp0 <- tapply(cl$xdat_inat[, "Human_pop", 1], cl$inat_cell100[row_g], mean)
hpb <- rep(NA_real_, cl$ncell100); hpb[as.integer(names(hp0))] <- as.numeric(hp0) - mean(hp0)

# ---- label check ----
c50 <- fread(file.path(PROJ, "prototype_spatial_trend/cell50_coordinates.csv"))
cells <- as.integer(rownames(b$inat_y))
geo <- data.table(cell50 = cells, region = E_OF[cl$inat_cell100])[c50, on = "cell50", nomatch = 0]
cat("region labels (ecoregion_levels) with mean lon/lat of footprint cells:\n")
print(geo[, .(name = LEV[region[1]], n_cells = .N, lon = round(mean(lon), 1), lat = round(mean(lat), 1)), keyby = region])

# ---- camera deployments ----
suppressMessages(library(terra))
g50 <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/grid50.tif")
oc <- b$real_data$occ.covs; y <- b$real_data$y
xy <- project(as.matrix(oc[, c("longitude", "latitude")]), from = "EPSG:4326", to = crs(g50))
dep <- data.table(cell50 = terra::extract(g50, xy)[, 1], c100 = cl$cell, region = E_OF[cl$cell],
                  year = as.integer(oc$year), project = as.character(oc$project_id),
                  det = rowSums(y > 0, na.rm = TRUE), occ = rowSums(!is.na(y)))
dep <- dep[!is.na(cell50) & occ > 0]; dep[, hp := hpb[c100]]; dep[, yr := year - 2016]
stopifnot(!anyNA(dep$hp))
# ---- iNat cell-years ----
Y <- b$inat_y; E <- b$inat_effort
inat <- data.table(cell50 = rep(cells, times = length(yrs)), year = rep(yrs, each = G), cnt = as.vector(Y), eff = as.vector(E),
                   hp = rep(hpb[cl$inat_cell100], times = length(yrs)), region = rep(E_OF[cl$inat_cell100], times = length(yrs)))
inat <- inat[!is.na(cnt) & !is.na(eff) & eff > 0]; inat[, yr := year - 2016]

ctrl <- glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
grad_re <- function(d, cam) {
  s_hp <- sd(d$hp); d[, `:=`(yd = yr / 10, hs = hp / s_hp)]
  f <- if (cam) cbind(det, occ - det) ~ yd * hs + (1 + yd | cell50) else cnt ~ yd * hs + offset(log(eff)) + (1 + yd | cell50)
  fam <- if (cam) binomial(link = "cloglog") else poisson
  fit <- function(agq) tryCatch(suppressWarnings(suppressMessages(glmer(f, data = d, control = ctrl, family = fam, nAGQ = agq))), error = function(e) e)
  m <- fit(1L); note <- ""
  if (inherits(m, "error")) { note <- paste("full fit failed (", conditionMessage(m), "); nAGQ=0 fit used"); m <- fit(0L) }
  if (inherits(m, "error")) return(list(est = NA_real_, se = NA_real_, mean_trend = NA_real_, note = paste("fit failed:", conditionMessage(m))))
  cf <- fixef(m); V <- as.matrix(vcov(m)); k <- 100 / (10 * s_hp)          # back to pp of %/yr per +1 hp unit
  est <- k * cf[["yd:hs"]]; se <- k * sqrt(V["yd:hs", "yd:hs"])
  # DEGENERATE-FIT GUARD (2026-10-01): elk's national camera fit returned +18.76 [+18.73, +18.78]; an interval
  # that narrow for an effect that size is a collapsed curvature estimate, not precision. Treat as not assessed.
  if (is.finite(est) && abs(est) > 2 && 1.96 * se < 0.02 * abs(est))
  { msg <- sprintf("!!! DEGENERATE FIT -- BUG SIGNAL, NOT AN ESTIMATE: %s %s stream: %+.2f +/- %.3f", SP, if (cam) "camera" else "iNat", est, 1.96 * se)
    message(msg); cat(msg, "\n"); ALARMS <<- c(ALARMS, msg)
    return(list(est = NA_real_, se = NA_real_, mean_trend = NA_real_, note = msg)) }
  list(est = est, se = se, mean_trend = pct(cf[["yd"]] / 10),
       note = trimws(paste(note, if (isSingular(m)) "singular fit" else "")))
}
tertile <- function(d, cam) {
  cut <- quantile(unique(d[, .(cell50, hp)])$hp, c(1/3, 2/3)); d <- d[hp <= cut[1] | hp > cut[2]]; d[, urban := as.integer(hp > cut[2])]
  f <- if (cam) cbind(det, occ - det) ~ yr * urban else cnt ~ yr * urban + offset(log(eff))
  m <- tryCatch(suppressWarnings(glm(f, data = d, family = if (cam) binomial(link = "cloglog") else poisson)), error = function(e) NULL)
  if (is.null(m)) return(NA_character_); cf <- coef(m); sprintf("%+.2f", pct(cf[["yr"]] + cf[["yr:urban"]]) - pct(cf[["yr"]]))
}
own_trend <- function(dd) { if (uniqueN(dd$year) < 2 || sum(dd$det) == 0) return(NA_real_)
  m <- tryCatch(suppressWarnings(glm(cbind(det, occ - det) ~ yr, family = binomial(link = "cloglog"), data = dd)), error = function(e) NULL)
  if (is.null(m)) NA_real_ else round(pct(coef(m)[["yr"]]), 1) }

rows <- list(); conc <- list(); ALARMS <- character(0)
for (st in c(0L, sort(unique(c(dep$region, inat$region))))) {
  nm <- if (st == 0) "NATIONAL" else LEV[st]
  d <- if (st == 0) dep else dep[region == st]; di <- if (st == 0) inat else inat[region == st]
  r <- list(species = SP, stratum = nm, cam_deployments = nrow(d), cam_detecting = sum(d$det > 0), cam_cells = uniqueN(d$cell50),
            cam_years = if (nrow(d)) sprintf("%d-%d", min(d$year), max(d$year)) else NA, inat_cells_with_records = uniqueN(di[cnt > 0, cell50]),
            inat_records = sum(di$cnt))
  reg <- st != 0
  cam_ok <- nrow(d) >= MIN_CAM_DEP && r$cam_detecting >= (if (reg) REG_MIN_CAM_DET else MIN_CAM_DET) &&
            r$cam_cells >= (if (reg) REG_MIN_CAM_CELLS else MIN_CAM_CELLS) && uniqueN(d$year) >= MIN_CAM_YEARS &&
            sd(unique(d[, .(cell50, hp)])$hp) >= MIN_HP_SD
  inat_ok <- r$inat_cells_with_records >= (if (reg) REG_MIN_INAT_CELLS else MIN_INAT_CELLS) &&
             nrow(di) > 0 && sd(unique(di[, .(cell50, hp)])$hp) >= MIN_HP_SD
  if (cam_ok) {
    g <- grad_re(copy(d), TRUE); gc_ <- g; r$cam_gradient <- if (is.na(g$est)) g$note else fmt(g$est, g$se); r$cam_excludes_zero <- abs(g$est) > 1.96 * g$se
    r$cam_ci_width <- round(2 * 1.96 * g$se, 2); r$cam_mean_trend_pct <- round(g$mean_trend, 2); r$cam_note <- g$note
    top <- d[, .(deployments = .N, detecting = sum(det > 0), years = sprintf("%d-%d", min(year), max(year)), n_years = uniqueN(year),
                 projects = uniqueN(project), hp = round(hp[1], 3)), by = cell50][order(-deployments)][1:min(10, uniqueN(d$cell50))]
    top[, share_of_deployments := round(cumsum(deployments) / nrow(d), 3)]
    top[, own_trend_pct := sapply(cell50, function(cc) own_trend(d[cell50 == cc]))]
    conc[[length(conc) + 1]] <- cbind(species = SP, stratum = nm, top)
    dl <- d[!cell50 %in% top$cell50]
    r$cam_gradient_without_top10 <- if (nrow(dl) >= MIN_CAM_DEP && sum(dl$det > 0) >= MIN_CAM_DET && uniqueN(dl$cell50) >= MIN_CAM_CELLS) {
      gl <- grad_re(copy(dl), TRUE); fmt(gl$est, gl$se) } else "NOT ASSESSED (too little data left)"
    # confounds
    pj <- d[, .(ny = uniqueN(year)), by = project]
    r$projects <- nrow(pj); r$projects_spanning_2plus_years <- pj[ny >= 2, .N]
    r$year_hp_spearman <- round(suppressWarnings(cor(d$year, d$hp, method = "spearman")), 3)
    q <- quantile(d$hp, c(1/3, 2/3))
    r$median_year_urban_vs_rural <- sprintf("%.0f vs %.0f", median(d[hp > q[2], year]), median(d[hp <= q[1], year]))
    r$top10_share_of_deployments <- tail(top$share_of_deployments, 1)
  } else r$cam_gradient <- "NOT ASSESSED (too little camera data)"
  if (inat_ok) {
    g <- grad_re(copy(di), FALSE); gi_ <- g; r$inat_gradient <- if (is.na(g$est)) g$note else fmt(g$est, g$se); r$inat_excludes_zero <- abs(g$est) > 1.96 * g$se
    r$inat_ci_width <- round(2 * 1.96 * g$se, 2); r$inat_mean_trend_pct <- round(g$mean_trend, 2); r$inat_note <- g$note
  } else r$inat_gradient <- "NOT ASSESSED (too little iNat data)"
  r$agreement <- "NOT ASSESSED (a stream lacks data)"
  if (cam_ok && inat_ok && exists("gc_") && exists("gi_") && !is.na(gc_$est) && !is.na(gi_$est)) {
    dd <- gi_$est - gc_$est; sd_ <- sqrt(gc_$se^2 + gi_$se^2)
    r$gradient_difference_inat_minus_cam <- fmt(dd, sd_)
    lo <- dd - 1.96 * sd_; hi <- dd + 1.96 * sd_
    r$agreement <- if (lo > 0 || hi < 0) { if (sign(gc_$est) != sign(gi_$est)) "OPPOSE" else "DIFFER (magnitude)" } else
      if (lo <= -10 || hi >= 10) "NOT ASSESSED (difference interval spans 0 and +/-10)" else "AGREE"
  }
  rm(list = intersect(c("gc_", "gi_"), ls()))
  r$fit_alarm <- if (any(grepl("DEGENERATE", c(r$cam_gradient, r$inat_gradient, r$cam_gradient_without_top10)))) "DEGENERATE FIT -- investigate" else ""
  rows[[length(rows) + 1]] <- as.data.table(r)
  cat("\n====", nm, "\n"); print(t(as.data.frame(r)))
}
res <- rbindlist(rows, fill = TRUE)
fwrite(res, file.path(OUT, sprintf("stream_screen_%s.csv", SP)))
if (length(conc)) fwrite(rbindlist(conc, fill = TRUE), file.path(OUT, sprintf("stream_screen_concentration_%s.csv", SP)))
if (length(ALARMS)) { cat("\n########## FIT ALARMS for", SP, "##########\n"); cat(ALARMS, sep = "\n") }
cat("\nSCREEN DONE", SP, "\n")
