#!/usr/bin/env Rscript
# (1) Bobcat urban-vs-rural camera trend: robustness of the pooled national contrast (+13.8 pp).
# (2) Coyote Eastern Temperate Forests: top-decile-by-population cells vs the rest, camera + iNat.
# Union (v2u) bundle data, no refit. 2026-10-01.
#
# Camera unit, urban/rural tertiles, model and %/yr reporting exactly as regional_urban_rural.R
# (national mode): cbind(det, occ - det) ~ yr * urban, cloglog; urban - rural in pp of %/yr.
# Checks on bobcat:
#   A baseline       reproduce the pooled contrast
#   B project FE     + factor(project_id): trend identified only within projects (ecology vs network turnover)
#   C leave-one-cell-out for the 10 urban cells with most deployments, plus dropping top-3 / top-10 together
#   D cell RE        glmer, (1 | cell50): random intercept
#   E cell RE+slope  glmer, (1 + yr | cell50): urban trend = mean of shrunken per-cell trends,
#                    so no single cell can dominate the slope
# CIs: A-C bootstrap over cells (B: 300 reps, project FE is slow); D-E Wald from glmer.
suppressMessages({library(data.table); library(terra)})
PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
H <- file.path(PROJ, "HPC"); OUT <- "/home/rwkays/isdm/fleet_v2u"
g50 <- rast("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/grid50.tif")
pct <- function(b) 100 * (exp(b) - 1)
fmt <- function(e, lo, hi) sprintf("%+.2f [%+.2f, %+.2f]", e, lo, hi)
has_lme4 <- requireNamespace("lme4", quietly = TRUE)
cat("lme4 available:", has_lme4, "\n")

load_sp <- function(sp) {
  ns <- readRDS(file.path(H, paste0(sp, "_v2u_national_scalar"), paste0("input_data_", sp, "_v2u_national_scalar.RDS")))
  eco <- readRDS(file.path(H, paste0(sp, "_v2u_ecoregion"), paste0("input_data_", sp, "_v2u_ecoregion.RDS")))
  E_OF <- eco$constants_list$ecoregion_of_cell100; rm(eco)
  cl <- ns$constants_list; oc <- ns$real_data$occ.covs; y <- ns$real_data$y
  xy <- project(as.matrix(oc[, c("longitude", "latitude")]), from = "EPSG:4326", to = crs(g50))
  dep <- data.table(cell50 = terra::extract(g50, xy)[, 1], region = E_OF[cl$cell], year = as.integer(oc$year),
                    project = as.character(oc$project_id), det = rowSums(y > 0, na.rm = TRUE),
                    occ = rowSums(!is.na(y)), hp = cl$occ_covars[, "Human_pop"])
  dep <- dep[!is.na(cell50) & occ > 0]; dep[, cell_hp := mean(hp), by = cell50]; dep[, yr := year - 2016]
  list(ns = ns, dep = dep, E_OF = E_OF)
}
fit_c <- function(d, rhs = "yr * urban") {
  m <- tryCatch(suppressWarnings(glm(as.formula(paste("cbind(det, occ - det) ~", rhs)), family = binomial(link = "cloglog"), data = d)),
                error = function(e) NULL)
  if (is.null(m) || !m$converged) return(c(NA_real_, NA_real_))
  cf <- coef(m); c(cf[["yr"]], cf[["yr"]] + cf[["yr:urban"]])
}
contrast <- function(d, rhs = "yr * urban", NB = 1000) {
  est <- fit_c(d, rhs); cr <- unique(d[urban == 0, cell50]); cu <- unique(d[urban == 1, cell50]); setkey(d, cell50)
  bs <- t(replicate(NB, fit_c(d[J(c(sample(cr, length(cr), TRUE), sample(cu, length(cu), TRUE))), allow.cartesian = TRUE], rhs)))
  bs <- bs[complete.cases(bs), , drop = FALSE]; q <- function(v) quantile(v, c(.025, .975))
  r <- pct(est[1]); u <- pct(est[2]); rb <- pct(bs[, 1]); ub <- pct(bs[, 2])
  data.table(rural = fmt(r, q(rb)[1], q(rb)[2]), urban = fmt(u, q(ub)[1], q(ub)[2]),
             urban_minus_rural = fmt(u - r, q(ub - rb)[1], q(ub - rb)[2]), P_urban_gt_rural = round(mean(ub > rb), 3), n_boot = nrow(bs))
}
point <- function(d) { e <- fit_c(d); data.table(rural = round(pct(e[1]), 2), urban = round(pct(e[2]), 2), urban_minus_rural = round(pct(e[2]) - pct(e[1]), 2)) }

# ======================= (1) BOBCAT =======================
set.seed(20261001)
B <- load_sp("bobcat"); d <- B$dep
cut <- quantile(d$cell_hp, c(1/3, 2/3)); d <- d[cell_hp <= cut[1] | cell_hp > cut[2]]; d[, urban := as.integer(cell_hp > cut[2])]
cat("\nBOBCAT national tertiles: rural", d[urban == 0, .N], "deployments /", d[urban == 0, uniqueN(cell50)], "cells; urban",
    d[urban == 1, .N], "/", d[urban == 1, uniqueN(cell50)], "cells\n")
res <- list()
res$A <- cbind(check = "A baseline (pooled, as reported)", contrast(copy(d)))
# B: project fixed effect
pj <- d[, .(n = .N, n_years = uniqueN(year), urban_share = mean(urban)), by = project]
cat("projects:", nrow(pj), "; spanning >= 2 years:", pj[n_years >= 2, .N], "holding", pj[n_years >= 2, sum(n)], "of", nrow(d), "deployments\n")
cat("urban deployments in projects spanning >= 2 years:", d[urban == 1 & project %in% pj[n_years >= 2, project], .N],
    "; rural:", d[urban == 0 & project %in% pj[n_years >= 2, project], .N], "\n")
res$B <- cbind(check = "B + project fixed effect (within-project trends)", contrast(copy(d), "yr * urban + factor(project)", NB = 300))
# C: leave-one-cell-out, top 10 urban cells by deployments
top <- d[urban == 1, .(n_dep = .N, n_det = sum(det > 0), years = sprintf("%d-%d", min(year), max(year)),
                       n_proj = uniqueN(project)), by = cell50][order(-n_dep)][1:10]
top[, share_of_urban_deps := round(cumsum(n_dep) / d[urban == 1, .N], 3)]
cell_slope <- function(cc) { dd <- d[cell50 == cc]; if (uniqueN(dd$year) < 2 || sum(dd$det) == 0) return(NA_real_)
  m <- tryCatch(suppressWarnings(glm(cbind(det, occ - det) ~ yr, family = binomial(link = "cloglog"), data = dd)), error = function(e) NULL)
  if (is.null(m)) NA_real_ else round(pct(coef(m)[["yr"]]), 1) }
top[, own_trend_pct := sapply(cell50, cell_slope)]
loo <- rbindlist(lapply(top$cell50, function(cc) cbind(dropped = as.character(cc), point(d[cell50 != cc]))))
loo <- rbind(cbind(dropped = "none", point(d)), loo,
             cbind(dropped = "top 3 together", point(d[!cell50 %in% top$cell50[1:3]])),
             cbind(dropped = "top 10 together", point(d[!cell50 %in% top$cell50])))
cat("\nTop 10 urban cells by deployments:\n"); print(top)
cat("\nLeave-out results (point estimates, %/yr):\n"); print(loo)
res$C3 <- cbind(check = "C drop top-3 urban cells", contrast(copy(d[!cell50 %in% top$cell50[1:3]])))
res$C10 <- cbind(check = "C drop top-10 urban cells", contrast(copy(d[!cell50 %in% top$cell50])))
# D/E: cell random effects
if (has_lme4) {
  re_fit <- function(form, lab) {
    m <- tryCatch(suppressWarnings(lme4::glmer(form, family = binomial(link = "cloglog"), data = d,
                  control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))), error = function(e) e)
    if (inherits(m, "error")) return(data.table(check = lab, urban_minus_rural = paste("FIT FAILED:", conditionMessage(m))))
    cf <- lme4::fixef(m); V <- as.matrix(vcov(m))
    br <- cf[["yr"]]; bu <- cf[["yr"]] + cf[["yr:urban"]]
    se_r <- sqrt(V["yr", "yr"]); se_u <- sqrt(V["yr", "yr"] + V["yr:urban", "yr:urban"] + 2 * V["yr", "yr:urban"])
    se_i <- sqrt(V["yr:urban", "yr:urban"])
    dif <- pct(bu) - pct(br)
    # delta-method-free reporting: interval of the log-scale interaction mapped to pp at the rural level
    lo <- pct(br + cf[["yr:urban"]] - 1.96 * se_i) - pct(br); hi <- pct(br + cf[["yr:urban"]] + 1.96 * se_i) - pct(br)
    data.table(check = lab, rural = fmt(pct(br), pct(br - 1.96 * se_r), pct(br + 1.96 * se_r)),
               urban = fmt(pct(bu), pct(bu - 1.96 * se_u), pct(bu + 1.96 * se_u)),
               urban_minus_rural = fmt(dif, lo, hi), P_urban_gt_rural = round(pnorm(cf[["yr:urban"]] / se_i), 3),
               singular = lme4::isSingular(m))
  }
  res$D <- re_fit(cbind(det, occ - det) ~ yr * urban + (1 | cell50), "D cell random intercept")
  res$E <- re_fit(cbind(det, occ - det) ~ yr * urban + (1 + yr | cell50), "E cell random intercept + random year-slope")
} else {
  res$D <- data.table(check = "D/E cell random effects", urban_minus_rural = "NOT RUN: lme4 not installed in this env")
}
bob <- rbindlist(res, fill = TRUE); options(width = 250); cat("\n"); print(bob)
fwrite(bob, file.path(OUT, "bobcat_urban_camera_robustness.csv"))
fwrite(top, file.path(OUT, "bobcat_urban_top10_cells.csv")); fwrite(loo, file.path(OUT, "bobcat_urban_leave_out.csv"))

# ======================= (2) COYOTE ETF top decile =======================
set.seed(20261002)
C <- load_sp("coyote"); dc <- C$dep[region == 1]
cutc <- quantile(unique(dc[, .(cell50, cell_hp)])$cell_hp, 0.9)
dc[, urban := as.integer(cell_hp > cutc)]
cam10 <- cbind(stream = "camera", contrast(copy(dc)),
               urban_deployments = dc[urban == 1, .N], urban_detecting = dc[urban == 1, sum(det > 0)], urban_cells = dc[urban == 1, uniqueN(cell50)],
               rest_deployments = dc[urban == 0, .N], rest_detecting = dc[urban == 0, sum(det > 0)], rest_cells = dc[urban == 0, uniqueN(cell50)],
               urban_years = dc[urban == 1, sprintf("%d-%d median %.0f", min(year), max(year), median(year))],
               rest_years = dc[urban == 0, sprintf("%d-%d median %.0f", min(year), max(year), median(year))])
ns <- C$ns; cl <- ns$constants_list; Y <- ns$inat_y; E <- ns$inat_effort; cells <- as.integer(rownames(Y)); yrs <- as.integer(colnames(Y))
G <- cl$ncell50; row_g <- integer(dim(cl$xdat_inat)[1]); for (g in 1:G) row_g[cl$inat_cell50_start[g]:cl$inat_cell50_end[g]] <- g
hp_cell <- as.vector(tapply(rowMeans(cl$xdat_inat[, "Human_pop", ]), row_g, mean))
inat <- data.table(cell50 = rep(cells, times = length(yrs)), year = rep(yrs, each = G), cnt = as.vector(Y), eff = as.vector(E),
                   hp = rep(hp_cell, times = length(yrs)), region = rep(C$E_OF[cl$inat_cell100], times = length(yrs)))
inat <- inat[region == 1 & !is.na(cnt) & !is.na(eff) & eff > 0]; inat[, yr := year - 2016]
cuti <- quantile(unique(inat[, .(cell50, hp)])$hp, 0.9); inat[, urban := as.integer(hp > cuti)]
fit_i <- function(dd, rhs = "yr * urban") { m <- tryCatch(suppressWarnings(glm(as.formula(paste("cnt ~", rhs, "+ offset(log(eff))")), family = poisson, data = dd)), error = function(e) NULL)
  if (is.null(m) || !m$converged) return(c(NA_real_, NA_real_)); cf <- coef(m)
  if ("yr:urban" %in% names(cf)) c(cf[["yr"]], cf[["yr"]] + cf[["yr:urban"]]) else c(cf[["yr"]], cf[["yr:hp"]]) }
ci_i <- function(dd) { est <- fit_i(dd); cr <- unique(dd[urban == 0, cell50]); cu <- unique(dd[urban == 1, cell50]); setkey(dd, cell50)
  bs <- t(replicate(1000, fit_i(dd[J(c(sample(cr, length(cr), TRUE), sample(cu, length(cu), TRUE))), allow.cartesian = TRUE])))
  bs <- bs[complete.cases(bs), , drop = FALSE]; q <- function(v) quantile(v, c(.025, .975))
  r <- pct(est[1]); u <- pct(est[2]); rb <- pct(bs[, 1]); ub <- pct(bs[, 2])
  data.table(rural = fmt(r, q(rb)[1], q(rb)[2]), urban = fmt(u, q(ub)[1], q(ub)[2]), urban_minus_rural = fmt(u - r, q(ub - rb)[1], q(ub - rb)[2]),
             P_urban_gt_rural = round(mean(ub > rb), 3), n_boot = nrow(bs)) }
inat10 <- cbind(stream = "iNat", ci_i(copy(inat)), urban_cells = inat[urban == 1 & cnt > 0, uniqueN(cell50)], rest_cells = inat[urban == 0 & cnt > 0, uniqueN(cell50)],
                urban_records = inat[urban == 1, sum(cnt)], rest_records = inat[urban == 0, sum(cnt)])
# continuous gradient with and without the top decile (does the gradient live in the tail?)
cg <- function(dd) { cs <- unique(dd$cell50); setkey(dd, cell50)
  one <- function(x) { m <- tryCatch(suppressWarnings(glm(cnt ~ yr * hp + offset(log(eff)), family = poisson, data = x)), error = function(e) NULL)
    if (is.null(m) || !m$converged) NA_real_ else 100 * coef(m)[["yr:hp"]] }
  e <- one(dd); bs <- replicate(1000, one(dd[J(sample(cs, length(cs), TRUE)), allow.cartesian = TRUE])); bs <- bs[!is.na(bs)]
  fmt(e, quantile(bs, .025), quantile(bs, .975)) }
grad <- data.table(stream = "iNat continuous gradient (pp per +1 unit Human_pop)",
                   all_cells = cg(copy(inat)), without_top_decile = cg(copy(inat[urban == 0])))
cat("\nCOYOTE ETF top decile vs rest (camera cut", round(cutc, 3), "; iNat cut", round(cuti, 3), ")\n")
print(cam10); print(inat10); print(grad)
fwrite(rbindlist(list(cam10, inat10), fill = TRUE), file.path(OUT, "coyote_etf_topdecile.csv"))
fwrite(grad, file.path(OUT, "coyote_etf_gradient_tail.csv"))
cat("ROBUSTNESS DONE\n")
