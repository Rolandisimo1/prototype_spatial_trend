## ess.R <species> — STEP 3: effective sample size + CI-width inflation
## Fits the hierarchical spec and the naive deployment-only spec on identical data,
## then compares posterior CI width of the mapped diel METRICS at matched locations.
source("cvlib.R")
SP <- commandArgs(TRUE)[1]; NTH <- as.integer(commandArgs(TRUE)[2]); if (is.na(NTH)) NTH <- 4
CAP <- 4; NDRAW <- 400; NLOC <- 60


d <- as.data.table(read_parquet(Sys.getenv("SUNP"))); x0 <- load_species(d, SP); rm(d); gc()
x <- thin_deps(x0, cap = CAP)

## ---- intra-array correlation of diel SHAPE (design effect) ----
## per-deployment normalized 48-bin curve (deployments with >= 20 events),
## reduced to first-harmonic coordinates; ICC from a one-way random-effects
## decomposition of each coordinate across arrays.
dep <- x0[, .(cnt = sum(count)), by = .(dep_key, array_id)][cnt >= 20]
cu  <- x0[dep_key %in% dep$dep_key]
cu[, rate := count / effort_h_aso]
sh  <- cu[, {
  p <- rate / sum(rate)
  .(h1c = sum(p * cos(bin_center_rad)), h1s = sum(p * sin(bin_center_rad)),
    h2c = sum(p * cos(2*bin_center_rad)), h2s = sum(p * sin(2*bin_center_rad)))
}, by = .(dep_key, array_id)]
icc_of <- function(v, g) {
  g <- droplevels(factor(g)); k <- table(g)
  if (length(k) < 2) return(NA_real_)
  gm <- tapply(v, g, mean); n <- length(v); a <- length(k)
  msb <- sum(k * (gm - mean(v))^2) / (a - 1)
  msw <- sum((v - gm[as.character(g)])^2) / (n - a)
  k0  <- (n - sum(k^2)/n) / (a - 1)
  s2b <- max((msb - msw) / k0, 0)
  s2b / (s2b + msw)
}
iccs <- sapply(c("h1c","h1s","h2c","h2s"), function(cc) icc_of(sh[[cc]], sh$array_id))
icc  <- mean(iccs, na.rm = TRUE)
mbar <- sh[, .N, by = array_id][, mean(N)]
deff <- 1 + (mbar - 1) * icc
n_dep_shape <- nrow(sh); n_arr_shape <- uniqueN(sh$array_id)

## ---- fits ----
t0 <- Sys.time(); mh <- fit_bam(f_spatial, x, nthreads = NTH)
t1 <- Sys.time(); mn <- fit_bam(f_naive,   x, nthreads = NTH)
cat(sprintf("[%s] fits done h=%.1f min n=%.1f min\n", SP,
    as.numeric(difftime(t1,t0,units="min")), as.numeric(difftime(Sys.time(),t1,units="min"))))

## ---- posterior CI width of mapped metrics at matched locations ----
set.seed(99)
locs <- unique(x[, .(array_id, x_km, y_km)])[, .(x_km = mean(x_km), y_km = mean(y_km)), by = array_id]
locs <- locs[sample(.N, min(NLOC, .N))]
grid <- CJ(row = seq_len(nrow(locs)), bin = 0:47)
grid[, `:=`(x_km = locs$x_km[row], y_km = locs$y_km[row], array_id = locs$array_id[row])]
grid[, bin_center_rad := (bin + 0.5) * 0.5 / 24 * 2*pi]
grid[, `:=`(s1 = sin(bin_center_rad), c1 = cos(bin_center_rad),
            s2 = sin(2*bin_center_rad), c2 = cos(2*bin_center_rad),
            log_effort_aso = 0)]
grid$dep_key <- factor(levels(x$dep_key)[1], levels = levels(x$dep_key))

## array-level random-effect SDs (for a NEW-array prediction interval).
## mgcv 're' smooths have identity penalty, so coef ~ N(0, scale/sp); nb() has scale=1.
re_sds <- function(m) {
  labs <- vapply(m$smooth, function(s) s$label, character(1))
  out  <- c(int = 0, s1 = 0, c1 = 0, s2 = 0, c2 = 0)
  emp  <- out
  for (i in seq_along(m$smooth)) {
    L <- labs[i]
    if (!grepl("array_id", L)) next
    ## bam() leaves smooth$sp.index empty; m$sp is NAMED by smooth label instead.
    spi <- m$sp[[L]]
    sd_ <- if (is.null(spi) || !is.finite(spi) || spi <= 0) NA_real_ else sqrt(1/spi)
    ## empirical SD of the fitted coefficients, as an independent check
    cf  <- coef(m)[m$smooth[[i]]$first.para:m$smooth[[i]]$last.para]
    key <- if (L == "s(array_id)") "int"
           else if (grepl("s1", L)) "s1" else if (grepl("c1", L)) "c1"
           else if (grepl("s2", L)) "s2" else if (grepl("c2", L)) "c2" else NA
    if (!is.na(key)) { out[key] <- sd_; emp[key] <- sd(cf) }
  }
  ## fall back to the empirical coefficient SD wherever the penalty-based value
  ## is unavailable (sp at its upper bound, i.e. variance shrunk toward zero)
  bad <- !is.finite(out)
  out[bad] <- emp[bad]
  out[!is.finite(out)] <- 0
  attr(out, "empirical") <- emp
  out
}

metric_draws <- function(m, gr, new_array = FALSE) {
  excl <- re_labels(m); assert_excluded(m, excl)
  Xp <- predict(m, gr, type = "lpmatrix", exclude = excl, discrete = FALSE)
  V  <- vcov(m, unconditional = TRUE)
  cf <- coef(m)
  ## zero the excluded (random-effect) columns so draws reflect the surface only
  keep <- which(colSums(abs(Xp)) > 0)
  L <- try(chol(V[keep, keep, drop = FALSE], pivot = TRUE), silent = TRUE)
  if (inherits(L, "try-error")) return(NULL)
  piv <- attr(L, "pivot"); rk <- attr(L, "rank")
  B <- matrix(rnorm(NDRAW * length(keep)), NDRAW)
  Bd <- matrix(0, NDRAW, length(keep))
  Bd[, piv[1:rk]] <- B[, 1:rk, drop = FALSE] %*% L[1:rk, 1:rk, drop = FALSE]
  cfd <- matrix(cf[keep], NDRAW, length(keep), byrow = TRUE) + Bd
  eta <- cfd %*% t(Xp[, keep, drop = FALSE])
  ## honest interval for an UNSEEN array: add a draw of that array's own
  ## diel-curve deviation (intercept + harmonic slopes), not just surface error
  if (new_array) {
    sds <- re_sds(m)
    stopifnot(all(is.finite(sds)))
    if (any(sds > 0)) {
      b0 <- rnorm(NDRAW, 0, sds[["int"]])
      b1 <- rnorm(NDRAW, 0, sds[["s1"]]); b2 <- rnorm(NDRAW, 0, sds[["c1"]])
      b3 <- rnorm(NDRAW, 0, sds[["s2"]]); b4 <- rnorm(NDRAW, 0, sds[["c2"]])
      eta <- eta + b0 +
        outer(b1, gr$s1) + outer(b2, gr$c1) + outer(b3, gr$s2) + outer(b4, gr$c2)
    }
  }
  mu  <- exp(eta)
  nr <- nrow(gr) / 48
  out <- vector("list", nr)
  for (i in seq_len(nr)) {
    idx <- ((i-1)*48 + 1):(i*48)
    M <- mu[, idx, drop = FALSE]
    P <- M / rowSums(M)
    noct <- 100 * rowSums(P[, is_night, drop = FALSE])
    crep <- 100 * rowSums(P[, is_crep, drop = FALSE])
    C <- P %*% cos(ang); S <- P %*% sin(ang)
    conc <- sqrt(C^2 + S^2)
    mh_  <- ((atan2(S, C)) %% (2*pi)) / (2*pi) * 24
    out[[i]] <- data.table(row = i,
      w_noct = diff(quantile(noct, c(.05,.95))), w_crep = diff(quantile(crep, c(.05,.95))),
      w_conc = diff(quantile(conc, c(.05,.95))),
      w_mean = { q <- quantile(((mh_ - median(mh_) + 12) %% 24), c(.05,.95)); as.numeric(diff(q)) })
  }
  rbindlist(out)
}
wh <- metric_draws(mh, grid)                      # hierarchical, surface only
wa <- metric_draws(mh, grid, new_array = TRUE)    # hierarchical, unseen-array PI
wn <- metric_draws(mn, grid)                      # naive deployment-only surface
cmp <- merge(merge(wh, wn, by = "row", suffixes = c("_h", "_n")),
             setnames(wa, setdiff(names(wa), "row"),
                      paste0(setdiff(names(wa), "row"), "_a")), by = "row")
rsd <- re_sds(mh); rsd_emp <- attr(rsd, "empirical")
cat("array RE sd (penalty):", paste(sprintf("%.3f", rsd), collapse=" "),
    "| (empirical):", paste(sprintf("%.3f", rsd_emp), collapse=" "), "\n")

res <- data.table(
  species = SP,
  n_dep_nominal_all = uniqueN(x0$dep_key), n_arrays = uniqueN(x0$array_id),
  n_dep_fitted = uniqueN(x$dep_key), n_events = sum(x0$count),
  mean_dep_per_array_shape = mbar, icc_shape = icc,
  icc_h1c = iccs[["h1c"]], icc_h1s = iccs[["h1s"]], icc_h2c = iccs[["h2c"]], icc_h2s = iccs[["h2s"]],
  n_dep_shape = n_dep_shape, n_arr_shape = n_arr_shape,
  design_effect = deff, n_eff_shape = n_dep_shape / deff,
  edf_hier = sum(mh$edf), edf_naive = sum(mn$edf),
  edf_spatial_hier = sum(mh$edf[mh$smooth[[1]]$first.para:mh$smooth[[1]]$last.para]),
  edf_spatial_naive = sum(mn$edf[mn$smooth[[1]]$first.para:mn$smooth[[1]]$last.para]),
  dev_expl_hier = summary(mh)$dev.expl, dev_expl_naive = summary(mn)$dev.expl,
  ciw_noct_naive = median(cmp$w_noct_n), ciw_noct_hier = median(cmp$w_noct_h),
  ciw_crep_naive = median(cmp$w_crep_n), ciw_crep_hier = median(cmp$w_crep_h),
  ciw_conc_naive = median(cmp$w_conc_n), ciw_conc_hier = median(cmp$w_conc_h),
  ciw_mean_naive = median(cmp$w_mean_n), ciw_mean_hier = median(cmp$w_mean_h),
  ciw_noct_newarray = median(cmp$w_noct_a), ciw_crep_newarray = median(cmp$w_crep_a),
  ciw_conc_newarray = median(cmp$w_conc_a), ciw_mean_newarray = median(cmp$w_mean_a),
  sd_array_int = rsd[["int"]], sd_array_s1 = rsd[["s1"]], sd_array_c1 = rsd[["c1"]],
  sd_array_s2 = rsd[["s2"]], sd_array_c2 = rsd[["c2"]],
  theta_hier = mh$family$getTheta(TRUE), theta_naive = mn$family$getTheta(TRUE))
res[, `:=`(infl_noct = ciw_noct_hier/ciw_noct_naive, infl_crep = ciw_crep_hier/ciw_crep_naive,
           infl_conc = ciw_conc_hier/ciw_conc_naive, infl_mean = ciw_mean_hier/ciw_mean_naive,
           infl_noct_newarray = ciw_noct_newarray/ciw_noct_naive,
           infl_conc_newarray = ciw_conc_newarray/ciw_conc_naive,
           infl_mean_newarray = ciw_mean_newarray/ciw_mean_naive,
           infl_noct_new_vs_surface = ciw_noct_newarray/ciw_noct_hier,
           infl_conc_new_vs_surface = ciw_conc_newarray/ciw_conc_hier)]
fwrite(res, sprintf("ess2_%s.csv", gsub(" ", "_", SP)))
fwrite(cmp, sprintf("ciw_%s.csv", gsub(" ", "_", SP)))
cat(sprintf("[%s] ICC=%.3f mbar=%.1f deff=%.2f n_eff=%.0f | ciw_noct naive=%.1f hier=%.1f newarr=%.1f | edf_sp h=%.0f n=%.0f | theta h=%.2f n=%.2f\n",
    SP, icc, mbar, deff, res$n_eff_shape, res$ciw_noct_naive, res$ciw_noct_hier,
    res$ciw_noct_newarray, res$edf_spatial_hier, res$edf_spatial_naive,
    res$theta_hier, res$theta_naive))
