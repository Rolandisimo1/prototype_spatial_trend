## step10_resid_moran.R -- two spatial-confounding diagnostics on the FINAL fitted models.
##
## (1) RESIDUAL MORAN'S I by distance band, per species, Set A and Set B, on the
##     harmonic coefficients AND on the reconstructed % nocturnal. Smooths should
##     remove more clustering than the linear stand-in did; whatever remains at
##     0-50 km is fine-scale structure our covariates do NOT capture.
## (2) COEFFICIENT MOVEMENT with vs without a spatial term: refit each covariate
##     model adding s(x_km,y_km) and measure how much the covariate's effect on
##     % nocturnal (p90 vs p10) changes. A large move means position was carrying it.
source("covlib.R")
BANDS <- list(c(0,50), c(50,200), c(200,800))

H <- fread("array_harmonics_raw.csv"); AC <- fread("array_covariates.csv")
CB <- fread("camera_setb_array.csv")
D  <- merge(H, merge(AC, CB, by = "array_id", all.x = TRUE), by = "array_id")
MODELS <- readRDS("final_models.rds")

## inverse-variance weighted Moran's I in a distance band (binary band weights)
moran_band <- function(x, y, z, w, lo, hi) {
  ok <- is.finite(z) & is.finite(w) & w > 0
  x <- x[ok]; y <- y[ok]; z <- z[ok]; w <- w[ok]
  n <- length(z); if (n < 20) return(c(NA_real_, NA_real_, NA_real_))
  zc <- z - sum(w*z)/sum(w)
  num <- 0; den <- sum(w*zc^2); S0 <- 0
  for (i in seq_len(n)) {
    d <- sqrt((x - x[i])^2 + (y - y[i])^2)
    nb <- which(d >= lo & d < hi & seq_len(n) != i)
    if (!length(nb)) next
    ww <- w[i]*w[nb]
    num <- num + sum(ww * zc[i] * zc[nb]); S0 <- S0 + sum(ww)
  }
  if (S0 <= 0 || den <= 0) return(c(NA_real_, NA_real_, NA_real_))
  I <- (n/S0) * (num/den)
  ## permutation p (200 shuffles) -- cheap and assumption-free
  set.seed(5); pm <- numeric(200)
  for (b in seq_len(200)) {
    zp <- zc[sample.int(n)]
    nu <- 0
    for (i in seq_len(n)) {
      d <- sqrt((x - x[i])^2 + (y - y[i])^2)
      nb <- which(d >= lo & d < hi & seq_len(n) != i)
      if (!length(nb)) next
      nu <- nu + sum(w[i]*w[nb] * zp[i] * zp[nb])
    }
    pm[b] <- (n/S0) * (nu/den)
  }
  c(I, mean(abs(pm) >= abs(I)), n)
}

rows <- list(); movement <- list()
for (SP in unique(D$species)) {
  for (vn in c("setA_lightsonly","setB_map")) {
    key <- paste(SP, vn); MM <- MODELS[[key]]
    if (is.null(MM)) next
    dv <- MM$dat; vars <- MM$vars
    ## observed and fitted % nocturnal at the training arrays
    pc <- predict_curves(MM$models, dv)
    res_noct <- dv$obs_noct - pc$M$noct
    w_noct   <- 1/pmax(dv$se_noct^2, 1e-6)
    for (bd in BANDS) {
      raw <- moran_band(dv$x_km, dv$y_km, dv$obs_noct, w_noct, bd[1], bd[2])
      rsd <- moran_band(dv$x_km, dv$y_km, res_noct,   w_noct, bd[1], bd[2])
      rows[[length(rows)+1]] <- data.table(species = SP, variant = vn,
        band = sprintf("%d-%d km", bd[1], bd[2]), n = raw[3],
        moran_raw = raw[1], p_raw = raw[2],
        moran_resid = rsd[1], p_resid = rsd[2],
        pct_clustering_removed = 100*(1 - rsd[1]/raw[1]))
    }
    ## --- coefficient movement with vs without a spatial term ---
    if (vn != "setA_lightsonly") next
    for (focal in intersect(c("ntl_1km","nlcd_1k_crop","nlcd_1k_pasture","tcc_1km"), vars)) {
      eff <- function(with_space) {
        mods <- lapply(setNames(HARM, HARM), function(h) {
          wgt <- 1/pmax(dv[[paste0("v_", h)]], 1e-8); wgt <- wgt/mean(wgt)
          d2 <- copy(dv); d2$.w <- wgt
          tt <- if (MM$bear) vars else sprintf("s(%s, k=5)", vars)
          if (with_space) tt <- c(tt, "s(x_km, y_km, k=30)")
          try(gam(as.formula(paste(h, "~", paste(tt, collapse = " + "))), data = d2,
                  weights = .w, method = "REML", select = TRUE), silent = TRUE)
        })
        if (any(sapply(mods, inherits, "try-error"))) return(NA_real_)
        q <- as.numeric(quantile(dv[[focal]], c(0.10, 0.90), na.rm = TRUE))
        nd <- dv[rep(1, 2)]
        for (v in vars) set(nd, j = v, value = rep(as.numeric(median(dv[[v]], na.rm = TRUE)), 2))
        set(nd, j = focal, value = q)
        if (with_space) { set(nd, j = "x_km", value = rep(median(dv$x_km), 2))
                          set(nd, j = "y_km", value = rep(median(dv$y_km), 2)) }
        M <- predict_curves(mods, nd)$M
        M$noct[2] - M$noct[1]
      }
      e_no <- eff(FALSE); e_sp <- eff(TRUE)
      movement[[length(movement)+1]] <- data.table(species = SP, covariate = focal,
        effect_no_space = e_no, effect_with_space = e_sp,
        abs_move_pp = abs(e_sp - e_no),
        sign_kept = is.finite(e_no) && is.finite(e_sp) && sign(e_no) == sign(e_sp))
    }
    cat("done", SP, vn, "\n"); flush.console()
  }
}
RM <- rbindlist(rows); MV <- rbindlist(movement)
fwrite(RM, "residual_moran.csv"); fwrite(MV, "coef_movement_space.csv")
print(RM[variant == "setA_lightsonly", .(species, band, moran_raw = round(moran_raw,3),
        moran_resid = round(moran_resid,3), pct_removed = round(pct_clustering_removed,1), p_resid)])
print(MV[, .(species, covariate, no_space = round(effect_no_space,2),
             with_space = round(effect_with_space,2), move = round(abs_move_pp,2), sign_kept)])
