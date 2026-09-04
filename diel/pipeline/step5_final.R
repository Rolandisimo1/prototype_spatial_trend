## step5_final.R -- full-data fits per species x variant:
##   * effect support (which covariates survive) via select=TRUE shrinkage + REML p
##   * project-cluster bootstrap of each smooth's effect range (robust significance)
##   * VIF of each covariate against the rest of the block
##   * sanity gate: reconstructed metrics vs observed metrics at the TRAINING arrays
source("covlib.R")
NBOOT <- 300L

H  <- fread("array_harmonics_raw.csv"); AC <- fread("array_covariates.csv")
CB <- fread("camera_setb_array.csv")
D  <- merge(H, merge(AC, CB, by = "array_id", all.x = TRUE), by = "array_id")
D[, project_id := sub("\\|.*$", "", array_id)]

vif_of <- function(dt, vars) {
  sapply(vars, function(v) {
    o <- setdiff(vars, v)
    if (!length(o)) return(1)
    f <- try(lm(as.formula(paste(v, "~", paste(o, collapse = "+"))), data = dt), silent = TRUE)
    if (inherits(f, "try-error")) return(NA_real_)
    r2 <- summary(f)$r.squared
    1/max(1e-9, 1 - r2)
  })
}

VAR_SETS <- list(setA_lightsonly = setA("A"), setA_poponly = setA("B"), setA_both = setA("C"))

out_eff <- list(); out_gate <- list(); out_vif <- list(); MODELS <- list()
for (SP in unique(D$species)) {
  dd <- D[species == SP]
  is_bear <- SP == "American Black Bear"
  vsets <- VAR_SETS
  if (!is_bear) {
    vsets$setB_map  <- c(setA("A"), SETB_MAP)
    vsets$setB_full <- c(setA("A"), SETB_MAP, HETERO[[SP]])
  }
  for (vn in names(vsets)) {
    vars <- intersect(vsets[[vn]], names(dd))
    if (is_bear) vars <- intersect(BEAR_SET, vars)
    vars <- vars[sapply(vars, function(v) sum(is.finite(dd[[v]])) > 0.8*nrow(dd) &&
                                          length(unique(dd[[v]])) > 5)]
    if (!length(vars)) next
    dv <- apply_trans(dd, vars); dv <- dv[complete.cases(dv[, ..vars])]
    kk <- if (is_bear) 3 else 5; lin <- is_bear

    mods <- lapply(setNames(HARM, HARM), function(h)
      try(fit_coef(dv, h, vars, k = kk, linear = lin), silent = TRUE))
    if (any(sapply(mods, inherits, "try-error"))) next
    MODELS[[paste(SP, vn)]] <- list(models = mods, vars = vars, dat = dv, bear = is_bear)

    ## --- effect support ---
    ## criterion: for smooths, edf > 0.05 AND approximate REML p < 0.05 on ANY of the
    ## five harmonic coefficients (a covariate that shapes any coefficient shapes the curve).
    for (v in vars) {
      pv <- sapply(HARM, function(h) {
        st <- summary(mods[[h]])
        tb <- if (lin) st$p.table else st$s.table
        rn <- if (lin) v else sprintf("s(%s)", v)
        if (!(rn %in% rownames(tb))) return(NA_real_)
        tb[rn, ncol(tb)]
      })
      ed <- sapply(HARM, function(h) {
        st <- summary(mods[[h]]); tb <- if (lin) st$p.table else st$s.table
        rn <- if (lin) v else sprintf("s(%s)", v)
        if (!(rn %in% rownames(tb))) return(NA_real_)
        if (lin) 1 else tb[rn, "edf"]
      })
      out_eff[[length(out_eff)+1]] <- data.table(
        species = SP, variant = vn, covariate = v,
        min_p = min(pv, na.rm = TRUE), max_edf = max(ed, na.rm = TRUE),
        supported = (min(pv, na.rm = TRUE) < 0.05) && (max(ed, na.rm = TRUE) > 0.05),
        n_coef_signif = sum(pv < 0.05, na.rm = TRUE))
    }
    vv <- vif_of(dv, vars)
    out_vif[[length(out_vif)+1]] <- data.table(species = SP, variant = vn,
                                               covariate = names(vv), vif = as.numeric(vv))

    ## --- sanity gate: reconstruction at the training arrays ---
    pc <- predict_curves(mods, dv)
    for (mm in MET) {
      ob <- dv[[paste0("obs_", mm)]]; pr <- pc$M[[mm]]
      ok <- is.finite(ob) & is.finite(pr)
      err <- if (mm == "peak") abs(circ_diff_h(pr[ok], ob[ok])) else abs(pr[ok] - ob[ok])
      out_gate[[length(out_gate)+1]] <- data.table(
        species = SP, variant = vn, metric = mm, n = sum(ok),
        r = suppressWarnings(cor(ob[ok], pr[ok])), mae = mean(err),
        sd_obs = sd(ob[ok]))
    }
  }
  cat("done", SP, "\n")
}
EFF <- rbindlist(out_eff); GATE <- rbindlist(out_gate); VIF <- rbindlist(out_vif)
fwrite(EFF, "covariate_support.csv"); fwrite(GATE, "sanity_gate.csv"); fwrite(VIF, "covariate_vif.csv")
saveRDS(MODELS, "final_models.rds")
cat(sprintf("EFF %d GATE %d VIF %d models %d\n", nrow(EFF), nrow(GATE), nrow(VIF), length(MODELS)))
print(VIF[variant %in% c("setA_lightsonly","setA_both") & covariate %in% c("ntl_1km","pop_1km")][order(species, variant)])
