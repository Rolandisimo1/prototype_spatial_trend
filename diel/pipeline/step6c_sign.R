## step6c_sign.R -- collinearity sign test on % nocturnal, JOINT vs ALONE.
## Point estimates from the same weighted GAMs used for the map.
## Uncertainty via POSTERIOR SIMULATION of the fitted coefficient models (fast, and
## the same machinery the maps use), plus a PROJECT-CLUSTER bootstrap at a modest
## number of reps for the sign-stability rate. Reported separately so neither is
## mistaken for the other.
source("covlib.R")
NBP <- 0L; NSIM <- 500L
H <- fread("array_harmonics_raw.csv"); AC <- fread("array_covariates.csv")
D <- merge(H, AC, by = "array_id"); D[, project_id := sub("\\|.*$", "", array_id)]

noct_delta <- function(mods, dv, vars, focal, nsim = 0) {
  q  <- as.numeric(quantile(dv[[focal]], c(0.10, 0.90), na.rm = TRUE))
  nd <- dv[rep(1, 2)]
  for (v in vars) set(nd, j = v, value = rep(as.numeric(median(dv[[v]], na.rm = TRUE)), 2))
  set(nd, j = focal, value = q)
  pc <- predict_curves(mods, nd, nsim = nsim, seed = 9)
  pt <- pc$M$noct[2] - pc$M$noct[1]
  if (!nsim) return(list(est = pt, lo = NA_real_, hi = NA_real_))
  sm <- sim_metrics(pc$sims)$noct          # 2 x nsim
  d  <- sm[2, ] - sm[1, ]
  list(est = pt, lo = quantile(d, 0.025, na.rm = TRUE), hi = quantile(d, 0.975, na.rm = TRUE))
}
fitmods <- function(dv, vars, kk, lin)
  lapply(setNames(HARM, HARM), function(h) try(fit_coef(dv, h, vars, k = kk, linear = lin), silent = TRUE))

rows <- list()
for (SP in unique(D$species)) {
  dd <- D[species == SP]; bear <- SP == "American Black Bear"
  kk <- if (bear) 3 else 5; lin <- bear
  for (focal in c("ntl_1km","pop_1km")) {
    for (mode in c("joint","alone")) {
      vars <- if (mode == "alone") focal else
              if (bear) unique(c(intersect(BEAR_SET, setA("C")), focal)) else setA("C")
      dv <- apply_trans(dd, vars); dv <- dv[complete.cases(dv[, ..vars])]
      if (nrow(dv) < 25) next
      mods <- fitmods(dv, vars, kk, lin)
      if (any(sapply(mods, inherits, "try-error"))) next
      r <- noct_delta(mods, dv, vars, focal, nsim = NSIM)
      ## project-cluster bootstrap for SIGN STABILITY only (few reps, sign is robust)
      pr <- unique(dv$project_id); set.seed(3); sg <- integer(0)
      for (b in seq_len(NBP)) {
        db <- dv[unlist(lapply(sample(pr, length(pr), TRUE), function(p) which(dv$project_id == p)))]
        if (uniqueN(db[[focal]]) < 8) next
        mb <- fitmods(db, vars, kk, lin)
        if (any(sapply(mb, inherits, "try-error"))) next
        v <- tryCatch(noct_delta(mb, db, vars, focal)$est, error = function(e) NA_real_)
        if (is.finite(v)) sg <- c(sg, sign(v))
      }
      rows[[length(rows)+1]] <- data.table(species = SP, covariate = focal, mode = mode,
        n_arrays = nrow(dv), effect_pp = r$est, sim_lo = r$lo, sim_hi = r$hi,
        pct_same_sign_projboot = if (length(sg)) mean(sg == sign(r$est)) else NA_real_,
        nboot = length(sg))
      cat(sprintf("  %-22s %-8s %-6s %+6.2f pp [%+.2f,%+.2f] signstab %.2f\n", SP, focal, mode,
                  r$est, r$lo, r$hi, if (length(sg)) mean(sg == sign(r$est)) else NA_real_)); flush.console()
    }
  }
}
R <- rbindlist(rows); fwrite(R, "human_pressure_effects.csv")
w <- dcast(R, species + covariate ~ mode,
           value.var = c("effect_pp","sim_lo","sim_hi","pct_same_sign_projboot"))
w[, sign_stable := sign(effect_pp_joint) == sign(effect_pp_alone)]
fwrite(w, "human_pressure_signtest.csv")
print(w[, .(species, covariate, joint = round(effect_pp_joint,2), alone = round(effect_pp_alone,2),
            sign_stable, stab_alone = round(pct_same_sign_projboot_alone,2))])
