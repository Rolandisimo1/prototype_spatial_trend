## step4_summary.R -- unit-weighted skill + cluster bootstrap over FOLDS (RULE 5),
## plus effect support / VIF diagnostics and the final-model fits.
source("covlib.R")
NB <- 4000L
CV <- fread("cv_folds_raw.csv")

## unit-weighted skill: weight each fold by its number of held-out arrays.
## Event-weighting let a few large arrays dominate in Phase 1 -- reported separately if needed.
uw <- function(w, m, mn) {
  ok <- is.finite(m) & is.finite(mn) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  1 - sum(w[ok]*m[ok])/sum(w[ok]*mn[ok])
}
boot_ci <- function(dt, wcol, mcol, ncol_, nb = NB, seed = 11) {
  set.seed(seed); n <- nrow(dt)
  if (n < 2) return(c(NA_real_, NA_real_))
  b <- replicate(nb, { i <- sample.int(n, n, TRUE)
                       uw(dt[[wcol]][i], dt[[mcol]][i], dt[[ncol_]][i]) })
  b <- b[is.finite(b)]
  if (!length(b)) return(c(NA_real_, NA_real_))
  as.numeric(quantile(b, c(0.025, 0.975)))
}

S <- CV[, {
  ci  <- boot_ci(.SD, "n_te", "mae_cov", "mae_null")
  cip <- boot_ci(.SD, "n_te", "mae_pos", "mae_null")
  cii <- boot_ci(.SD[is.finite(mae_cov_inside)], "n_inside_eval", "mae_cov_inside", "mae_null_inside")
  .(cov_skill = uw(n_te, mae_cov, mae_null),
    ci_lo = ci[1], ci_hi = ci[2],
    pos_skill = uw(n_te, mae_pos, mae_null),
    pos_ci_lo = cip[1], pos_ci_hi = cip[2],
    cov_skill_inside = uw(n_inside_eval, mae_cov_inside, mae_null_inside),
    in_ci_lo = cii[1], in_ci_hi = cii[2],
    cov_skill_outside = uw(n_outside_eval, mae_cov_outside, mae_null_outside),
    folds_pos = sum(skill_cov > 0, na.rm = TRUE), n_folds = .N,
    n_te_tot = sum(n_te), n_inside_tot = sum(n_inside_eval, na.rm = TRUE),
    n_vars = n_vars[1])
}, by = .(species, variant, metric)]
S[, `:=`(beats_null = is.finite(ci_lo) & ci_lo > 0,
         beats_position = is.finite(cov_skill) & is.finite(pos_skill) & cov_skill > pos_skill,
         beats_null_inside = is.finite(in_ci_lo) & in_ci_lo > 0,
         ## the harmonic-framework position competitor, bootstrapped the SAME way.
         ## Distinct from the PUBLISHED bam position benchmark (skill_bootstrap_cis.csv,
         ## 0/20 CIs excluding zero) -- both are reported, they are not the same model.
         pos_beats_null = is.finite(pos_ci_lo) & pos_ci_lo > 0)]
fwrite(S, "covariate_model_cv.csv")
cat(sprintf("rows %d | CI>0: %d | beats position: %d\n", nrow(S), sum(S$beats_null), sum(S$beats_position)))
print(S[beats_null == TRUE, .(species, variant, metric, cov_skill = round(cov_skill,3),
       ci_lo = round(ci_lo,3), ci_hi = round(ci_hi,3), pos = round(pos_skill,3), folds_pos)])
