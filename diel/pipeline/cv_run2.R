## cv_run2.R <species> <nthreads> — leave-array-out + spatial-block CV, one species.
## Differences from v1, all for tractability/robustness (documented in the report):
##  * 5 folds per scheme instead of 10/8 — every array is still held out exactly
##    once, so unit coverage is identical; only the number of refits halves.
##  * smoothing parameters estimated on fold 1's TRAINING data and reused for the
##    remaining folds (no full-data leak; saves the fREML outer loop each time).
##  * per-fold incremental CSV append, so a restart never loses completed folds.
source("cvlib.R")
args <- commandArgs(TRUE); SP <- args[1]; NTH <- as.integer(args[2])
if (is.na(NTH)) NTH <- 6
BLOCK_KM <- 400; NFOLD <- 5
MIN_EV_ARRAY <- 30; MIN_EV_BLOCK <- 100; CAP <- 4
TAG  <- gsub(" ", "_", SP)
OUTF <- sprintf("cvf_%s.csv", TAG)

d <- as.data.table(read_parquet(Sys.getenv("SUNP")))
x <- load_species(d, SP); rm(d); gc()
cat(sprintf("[%s] rows=%d deps=%d arrays=%d events=%d\n", SP, nrow(x),
            uniqueN(x$dep_key), uniqueN(x$array_id), sum(x$count))); flush.console()

arr <- x[, .(x_km = mean(x_km), y_km = mean(y_km), ev = sum(count),
             ndep = uniqueN(dep_key)), by = array_id]
arr[, block := paste0("B", floor(x_km/BLOCK_KM), "_", floor(y_km/BLOCK_KM))]
set.seed(20260814); arr[, fold_array := sample(rep_len(1:NFOLD, .N))]
blk <- arr[, .(x_km = mean(x_km), y_km = mean(y_km)), by = block]
set.seed(7)
km <- kmeans(as.matrix(blk[, .(x_km, y_km)]), centers = min(NFOLD, nrow(blk)),
             nstart = 25, iter.max = 100)
blk[, fold_block := km$cluster]
arr <- merge(arr, blk[, .(block, fold_block)], by = "block")
x   <- merge(x, arr[, .(array_id, block, fold_array, fold_block)], by = "array_id")

## ---- per-unit observed/predicted curve metrics ----
eval_unit <- function(te, pm, pn, by_col, min_ev) {
  te <- copy(te); te[, `:=`(pm = pm, pn = pn)]
  agg <- te[, .(cnt = sum(count), eff = sum(effort_h_aso), pm = sum(pm), pn = sum(pn)),
            by = c(by_col, "bin")]
  keep <- agg[, .(ev = sum(cnt)), by = c(by_col)][ev >= min_ev][[by_col]]
  agg <- agg[get(by_col) %in% keep]; setorderv(agg, c(by_col, "bin"))
  agg[, `:=`(o_rate = cnt/eff, m_rate = pm/eff, n_rate = pn/eff)]
  out <- agg[, {
    po <- o_rate/sum(o_rate); pm_ <- m_rate/sum(m_rate); pn_ <- n_rate/sum(n_rate)
    mo <- curve_metrics(o_rate); mm <- curve_metrics(m_rate); mn <- curve_metrics(n_rate)
    .(n_events = sum(cnt),
      r_shape = cor(po, pm_), r_shape_null = cor(po, pn_),
      obs_noct = mo$pct_noct, pred_noct = mm$pct_noct, null_noct = mn$pct_noct,
      obs_crep = mo$pct_crep, pred_crep = mm$pct_crep, null_crep = mn$pct_crep,
      obs_conc = mo$conc,     pred_conc = mm$conc,     null_conc = mn$conc,
      obs_peak = mo$peak_h,   pred_peak = mm$peak_h,   null_peak = mn$peak_h,
      obs_mean = mo$mean_h,   pred_mean = mm$mean_h,   null_mean = mn$mean_h,
      mnll_model = sum(cnt*log(pmax(pm_,1e-12)))/sum(cnt),
      mnll_null  = sum(cnt*log(pmax(pn_,1e-12)))/sum(cnt),
      mnll_sat   = sum(cnt*log(pmax(po, 1e-12)))/sum(cnt))
  }, by = c(by_col)]
  out[, `:=`(ae_noct = abs(pred_noct-obs_noct), ae_noct_null = abs(null_noct-obs_noct),
             ae_crep = abs(pred_crep-obs_crep), ae_crep_null = abs(null_crep-obs_crep),
             ae_conc = abs(pred_conc-obs_conc), ae_conc_null = abs(null_conc-obs_conc),
             ae_peak = abs(circ_diff_h(pred_peak,obs_peak)),
             ae_peak_null = abs(circ_diff_h(null_peak,obs_peak)),
             ae_mean = abs(circ_diff_h(pred_mean,obs_mean)),
             ae_mean_null = abs(circ_diff_h(null_mean,obs_mean)))]
  out[]
}

SP_CACHE <- list(spatial = NULL, null = NULL)

run_fold <- function(scheme, fold_col, by_col, min_ev, k) {
  ho <- arr[get(fold_col) == k, as.character(array_id)]
  if (!length(ho)) return(NULL)
  tr <- thin_deps(x[!array_id %in% ho], cap = CAP)
  te <- x[array_id %in% ho]
  t0 <- Sys.time()
  m  <- try(fit_bam(f_spatial, tr, sp = SP_CACHE$spatial, nthreads = NTH), silent = TRUE)
  if (inherits(m, "try-error")) { cat("FIT FAIL", scheme, k, "\n"); return(NULL) }
  mn <- try(fit_bam(f_null, tr, sp = SP_CACHE$null, nthreads = NTH), silent = TRUE)
  if (inherits(mn, "try-error")) { cat("NULL FAIL", scheme, k, "\n"); return(NULL) }
  if (is.null(SP_CACHE$spatial)) {
    SP_CACHE$spatial <<- m$sp; SP_CACHE$null <<- mn$sp
    cat("  [sp cached from first fit]\n")
  }
  pm <- pred_surface(m,  te, tr)
  pn <- pred_surface(mn, te, tr)
  r  <- eval_unit(te, pm, pn, by_col, min_ev)
  if (!nrow(r)) return(NULL)
  trA <- arr[!as.character(array_id) %in% ho]; hoA <- arr[as.character(array_id) %in% ho]
  D <- as.matrix(dist(rbind(as.matrix(hoA[, .(x_km,y_km)]), as.matrix(trA[, .(x_km,y_km)]))))
  nh <- nrow(hoA)
  hoA[, d_nn_train_km := apply(D[1:nh, (nh+1):ncol(D), drop = FALSE], 1, min)]
  if (by_col == "array_id") {
    r <- merge(r, hoA[, .(array_id, d_nn_train_km, ndep)], by = "array_id", all.x = TRUE)
  } else {
    r <- merge(r, hoA[, .(d_nn_train_km = median(d_nn_train_km), ndep = sum(ndep)),
                      by = block], by = "block", all.x = TRUE)
  }
  r[, `:=`(species = SP, scheme = scheme, fold = k,
           n_train_arrays = uniqueN(tr$array_id), n_train_dep = uniqueN(tr$dep_key),
           n_ho_arrays = length(ho), dev_expl = summary(m)$dev.expl,
           theta = m$family$getTheta(TRUE),
           fit_min = as.numeric(difftime(Sys.time(), t0, units = "min")))]
  r[, unit := as.character(get(by_col))]
  fwrite(r, OUTF, append = file.exists(OUTF))
  cat(sprintf("[%s] %s fold %d/%d units=%d r=%.3f (null %.3f) mnll %.4f vs %.4f | %.1f min\n",
      SP, scheme, k, NFOLD, nrow(r), mean(r$r_shape, na.rm=TRUE),
      mean(r$r_shape_null, na.rm=TRUE), mean(r$mnll_model), mean(r$mnll_null),
      r$fit_min[1])); flush.console()
  rm(m, mn, tr, te); gc(); invisible(NULL)
}

for (k in 1:NFOLD) run_fold("array", "fold_array", "array_id", MIN_EV_ARRAY, k)
for (k in 1:NFOLD) run_fold("block", "fold_block", "block",    MIN_EV_BLOCK, k)
cat(sprintf("[%s] ALLDONE\n", SP))
