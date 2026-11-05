#!/usr/bin/env Rscript
# compare_to_fitted.R <src_run_dir> <rebuilt_pre_reduce.RDS> <rebuilt_post_reduce.RDS> <mode>
# Regression gate for the windowed builder (2026-10-17): a rebuild through run_fleet_prep.R must reproduce
# the FITTED 18-yr bundle element for element. Pre-reduce is compared with <src>/input_data_<src>.RDS.bak_10cov,
# post-reduce with the fitted input_data_<src>.RDS. Ignored by design: v2_notes (carries the build time),
# window (attached only by the windowed path), site_reorder + the permutation it applied (fisher only: the
# fitted bundle was reordered after its build, so for a reordered source the post-reduce comparison is
# done after applying the recorded permutation to the rebuild). S4/external-pointer objects (terra SpatRaster
# in real_data$grid100) are compared by class only.
a <- commandArgs(trailingOnly = TRUE); src <- a[1]; pre <- a[2]; post <- a[3]; mode <- a[4]
H <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/HPC"
IGN <- c("v2_notes", "window", "site_reorder")
nd <- 0
cmp <- function(x, y, path = "") {
  if (isS4(x) || isS4(y)) { if (!identical(class(x), class(y))) { nd <<- nd + 1; cat("DIFF class ", path, "\n") }; return(invisible()) }
  if (is.list(x) && is.list(y) && !is.data.frame(x)) {
    for (k in union(names(x), names(y))) {
      if (path == "" && k %in% IGN) next
      if (!k %in% names(x)) { nd <<- nd + 1; cat("ONLY in fitted:  ", path, "$", k, "\n", sep = ""); next }
      if (!k %in% names(y)) { nd <<- nd + 1; cat("ONLY in rebuild: ", path, "$", k, "\n", sep = ""); next }
      cmp(x[[k]], y[[k]], paste0(path, "$", k))
    }
  } else if (!identical(x, y)) {
    nd <<- nd + 1
    cat("DIFF ", path, " class ", class(x)[1], "/", class(y)[1], " dim ", paste(dim(x), collapse = "x"), "/",
        paste(dim(y), collapse = "x"), " all.equal: ", isTRUE(all.equal(x, y, check.attributes = FALSE)), "\n", sep = "")
  }
}
permute_sites <- function(b, idx) {
  ns <- length(idx)
  f <- function(x) {
    if (isS4(x) || is.function(x) || is.environment(x)) return(x)
    if (is.data.frame(x)) { if (nrow(x) == ns) { x <- x[idx, , drop = FALSE]; rownames(x) <- NULL }; return(x) }
    if (is.list(x)) { for (nm in names(x)) x[[nm]] <- f(x[[nm]]); return(x) }
    if (!is.null(dim(x)) && dim(x)[1] == ns)
      return(switch(length(dim(x)), x[idx, drop = FALSE], x[idx, , drop = FALSE], x[idx, , , drop = FALSE]))
    if (is.null(dim(x)) && length(x) == ns && ns > 1) return(x[idx])
    x
  }
  for (nm in names(b)) b[[nm]] <- f(b[[nm]]); b
}
fit_post <- readRDS(file.path(H, src, paste0("input_data_", src, ".RDS")))
reb_post <- readRDS(post)
if (!is.null(fit_post$site_reorder)) { cat("source carries site_reorder: applying its permutation to the rebuild\n")
  reb_post <- permute_sites(reb_post, fit_post$site_reorder$permutation) }
cat("== post-reduce vs fitted", src, "(", mode, ")\n"); nd <- 0; cmp(fit_post, reb_post); n_post <- nd
rm(fit_post, reb_post); invisible(gc())
f10 <- file.path(H, src, paste0("input_data_", src, ".RDS.bak_10cov"))
n_pre <- NA
if (file.exists(f10)) { cat("== pre-reduce vs fitted .bak_10cov\n"); nd <- 0; cmp(readRDS(f10), readRDS(pre)); n_pre <- nd
} else cat("no .bak_10cov for", src, "-- pre-reduce comparison skipped\n")
cat(sprintf("REGRESSION %s %s: post-reduce diffs %d, pre-reduce diffs %s -> %s\n", mode, src, n_post,
            as.character(n_pre), if (n_post == 0 && (is.na(n_pre) || n_pre == 0)) "PASS" else "FAIL"))
