#!/usr/bin/env Rscript
# Re-smoke check (2026-10-17): yb_rhat_split and yb_ess present, finite and plausible in every row.
fs <- list.files(c("tier1_smoke_rhat_it1400_RW", "tier1_smoke_rhat_it1400_VOL"), "^row_\\d+\\.rds$", full.names = TRUE)
r <- do.call(rbind, lapply(fs, readRDS))
print(r[, c("block", "row_id", "status", "range_cells", "n_site", "effect", "yb_true", "yb_mean", "yb_rhat_split", "yb_ess", "elapsed_sec")], row.names = FALSE)
stopifnot(nrow(r) == 3, all(r$status == "OK"), all(c("yb_rhat_split", "yb_ess") %in% names(r)),
          all(is.finite(r$yb_rhat_split)), all(r$yb_rhat_split > 0.9), all(is.finite(r$yb_ess)), all(r$yb_ess > 0))
cat("RHAT SMOKE PASSED\n")
