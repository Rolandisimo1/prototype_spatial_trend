#!/usr/bin/env Rscript
# Tier 1 smoke-test checker: is every axis live, are yb_* populated, does
# yb_true move with the effect level? Stops on the first failed check.
dirs <- Sys.glob("tier1_smoke_it*_*")
r <- do.call(rbind, lapply(dirs, function(d) {
  it <- as.integer(sub("tier1_smoke_it(\\d+)_.*", "\\1", d))
  do.call(rbind, lapply(list.files(d, "^row_.*rds$", full.names = TRUE), function(f) {
    x <- readRDS(f); x$smoke_iter <- it; x[setdiff(names(x), "error_msg")] }))
}))
options(width = 250)
show <- c("block","row_id","range_cells","window_frac","n_site","inat_records","effect","smoke_iter","status",
          "d_inat_cells","d_cam_cells","d_ncell50","d_inat_mult","d_nsite","d_nsite_unique","d_site_pool","d_cam_year_min","d_cam_nyears","d_s_occ",
          "r_cam_det_sites","r_cam_naive_occ","r_inat_records","d_inat_expected_null","cam_pct_true","yb_true","yv_true","yb_mean","yb_ci_lb","yb_ci_ub",
          "yb_covered","yb_detected","tvb_true","tvb_mean","elapsed_sec")
print(r[order(r$block, r$row_id, r$smoke_iter), intersect(show, names(r))], row.names = FALSE, digits = 4)
chk <- function(ok, msg) { cat(sprintf("[%s] %s\n", if (isTRUE(ok)) "PASS" else "FAIL", msg)); if (!isTRUE(ok)) stop(msg) }
chk(all(r$status == "OK"), "every smoke row status OK")
yb <- c("yb_true","yb_mean","yb_ci_lb","yb_ci_ub","yb_bias","yb_ci_width","yb_covered","yb_detected")
chk(all(yb %in% names(r)) && !anyNA(r[, yb]), "all yb_* columns present and populated")
chk(all(r$d_range_cells == r$range_cells), "realized range size == range label, every row")
rw <- r[r$block == "RW", ]
rr <- unique(rw[order(rw$range_cells), c("range_cells","d_inat_cells","d_cam_cells")])
chk(nrow(rr) == 5 && all(diff(rr$d_inat_cells) > 0) && all(diff(rr$d_cam_cells) >= 0) && all(rr$d_inat_cells <= rr$range_cells),
    "range axis live: cells carrying iNat data rise across the 5 ranges, cameras never leave the range")
chk(all(rw$d_cam_nyears == round(pmax(2, rw$window_frac * 18))) && length(unique(rw$d_cam_nyears)) == 4 &&
    all(rw$d_cam_year_min + rw$d_cam_nyears - 1 == 2025), "window axis live: 18/11/7/4 camera years ending 2025")
chk(length(unique(round(rw$d_s_occ, 6))) == 4, "window axis reaches year_occ scaling (4 distinct s_occ)")
vol <- r[r$block == "VOL", ]
chk(all(vol$d_nsite == vol$n_site) && length(unique(vol$d_nsite)) == 4, "camera axis live: 4 realized site counts")
v1 <- vol[vol$n_site == min(vol$n_site) & vol$effect == -0.05 & vol$smoke_iter == min(vol$smoke_iter), ]
v1 <- v1[order(v1$inat_records), ]
chk(nrow(v1) == 4 && isTRUE(all.equal(v1$d_inat_expected_null, v1$inat_records)) && all(diff(v1$r_inat_records) > 0) &&
    all(v1$r_inat_records / v1$inat_records > 0.3 & v1$r_inat_records / v1$inat_records < 3),
    "iNat axis live: expected records == label, realized records rise with it and land near the target")
v2 <- vol[vol$inat_records == 20000 & vol$effect == -0.05 & vol$smoke_iter == min(vol$smoke_iter), ]; v2 <- v2[order(v2$n_site), ]
chk(all(diff(v2$r_cam_det_sites) > 0), "camera axis reaches the data: detecting sites rise with camera count")
e <- rw[rw$range_cells == 300 & rw$window_frac == 0.4, ]; e <- e[order(e$effect), ]
chk(length(unique(e$yb_true)) == 4 && sum(e$yb_true == 0) == 1 && all(e$yv_true[e$effect == 0] == 0) &&
    all(e$tvb_true[e$effect == 0] == 0), "yb_true differs across the 4 effect levels; null has year_beta = year_var = 0")
chk(all(abs(r$cam_pct_true - 100 * r$effect) < 1e-6), "camera %/yr truth equals the effect label, every row")
# timing: per-iteration cost from rows run at two N_ITER values
two <- merge(r[r$smoke_iter == min(r$smoke_iter), c("block","row_id","elapsed_sec")],
             r[r$smoke_iter == max(r$smoke_iter), c("block","row_id","elapsed_sec","smoke_iter")], by = c("block","row_id"))
if (nrow(two) && max(r$smoke_iter) > min(r$smoke_iter)) {
  two$sec_per_iter <- (two$elapsed_sec.y - two$elapsed_sec.x) / (max(r$smoke_iter) - min(r$smoke_iter))
  two$fixed_sec <- two$elapsed_sec.x - two$sec_per_iter * (min(r$smoke_iter) + 100)
  two$full_10k_hours <- (two$fixed_sec + two$sec_per_iter * 10000) / 3600
  print(merge(two, unique(r[, c("block","row_id","n_site","range_cells")])), row.names = FALSE, digits = 3)
}
write.csv(r, "tier1_smoke_rows.csv", row.names = FALSE)
cat("SMOKE CHECK PASSED\n")
