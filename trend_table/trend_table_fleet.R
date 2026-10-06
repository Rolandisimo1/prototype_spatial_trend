#!/usr/bin/env Rscript
# Fleet trend table (2026-10-05). Same quantities as trend_decomp.R (trend_decomp_union_vs_presence.csv),
# extended to all 18 v2u national_scalar species (PART=ns) and to the 5 ecoregion fits (PART=eco).
#   camera %/yr = 100*(exp(year_beta*s_occ)-1)          iNat %/yr = 100*(exp(total_var_beta*s_vals)-1)
#   gap_pp = iNat%/yr - camera%/yr, per draw.            s_* = median step of the covariate's distinct values.
# R-hat is computed here on the trend parameters (3 chains, gelman.diag, no autoburnin) so non-convergence
# shows in the table itself. Flags:
#   cam_direction: "increase"/"decline" if camera 95% CrI excludes 0; "near zero" if it includes 0 but both
#     bounds lie within +/-2 %/yr (a precise null); "NOT RESOLVABLE" if it includes 0 and is wider than that.
#   converged: "yes" rhat_trend_max <= 1.05; "marginal" <= 1.1; "NO" > 1.1.
# Ecoregion: union ecoregion puts year_region on iNat ONLY (camera = national year_beta);
#   shared ecoregion puts it on both: camera_r = year_beta + year_region[r], iNat_r = total_var_beta + year_region[r].
suppressPackageStartupMessages(library(coda))
H <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/HPC"
OUT <- "/home/rwkays/isdm/fleet_v2u/trend_table"; dir.create(OUT, showWarnings = FALSE)
getv <- function(b, v) { x <- NULL; for (nm in c("constants_list", "data_list")) if (!is.null(b[[nm]][[v]])) x <- b[[nm]][[v]]; if (is.null(x)) x <- b[[v]]; x }
stepv <- function(x) { u <- sort(unique(round(x, 7))); median(diff(u)) }
pct <- function(x, s) 100 * (exp(x * s) - 1)
q   <- function(x, d = 2) sprintf(paste0("%+.", d, "f [%+.", d, "f, %+.", d, "f]"), mean(x), quantile(x, .025), quantile(x, .975))
dirflag <- function(x) { lo <- quantile(x, .025); hi <- quantile(x, .975)
  if (lo > 0) "increase" else if (hi < 0) "decline" else if (lo >= -2 && hi <= 2) "near zero" else "NOT RESOLVABLE" }
convflag <- function(r) if (r <= 1.05) "yes" else if (r <= 1.1) "marginal" else "NO"
load3 <- function(tok, pars) {
  f <- file.path(H, tok, sprintf("chain_%s_%d_single.RDS", tok, 1:3)); stopifnot(all(file.exists(f)))
  lapply(f, function(x) { s <- readRDS(x)$samples; p <- if (is.function(pars)) pars(colnames(s)) else pars
                          s[, p, drop = FALSE] })
}
rhat <- function(ch) gelman.diag(mcmc.list(lapply(ch, mcmc)), autoburnin = FALSE, multivariate = FALSE)$psrf[, 1]
steps <- function(tok) { b <- readRDS(file.path(H, tok, paste0("input_data_", tok, ".RDS")))
  list(so = stepv(getv(b, "year_occ")), sv = stepv(getv(b, "year_vals")), b = b) }
row_for <- function(a, rh, so, sv) {
  cam <- pct(a[, "year_beta"], so); ina <- pct(a[, "total_var_beta"], sv); yv <- pct(a[, "year_var"], sv)
  data.frame(s_occ = round(so, 4), s_vals = round(sv, 4),
    year_beta = q(a[, "year_beta"], 4), year_var = q(a[, "year_var"], 4), total_var_beta = q(a[, "total_var_beta"], 4),
    rhat_year_beta = round(rh[["year_beta"]], 4), rhat_year_var = round(rh[["year_var"]], 4),
    rhat_total_var_beta = round(rh[["total_var_beta"]], 4),
    camera_pct = q(cam), inat_pct = q(ina), year_var_pct = q(yv), gap_pp = q(ina - cam),
    P_cam_gt0 = round(mean(cam > 0), 3), P_inat_lt0 = round(mean(ina < 0), 3),
    P_opposite = round(mean(cam > 0 & ina < 0) + mean(cam < 0 & ina > 0), 3),
    cam_ci_width_pp = round(diff(quantile(cam, c(.025, .975))), 2),
    cam_direction = dirflag(cam), rhat_trend_max = round(max(rh), 4), converged = convflag(max(rh)))
}
TP <- c("year_beta", "year_var", "total_var_beta")
part <- Sys.getenv("PART", "ns"); options(width = 300)

if (part == "ns") {
  spp <- c("moose", "bobcat", "white-tailed_deer", "coyote",
           "american_black_bear", "black-tailed_jackrabbit", "eastern_chipmunk", "eastern_fox_squirrel",
           "eastern_gray_squirrel", "elk", "grey_fox", "mule_deer", "nine-banded_armadillo", "northern_raccoon",
           "red_fox", "red_squirrel", "striped_skunk", "virginia_opossum")
  rows <- list()
  for (sp in spp) {
    tok <- paste0(sp, "_v2u_national_scalar"); st <- steps(tok); st$b <- NULL
    ch <- load3(tok, TP); rh <- rhat(ch); a <- do.call(rbind, ch)
    rows[[sp]] <- cbind(species = sp, mask = "v2u", row_for(a, rh, st$so, st$sv))
    cat(sp, "done; rhat_trend_max", max(rh), "\n"); rm(ch, a); invisible(gc())
  }
  out <- do.call(rbind, rows)
  write.csv(out, file.path(OUT, "trend_table_v2u_national_scalar_18spp.csv"), row.names = FALSE)
  print(out[, c("species", "s_occ", "camera_pct", "inat_pct", "gap_pp", "P_cam_gt0", "cam_direction", "rhat_trend_max", "converged")], row.names = FALSE)
  # cross-check the four already-reported species against trend_decomp_union_vs_presence.csv
  old <- read.csv("/home/rwkays/isdm/fleet_v2u/trend_decomp_union_vs_presence.csv"); old <- old[old$mask == "v2u", ]
  m <- merge(old, out, by = "species", suffixes = c(".old", ".new"))
  cat("\nCROSS-CHECK vs trend_decomp_union_vs_presence.csv (should be identical):\n")
  print(m[, c("species", "camera_pct.old", "camera_pct.new", "gap_pp.old", "gap_pp.new", "rhat_max", "rhat_trend_max")], row.names = FALSE)
  cat("all identical:", all(m$camera_pct.old == m$camera_pct.new & m$gap_pp.old == m$gap_pp.new & m$inat_pct.old == m$inat_pct.new), "\n")
}

if (part == "eco") {
  prepped <- readRDS(file.path(dirname(H), "prototype_spatial_trend", "prepped_sim_inputs.RDS"))
  lev <- prepped$ecoregion_levels
  fits <- c("moose_v2u_ecoregion", "bobcat_v2u_ecoregion", "white-tailed_deer_v2u_ecoregion",
            "white-tailed_deer_v2u_ecoregion_shared", "bobcat_v2u_ecoregion_shared")
  nat <- list(); reg <- list()
  for (tok in fits) {
    shared <- grepl("_shared$", tok); sp <- sub("_v2u_ecoregion.*", "", tok)
    st <- steps(tok); cl <- st$b$constants_list
    cell_reg <- cl$ecoregion_of_cell100[cl$cell]; det <- rowSums(st$b$real_data$y > 0, na.rm = TRUE) > 0
    inat_reg <- cl$ecoregion_of_cell100[cl$inat_cell100]; st$b <- NULL
    ch <- load3(tok, function(cn) c(TP, "sigma_region", grep("^year_region\\[", cn, value = TRUE)))
    rh <- rhat(ch); a <- do.call(rbind, ch); R <- sum(grepl("^year_region\\[", colnames(a)))
    nat[[tok]] <- cbind(fit = tok, species = sp, model = if (shared) "shared ecoregion (camera+iNat)" else "union ecoregion (iNat-only regional)",
                        row_for(a, rh[TP], st$so, st$sv), sigma_region = q(a[, "sigma_region"], 4),
                        rhat_sigma_region = round(rh[["sigma_region"]], 4),
                        rhat_year_region_max = round(max(rh[grep("^year_region", names(rh))]), 4))
    for (r in 1:R) {
      yr <- a[, sprintf("year_region[%d]", r)]
      cam <- pct(a[, "year_beta"] + if (shared) yr else 0, st$so); ina <- pct(a[, "total_var_beta"] + yr, st$sv)
      rr <- c(rh[TP], rh[[sprintf("year_region[%d]", r)]])
      reg[[paste(tok, r)]] <- data.frame(fit = tok, species = sp, shared = shared, region = r,
        region_name = if (!is.null(lev)) lev[r] else NA, camera_sites = sum(cell_reg == r), camera_detecting_sites = sum(cell_reg == r & det),
        inat_cell50s = sum(inat_reg == r), year_region = q(yr, 4),
        camera_pct = if (shared) q(cam) else paste(q(cam), "(national; no regional camera term)"),
        inat_pct = q(ina), gap_pp = q(ina - cam), P_cam_gt0 = round(mean(cam > 0), 3), P_inat_gt0 = round(mean(ina > 0), 3),
        cam_direction = if (shared) dirflag(cam) else "n/a (national)", inat_direction = dirflag(ina),
        data_note = if (sum(cell_reg == r & det) == 0 && sum(inat_reg == r) == 0) "PRIOR ONLY: no detecting cameras, no iNat cells"
                    else if (sum(cell_reg == r & det) == 0) "no detecting camera sites" else "",
        rhat_max = round(max(rr), 4), converged = convflag(max(rr)))
    }
    cat(tok, "done; trend rhat max", max(rh[TP]), " year_region rhat max", max(rh[grep("^year_region", names(rh))]), "\n")
    rm(ch, a); invisible(gc())
  }
  n <- do.call(rbind, nat); g <- do.call(rbind, reg)
  write.csv(n, file.path(OUT, "trend_table_ecoregion_national.csv"), row.names = FALSE)
  write.csv(g, file.path(OUT, "trend_table_ecoregion_regional.csv"), row.names = FALSE)
  print(n[, c("fit", "camera_pct", "inat_pct", "gap_pp", "cam_direction", "sigma_region", "rhat_trend_max", "rhat_year_region_max")], row.names = FALSE)
  print(g[, c("fit", "region", "region_name", "camera_detecting_sites", "inat_cell50s", "camera_pct", "inat_pct", "gap_pp", "cam_direction", "inat_direction", "data_note", "rhat_max")], row.names = FALSE)
}
cat("DONE\n")
