## step1_harmonics.R -- array-level diel curves as harmonic coefficients.
## Unit of analysis = camera ARRAY (finest scale with replication).
## Rate model: log E[count_b] = log(effort_b) + b0 + b1 sin(t) + b2 cos(t) + b3 sin(2t) + b4 cos(2t)
## quasipoisson: point estimates = Poisson, SEs scaled by estimated dispersion (clamped >= 1).
suppressPackageStartupMessages({library(data.table); library(arrow); library(sf)})
source("cvlib.R")                       # REF_SR/REF_SS, AEA, curve_metrics, NBIN

SUNP  <- Sys.getenv("SUNP"); COVP <- Sys.getenv("COVP")
GRIDP <- Sys.getenv("GRIDP"); RMP  <- Sys.getenv("RMP"); L48P <- Sys.getenv("L48P")
stopifnot(file.exists(SUNP), file.exists(COVP))
NSIM <- 1000L
SPECIES <- c("White-tailed Deer","Northern Raccoon","Eastern Gray Squirrel",
             "Coyote","American Black Bear")
SCI <- c("White-tailed Deer"="Odocoileus virginianus","Northern Raccoon"="Procyon lotor",
         "Eastern Gray Squirrel"="Sciurus carolinensis","Coyote"="Canis latrans",
         "American Black Bear"="Ursus americanus")

## ---- design matrix on the 48 sun-time bins (shared by every array) ----
bcen <- (seq_len(NBIN) - 0.5) * 0.5          # bin centres on the 24 h sun-time scale
trad <- bcen / 24 * 2*pi
Xd   <- cbind(1, sin(trad), cos(trad), sin(2*trad), cos(2*trad))
colnames(Xd) <- c("b0","s1","c1","s2","c2")

## Rowcliffe activity level = mean/max of the fitted rate curve
act_level <- function(v) { v <- pmax(v, 0); m <- max(v); if (!is.finite(m) || m <= 0) NA_real_ else mean(v)/m }
all_metrics <- function(v) { m <- curve_metrics(v); m$act <- act_level(v); m }
MET <- c("noct","crep","conc","peak","mean","act")

## ---- load: aso window only, finite effort, lower 48 (deployments present in cov table) ----
d <- as.data.table(read_parquet(SUNP))
cv <- fread(COVP)
cv[, dep_key := paste0(project_id, "|", deployment_id)]          # RULE 2: composite key
l48_dep <- unique(cv$dep_key)
d <- d[aso_window == TRUE & effort_h_aso > 0 & is.finite(log_effort_aso) & dep_key %in% l48_dep]
cat(sprintf("rows after aso+L48 filter: %d\n", nrow(d)))

## ---- array -> 100 km cell, for the range mask (point in polygon on the viewer grid) ----
gr  <- st_read(GRIDP, quiet = TRUE)
rmk <- fread(RMP)
arr_ll <- d[, .(lon = mean(longitude), lat = mean(latitude)), by = array_id]
ap  <- st_as_sf(arr_ll, coords = c("lon","lat"), crs = 4326)
ix  <- as.integer(st_within(ap, st_make_valid(gr), sparse = TRUE))
arr_ll[, cell_id := ifelse(is.na(ix), NA_character_, gr$cell_id[ix])]
cat(sprintf("arrays with a 100km cell: %d / %d\n", sum(!is.na(arr_ll$cell_id)), nrow(arr_ll)))

## AEA km coords per array (same projection the position models used)
co <- st_coordinates(st_transform(ap, AEA))
arr_ll[, `:=`(x_km = co[,1], y_km = co[,2])]

out <- list()
for (SP in SPECIES) {
  x <- d[species == SP]
  if (!nrow(x)) next
  ## RULE 4: mask to this species' UNION range mask (100 km cells)
  keep_cells <- rmk[species == SCI[[SP]], cell_id]
  ok_arr <- arr_ll[cell_id %in% keep_cells, array_id]
  x <- x[array_id %in% ok_arr]

  ## aggregate cameras -> ARRAY x bin. RULE 1: effort summed, rate = count/effort.
  ag <- x[, .(count = sum(count), eff = sum(effort_h_aso), ndep = uniqueN(dep_key)),
          by = .(array_id, bin)]
  setorder(ag, array_id, bin)
  ev <- ag[, .(events = sum(count), nbin = .N, ndep = max(ndep), eff_tot = sum(eff)), by = array_id]

  for (thr in c(30L, 50L, 100L)) cat(sprintf("  [%s] thr=%3d -> %d arrays\n", SP, thr, sum(ev$events >= thr)))
  cat(sprintf("[%s] arrays total %d, events %d\n", SP, nrow(ev), sum(ev$events)))
  out[[SP]] <- list(ag = ag, ev = ev)
}
saveRDS(list(out = out, arr_ll = arr_ll, Xd = Xd, bcen = bcen), "s1_agg.rds")
cat("STEP1-AGG DONE\n")
