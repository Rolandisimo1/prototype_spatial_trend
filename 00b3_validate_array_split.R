#!/usr/bin/env Rscript
# =============================================================================
# validate_array_split.R -- STAGE A. READ-ONLY.
# Re-derive the 5 km complete-linkage array-year split on the FULL raw data and
# check it reproduces the already-agreed totals from array_units_r5km.csv /
# SEED Sec 4:  4,002 units | 26,748 cameras | 963,494 trap-nights | median 4.
# If it does not reproduce them, STOP -- do not build site_array on a different
# structure.
# =============================================================================
suppressMessages({ library(tidyverse) })

RAW  <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data"
SPT  <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final/prototype_spatial_trend"
OUT  <- Sys.getenv("OUTDIR", "/home/rwkays/isdm/diag/out_arraysplit")
RADIUS <- as.numeric(Sys.getenv("RADIUS_KM", "5"))
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

dep <- suppressMessages(read_csv(file.path(RAW, "combined_deployments_all.csv"),
                                show_col_types = FALSE))
cat("deployments:", nrow(dep), "\n")

# camera_trap_array lives in the SEQUENCES file (col 15), one row per sequence-
# class; reduce to distinct (deployment_id, camera_trap_array).
# MUST use a real CSV parser, not awk -F',': deployment_id is a QUOTED field
# that contains commas, e.g.
#   "OR_Wetland_IslandCWS_Water Garden, Lodge Tree 08/26/2024"
# which shifts every naive comma-split column past it and corrupts both the id
# and camera_trap_array's position. That defect produced an under-populated
# camera_trap_array map and spurious subproject_name fallbacks on a first pass.
cta_file <- file.path(OUT, "dep_to_cta.csv")
if (!file.exists(cta_file)) {
  cat("extracting deployment_id -> camera_trap_array from sequences (CSV-aware) ...\n")
  seqs <- suppressMessages(read_csv(file.path(RAW, "combined_sequences_all.csv"),
                                    col_select = c(deployment_id, camera_trap_array),
                                    show_col_types = FALSE, progress = FALSE))
  cat("  sequence rows read:", nrow(seqs), "\n")
  seqs <- distinct(seqs, deployment_id, camera_trap_array)
  write_csv(seqs, cta_file)
  rm(seqs); gc()
}
cta <- suppressMessages(read_csv(cta_file, show_col_types = FALSE))
cta <- cta %>% filter(!is.na(camera_trap_array), camera_trap_array != "") %>%
  distinct(deployment_id, .keep_all = TRUE)
cat("deployments with a camera_trap_array:", nrow(cta), "\n")

d <- dep %>%
  mutate(deployment_id = as.character(deployment_id)) %>%
  left_join(cta %>% mutate(deployment_id = as.character(deployment_id)),
            by = "deployment_id") %>%
  mutate(array_field = ifelse(!is.na(camera_trap_array) & camera_trap_array != "",
                              camera_trap_array, subproject_name),
         src = ifelse(!is.na(camera_trap_array) & camera_trap_array != "",
                      "camera_trap_array", "subproject_name"))
cat("src tally (all deployments):\n"); print(table(d$src, useNA = "ifany"))

# usable rows: need coords, year, an array label
d <- d %>% filter(!is.na(longitude), !is.na(latitude), !is.na(year), !is.na(array_field))
cat("usable rows:", nrow(d), "\n")

# ---- complete-linkage split at RADIUS km, within (array_field, year) --------
# Great-circle distances; complete linkage guarantees the diameter bound
# (single linkage would chain and can exceed it -- SEED Sec 4).
hav <- function(lon, lat) {
  n <- length(lon); R <- 6371
  la <- lat * pi/180; lo <- lon * pi/180
  m <- matrix(0, n, n)
  for (i in seq_len(n)) {
    dlat <- la - la[i]; dlon <- lo - lo[i]
    a <- sin(dlat/2)^2 + cos(la[i]) * cos(la) * sin(dlon/2)^2
    m[i, ] <- 2 * R * asin(pmin(1, sqrt(a)))
  }
  as.dist(m)
}

d$grp <- paste(d$array_field, d$year, sep = "@@")
split_one <- function(ix) {
  if (length(ix) == 1L) return("0")
  dd <- hav(d$longitude[ix], d$latitude[ix])
  if (max(dd) <= RADIUS) return(rep("0", length(ix)))
  as.character(cutree(hclust(dd, method = "complete"), h = RADIUS) - 1L)
}
idx <- split(seq_len(nrow(d)), d$grp)
sub <- character(nrow(d))
for (g in idx) sub[g] <- split_one(g)
d$unit_5 <- paste0(d$grp, "#", sub)

u <- d %>% group_by(unit_5) %>%
  summarise(n_cameras = n(),
            trap_nights = sum(survey_nights, na.rm = TRUE),
            year = first(year), src = first(src), .groups = "drop")

cat("\n===== RE-DERIVED (radius", RADIUS, "km) =====\n")
cat("  units:        ", nrow(u), "\n")
cat("  cameras:      ", sum(u$n_cameras), "\n")
cat("  trap_nights:  ", round(sum(u$trap_nights)), "\n")
cat("  median cam/unit:", median(u$n_cameras), "\n")
cat("  singletons:   ", sum(u$n_cameras == 1), "\n")
cat("  src tally:\n"); print(table(u$src))

ref <- suppressMessages(read_csv(file.path(SPT, "array_units_r5km.csv"), show_col_types = FALSE))
cat("\n===== REFERENCE array_units_r5km.csv =====\n")
cat("  units:        ", nrow(ref), "\n")
cat("  cameras:      ", sum(ref$n_cameras), "\n")
cat("  trap_nights:  ", round(sum(ref$trap_nights)), "\n")
cat("  median cam/unit:", median(ref$n_cameras), "\n")
cat("  singletons:   ", sum(ref$n_cameras == 1), "\n")
cat("  src tally:\n"); print(table(ref$src))

ok <- (nrow(u) == nrow(ref)) && (sum(u$n_cameras) == sum(ref$n_cameras)) &&
      (median(u$n_cameras) == median(ref$n_cameras))
cat("\n===== MATCH:", ok, "=====\n")
if (!ok) cat("*** counts differ -- STOP and report before building site_array ***\n")
write.csv(u, file.path(OUT, "rederived_units.csv"), row.names = FALSE)
saveRDS(d[, c("deployment_id","year","array_field","src","unit_5","longitude","latitude")],
        file.path(OUT, "dep_unit_map.RDS"))
cat("wrote dep_unit_map.RDS (per-deployment unit label)\n")
cat("DONE\n")
