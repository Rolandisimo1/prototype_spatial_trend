#!/usr/bin/env Rscript
# =============================================================================
# rebuild_detections.R -- downstream correction of umflist.RDS for camera names lost upstream (2026-10-08).
#
# umflist.RDS is built by codebase/Data_prep/01_Prep_UMFs.R (B. Wagler, 3/16/2026), which keeps ONE
# scientific name per species (line 150: sequences[sciname %in% species_list]). Wildlife Insights renamed
# collared peccary Pecari tajacu -> Dicotyles tajacu at 2023/24, so 47.9% of its record (all 2024-2025)
# was dropped; 855 elk sequences labelled Cervus elaphus were dropped too. Her script cannot run on Hazel
# (CONUS boundary, canopy raster and GRIP roads are not here), so this rebuilds ONLY the detection matrix
# for the affected species and keeps her deployments, rows, siteCovs and obsCovs untouched.
#
# GATE (rwkays: "the gate is the whole thing"): for EVERY species in umflist, the matrix is first rebuilt from
# its single original name with her exact window logic and must reproduce her y EXACTLY (same NA pattern,
# same 0/1 values, same dimnames). Any difference, for any species, stops the script before anything is
# written. A close match does not pass.
#
# Then, for SCOPE species only (rwkays: peccary and elk), y is rebuilt from the crosswalk's full camera name
# list (camera_name_map.csv: MDD renames + recorded decisions) and written to umflist_v2.RDS beside hers.
# Her umflist.RDS is never modified. attr(umflist_v2, "wi_names_used") records the names behind every y.
#
# Also the standing detecting-site ratio check (rwkays 2026-10-08): per fleet species, deployments with >= 1
# detection in the camera record vs raw detecting deployments (distinct project___deployment uids with >= 1
# sequence under any of the species' camera names -- keyed on MDD-resolved scientific names, never on
# common name). Flags ratio < RATIO_MIN. Written for umflist and umflist_v2.
# =============================================================================
suppressPackageStartupMessages(library(data.table))
D <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data"
W <- "/home/rwkays/isdm/fleet_v2u"
OUT_UMF <- file.path(D, "umflist_v2.RDS")
SCOPE <- strsplit(Sys.getenv("SCOPE", "collared_peccary,elk"), ",")[[1]]
RATIO_MIN <- 0.75
options(width = 220)
if (file.exists(OUT_UMF)) stop("refusing to overwrite ", OUT_UMF)

u <- readRDS(file.path(D, "umflist.RDS"))
ref <- lapply(u, function(x) attr(x, "y"))
uids <- rownames(ref[[1]])
stopifnot(!is.null(uids), !anyDuplicated(uids), all(vapply(ref, function(y) identical(rownames(y), uids), TRUE)))
cat("umflist:", length(u), "species x", length(uids), "deployments x", ncol(ref[[1]]), "windows\n")

# ---- deployments, exactly as 01_Prep_UMFs.R (read.csv, !is.na(start_date), distinct, as.Date) ----------
all_dep <- read.csv(file.path(D, "combined_deployments_all.csv"))
all_dep <- all_dep[!is.na(all_dep$start_date), ]
all_dep <- unique(all_dep)                                   # dplyr::distinct() on all columns
all_dep$start_date <- as.Date(all_dep$start_date); all_dep$end_date <- as.Date(all_dep$end_date)
dep <- as.data.table(all_dep)[, uid := paste(project_id, deployment_id, sep = "___")]
# Her CONUS filter needs a boundary file that is not on Hazel; stand-in = her own final deployment set.
dep <- dep[uid %in% uids]
dep[, inner_days := as.integer(end_date - start_date - 1)]
dep <- dep[inner_days >= 10 & inner_days <= 100]
if (anyDuplicated(dep$uid)) stop("GATE: duplicate deployment uids after her filters: ", sum(duplicated(dep$uid)))
if (!identical(dep$uid, uids)) stop("GATE: deployment set/order differs from umflist rows (", nrow(dep), " vs ", length(uids), ")")
cat("deployments reproduce umflist rows exactly:", nrow(dep), "\n")

valid_days <- dep[, .(date = seq(start_date + 1, end_date - 1, by = "1 day")), by = .(uid)]
valid_days[, day_index := seq_len(.N), by = uid]
valid_days[, window := (day_index - 1) %/% 10 + 1]
vw <- valid_days[, .N, by = .(uid, window)][N == 10]
valid_days <- merge(valid_days, vw[, .(uid, window)], by = c("uid", "window"))
stopifnot(max(valid_days$window) == ncol(ref[[1]]))
master <- unique(valid_days[, .(uid, window)])
master[, row := match(uid, uids)]

# ---- sequences (read.csv on the needed columns, so id columns convert exactly as hers did) -------------
hdr <- names(read.csv(file.path(D, "combined_sequences_all.csv"), nrows = 1))
need <- c("project_id", "deployment_id", "genus", "species", "start_time")
cc <- ifelse(hdr %in% need, NA, "NULL")
sq <- as.data.table(read.csv(file.path(D, "combined_sequences_all.csv"), colClasses = cc))
sq[, sciname := paste(genus, species)]
sq[, uid := paste(project_id, deployment_id, sep = "___")]
sq[, date := as.Date(start_time)]
det_all <- unique(sq[, .(uid, sciname, date)])
det_win <- unique(merge(det_all, valid_days[, .(uid, date, window)], by = c("uid", "date"))[, .(uid, window, sciname)])
cat("sequences:", nrow(sq), " detection-days:", nrow(det_all), " detection-windows:", nrow(det_win), "\n")

build_y <- function(names) {
  y <- matrix(NA_integer_, nrow = length(uids), ncol = ncol(ref[[1]]), dimnames = dimnames(ref[[1]]))
  y[cbind(master$row, master$window)] <- 0L
  h <- unique(det_win[sciname %in% names, .(uid, window)])
  if (nrow(h)) y[cbind(match(h$uid, uids), h$window)] <- 1L
  y
}
same_y <- function(a, b) identical(dim(a), dim(b)) && identical(dimnames(a), dimnames(b)) &&
  identical(is.na(a), is.na(b)) && all(as.integer(a[!is.na(a)]) == as.integer(b[!is.na(b)]))

# ---- GATE: every species, single original name -----------------------------------------------------------
gate <- rbindlist(lapply(names(u), function(k) {
  y <- build_y(k); a <- ref[[k]]
  data.table(name = k, identical = same_y(y, a), n_det_hers = sum(rowSums(a, na.rm = TRUE) > 0),
             n_det_rebuilt = sum(rowSums(y, na.rm = TRUE) > 0), cells_differ = sum(xor(is.na(a), is.na(y))) +
               sum(a != y, na.rm = TRUE))
}))
print(gate)
fwrite(gate, file.path(W, "rebuild_gate_umflist.csv"))
if (!all(gate$identical)) stop("GATE FAILED: single-name rebuild does not reproduce umflist exactly for: ",
                               paste(gate[identical == FALSE, name], collapse = ", "), " -- nothing written")
cat("\nGATE PASSED: single-name rebuild reproduces umflist exactly for all", nrow(gate), "species\n")

# ---- v2: SCOPE species from the crosswalk's camera name list ---------------------------------------------
cw <- fread(file.path(W, "species_name_crosswalk.csv")); map <- fread(file.path(W, "camera_name_map.csv"))
stopifnot("wi_names" %in% names(cw), all(SCOPE %in% cw$species))
u2 <- u; used <- setNames(as.list(names(u)), names(u))
yr <- as.integer(attr(u[[1]], "siteCovs")$year)
cmp <- list()
for (sp in SCOPE) {
  k <- cw[species == sp, wi_name]; nms <- map[species == sp, sciname]
  stopifnot(k %in% names(u), k %in% nms)
  y2 <- build_y(nms); y1 <- ref[[k]]
  stopifnot(identical(dimnames(y2), dimnames(y1)), identical(is.na(y2), is.na(y1)), all(y2 >= y1, na.rm = TRUE))
  attr(u2[[k]], "y") <- y2; used[[k]] <- nms
  d1 <- rowSums(y1, na.rm = TRUE) > 0; d2 <- rowSums(y2, na.rm = TRUE) > 0
  cmp[[sp]] <- data.table(species = sp, year = sort(unique(yr)))[, `:=`(
    deployments = sapply(year, function(t) sum(yr == t)),
    detecting_umflist = sapply(year, function(t) sum(d1 & yr == t)),
    detecting_v2 = sapply(year, function(t) sum(d2 & yr == t)))]
  cat(sprintf("%s: %s -> names [%s]; detecting deployments %d -> %d\n", sp, k, paste(nms, collapse = " | "), sum(d1), sum(d2)))
}
cmp <- rbindlist(cmp); print(cmp[detecting_umflist > 0 | detecting_v2 > 0]); fwrite(cmp, file.path(W, "rebuild_v2_by_year.csv"))
stopifnot(identical(attr(u2[[1]], "siteCovs"), attr(u[[1]], "siteCovs")))
for (k in setdiff(names(u), cw[species %in% SCOPE, wi_name])) stopifnot(identical(attr(u2[[k]], "y"), ref[[k]]))
attr(u2, "wi_names_used") <- used
attr(u2, "provenance") <- sprintf("rebuild_detections.R %s: y rebuilt for %s from camera_name_map.csv; all other y identical to umflist.RDS; gate passed for %d species",
                                  format(Sys.time()), paste(SCOPE, collapse = ","), nrow(gate))
saveRDS(u2, OUT_UMF); cat("wrote", OUT_UMF, "\n")

# ---- standing check: detecting-site ratio, keyed on MDD-resolved names ------------------------------------
raw_det <- function(nms) uniqueN(sq[sciname %in% nms, uid])
rat <- rbindlist(lapply(seq_len(nrow(cw)), function(i) {
  sp <- cw$species[i]; k <- cw$wi_name[i]; nms <- map[species == sp, sciname]
  if (!k %in% names(u)) return(NULL)
  r <- raw_det(nms)
  data.table(species = sp, raw_detecting_deployments = r,
             umflist_detecting = sum(rowSums(ref[[k]], na.rm = TRUE) > 0),
             v2_detecting = sum(rowSums(attr(u2[[k]], "y"), na.rm = TRUE) > 0))
}))
rat[, `:=`(ratio_umflist = round(umflist_detecting / raw_detecting_deployments, 3),
           ratio_v2 = round(v2_detecting / raw_detecting_deployments, 3))]
rat[, flag_umflist := ifelse(ratio_umflist < RATIO_MIN, "LOW", "")][, flag_v2 := ifelse(ratio_v2 < RATIO_MIN, "LOW", "")]
cat("\nDETECTING-SITE RATIO CHECK (flag <", RATIO_MIN, "); median umflist", median(rat$ratio_umflist), "\n")
print(rat[order(ratio_umflist)])
fwrite(rat, file.path(W, "camera_detect_ratio_check.csv"))
cat("DONE\n")
