#!/usr/bin/env Rscript
# =============================================================================
# build_site_array.R -- STAGE B.
# (1) Build inputs$site_array: the 5 km array-year label per camera, in the
#     SAME order as constants_list$cell (i.e. the bobcat bundle's site order).
# (2) Report the realized bobcat array structure post-range-mask.
# (3) DECISIVE FEASIBILITY CHECK: build_reduced_constants() keeps only
#     n_site_keep = 700 of 20,531 sites. Measure what the array structure
#     looks like AFTER that reduction and after drop_singletons = TRUE, since
#     that is what the sweep actually fits.
# Writes prepped_sim_inputs_with_array.RDS (a COPY; original untouched).
# =============================================================================
suppressMessages({ library(tidyverse) })
PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
SPT  <- file.path(PROJ, "prototype_spatial_trend")
DIAG <- "/home/rwkays/isdm/diag"
source(file.path(SPT, "sim_helpers.R"))
source(file.path(SPT, "sim_helpers_array.R"))

cat("=== load ===\n")
prepped <- readRDS(file.path(SPT, "prepped_sim_inputs.RDS"))
cl <- prepped$constants_list
cat("nsite:", cl$nsite, "\n")

b  <- readRDS(file.path(PROJ, "HPC", "bobcat", "input_data_bobcat.RDS"))
oc <- b$real_data$occ.covs
cat("occ.covs rows:", nrow(oc), " (must equal nsite)\n")
stopifnot(nrow(oc) == cl$nsite)

map <- readRDS(file.path(DIAG, "out_arraysplit", "dep_unit_map.RDS"))
cat("dep_unit_map rows:", nrow(map), "\n")

# ---- recover deployment_id LABELS ------------------------------------------
# occ.covs$deployment_id is an integer FACTOR CODE (range 1..24623 == nlevels
# of umflist siteCovs$deployment_id), NOT a label -- terra/unmarked dropped the
# levels when the siteCovs data.frame was carried through. Joining on the codes
# matches nothing. Recover the labels by indexing the original factor levels;
# verified 100% match against dep_unit_map on (label, year).
# Coordinates+year also match 100% but have 2,065 duplicate keys in the map,
# so the label key is the unambiguous one and is used here.
umf <- readRDS("/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/umflist.RDS")
lev <- levels(umf[["Lynx rufus"]]@siteCovs$deployment_id)
rm(umf); gc()
code <- as.integer(oc$deployment_id)
stopifnot(all(code >= 1 & code <= length(lev), na.rm = TRUE))
oc$deployment_id <- lev[code]
cat("recovered deployment_id labels, e.g.:", paste(head(oc$deployment_id, 2), collapse = ", "), "\n")

map$deployment_id <- as.character(map$deployment_id)
key_oc  <- paste(oc$deployment_id, oc$year, sep = "@@")
key_map <- paste(map$deployment_id, map$year, sep = "@@")
cat("duplicate keys in map:", sum(duplicated(key_map)), "\n")
m <- match(key_oc, key_map)
site_array <- map$unit_5[m]
cat("matched:", sum(!is.na(m)), "of", length(m),
    sprintf("(%.2f%%)", 100 * mean(!is.na(m))), "\n")
if (any(is.na(m))) {
  cat("UNMATCHED sample:\n"); print(utils::head(oc[is.na(m), c("deployment_id","year","subproject_name")], 5))
}

# ---- (2) full bobcat set structure -----------------------------------------
cat("\n=== bobcat FULL modelled set (nsite = 20,531), 5 km array-years ===\n")
tb <- table(site_array[!is.na(site_array)])
cat("  array-years represented:", length(tb), "\n")
cat("  cameras assigned:", sum(tb), "\n")
cat("  cameras/unit min/median/mean/max:", min(tb), "/", median(tb), "/",
    round(mean(tb), 2), "/", max(tb), "\n")
cat("  singleton units:", sum(tb == 1), sprintf("(%.1f%%)", 100*mean(tb == 1)), "\n")
cat("  cameras in units >= 2:", sum(tb[tb >= 2]),
    sprintf("(%.1f%% retained)", 100*sum(tb[tb>=2])/sum(tb)), "\n")
print(table(cut(as.integer(tb), c(0,1,2,5,10,20,50,Inf),
                labels = c("1","2","3-5","6-10","11-20","21-50",">50"))))

# ---- (3) THE FEASIBILITY CHECK: after build_reduced_constants ---------------
cat("\n=== AFTER build_reduced_constants (n_site_keep = 700) ===\n")
DESIGN_SEED <- 20260823L
reduced <- build_reduced_constants(
  cl = prepped$constants_list, inat_effort_real = prepped$inat_effort_real,
  real_y_template = prepped$real_y_template, cell100_geo = prepped$cell100_geo,
  seed = DESIGN_SEED)
sk <- reduced$site_keep
cat("  site_keep length:", length(sk), " unique:", length(unique(sk)), "\n")
sa_red <- site_array[sk]
cat("  non-NA labels among kept sites:", sum(!is.na(sa_red)), "\n")

aid <- assign_arrays_from_field(array_field = sa_red,
                               site_year = reduced$constants$year_occ,
                               drop_singletons = TRUE)
cat("  sites surviving drop_singletons:", sum(!is.na(aid)),
    sprintf("(%.1f%% of %d kept sites)", 100*mean(!is.na(aid)), length(aid)), "\n")
cat("  arrays formed:", length(unique(aid[!is.na(aid)])), "\n")
if (sum(!is.na(aid)) > 0) {
  tb2 <- table(aid[!is.na(aid)])
  cat("  cameras/array min/median/mean/max:", min(tb2), "/", median(tb2), "/",
      round(mean(tb2), 2), "/", max(tb2), "\n")
  print(table(cut(as.integer(tb2), c(0,1,2,5,10,20,Inf),
                  labels = c("1","2","3-5","6-10","11-20",">20"))))
} else cat("  *** NO ARRAYS SURVIVE -- array arms cannot be fit at this reduction ***\n")

cat("\n  >>> camera-level arm fits", length(sk), "sites;",
    "array arms would fit", sum(!is.na(aid)), "cameras in",
    length(unique(aid[!is.na(aid)])), "arrays <<<\n")

# ---- write the augmented inputs (COPY, original untouched) ------------------
prepped$site_array <- site_array
outf <- file.path(SPT, "prepped_sim_inputs_with_array.RDS")
saveRDS(prepped, outf)
cat("\nwrote", outf, "\n")
saveRDS(site_array, file.path(DIAG, "out_arraysplit", "site_array_bobcat.RDS"))
cat("DONE\n")
