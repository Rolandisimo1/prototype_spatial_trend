#!/usr/bin/env Rscript
# =============================================================================
# crosswalk_camera_names.R -- STAGE 2 of the species-name crosswalk (2026-10-08).
# Stage 1 is build_name_crosswalk.R (one name per source). This stage fixes the schema hole that
# lost 47.9% of the peccary camera record: one wi_name slot cannot represent a source that
# changed names mid-series (Wildlife Insights: Pecari tajacu through 2023, Dicotyles tajacu from 2024).
#
# Adds to species_name_crosswalk.csv (backup kept), without touching existing columns:
#   mdd_id, mdd_version   -- MDD id of mdd_binomial
#   wi_names              -- EVERY camera name for the species, with the years it was used:
#                            "Pecari tajacu[2015-2023]; Dicotyles tajacu[2024-2025]". A name is listed if
#                            it resolves to the species' MDD id (rename), or if a recorded decision in
#                            manual_taxonomy_decisions.csv reassigns it (e.g. Cervus elaphus -> elk).
#   n_wi_names
# and writes camera_name_map.csv (one row per species x camera name) and camera_name_checks.csv.
#
# BUILD CHECKS (each stops the build):
#   A. every mdd_binomial is an MDD accepted name (the peccary row held a synonym, Pecari tajacu).
#   B. the stage-1 wi_name is among the species' wi_names.
#   C. any camera common-name label that carries a fleet species' name AND another scientific name
#      fails unless that other name is listed in the species' wi_names or covered by a recorded
#      decision (reassign / no_change_confirmed / accept_single_name). Labels are matched
#      case-insensitively only to FIND candidates; nothing is ever keyed or merged on common name.
# The consumer-side check (a species with n_wi_names > 1 must use a camera record built from all of
# them) lives in prep_camera_sites.R.
# =============================================================================
suppressPackageStartupMessages(library(data.table))
W   <- "/home/rwkays/isdm/fleet_v2u"
source(file.path(W, "mdd_lib.R"))
CW  <- Sys.getenv("CROSSWALK", file.path(W, "species_name_crosswalk.csv"))
DEC <- Sys.getenv("TAXONOMY_DECISIONS", file.path(W, "manual_taxonomy_decisions.csv"))
SEQ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/data/combined_sequences_all.csv"
options(width = 250)

mdd <- load_mdd(); cat("MDD", mdd$version, ":", length(mdd$acc), "accepted species\n")
cw  <- fread(CW)
dec <- fread(DEC); dec <- dec[scope == "camera_sequences"]
cat("decisions read:", nrow(dec), "from", DEC, "\n")

# ---- A. accepted-name assertion -------------------------------------------------------------
badA <- cw[!mdd_binomial %in% names(mdd$acc)]
if (nrow(badA)) {
  r <- mdd_resolve(badA$mdd_binomial, mdd)
  print(merge(badA[, .(species, mdd_binomial)], r, by.x = "mdd_binomial", by.y = "verbatim"))
  stop("CHECK A FAILED -- mdd_binomial is not an MDD ", mdd$version, " accepted name for: ",
       paste(sprintf("%s (%s -> should be %s)", badA$species, badA$mdd_binomial,
                     r$resolved[match(badA$mdd_binomial, r$verbatim)]), collapse = ", "))
}
cw[, mdd_id := unname(mdd$acc[mdd_binomial])]; cw[, mdd_version := mdd$version]
cat("CHECK A passed: all", nrow(cw), "mdd_binomial are MDD", mdd$version, "accepted names\n")

# ---- raw camera names ----------------------------------------------------------------------------
sq <- fread(SEQ, select = c("project_id", "deployment_id", "genus", "species", "common_name", "year"))
sq[, sciname := trimws(paste(genus, species))]
sq[, uid := paste(project_id, deployment_id, sep = "___")]
nm <- sq[, .(n_sequences = .N, n_deployments = uniqueN(uid), first_year = min(year), last_year = max(year)), by = sciname]
res <- mdd_resolve(nm$sciname, mdd)
nm <- merge(nm, res, by.x = "sciname", by.y = "verbatim", all.x = TRUE)

# ---- per-species name list -----------------------------------------------------------------------
map <- rbindlist(lapply(seq_len(nrow(cw)), function(k) {
  sp <- cw$species[k]; id <- cw$mdd_id[k]
  a <- nm[!is.na(mdd_id) & mdd_id == id, .(sciname, n_sequences, n_deployments, first_year, last_year, match_type,
                                           basis = ifelse(match_type == "accepted_exact", "accepted name", "MDD synonym (rename)"))]
  d <- dec[action == "reassign" & to_mdd_id == id]
  if (nrow(d)) {
    b <- nm[sciname %in% d$from_name, .(sciname, n_sequences, n_deployments, first_year, last_year, match_type)]
    b[, basis := paste0("REASSIGNED by decision (", d$authority[match(sciname, d$from_name)], ", ", d$date[match(sciname, d$from_name)], ")")]
    a <- rbind(a, b)
  }
  if (!nrow(a)) stop("no camera name resolves to the MDD id of ", sp, " (", cw$mdd_binomial[k], ")")
  cbind(species = sp, mdd_id = id, a[order(first_year, sciname)])
}))
map[, mdd_version := mdd$version]
ws <- map[, .(wi_names = paste(sprintf("%s[%d-%d]", sciname, first_year, last_year), collapse = "; "),
              n_wi_names = .N, wi_name_list = list(sciname)), by = species]
cw <- merge(cw, ws[, .(species, wi_names, n_wi_names)], by = "species", sort = FALSE)

# ---- B. stage-1 wi_name must be among the list ---------------------------------------------------
badB <- cw[!mapply(function(w, s) w %in% ws[species == s, wi_name_list[[1]]], wi_name, species)]
if (nrow(badB)) stop("CHECK B FAILED -- stage-1 wi_name not among the MDD-resolved camera names for: ",
                     paste(badB$species, collapse = ", "))
cat("CHECK B passed\n")

# ---- C. common-name labels carrying more than one scientific name --------------------------------
lab <- unique(sq[, .(label = tolower(trimws(common_name)), sciname)])
lab <- lab[, if (uniqueN(sciname) > 1) .SD, by = label]
covered <- c(dec$from_name)                                  # any recorded decision covers its from_name
chk <- rbindlist(lapply(seq_len(nrow(ws)), function(k) {
  own <- ws$wi_name_list[[k]]
  hits <- lab[label %in% lab[sciname %in% own, label] & !sciname %in% own]
  if (!nrow(hits)) return(NULL)
  hits[, .(species = ws$species[k], label, other_name = sciname,
           status = ifelse(sciname %in% covered, paste("covered by decision:", dec$action[match(sciname, dec$from_name)]), "UNLISTED"))]
}))
if (nrow(chk)) { cat("\nCommon-name labels shared with another scientific name:\n"); print(chk) }
fwrite(if (nrow(chk)) chk else data.table(species = character(), label = character(), other_name = character(), status = character()),
       file.path(W, "camera_name_checks.csv"))
badC <- if (nrow(chk)) chk[status == "UNLISTED"] else chk
if (nrow(badC)) stop("CHECK C FAILED -- a fleet species' camera label also carries an unlisted scientific name: ",
                     paste(unique(sprintf("%s: '%s' also %s", badC$species, badC$label, badC$other_name)), collapse = "; "),
                     ". List it in wi_names (via MDD) or record a decision in ", DEC)
cat("CHECK C passed\n")

# ---- write ---------------------------------------------------------------------------------------
setcolorder(cw, c(setdiff(names(cw), c("mdd_id", "mdd_version", "wi_names", "n_wi_names", "note")),
                  "mdd_id", "mdd_version", "wi_names", "n_wi_names", "note"))
bak <- paste0(CW, ".bak_", format(Sys.time(), "%Y%m%d_%H%M%S")); file.copy(CW, bak)
fwrite(cw, CW); fwrite(map, file.path(W, "camera_name_map.csv"))
cat("\nwrote", CW, "(backup", bak, ") and camera_name_map.csv\n")
print(cw[n_wi_names > 1, .(species, wi_name, wi_names)])
# Informational (rwkays 2026-10-08): the common-name field is not clean ("gray fox" 17 vs "grey fox" 14,496 under one
# scientific name), which is why nothing keys on it. Record every fleet name's label variants.
lv <- sq[sciname %in% map$sciname, .N, by = .(sciname, common_name)][, if (.N > 1) .SD, by = sciname][order(sciname, -N)]
fwrite(lv, file.path(W, "camera_label_variants.csv")); cat("fleet camera names with >1 common-name label:", uniqueN(lv$sciname), "\n")
cat("unresolved camera names (all, not only fleet):", nm[is.na(mdd_id), .N], "\n")
fwrite(nm[is.na(mdd_id)][order(-n_sequences)], file.path(W, "camera_names_unresolved.csv"))
cat("DONE\n")
