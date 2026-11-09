#!/usr/bin/env Rscript
# apply_campaign_exclusion.R -- decision M2: exclude organised-campaign footprints from the iNaturalist stream.
#
# For bobcat, coyote, California ground squirrel and northern raccoon, removes the iNaturalist data in the
# campaign cells listed in table_campaign_strong_effects.csv (column cell50_list), in every year. Camera data
# in those cells are untouched, and the cells stay on the model grid, so abundance is still estimated there.
#
# How a cell is removed: inat_y is set to NA and inat_effort to 0 for that cell, and n_cells_year /
# inat_cells_by_year are rebuilt exactly as run_fleet_prep.R builds them (cells with effort > 0). The
# iNaturalist likelihood only visits cells in inat_cells_by_year, so an excluded cell contributes nothing
# and is not read as an observed zero count. NA + effort 0 is the state the builder gives any no-effort cell.
#
# Source bundles are read only. Output: <OUT_DIR>/<run>/input_data_<run>.RDS, refused if it already exists.
# Env: CAMPAIGN_CSV, OUT_DIR (required); HPC_DIR, UNION_DIR (default: project locations on Hazel).
# Env ONLY_SPECIES (optional): build one species label and merge its rows into the existing output tables.
# Red fox (Roland, 2026-11-09): excluded for the five cells of the Twin Cities Coyote and Fox Project, where its
# selectivity is 0.54 (0.44-0.66). The cells are read from that project's row, listed under Coyote in the table.

env <- function(k, d = NA) { v <- Sys.getenv(k, ""); if (nzchar(v)) v else if (is.na(d)) stop("set ", k) else d }
PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
CSV  <- env("CAMPAIGN_CSV"); OUT <- env("OUT_DIR")
HPC  <- env("HPC_DIR", file.path(PROJ, "HPC"))
UG   <- env("UNION_DIR", file.path(PROJ, "output", "inat_grids_union_fleet34"))

# label in the campaign table -> token, expected cell count, wide lon/lat box the cells must fall inside
SPP <- list(
  "Bobcat"                     = list(tok = "bobcat",                     n = 47L, lon = c(-124, -112),   lat = c(30, 38)),
  "Coyote"                     = list(tok = "coyote",                     n = 5L,  lon = c(-95.5, -91.5), lat = c(43.5, 46.5)),
  "California ground squirrel" = list(tok = "california_ground_squirrel", n = 3L,  lon = c(-123.5, -120.5), lat = c(36.5, 39)),
  "Northern raccoon"           = list(tok = "northern_raccoon",           n = 3L,  lon = c(-124, -120.5), lat = c(45.5, 48.5)),
  "Red fox"                    = list(tok = "red_fox",                    n = 5L,  lon = c(-95.5, -91.5), lat = c(43.5, 46.5),
                                      from_species = "Coyote", from_campaign = "Twin Cities Coyote and Fox Project"))
WINDOWS <- c("", "_win10", "_win5")

tab <- read.csv(CSV, stringsAsFactors = FALSE, check.names = FALSE)
stopifnot(setequal(unique(tab$species), setdiff(names(SPP), "Red fox")))
pairs_out <- list(); counts_out <- list()
ONLY <- Sys.getenv("ONLY_SPECIES", "")
if (nzchar(ONLY)) stopifnot(ONLY %in% names(SPP))

for (label in if (nzchar(ONLY)) ONLY else names(SPP)) {
  s <- SPP[[label]]
  rows <- if (is.null(s$from_species)) tab[tab$species == label, ] else
    tab[tab$species == s$from_species & tab$campaign == s$from_campaign, ]
  stopifnot(nrow(rows) >= 1)
  cells <- sort(unique(as.integer(unlist(strsplit(rows$cell50_list, ";", fixed = TRUE)))))
  if (length(cells) != s$n || anyNA(cells)) stop(label, ": expected ", s$n, " cells, table gives ", length(cells))

  geo <- read.csv(file.path(UG, paste0(s$tok, "_union_cells.csv")))
  g <- geo[match(cells, geo$cell50), ]
  if (anyNA(g$lon)) stop(label, ": campaign cell not on the union grid: ", paste(cells[is.na(g$lon)], collapse = ","))
  if (any(g$lon < s$lon[1] | g$lon > s$lon[2] | g$lat < s$lat[1] | g$lat > s$lat[2]))
    stop(label, ": campaign cells fall outside the expected region; cell50 numbering may not match")
  cat(sprintf("%s: %d cells, lon %.2f to %.2f, lat %.2f to %.2f\n", label, length(cells),
              min(g$lon), max(g$lon), min(g$lat), max(g$lat)))

  for (w in WINDOWS) {
    run <- paste0(s$tok, "_v2u_national_scalar", w)
    src <- file.path(HPC, run, paste0("input_data_", run, ".RDS"))
    dst_dir <- file.path(OUT, run); dst <- file.path(dst_dir, paste0("input_data_", run, ".RDS"))
    if (file.exists(dst)) stop("refusing to overwrite ", dst)
    a <- readRDS(src); b <- a; cl <- b$constants_list
    Y <- b$inat_y; E <- b$inat_effort
    stopifnot(!is.null(rownames(Y)), identical(rownames(Y), rownames(E)), identical(colnames(Y), colnames(E)),
              identical(dim(Y), dim(E)), nrow(Y) == cl$ncell50)

    # the bundle must already follow the builder's rule, or the rebuild below would change other cells
    nyr <- ncol(Y)
    index_of <- function(eff) lapply(seq_len(nyr), function(t) which(eff[, t] > 0))
    old_idx <- index_of(E)
    for (t in seq_len(nyr)) stopifnot(cl$n_cells_year[t] == length(old_idx[[t]]),
                                      identical(as.integer(cl$inat_cells_by_year[seq_len(cl$n_cells_year[t]), t]),
                                                as.integer(old_idx[[t]])))
    in_lik <- function(M, idx) sum(vapply(seq_len(nyr), function(t) sum(M[idx[[t]], t]), 0))

    row_ids <- as.integer(rownames(Y)); drop <- which(row_ids %in% cells)
    # a campaign cell can legitimately be absent from a shorter window's grid; it is recorded as such
    before <- in_lik(Y, old_idx); cy_before <- sum(lengths(old_idx))
    stopifnot(!anyNA(E), !any(is.na(Y) & E > 0))
    rec_by_cell <- setNames(rowSums(Y[drop, , drop = FALSE], na.rm = TRUE), row_ids[drop])
    cy_by_cell  <- setNames(rowSums(E[drop, , drop = FALSE] > 0), row_ids[drop])

    Y[drop, ] <- NA; E[drop, ] <- 0
    new_idx <- index_of(E)
    n_new <- lengths(new_idx); m <- matrix(NA_integer_, max(n_new), nyr)
    for (t in seq_len(nyr)) m[seq_len(n_new[t]), t] <- new_idx[[t]]
    storage.mode(m) <- storage.mode(cl$inat_cells_by_year)

    b$inat_y <- Y; b$inat_effort <- E
    b$constants_list$n_cells_year <- as.integer(n_new)
    b$constants_list$inat_cells_by_year <- m
    after <- in_lik(Y, new_idx); cy_after <- sum(n_new)
    if (before - after != sum(rec_by_cell)) stop(run, ": removed records do not reconcile")

    b$campaign_exclusion <- list(decision = "M2", campaigns = rows$campaign, cell50_removed = row_ids[drop],
                                 inat_records_before = before, inat_records_after = after,
                                 source_bundle = src, built = format(Sys.time()))

    # everything except the four iNaturalist fields must be unchanged, cameras above all
    same <- setdiff(names(a), c("inat_y", "inat_effort", "constants_list"))
    for (n in same) stopifnot(identical(a[[n]], b[[n]]))
    cn <- setdiff(names(a$constants_list), c("n_cells_year", "inat_cells_by_year"))
    for (n in cn) stopifnot(identical(a$constants_list[[n]], b$constants_list[[n]]))
    keep <- setdiff(seq_len(nrow(Y)), drop)
    stopifnot(identical(a$inat_y[keep, ], Y[keep, ]), identical(a$inat_effort[keep, ], E[keep, ]),
              all(is.na(Y[drop, ])), all(E[drop, ] == 0))

    dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(b, dst)

    pairs_out[[run]] <- data.frame(
      run = run, species = label, cell50 = cells, lon = round(g$lon, 3), lat = round(g$lat, 3),
      in_bundle_grid = cells %in% row_ids,
      inat_records_removed = as.numeric(ifelse(cells %in% row_ids, rec_by_cell[as.character(cells)], 0)),
      cell_years_removed = as.integer(ifelse(cells %in% row_ids, cy_by_cell[as.character(cells)], 0)))
    counts_out[[run]] <- data.frame(
      run = run, species = label, window = if (w == "") "full_record" else sub("_", "", w),
      campaign_cells_listed = length(cells), campaign_cells_in_bundle_grid = length(drop),
      inat_records_before = before, inat_records_after = after, inat_records_removed = before - after,
      pct_removed = round(100 * (before - after) / before, 2),
      cell_years_before = cy_before, cell_years_after = cy_after)
    cat(sprintf("  %-55s records %d -> %d (-%.2f%%), cells in grid %d/%d\n", run, before, after,
                100 * (before - after) / before, length(drop), length(cells)))
  }
}
save_table <- function(new, name) {
  f <- file.path(OUT, name)
  if (nzchar(ONLY) && file.exists(f)) {
    old <- read.csv(f, stringsAsFactors = FALSE, check.names = FALSE)
    new <- rbind(old[!old$run %in% new$run, ], new)
  }
  write.csv(new, f, row.names = FALSE)
}
save_table(do.call(rbind, pairs_out), "table_campaign_cells_removed.csv")
save_table(do.call(rbind, counts_out), "table_campaign_exclusion_counts.csv")
cat("DONE\n")
