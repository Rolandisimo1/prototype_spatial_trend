#!/usr/bin/env Rscript
# =============================================================================
# diagnose_bobcat_v2b_ecoregion_car.R
# -----------------------------------------------------------------------------
# Why do 12 contiguous link_occ_intercept cells (854-855, 873-875, 883-885,
# 891-894) plus MCMT_tau fail the R-hat gate in bobcat_v2b_ecoregion, when the
# trend parameters all converge at R-hat <= 1.004?
#
# Two hypotheses to separate, per the 8/27 pull note:
#   (A) DATA-SPARSE  -- these cells carry little or no camera/iNat information,
#       so the CAR prior is doing nearly all the work and the posterior is
#       weakly identified.
#   (B) MASK EDGE    -- these cells sit on the boundary of the retained cell
#       set, so they have few CAR neighbours and inherit little smoothing.
# The two are separable: (A) is about data attached to the cell, (B) is about
# the cell's degree in the CAR adjacency graph. They can also both hold.
#
# The comparison that matters is against the REST OF THE FIELD, not against
# zero -- plenty of cells are data-poor without failing R-hat. So every
# quantity is reported as the offenders' value AND its percentile in the
# all-cell distribution.
#
# Read-only. Loads input_data (small) and the R-hat CSV; does NOT load chains.
# =============================================================================

PROJ <- "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final"
SPECIES <- "bobcat_v2b_ecoregion"
D <- file.path(PROJ, "HPC", SPECIES)

OFFENDERS <- c(854, 855, 873, 874, 875, 883, 884, 885, 891, 892, 893, 894)

cat("=== loading input_data ===\n")
input_data <- readRDS(file.path(D, paste0("input_data_", SPECIES, ".RDS")))
cl <- input_data$constants_list
cat("fields:", paste(names(cl), collapse = ", "), "\n\n")

ncell100 <- cl$ncell100
stopifnot(all(OFFENDERS <= ncell100))

# ---- (A) data support per cell100 ------------------------------------------
# camera sites attached to each cell100
n_cam <- as.integer(table(factor(cl$cell, levels = seq_len(ncell100))))

# camera DETECTIONS per cell100 -- a cell can hold sites that never detected,
# which is a different kind of sparse than holding no sites at all.
y <- input_data$data$y
if (is.null(y)) y <- input_data$y
n_det <- rep(NA_integer_, ncell100)
if (!is.null(y)) {
  det_per_site <- rowSums(y, na.rm = TRUE)
  n_det <- as.integer(tapply(det_per_site, factor(cl$cell, levels = seq_len(ncell100)),
                             sum, default = 0))
  n_det[is.na(n_det)] <- 0L
}

# iNat monitored cell50s mapping up to each cell100
monitored_g  <- unique(na.omit(as.vector(cl$inat_cells_by_year)))
n_inat <- as.integer(table(factor(cl$inat_cell100[monitored_g],
                                  levels = seq_len(ncell100))))

# ---- (B) CAR graph degree ---------------------------------------------------
deg <- cl$num
stopifnot(length(deg) == ncell100)
# neighbour indices per cell, from the flattened adj/num CAR representation
adj_start <- cumsum(c(1L, deg))[seq_len(ncell100)]
nbrs <- function(i) {
  if (deg[i] == 0) return(integer(0))
  cl$adj[adj_start[i]:(adj_start[i] + deg[i] - 1L)]
}

# a cell is "interior" if it and all its neighbours have data; an edge cell
# has neighbours that are themselves data-empty, so smoothing borrows nothing.
nbr_mean_cam <- vapply(seq_len(ncell100),
                       function(i) { nb <- nbrs(i); if (!length(nb)) NA_real_ else mean(n_cam[nb]) },
                       numeric(1))
nbr_frac_empty <- vapply(seq_len(ncell100),
                         function(i) { nb <- nbrs(i); if (!length(nb)) NA_real_ else mean(n_cam[nb] == 0) },
                         numeric(1))

pct <- function(v, x) round(100 * mean(v <= x, na.rm = TRUE), 1)

# ---- per-offender table -----------------------------------------------------
tab <- data.frame(
  cell100        = OFFENDERS,
  ecoregion      = cl$ecoregion_of_cell100[OFFENDERS],
  n_cam_sites    = n_cam[OFFENDERS],
  n_cam_det      = n_det[OFFENDERS],
  n_inat_cell50  = n_inat[OFFENDERS],
  car_degree     = deg[OFFENDERS],
  nbr_mean_cam   = round(nbr_mean_cam[OFFENDERS], 2),
  nbr_frac_empty = round(nbr_frac_empty[OFFENDERS], 2),
  pct_cam        = vapply(OFFENDERS, function(i) pct(n_cam, n_cam[i]), numeric(1)),
  pct_degree     = vapply(OFFENDERS, function(i) pct(deg, deg[i]), numeric(1)),
  stringsAsFactors = FALSE
)
cat("=== the 12 offending CAR cells ===\n")
print(tab, row.names = FALSE)

# ---- offenders vs the rest of the field ------------------------------------
rest <- setdiff(seq_len(ncell100), OFFENDERS)
cmp <- function(nm, v) data.frame(
  quantity      = nm,
  offenders_med = round(median(v[OFFENDERS], na.rm = TRUE), 3),
  rest_med      = round(median(v[rest], na.rm = TRUE), 3),
  offenders_mean= round(mean(v[OFFENDERS], na.rm = TRUE), 3),
  rest_mean     = round(mean(v[rest], na.rm = TRUE), 3),
  stringsAsFactors = FALSE)

cat("\n=== offenders vs rest of field ===\n")
print(rbind(
  cmp("n_cam_sites",    as.numeric(n_cam)),
  cmp("n_cam_det",      as.numeric(n_det)),
  cmp("n_inat_cell50",  as.numeric(n_inat)),
  cmp("car_degree",     as.numeric(deg)),
  cmp("nbr_mean_cam",   nbr_mean_cam),
  cmp("nbr_frac_empty", nbr_frac_empty)
), row.names = FALSE)

# ---- how special are they, really? -----------------------------------------
# If N cells are as data-poor / as low-degree as the offenders but converge
# fine, then that property alone is NOT the explanation.
cat("\n=== discriminating power of each hypothesis ===\n")
thr_cam <- max(n_cam[OFFENDERS]); thr_deg <- max(deg[OFFENDERS])
n_like_cam <- sum(n_cam[rest] <= thr_cam)
n_like_deg <- sum(deg[rest]   <= thr_deg)
cat(sprintf("cells with n_cam_sites <= %d (offender max): %d of %d others (%.1f%%) converge fine\n",
            thr_cam, n_like_cam, length(rest), 100 * n_like_cam / length(rest)))
cat(sprintf("cells with car_degree  <= %d (offender max): %d of %d others (%.1f%%) converge fine\n",
            thr_deg, n_like_deg, length(rest), 100 * n_like_deg / length(rest)))
cat("\nIf those percentages are large, the property is common among CONVERGED\n",
    "cells and cannot by itself explain the 12 -- look to contiguity (one\n",
    "spatially-connected weakly-identified patch) plus the MCMT_tau coupling.\n", sep = "")

# ---- are the 12 actually contiguous in the CAR graph? ----------------------
cat("\n=== contiguity of the offender set ===\n")
adjmat <- lapply(OFFENDERS, function(i) intersect(nbrs(i), OFFENDERS))
names(adjmat) <- OFFENDERS
for (k in names(adjmat))
  cat(sprintf("  cell %s: %d of its %d CAR neighbours are also offenders (%s)\n",
              k, length(adjmat[[k]]), deg[as.integer(k)],
              paste(adjmat[[k]], collapse = ",")))
n_internal <- sum(vapply(adjmat, length, integer(1)))
cat(sprintf("\ninternal offender-offender edges: %d (of %d total offender edges)\n",
            n_internal, sum(deg[OFFENDERS])))

# ---- MCMT_tau coupling ------------------------------------------------------
# MCMT_tau is the precision of the MCMT SVC CAR field. If the offending cells
# are also where MCMT varies most sharply, the intercept field and the MCMT
# field are competing to explain the same local signal -- which would make
# MCMT_tau's own R-hat (1.182, the worst in the fit) the driver rather than a
# coincidence.
MCMT <- cl$MCMT
if (!is.null(MCMT) && length(MCMT) == length(cl$cell)) {
  mcmt_by_cell <- tapply(MCMT, factor(cl$cell, levels = seq_len(ncell100)),
                         function(v) if (length(v)) sd(v) else NA_real_)
  cat("\n=== MCMT within-cell SD (proxy for local MCMT gradient) ===\n")
  cat("offenders median:", round(median(mcmt_by_cell[OFFENDERS], na.rm = TRUE), 4),
      "  rest median:",   round(median(mcmt_by_cell[rest], na.rm = TRUE), 4), "\n")
} else {
  cat("\n(MCMT not available per-site in constants_list -- skipping the MCMT coupling check)\n")
}

cat("\nDIAGNOSTIC DONE\n")
