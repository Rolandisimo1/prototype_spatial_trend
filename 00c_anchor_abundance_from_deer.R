# 00c_anchor_abundance_from_deer.R
# ---------------------------------------------------------------------------
# OPTIONAL, RUN ON HAZEL (login node or a compute node that can read $PROJ).
# Anchors the TOP of the abundance ladder to a REAL common species instead of
# a guessed multiplier. White-tailed deer is our most-abundant fleet species.
#
# IMPORTANT: deer on Hazel was DATA-PREP ONLY -- it was never fit, so there is
# NO deer posterior. We therefore anchor using OBSERVED DATA ratios (raw iNat
# counts and camera detection rates) from the deer PREP bundle vs bobcat's,
# NOT fitted parameters. Observed ratios are exactly what we need: they say how
# much more raw information a common species delivers per cell/site.
#
# Reads (READ-ONLY; never writes into Arielle's dirs):
#   $PROJ/HPC/white-tailed_deer/input_data_white-tailed_deer.RDS   (deer prep)
#   $PROJ/HPC/bobcat/input_data_bobcat.RDS                          (bobcat prep)
# Writes:
#   prototype_spatial_trend/abundance_anchor.RDS  (data.frame level/occ_shift/count_log_mult)
#
# If the deer bundle is absent/unreadable, this script exits cleanly and the
# sweep falls back to abundance_levels_default().
# ---------------------------------------------------------------------------

suppressMessages({ library(dplyr) })

PROJ <- Sys.getenv("PROJ",
  "/rsstu/users/j/jkpacifi/NSFiSDMs/Arielle_iSDM_temporal/integrated_code/Temporal_trend_final")
PROTO <- getwd()  # assume run from prototype_spatial_trend/

deer_path   <- file.path(PROJ, "HPC", "white-tailed_deer", "input_data_white-tailed_deer.RDS")
bobcat_path <- file.path(PROJ, "HPC", "bobcat", "input_data_bobcat.RDS")

if (!file.exists(deer_path)) {
  cat("Deer prep bundle NOT found at:\n  ", deer_path,
      "\n-> Skipping deer anchor; sweep will use abundance_levels_default().\n")
  quit(save = "no", status = 0)
}

# --- helper: pull raw observed information from a prep bundle ----------------
# The prep bundle is the input_data list handed to the model. We extract the
# iNat count matrix (cell50 x year) and the camera detection array (site x
# occasion) if present, and reduce to interpretable per-unit rates.
# CONFIRMED field names (input_data_bobcat.RDS, verified directly against the
# real prep bundle earlier in this project): top-level $inat_y (NOT $y_inat),
# camera detections at $real_data$y (NOT top-level $y or $data$y) --
# $constants_list (not $constants) has no y-data at all, just model
# constants. The tryCatch fallbacks below are kept in case the deer bundle
# (data-prep only, never fit) drifted from this structure, but the PRIMARY
# lookup now matches the real, verified field names rather than guessed ones.
observed_information <- function(path) {
  dat <- readRDS(path)
  # iNat counts per cell50 x year (the abundance signal for iNat)
  yi <- tryCatch(dat$inat_y, error = function(e) NULL)
  if (is.null(yi)) yi <- tryCatch(dat$y_inat, error = function(e) NULL)
  if (is.null(yi)) yi <- tryCatch(dat$constants_list$inat_y, error = function(e) NULL)
  # camera detection history site x occasion (0/1)
  yc <- tryCatch(dat$real_data$y, error = function(e) NULL)
  if (is.null(yc)) yc <- tryCatch(dat$y, error = function(e) NULL)
  if (is.null(yc)) yc <- tryCatch(dat$data$y, error = function(e) NULL)

  inat_per_cell <- if (!is.null(yi)) sum(yi, na.rm = TRUE) / max(nrow(as.matrix(yi)), 1) else NA
  cam_det_rate  <- if (!is.null(yc)) mean(yc > 0, na.rm = TRUE) else NA  # naive detection rate
  list(inat_per_cell = inat_per_cell, cam_det_rate = cam_det_rate,
       has_inat = !is.null(yi), has_cam = !is.null(yc))
}

deer <- observed_information(deer_path)
cat("Deer observed information:\n"); str(deer)

if (file.exists(bobcat_path)) {
  bob <- observed_information(bobcat_path)
  cat("Bobcat observed information:\n"); str(bob)
} else {
  cat("Bobcat prep bundle not found; cannot form ratio. Using default ladder.\n")
  quit(save = "no", status = 0)
}

# --- form empirical multipliers deer/bobcat ---------------------------------
# iNat count multiplier -> count_log_mult = log(ratio).
# Occupancy lift -> convert naive camera detection rates to cloglog and take the
# difference (additive shift on the cloglog intercept scale). Guard against
# saturation (rate==1) and empties.
safe_ratio <- function(a, b) if (is.finite(a) && is.finite(b) && b > 0) a / b else NA
cloglog <- function(p) log(-log(1 - pmin(pmax(p, 1e-4), 1 - 1e-4)))

count_mult  <- safe_ratio(deer$inat_per_cell, bob$inat_per_cell)
occ_shift_deer <- if (is.finite(deer$cam_det_rate) && is.finite(bob$cam_det_rate))
  cloglog(deer$cam_det_rate) - cloglog(bob$cam_det_rate) else NA

cat(sprintf("\nEmpirical deer/bobcat: iNat count mult = %.2f, cloglog occ shift = %.2f\n",
            count_mult, occ_shift_deer))

# --- build a 3-level ladder anchored so 'common' ~ deer ---------------------
# bobcat baseline = 0; 'common' = deer-derived; 'moderate' = halfway (geometric
# for the count multiplier, linear for the cloglog shift). Fall back to defaults
# for any piece that couldn't be computed.
occ_common   <- if (is.finite(occ_shift_deer)) occ_shift_deer else 1.5
count_common <- if (is.finite(count_mult) && count_mult > 0) log(count_mult) else log(8)

ladder <- data.frame(
  level          = c("bobcat_baseline", "moderate", "common_deerlike"),
  occ_shift      = c(0, occ_common / 2, occ_common),
  count_log_mult = c(0, count_common / 2, count_common),
  stringsAsFactors = FALSE
)
cat("\nAbundance ladder (deer-anchored):\n"); print(ladder)
saveRDS(ladder, file.path(PROTO, "abundance_anchor.RDS"))
cat("\nWrote abundance_anchor.RDS\n")
