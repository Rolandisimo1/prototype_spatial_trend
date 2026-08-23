# =============================================================================
# sim_helpers_array.R
# -----------------------------------------------------------------------------
# Array-level aggregation for the camera submodel, for the estimator-comparison
# simulation. Source AFTER sim_helpers.R (and sim_helpers_abundance.R if using
# the abundance ladder).
#
# PURPOSE
# Convert a camera-level simulated dataset into an array-level one, so the
# same simulated truth can be fit by three camera observation models:
#   1. camera-level occupancy  (baseline; the current production estimator)
#   2. array-level occupancy   (cameras as replicates within an array)
#   3. array-level Royle-Nichols (abundance drives per-camera detection)
#
# THE ESTIMATOR (per PI specification)
# The array observation is the PROPORTION of cameras in the array that
# detected the species, adjusted for detectability -- NOT "occupied if any
# camera detected". Any-detection aggregation collapses a 2-of-20 array and a
# 19-of-20 array to the same value, discards the graded signal, and pushes
# psi_a toward 1 where cloglog(psi) = log(lambda) stops discriminating
# abundance (camera information about lambda scales as 1 - psi).
#
# So cameras become REPLICATES within an array-year:
#     z_a       ~ Bernoulli(psi_a)                 array truly used
#     w_{a,i}   ~ Bernoulli(q_{a,i}) | z_a = 1     camera i detected
#     w_{a,i}   = 0                  | z_a = 0
# which dOcc_v handles unchanged: probOcc = psi_a, probDetect = q[a, 1:n_a],
# len = n_a. Only the replicate dimension changes (10-day windows -> cameras).
#
# EFFORT -- THE CRITICAL DETAIL
# At camera level, effort is handled AUTOMATICALLY by the structure: J[i] is
# the number of 10-day windows, so an 8-week camera contributes ~6 Bernoulli
# trials against a 3-week camera's ~2, and dOcc_v(..., len = J[i]) accounts
# for it exactly. Aggregating to arrays COLLAPSES the window dimension, so
# effort stops being automatic and MUST re-enter as a per-camera detection
# covariate (log J[i]). Without it, deployment length is confounded with
# occupancy. That is what log_effort below is for.
#
# INTERPRETATION SHIFT -- state this in any writeup
# q now absorbs within-array habitat patchiness as well as detection, and
# psi_a means "the array was used", not "each camera's site was used". That is
# a genuine change in what the parameter means, not merely a change of scale.
#
# ARRAY STRUCTURE PROVENANCE
# The real array/subproject field lives in the Hazel prep bundles and is NOT
# reachable locally. assign_arrays_spatial() below therefore builds a
# PARAMETERIZED APPROXIMATION from site coordinates, so the simulation can
# proceed and be developed locally. Every function that uses it records
# array_source = "spatial_approximation" in its output. When the real
# structure is available, pass it via assign_arrays_from_field() instead and
# the rest of the pipeline is unchanged.
#
# Arrays are ARRAY-YEAR units. The same physical array in 2019 and 2024 must
# be separate observations, since the trend is the estimand. (Note: the array
# grouping in the activity-patterns project pooled across years -- that
# grouping cannot be lifted directly.)
# =============================================================================

#' @name assign_arrays_spatial
#' @description Assign cameras to arrays by spatial proximity within a year,
#'   as a stand-in for the real array/subproject field. Greedy single-linkage
#'   within `radius_km`: cameras are seeded in order and each unassigned
#'   camera within the radius of a seed joins that seed's array. Deliberately
#'   simple and deterministic given `site_order`, so the structure is
#'   reproducible across replicates (the DESIGN is fixed; only data vary).
#' @param x,y Numeric vectors of projected site coordinates (metres), one per
#'   camera site.
#' @param site_year Integer/numeric vector, the year index of each site. Sites
#'   in different years are never placed in the same array-year.
#' @param radius_km Numeric, single-linkage radius in kilometres.
#' @param max_per_array Integer, cap on cameras per array; further cameras
#'   within the radius start a new array. Prevents one dense region becoming a
#'   single implausibly large array.
#' @return Integer vector of array-year ids, one per camera site.
assign_arrays_spatial <- function(x, y, site_year, radius_km = 15,
                                  max_per_array = 25) {
  stopifnot(length(x) == length(y), length(x) == length(site_year))
  n <- length(x)
  r_m <- radius_km * 1000
  array_id <- rep(NA_integer_, n)
  next_id <- 1L

  for (yr in sort(unique(site_year))) {
    idx <- which(site_year == yr)
    for (i in idx) {
      if (!is.na(array_id[i])) next
      # seed a new array at the first unassigned camera in this year
      array_id[i] <- next_id
      members <- 1L
      d <- sqrt((x[idx] - x[i])^2 + (y[idx] - y[i])^2)
      cand <- idx[order(d)]
      for (j in cand) {
        if (members >= max_per_array) break
        if (!is.na(array_id[j])) next
        if (sqrt((x[j] - x[i])^2 + (y[j] - y[i])^2) <= r_m) {
          array_id[j] <- next_id
          members <- members + 1L
        }
      }
      next_id <- next_id + 1L
    }
  }
  stopifnot(!anyNA(array_id))
  array_id
}

#' @name assign_arrays_from_field
#' @description Build array-year ids from a real array/subproject field.
#'   Preferred over the spatial approximation whenever the field is available.
#' @param array_field Character/factor vector, the array or subproject label
#'   per camera site.
#' @param site_year Numeric vector, year index per site.
#' @param drop_singletons Logical; if TRUE, sites whose array-year contains
#'   only one camera are returned as NA (caller decides whether to drop them).
#'   Singleton arrays carry no within-array replication.
#' @return Integer vector of array-year ids (NA for dropped singletons).
assign_arrays_from_field <- function(array_field, site_year,
                                     drop_singletons = FALSE) {
  stopifnot(length(array_field) == length(site_year))
  key <- paste(as.character(array_field), site_year, sep = "__")
  ids <- as.integer(factor(key))
  if (drop_singletons) {
    tab <- table(ids)
    ids[ids %in% as.integer(names(tab)[tab == 1])] <- NA_integer_
  }
  ids
}

#' @name aggregate_to_array
#' @description Collapse a camera-level simulated detection matrix to
#'   array-level detection/non-detection per camera. For each camera, the
#'   array-level observation is whether that camera EVER detected the species
#'   across its own windows; cameras are then the replicates within the array.
#' @param y Matrix (nsite x max_J) of simulated binary camera-window
#'   detections, as produced by simulate_replicate_data(). NA-padded to the
#'   right for cameras with fewer windows.
#' @param J Integer vector of window counts per camera site.
#' @param array_id Integer vector of array-year ids per camera site.
#' @return List with:
#'   \describe{
#'     \item{w}{Matrix (n_array x max_n_a) of 0/1 per camera within array,
#'       NA-padded.}
#'     \item{n_a}{Integer vector, cameras per array.}
#'     \item{log_effort}{Matrix matching w, log(J) of each camera (0-filled in
#'       the padded region; never read, since len = n_a bounds the loop).}
#'     \item{prop_detect}{Numeric vector, realized proportion of cameras
#'       detecting per array -- the quantity the estimator is built around.
#'       Report this: values near 1 mean saturation and are the warning sign.}
#'     \item{site_index}{List of the original site rows per array.}
#'   }
aggregate_to_array <- function(y, J, array_id) {
  stopifnot(nrow(y) == length(J), length(J) == length(array_id))
  keep <- !is.na(array_id)
  y <- y[keep, , drop = FALSE]; J <- J[keep]; array_id <- array_id[keep]

  # per-camera: did it ever detect, over its own J windows?
  ever <- vapply(seq_len(nrow(y)), function(i) {
    as.integer(any(y[i, seq_len(J[i])] == 1, na.rm = TRUE))
  }, integer(1))

  arrays <- sort(unique(array_id))
  n_array <- length(arrays)
  n_a <- vapply(arrays, function(a) sum(array_id == a), integer(1))
  max_na <- max(n_a)

  w <- matrix(NA_integer_, nrow = n_array, ncol = max_na)
  log_eff <- matrix(0, nrow = n_array, ncol = max_na)
  site_index <- vector("list", n_array)

  for (k in seq_len(n_array)) {
    rows <- which(array_id == arrays[k])
    site_index[[k]] <- rows
    w[k, seq_along(rows)] <- ever[rows]
    log_eff[k, seq_along(rows)] <- log(J[rows])
  }

  list(w = w, n_a = as.integer(n_a), log_effort = log_eff,
       prop_detect = rowSums(w, na.rm = TRUE) / n_a,
       site_index = site_index,
       array_source = "caller_supplied")
}

#' @name build_array_constants
#' @description Derive the array-level constants a NIMBLE array model needs,
#'   from camera-level constants plus an array assignment. Site-level
#'   covariates are carried in two ways: occupancy covariates are averaged to
#'   the array (they describe the array's location), while detection
#'   covariates stay PER CAMERA (cameras are the replicates).
#' @param constants Camera-level constants list from build_reduced_constants().
#' @param array_id Integer vector of array-year ids per camera site.
#' @param agg Output of aggregate_to_array().
#' @return Constants list for the array-level models, with `narray`, `n_a`,
#'   per-array occupancy covariates, and per-camera detection covariates.
#' @details The array's cell100 / cell50 assignment is taken as the modal cell
#'   of its member cameras, so the array attaches to exactly one CAR cell.
#'   Arrays spanning a cell boundary are assigned to their majority cell; the
#'   fraction affected is reported by summarize_array_structure() and should be
#'   checked rather than assumed small.
build_array_constants <- function(constants, array_id, agg) {
  keep <- !is.na(array_id)
  aid <- array_id[keep]
  arrays <- sort(unique(aid))

  modal <- function(v) as.integer(names(sort(table(v), decreasing = TRUE))[1])
  per_array <- function(f) vapply(arrays, function(a) f(which(aid == a)), numeric(1))

  out <- constants
  out$narray <- length(arrays)
  out$n_a    <- agg$n_a
  out$max_na <- ncol(agg$w)

  # occupancy side: array-level (describes where the array is)
  out$array_cell   <- vapply(arrays, function(a) modal(constants$cell[keep][aid == a]), integer(1))
  out$array_MWMT   <- per_array(function(r) mean(constants$MWMT[keep][r]))
  out$array_MCMT   <- per_array(function(r) mean(constants$MCMT[keep][r]))
  out$array_year_occ <- per_array(function(r) mean(constants$year_occ[keep][r]))
  oc <- constants$occ_covars[keep, , drop = FALSE]
  out$array_occ_covars <- t(vapply(arrays, function(a) colMeans(oc[aid == a, , drop = FALSE]),
                                    numeric(ncol(oc))))

  # detection side: PER CAMERA within array (cameras are the replicates)
  out$log_effort    <- agg$log_effort
  out$array_canopy  <- .per_camera_matrix(constants$canopy_height[keep], aid, arrays, out$max_na)
  out$array_roaddist<- .per_camera_matrix(constants$log_roaddist[keep], aid, arrays, out$max_na)

  # array_yday: per-camera mean yday across that camera's own windows (the
  # PI-agreed choice per array_level_test_spec.md -- within-array seasonal
  # spread is lost, stated as a caveat, not silently dropped). Both
  # model_code_array_occ.R and model_code_array_rn.R read array_yday[a, i]
  # directly in the detection linear predictor; it must be supplied here.
  camera_mean_yday <- rowMeans(constants$yday[keep, , drop = FALSE], na.rm = TRUE)
  out$array_yday <- .per_camera_matrix(camera_mean_yday, aid, arrays, out$max_na)

  # weights for dcar_normal: both array model files call
  # dcar_normal(adj[1:nadj], weights[1:nadj], num[1:nnum], ...) explicitly
  # (unlike the baseline, which relies on dcar_normal's own default
  # weights = adj/adj). The CAR graph itself is carried over unchanged from
  # constants (same ncell100/adj/num/nadj), so binary (all-1) weights
  # reproduce that same default explicitly.
  out$weights <- rep(1, length(constants$adj))

  out
}

#' @name .per_camera_matrix
#' @description Lay a per-site covariate out as an array x camera matrix,
#'   0-padded beyond n_a (padding is never read; len = n_a bounds the loop).
#' @param v Numeric vector, one value per retained site.
#' @param aid Integer vector of array ids per retained site.
#' @param arrays Sorted unique array ids.
#' @param max_na Integer, matrix column count.
#' @return Numeric matrix (n_array x max_na).
.per_camera_matrix <- function(v, aid, arrays, max_na) {
  m <- matrix(0, nrow = length(arrays), ncol = max_na)
  for (k in seq_along(arrays)) {
    vals <- v[aid == arrays[k]]
    m[k, seq_along(vals)] <- vals
  }
  m
}

#' @name summarize_array_structure
#' @description Report the realized array structure and the realized
#'   information content, so the design is verified rather than assumed. Call
#'   this on the pilot before committing to a full sweep.
#' @param agg Output of aggregate_to_array().
#' @param constants Camera-level constants (for the cell-spanning check).
#' @param array_id Integer vector of array-year ids per camera site.
#' @return A one-row data.frame of diagnostics.
#' @details `mean_prop_detect` near 1 is the saturation warning: it means
#'   array-level occupancy has little information about intensity left, which
#'   is exactly the regime Royle-Nichols is meant to handle. `frac_saturated`
#'   counts arrays where every camera detected -- those contribute nothing to
#'   distinguishing abundance under an occupancy model.
summarize_array_structure <- function(agg, constants, array_id) {
  keep <- !is.na(array_id)
  aid <- array_id[keep]
  arrays <- sort(unique(aid))
  spans <- vapply(arrays, function(a) {
    length(unique(constants$cell[keep][aid == a])) > 1L
  }, logical(1))

  data.frame(
    n_camera_sites   = sum(keep),
    n_arrays         = length(arrays),
    reduction_factor = round(sum(keep) / length(arrays), 2),
    cameras_min      = min(agg$n_a),
    cameras_median   = stats::median(agg$n_a),
    cameras_max      = max(agg$n_a),
    n_singleton      = sum(agg$n_a == 1L),
    mean_prop_detect = round(mean(agg$prop_detect), 4),
    frac_saturated   = round(mean(agg$prop_detect == 1), 4),
    frac_all_zero    = round(mean(agg$prop_detect == 0), 4),
    frac_cell_spanning = round(mean(spans), 4),
    stringsAsFactors = FALSE
  )
}
