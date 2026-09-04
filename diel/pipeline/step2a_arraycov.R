## step2a_arraycov.R -- aggregate camera covariates to ARRAY level (the unit of analysis)
## and compute the two camera-derived Set B axes (heterospecific rate, day length).
suppressPackageStartupMessages({library(data.table); library(arrow)})
SUNP <- Sys.getenv("SUNP"); COVP <- Sys.getenv("COVP")

H  <- fread("array_harmonics_raw.csv")
cv <- fread(COVP)
cv[, dep_key := paste0(project_id, "|", deployment_id)]

## map deployments -> arrays using the diel table (the array grouping lives there)
d <- as.data.table(read_parquet(SUNP))
d <- d[aso_window == TRUE & effort_h_aso > 0]
map <- unique(d[, .(dep_key, array_id)])
stopifnot(!anyDuplicated(map$dep_key))
cv <- merge(cv, map, by = "dep_key")

## SET A covariates: 1 km buffer means at each camera -> ARRAY MEAN.
## tmax_window is DROPPED (survey-window quantity, undefined at an unsampled cell);
## replaced by WorldClim normals t_warmmonth / t_coldmonth, evaluable everywhere.
SETA <- c("ntl_1km","pop_1km","nlcd_1k_crop","nlcd_1k_pasture","tcc_1km",
          "t_warmmonth","t_coldmonth","rug_1km")
stopifnot(all(SETA %in% names(cv)))
A <- cv[, c(lapply(.SD, mean, na.rm = TRUE), .(n_dep_cov = .N)),
        by = array_id, .SDcols = SETA]

## ---- Set B, camera-derived axis 1: DAY LENGTH over the survey window ----
## CBM model; sun-time anchoring fixes sunrise/sunset positions, so day LENGTH
## is not otherwise in the model. Uses the array's own active dates.
daylen <- function(lat, doy) {
  th <- 0.2163108 + 2*atan(0.9671396*tan(0.00860*(doy-186)))
  ph <- asin(0.39795*cos(th))
  x  <- (sin(0.8333/180*pi) + sin(lat/180*pi)*sin(ph)) / (cos(lat/180*pi)*cos(ph))
  24 - (24/pi)*acos(pmin(pmax(x, -1), 1))
}
dep_dates <- cv[, .(dep_key, array_id, latitude, start_date, end_date)]
dep_dates[, `:=`(sd_ = as.IDate(start_date), ed_ = as.IDate(end_date))]
dep_dates[, mid := as.IDate((as.numeric(sd_) + as.numeric(ed_))/2, origin = "1970-01-01")]
dep_dates[, doy := as.integer(format(mid, "%j"))]
dep_dates[, dl := daylen(latitude, doy)]
DL <- dep_dates[is.finite(dl), .(daylength_h = mean(dl)), by = array_id]

## ---- Set B, camera-derived axis 2: HETEROSPECIFIC DETECTION RATE ----
## effort-corrected detections/100 h at the SAME array, log1p scale.
## CAVEAT (carried into the report): conflates abundance with detectability.
## Never a species' own rate as its own covariate; NOT evaluable at an unsampled
## cell, so it can enter the CV signal test but NOT the continental map.
rate <- d[, .(ev = sum(count), eff = sum(effort_h_aso)), by = .(species, array_id)]
rate <- rate[eff > 0][, r100 := 100*ev/eff]
wide <- dcast(rate, array_id ~ species, value.var = "r100", fill = 0)
setnames(wide, old = setdiff(names(wide), "array_id"),
         new = paste0("rate_", gsub("[^A-Za-z]", "", setdiff(names(wide), "array_id"))))
for (cc in setdiff(names(wide), "array_id")) set(wide, j = cc, value = log1p(wide[[cc]]))

ARR <- Reduce(function(a,b) merge(a,b,by="array_id",all.x=TRUE), list(A, DL, wide))
fwrite(ARR, "array_covariates.csv")
cat(sprintf("arrays with covariates: %d\n", nrow(ARR)))
cat("Set A cols:", paste(SETA, collapse=", "), "\n")
cat("hetero cols:", paste(setdiff(names(wide),"array_id"), collapse=", "), "\n")
print(ARR[, lapply(.SD, function(z) round(c(mean(z,na.rm=TRUE)),3)), .SDcols = c(SETA,"daylength_h")])
