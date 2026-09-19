#!/usr/bin/env Rscript
# Prints one JSON line {"ncell50":..,"nsite":..} for a fitted bundle.
# Used by submit_chains.py to size chains; cached per bundle in bundle_meta.json.
a <- commandArgs(trailingOnly = TRUE)
b <- readRDS(a[1])
pick <- function(nm) {
  for (s in c("constants_list", "data_list")) if (!is.null(b[[s]][[nm]])) return(b[[s]][[nm]])
  b[[nm]]
}
ncell50 <- pick("ncell50")
nsite   <- pick("nsite")
stopifnot(length(ncell50) == 1, length(nsite) == 1, ncell50 > 0, nsite > 0)
cat(sprintf('{"ncell50": %d, "nsite": %d}\n', as.integer(ncell50), as.integer(nsite)))
