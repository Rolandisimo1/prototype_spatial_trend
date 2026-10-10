# =============================================================================
# mdd_lib.R -- resolve scientific names to Mammal Diversity Database species ids (2026-10-08).
#
# The R crosswalk is authoritative; reportlib/mdd_audit.py (local) is the independent audit
# that checks it. Same resolution rule in both, so they can be compared:
#   1. name is an MDD accepted sciName            -> its id          (match_type accepted_exact)
#   2. name is an MDD_original_combination        -> that row's MDD_species -> accepted id (synonym)
#   3. otherwise                                  -> NA, never invented (unresolved)
# Comparison is on ID, never on name strings: two names with one id are a rename to union;
# two names with two ids are two species and are NEVER unioned by lookup (that needs a recorded
# decision in manual_taxonomy_decisions.csv).
# Pinned to MDD_VERSION; every output carries the version.
# =============================================================================
suppressPackageStartupMessages(library(data.table))
MDD_DIR <- Sys.getenv("MDD_DIR", "/home/rwkays/isdm/fleet_v2u/mdd_v2.5")
MDD_VERSION <- "v2.5"

load_mdd <- function(dir = MDD_DIR, version = MDD_VERSION) {
  spf <- list.files(dir, pattern = paste0("^MDD_", version, ".*\\.csv$"), full.names = TRUE)
  synf <- list.files(dir, pattern = paste0("^Species_Syn_Current_", version, "\\.csv$"), full.names = TRUE)
  if (length(spf) != 1 || length(synf) != 1) stop("MDD ", version, " species/synonym CSVs not found (exactly 1 each) in ", dir)
  sp <- fread(spf, select = c("sciName", "id"))
  syn <- fread(synf, select = c("MDD_species", "MDD_original_combination"))
  if (!nrow(sp) || !nrow(syn)) stop("MDD tables empty")
  sp[, name := trimws(gsub("_", " ", sciName))]
  syn[, orig := trimws(gsub("_", " ", MDD_original_combination))]
  syn[, sp_name := trimws(gsub("_", " ", MDD_species))]
  syn <- unique(syn[nzchar(orig) & nzchar(sp_name), .(orig, sp_name)])
  list(acc = setNames(sp$id, sp$name), syn = syn, version = version)
}

mdd_resolve <- function(names, mdd = load_mdd()) {
  u <- unique(trimws(as.character(names[!is.na(names)])))
  rbindlist(lapply(u, function(n) {
    if (n %in% names(mdd$acc))
      return(data.table(verbatim = n, resolved = n, mdd_id = unname(mdd$acc[n]), match_type = "accepted_exact"))
    s <- unique(mdd$syn[orig == n, sp_name]); s <- s[s %in% names(mdd$acc)]
    if (length(s) == 1)
      return(data.table(verbatim = n, resolved = s, mdd_id = unname(mdd$acc[s]), match_type = "synonym"))
    data.table(verbatim = n, resolved = NA_character_, mdd_id = NA_real_,
               match_type = if (length(s) > 1) "ambiguous_synonym" else "unresolved")
  }))[, mdd_version := mdd$version][]
}
