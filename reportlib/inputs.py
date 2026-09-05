"""Inputs resolved by filename, not by hand-staged path.

The sandbox workspace is wiped between sessions; the repo and the artifact
store are not. Everything here resolves against the repo directory or the
artifact store so the pipeline runs from a clean checkout.
"""
import os
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PULL = os.path.join(os.path.dirname(REPO), "hazel_pull_20260827", "posteriors_20260828")

# Simulation sweeps. Tagged pairs: each generation has its own rows+summary.
# n30_abund  = ladder-fixed 4-arm sweep, the bias/precision/coverage source
# n30_nullfix = corrected-null sweep, the false-positive source
# n90_rnnull = RN-only replicate bump at n=90/cell, which cleared the RN arms
SIM = {
    "abund_summary": "estimator_sweep_n30_abund_summary.csv",
    "nullfix_rows":  "estimator_sweep_n30_nullfix_rows.csv",
    "rnnull_rows":   "estimator_sweep_n90_rnnull_rows.csv",
}

# Real fits. All eight are single continuous chains (no resume boundary), so
# none can carry the chunked-checkpoint burn-in defect.
POSTERIOR = {
    "bobcat_ns":  (PULL, "bobcat_v2b_national_scalar_posterior.csv"),
    "bobcat_eco": (PULL, "bobcat_v2b_ecoregion_global.csv"),
    "bobcat_reg": (PULL, "bobcat_v2b_ecoregion_ecoregion_posterior.csv"),
    "wtd_ns":     (PULL, "white-tailed_deer_v2b_national_scalar_posterior.csv"),
    "wtd_eco":    (PULL, "white-tailed_deer_v2b_ecoregion_global.csv"),
    "wtd_reg":    (PULL, "white-tailed_deer_v2b_ecoregion_ecoregion_posterior.csv"),
    "moose_ns":   (REPO, "moose_v2b_national_scalar_posterior_single.csv"),
    "moose_v1fix9_ns": (PULL, "moose_v1fix9_national_scalar_posterior.csv"),
    "occbeta":    (PULL, "occbeta_posterior_all_models.csv"),
    "theta":      (PULL, "theta_posterior_all_models.csv"),
}

MU = {
    "Bobcat":            "bobcat_v2b_national_scalar_mu_snapshots_single.csv",
    "White-tailed deer": "white-tailed_deer_v2b_national_scalar_mu_snapshots_single.csv",
    "Moose":             "moose_v2b_national_scalar_mu_snapshots_single.csv",
}

AGENCY_MASTER = "deer_moose_trends_master.csv"   # WTD + moose only; no bobcat data exists
ECOREGIONS = "/Users/rwkays/claude_code/geospatialdata/epa_na_ecoregions_level1/NA_CEC_Eco_Level1.shp"
IUCN_RANGES = "/Users/rwkays/claude_code/geospatialdata/mammalranges/data_0.shp"

# Snapshot years: -10, -5, current. Inside the fitted 2008-2025 span of the
# 18-year fits. NOTE both fall outside a 5-year (2021-2025) window and 2015
# falls outside the 10-year (2016-2025) window -- see conventions.verify_snapshot_years.
SNAPSHOT_YEARS = [2015, 2020, 2025]


def sim(key):
    return pd.read_csv(os.path.join(REPO, SIM[key]))


def posterior(key):
    d, f = POSTERIOR[key]
    return pd.read_csv(os.path.join(d, f))


def mu(species):
    return pd.read_csv(os.path.join(REPO, MU[species]))


def states_conus(cache="states"):
    """Census state boundaries, lower 48 only, Albers equal-area (EPSG:5070).

    Downloads on a clean checkout. AK/HI/territories dropped: every map in this
    report is lower-48, and leaving them in silently expands the map extent so
    the CONUS data renders as a small cluster.
    """
    import geopandas as gpd, urllib.request, zipfile
    shp = os.path.join(cache, "cb_2022_us_state_20m.shp")
    if not os.path.exists(shp):
        os.makedirs(cache, exist_ok=True)
        z = os.path.join(cache, "cb.zip")
        urllib.request.urlretrieve(
            "https://www2.census.gov/geo/tiger/GENZ2022/shp/cb_2022_us_state_20m.zip", z)
        with zipfile.ZipFile(z) as zf:
            zf.extractall(cache)
    g = gpd.read_file(shp).to_crs(5070)
    return g[~g["STUSPS"].isin({"HI", "PR", "VI", "GU", "MP", "AS", "AK"})].copy()
