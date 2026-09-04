"""step7_gridassemble.py -- assemble the EE layer caches into grid_covariates_25km.csv.

Columns mirror deployment_covariates.csv names EXACTLY so the model can be applied
without renaming. Support is documented per column in grid_covariate_support.csv.
"""
import json, os, sys
import numpy as np, pandas as pd

WHICH = sys.argv[1] if len(sys.argv) > 1 else "1kb"
grid = pd.read_csv("grid25_full.csv")


def load(lab):
    d = {}
    for p in (f"cache_{WHICH}_{lab}.json", f"cache_new_{WHICH}_{lab}.json"):
        if os.path.exists(p):
            d.update(json.load(open(p)))
    return d


def col(lab, key="mean"):
    d = load(lab)
    return grid.cell25.map(lambda c: (d.get(c) or {}).get(key))


NLCD_GROUPS = {"nlcd_1k_crop": [82], "nlcd_1k_pasture": [81],
               "nlcd_1k_wetland": [90, 95], "nlcd_1k_shrub_herb": [52, 71],
               "nlcd_1k_forest": [41, 42, 43]}

out = pd.DataFrame({"cell25": grid.cell25, "lon": grid.lon, "lat": grid.lat,
                    "x_km": grid.x_km, "y_km": grid.y_km,
                    "X_5070_km": grid.X_5070_km, "Y_5070_km": grid.Y_5070_km,
                    "cell_id": grid.cell_id})
out["ntl_1km"] = pd.to_numeric(col("ntl"), errors="coerce")
out["tcc_1km"] = pd.to_numeric(col("tcc"), errors="coerce")
out["rug_1km"] = pd.to_numeric(col("rug", "stdDev"), errors="coerce")
out["elev"]    = pd.to_numeric(col("elev"), errors="coerce")
out["dist_water_m"] = pd.to_numeric(col("dw"), errors="coerce")
out["snowfrac"] = pd.to_numeric(col("swe"), errors="coerce")

# GHSL population: camera column pop_1km is the 1 km buffer MEAN of per-pixel
# population COUNT (100 m pixels). Same reducer/scale here, so directly comparable.
out["pop_1km"] = pd.to_numeric(col("pop"), errors="coerce")

# NLCD composition -> class percentages from the frequency histogram
nl = load("nlcd")
def nlcd_pct(cid, codes):
    h = (nl.get(cid) or {}).get("histogram")
    if isinstance(h, str):
        h = json.loads(h)
    if not h:
        return np.nan
    tot = sum(h.values())
    return 100.0 * sum(float(h.get(str(c), 0)) for c in codes) / tot if tot else np.nan
for name, codes in NLCD_GROUPS.items():
    out[name] = [nlcd_pct(c, codes) for c in grid.cell25]

# WorldClim monthly tavg (0.1 degC) -> warmest / coldest month means in degC,
# matching the camera columns t_warmmonth / t_coldmonth.
tv = load("tavg")
bands = [f"{m:02d}_tavg" for m in range(1, 13)]
def tstats(cid):
    p = tv.get(cid) or {}
    v = [p.get(b) for b in bands]
    v = [x for x in v if x is not None]
    if not v:
        return (np.nan, np.nan)
    v = np.array(v, float) * 0.1
    return (float(v.max()), float(v.min()))
tw, tc = zip(*[tstats(c) for c in grid.cell25])
out["t_warmmonth"], out["t_coldmonth"] = tw, tc

# Day length: mean over the Aug-Oct analysis window (the aso_window the counts use),
# computed from latitude -- identical formula to the camera-side value.
def daylen(lat, doy):
    th = 0.2163108 + 2 * np.arctan(0.9671396 * np.tan(0.00860 * (doy - 186)))
    ph = np.arcsin(0.39795 * np.cos(th))
    x = (np.sin(0.8333 / 180 * np.pi) + np.sin(lat / 180 * np.pi) * np.sin(ph)) / (
         np.cos(lat / 180 * np.pi) * np.cos(ph))
    return 24 - (24 / np.pi) * np.arccos(np.clip(x, -1, 1))
out["daylength_h"] = np.mean([daylen(out.lat.values, d) for d in range(213, 305)], axis=0)

out.to_csv(f"grid_covariates_25km_{WHICH}.csv", index=False)
cov = [c for c in out.columns if c not in
       ("cell25", "lon", "lat", "x_km", "y_km", "X_5070_km", "Y_5070_km", "cell_id")]
print(WHICH, "cells", len(out))
print(out[cov].notna().mean().round(3).to_dict())
