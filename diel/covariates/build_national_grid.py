"""build_national_grid.py -- extract the 4 grid-available manuscript covariates
(pop_1km, ag_5km, tcc_1km, rug_5km) on the 25 km CONUS prediction grid, at the
SAME buffer support used at the arrays (1 km buffer mean for pop_1km/tcc_1km,
5 km buffer mean for ag_5km, 5 km buffer SD of elevation for rug_5km).

tmax_hottest_month is DELIBERATELY NOT extracted here: in the manuscript's own
deployment_covariates_complete_final.parquet it is a per-deployment,
survey-window-matched Daymet value (median within-array SD ~2.7 degC, up to
11.6 degC for arrays surveyed across seasons/years) -- it has no defined value
at an unsampled grid cell and cannot be honestly mapped without substituting a
different climate product (e.g. WorldClim normals), which was out of scope for
this rerun. Models containing tmax_hottest_month are reported in the CV table
but never projected onto the national grid.

Requires an existing 25 km grid of (cell25, lon, lat) points AND the reference
grid_covariates_25km*.csv (any variant) for pop_1km / tcc_1km, which are reused
unchanged (their extraction is already at the correct 1 km-buffer support and
was verified against the array-level values before reuse, r=1.0). Only ag_5km
(NLCD crop+pasture % in a 5 km buffer) and rug_5km (elevation SD in a 5 km
buffer) are extracted fresh here.

Run: python build_national_grid.py <path/to/existing_grid_covariates_25km.csv>
Output: national_grid_covariates_4cov_25km.csv
"""
import ee, json, os, sys, time
import numpy as np, pandas as pd
from google.oauth2 import service_account
from concurrent.futures import ThreadPoolExecutor

KEY = os.environ.get("GEE_SERVICE_ACCOUNT_KEY", "keys/service_account.json")

def init_ee():
    info = json.load(open(KEY))
    creds = service_account.Credentials.from_service_account_info(
        info, scopes=["https://www.googleapis.com/auth/earthengine",
                       "https://www.googleapis.com/auth/cloud-platform"])
    ee.Initialize(creds, project=info["project_id"])
    print("EE initialized for", info.get("client_email"))

NLCD = None; DEM = None

def build_images():
    global NLCD, DEM
    NLCD = ee.Image(ee.ImageCollection("USGS/NLCD_RELEASES/2021_REL/NLCD")
                     .filter(ee.Filter.eq("system:index", "2021")).first()).select("landcover")
    DEM = ee.ImageCollection("COPERNICUS/DEM/GLO30").select("DEM").mosaic()

def fc(sub):
    return ee.FeatureCollection([
        ee.Feature(ee.Geometry.Point([float(r.lon), float(r.lat)]), {"cid": r.cell25})
        for r in sub.itertuples()])

def run(img, sub, scale, reducer, radius):
    f = fc(sub).map(lambda ft: ft.buffer(radius))
    return img.reduceRegions(collection=f, reducer=reducer, scale=scale).getInfo()

def batched(sub, fn, nb=20, label="", workers=10, cache=None):
    out = {}
    if cache and os.path.exists(cache):
        out = json.load(open(cache))
        sub = sub[~sub.cell25.isin(out.keys())]
        print(f"  {label} cache hit {len(out)}, remaining {len(sub)}", flush=True)
        if not len(sub):
            return out
    chunks = [sub.iloc[i:i+nb] for i in range(0, len(sub), nb)]
    def one(ch):
        for attempt in range(5):
            try:
                return fn(ch)
            except Exception as e:
                if attempt == 4:
                    print("FAIL", label, str(e)[:150], flush=True)
                    return {"features": []}
                time.sleep(2*(attempt+1))
        return {"features": []}
    done = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(one, chunks):
            for ft in r["features"]:
                out[ft["properties"]["cid"]] = ft["properties"]
            done += 1
            if done % 40 == 0:
                print(f"  {label} {done*nb}/{len(sub)}", flush=True)
                if cache: json.dump(out, open(cache, "w"))
    if cache: json.dump(out, open(cache, "w"))
    return out

def ag_pct(cid, ag_cache, codes=(81, 82)):
    h = (ag_cache.get(cid) or {}).get("histogram")
    if isinstance(h, str): h = json.loads(h)
    if not h: return np.nan
    tot = sum(h.values())
    return 100.0*sum(float(h.get(str(c), 0)) for c in codes)/tot if tot else np.nan

def main(existing_grid_csv):
    init_ee(); build_images()
    grid = pd.read_csv(existing_grid_csv)
    os.makedirs("gee_cache", exist_ok=True)

    ag_res  = batched(grid[["cell25","lon","lat"]],
                       lambda s: run(NLCD, s, 30, ee.Reducer.frequencyHistogram(), 5000),
                       label="ag5k", cache="gee_cache/ag5k.json")
    rug_res = batched(grid[["cell25","lon","lat"]],
                       lambda s: run(DEM, s, 30, ee.Reducer.stdDev(), 5000),
                       label="rug5k", cache="gee_cache/rug5k.json")

    out = grid[["cell25","lon","lat","x_km","y_km","X_5070_km","Y_5070_km","cell_id",
                "pop_1km","tcc_1km"]].copy()
    out["ag_5km"]  = out.cell25.map(lambda c: ag_pct(c, ag_res))
    out["rug_5km"] = out.cell25.map(lambda c: (rug_res.get(c) or {}).get("stdDev")).astype(float)

    # drop GHSL/JRC no-data sentinel (-200) and coastal buffer cells that blend it in
    out.loc[out.pop_1km == -200, "pop_1km"] = np.nan
    clean = out.dropna(subset=["pop_1km","tcc_1km","ag_5km","rug_5km"])
    clean = clean[clean.pop_1km >= 0]
    clean.to_csv("national_grid_covariates_4cov_25km.csv", index=False)
    print(f"wrote {len(clean)} / {len(grid)} cells ({100*len(clean)/len(grid):.1f}%)")

if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "grid_covariates_25km_seamfixed.csv")
