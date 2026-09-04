"""step2b_gridcov.py -- extract the Set A / Set B covariates on the 25 km prediction grid.

SUPPORT MATCHING (the silent-bias risk in this step):
  Camera covariates in deployment_covariates.csv are 1 km BUFFER MEANS around each
  camera point, then averaged over the array's cameras.
  Grid values are extracted at TWO supports:
    *_1kb  : 1 km buffer mean at the 25 km cell CENTROID  <- MATCHED to cameras, used for prediction
    *_cell : 25 km cell-area mean                          <- comparison only, quantifies support effect
Layers/bands are identical to the camera pass (verified with check_asset).
"""
import ee, json, math, sys, time
import pandas as pd, numpy as np
from google.oauth2 import service_account

KEY = os.environ.get("GEE_SERVICE_ACCOUNT_KEY", "keys/service_account.json")
_ki = json.load(open(KEY))
ee.Initialize(service_account.Credentials.from_service_account_file(
    KEY, scopes=["https://www.googleapis.com/auth/earthengine",
                 "https://www.googleapis.com/auth/cloud-platform"]),
    project=_ki["project_id"])

grid = pd.read_csv("grid25_l48.csv")
print("cells", len(grid), flush=True)

NLCD  = ee.Image(ee.ImageCollection("USGS/NLCD_RELEASES/2021_REL/NLCD")
                 .filter(ee.Filter.eq("system:index", "2021")).first()).select("landcover")
# TCC is tiled by REGION and YEAR (e.g. TCC_v2021-4_CONUS_2021); a plain
# filterDate().first() returns a null image, so mosaic the CONUS 2021 tiles.
TCC   = (ee.ImageCollection("USGS/NLCD_RELEASES/2021_REL/TCC/v2021-4")
         .filter(ee.Filter.stringContains("system:index", "CONUS_2021"))
         .select("Science_Percent_Tree_Canopy_Cover").mosaic().unmask(0))
NTL   = ee.ImageCollection("NOAA/VIIRS/DNB/ANNUAL_V22").filterDate(
            "2023-01-01", "2024-01-01").first().select("average")
POP   = ee.ImageCollection("JRC/GHSL/P2023A/GHS_POP").filterDate(
            "2020-01-01", "2021-01-01").first().select("population_count")
DEM   = ee.ImageCollection("COPERNICUS/DEM/GLO30").select("DEM").mosaic()
WCm   = ee.ImageCollection("WORLDCLIM/V1/MONTHLY")
TAVG  = WCm.select("tavg").toBands()                       # 12 bands, 0.1 degC
# Set B winter severity: MODIS 8-day maximum snow extent, fraction of Nov-Mar
# composites mapped as snow (2019-2023). Chosen over Daymet SWE purely on
# server cost (5x faster for an equivalent axis); it is snow COVER DURATION,
# which is what the winter-severity axis was asked to capture.
SWE   = (ee.ImageCollection("MODIS/061/MOD10A2").select("Maximum_Snow_Extent")
         .filterDate("2019-11-01", "2023-04-01")
         .filter(ee.Filter.calendarRange(11, 3, "month"))
         .map(lambda i: i.eq(200)).mean().rename("snowfrac"))
# 'days below freezing' is NOT extracted: WorldClim bio06 correlates r=0.984 with
# t_coldmonth (already in Set A), so it would add a dimension to the
# extrapolation envelope while carrying no independent information.
FRZ   = None
# Set B: distance to permanent water (JRC GSW occurrence >= 50%), metres
WAT   = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence").gte(50).unmask(0)
DIST  = WAT.fastDistanceTransform(256).sqrt().multiply(
            ee.Image.pixelArea().sqrt()).rename("dist_water")

def fc(sub):
    return ee.FeatureCollection([
        ee.Feature(ee.Geometry.Point([float(r.lon), float(r.lat)]), {"cid": r.cell25})
        for r in sub.itertuples()])

def run(img, sub, scale, reducer, radius=None, tile=None):
    f = fc(sub)
    if radius:
        f = f.map(lambda ft: ft.buffer(radius))
    elif tile:
        f = f.map(lambda ft: ft.buffer(tile, 1).bounds())
    return img.reduceRegions(collection=f, reducer=reducer, scale=scale).getInfo()

from concurrent.futures import ThreadPoolExecutor
import os

def batched(sub, fn, nb=25, label="", workers=12, cache=None):
    """Parallel batched reduceRegions with an on-disk cache per layer, so an
    interrupted run resumes instead of restarting."""
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
                    print("FAIL", label, str(e)[:100], flush=True)
                    return {"features": []}
                time.sleep(2 * (attempt + 1))

    done = 0
    with ThreadPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(one, chunks):
            for ft in r["features"]:
                out[ft["properties"]["cid"]] = ft["properties"]
            done += 1
            if done % 40 == 0:
                print(f"  {label} {done*nb}/{len(sub)}", flush=True)
                if cache:
                    json.dump(out, open(cache, "w"))
    if cache:
        json.dump(out, open(cache, "w"))
    return out

MEAN = ee.Reducer.mean()
HIST = ee.Reducer.frequencyHistogram()
SD   = ee.Reducer.stdDev()

def pct_from_hist(props, codes):
    h = props.get("histogram") or props.get("landcover") or {}
    if isinstance(h, str):
        h = json.loads(h)
    tot = sum(h.values()) if h else 0
    if not tot:
        return np.nan
    return 100.0 * sum(float(h.get(str(c), 0)) for c in codes) / tot

if __name__ == "__main__":
    which = sys.argv[1]                     # "1kb" or "cell"
    R, TILE = (1000, None) if which == "1kb" else (None, 12500)
    kw = dict(radius=R, tile=TILE)
    g = grid
    res = {}
    jobs = [("ntl", NTL, 500, MEAN), ("pop", POP, 100, MEAN), ("tcc", TCC, 30, MEAN),
            ("nlcd", NLCD, 30, HIST), ("rug", DEM, 30, SD), ("elev", DEM, 30, MEAN),
            ("tavg", TAVG, 1000, MEAN), ("swe", SWE, 500, MEAN), ("dw", DIST, 300, MEAN)]
    for lab, img, sc, red in jobs:
        res[lab] = batched(g, lambda s, i=img, c=sc, r=red: run(i, s, c, r, **kw),
                           label=lab, cache=f"cache_{which}_{lab}.json")
    with open(f"grid_ee_{which}.json", "w") as fh:
        json.dump(res, fh)
    print("DONE", which, {k: len(v) for k, v in res.items()}, flush=True)
