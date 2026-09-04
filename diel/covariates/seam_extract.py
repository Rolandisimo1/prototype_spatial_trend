"""seam_extract.py -- reproduce the ORIGINAL step2b_gridcov.py 1 km-buffer support exactly.

Support (recovered from step2b_gridcov.py + step7_gridassemble.py):
  geometry = ee.Geometry.Point([lon,lat]).buffer(1000)   # 1 km buffer at cell CENTROID
  reducers/scales per layer as in JOBS below.
Nothing is unmasked or cleaned post-hoc (GHSL -200 water nodata is carried through
as the original did), so new values are drawn from the same distribution as stored ones.
"""
import ee, json, os, time
import numpy as np, pandas as pd
from google.oauth2 import service_account
from concurrent.futures import ThreadPoolExecutor

KEY = "/Users/rwkays/claude_code/keys/snapshotusa-fd9e954c98b9.json"


def init():
    ki = json.load(open(KEY))
    ee.Initialize(service_account.Credentials.from_service_account_file(
        KEY, scopes=["https://www.googleapis.com/auth/earthengine",
                     "https://www.googleapis.com/auth/cloud-platform"]),
        project=ki["project_id"])
    return ki["project_id"]


def layers():
    NLCD = ee.Image(ee.ImageCollection("USGS/NLCD_RELEASES/2021_REL/NLCD")
                    .filter(ee.Filter.eq("system:index", "2021")).first()).select("landcover")
    TCC = (ee.ImageCollection("USGS/NLCD_RELEASES/2021_REL/TCC/v2021-4")
           .filter(ee.Filter.stringContains("system:index", "CONUS_2021"))
           .select("Science_Percent_Tree_Canopy_Cover").mosaic().unmask(0))
    NTL = ee.ImageCollection("NOAA/VIIRS/DNB/ANNUAL_V22").filterDate(
        "2023-01-01", "2024-01-01").first().select("average")
    POP = ee.ImageCollection("JRC/GHSL/P2023A/GHS_POP").filterDate(
        "2020-01-01", "2021-01-01").first().select("population_count")
    DEM = ee.ImageCollection("COPERNICUS/DEM/GLO30").select("DEM").mosaic()
    TAVG = ee.ImageCollection("WORLDCLIM/V1/MONTHLY").select("tavg").toBands()
    SWE = (ee.ImageCollection("MODIS/061/MOD10A2").select("Maximum_Snow_Extent")
           .filterDate("2019-11-01", "2023-04-01")
           .filter(ee.Filter.calendarRange(11, 3, "month"))
           .map(lambda i: i.eq(200)).mean().rename("snowfrac"))
    WAT = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").select("occurrence").gte(50).unmask(0)
    DIST = WAT.fastDistanceTransform(256).sqrt().multiply(
        ee.Image.pixelArea().sqrt()).rename("dist_water")
    MEAN, HIST, SD = ee.Reducer.mean(), ee.Reducer.frequencyHistogram(), ee.Reducer.stdDev()
    return [("ntl", NTL, 500, MEAN), ("pop", POP, 100, MEAN), ("tcc", TCC, 30, MEAN),
            ("nlcd", NLCD, 30, HIST), ("rug", DEM, 30, SD), ("elev", DEM, 30, MEAN),
            ("tavg", TAVG, 1000, MEAN), ("swe", SWE, 500, MEAN), ("dw", DIST, 300, MEAN)]


def _fc(sub):
    return ee.FeatureCollection([
        ee.Feature(ee.Geometry.Point([float(r.lon), float(r.lat)]).buffer(1000),
                   {"cid": r.cell25}) for r in sub.itertuples()])


def _run(img, sub, scale, reducer):
    return img.reduceRegions(collection=_fc(sub), reducer=reducer, scale=scale).getInfo()


def batched(sub, img, scale, reducer, nb=20, label="", workers=4, cache=None):
    out = {}
    if cache and os.path.exists(cache):
        out = json.load(open(cache))
        sub = sub[~sub.cell25.isin(out.keys())]
        if not len(sub):
            return out
    chunks = [sub.iloc[i:i + nb] for i in range(0, len(sub), nb)]

    def one(ch):
        for attempt in range(6):
            try:
                return _run(img, ch, scale, reducer)
            except Exception as e:
                if attempt == 5:
                    print("FAIL", label, str(e)[:120], flush=True)
                    return {"features": []}
                time.sleep(3 * (attempt + 1))   # graceful backoff on throttling

    with ThreadPoolExecutor(max_workers=workers) as ex:
        for r in ex.map(one, chunks):
            for ft in r["features"]:
                out[ft["properties"]["cid"]] = ft["properties"]
    if cache:
        json.dump(out, open(cache, "w"))
    return out


NLCD_GROUPS = {"nlcd_1k_crop": [82], "nlcd_1k_pasture": [81],
               "nlcd_1k_wetland": [90, 95], "nlcd_1k_shrub_herb": [52, 71],
               "nlcd_1k_forest": [41, 42, 43]}
BANDS = [f"{m:02d}_tavg" for m in range(1, 13)]


def daylen(lat, doy):
    th = 0.2163108 + 2 * np.arctan(0.9671396 * np.tan(0.00860 * (doy - 186)))
    ph = np.arcsin(0.39795 * np.cos(th))
    x = (np.sin(0.8333 / 180 * np.pi) + np.sin(lat / 180 * np.pi) * np.sin(ph)) / (
        np.cos(lat / 180 * np.pi) * np.cos(ph))
    return 24 - (24 / np.pi) * np.arccos(np.clip(x, -1, 1))


def assemble(df, res):
    def col(lab, key="mean"):
        d = res[lab]
        return pd.to_numeric(df.cell25.map(lambda c: (d.get(c) or {}).get(key)), errors="coerce")

    out = pd.DataFrame({"cell25": df.cell25.values})
    out["ntl_1km"] = col("ntl").values
    out["tcc_1km"] = col("tcc").values
    out["rug_1km"] = col("rug", "stdDev").values
    out["elev"] = col("elev").values
    out["dist_water_m"] = col("dw").values
    out["snowfrac"] = col("swe").values
    out["pop_1km"] = col("pop").values

    nl = res["nlcd"]

    def nlcd_pct(cid, codes):
        h = (nl.get(cid) or {}).get("histogram")
        if isinstance(h, str):
            h = json.loads(h)
        if not h:
            return np.nan
        tot = sum(h.values())
        return 100.0 * sum(float(h.get(str(c), 0)) for c in codes) / tot if tot else np.nan
    for name, codes in NLCD_GROUPS.items():
        out[name] = [nlcd_pct(c, codes) for c in df.cell25]

    tv = res["tavg"]

    def tstats(cid):
        p = tv.get(cid) or {}
        v = [x for x in (p.get(b) for b in BANDS) if x is not None]
        if not v:
            return (np.nan, np.nan)
        v = np.array(v, float) * 0.1
        return (float(v.max()), float(v.min()))
    tw, tc = zip(*[tstats(c) for c in df.cell25])
    out["t_warmmonth"], out["t_coldmonth"] = tw, tc
    out["daylength_h"] = np.mean([daylen(df.lat.values, d) for d in range(213, 305)], axis=0)
    return out


def extract_all(df, tag, workers=4):
    res = {}
    for lab, img, sc, red in layers():
        res[lab] = batched(df, img, sc, red, label=f"{tag}:{lab}",
                           workers=workers, cache=f"cache_{tag}_{lab}.json")
        print(f"  {tag}:{lab} {len(res[lab])}/{len(df)}", flush=True)
    return assemble(df, res)
