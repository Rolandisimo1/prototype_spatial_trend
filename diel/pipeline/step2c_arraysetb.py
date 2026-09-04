"""step2c_arraysetb.py -- Set B raster covariates (SWE, freeze days, distance to water)
at CAMERA support: 1 km buffer mean at each camera point, then averaged to the array.
Identical layers/reducers/support to the grid pass, so camera and grid values are comparable.
"""
import ee, json, time
import pandas as pd, numpy as np

src = open("step2b_gridcov.py").read().split('if __name__ == "__main__":')[0]
src = src.replace('grid = pd.read_csv("grid25_l48.csv")', 'grid = None')
src = src.replace('print("cells", len(grid), flush=True)', '')
ns = {}
exec(src, ns)

cov = pd.read_csv(COVP := __import__("os").environ["COVP"])
cov["dep_key"] = cov.project_id.astype(str) + "|" + cov.deployment_id.astype(str)
amap = pd.read_csv("dep_array_map.csv")            # dep_key -> array_id
cov = cov.merge(amap, on="dep_key", how="inner")
pts = cov[["dep_key", "array_id", "latitude", "longitude"]].copy()
pts["key"] = pts.latitude.round(5).astype(str) + "_" + pts.longitude.round(5).astype(str)

# CAP cameras per array. Within-array 1 km buffers are near-copies (median camera
# spacing 51 m, 94% within 1 km), so an array mean over <=6 sampled cameras is an
# accurate stand-in for the mean over all of them. Cuts EE calls ~4x.
CAP = 6
_uc = pts.drop_duplicates("key")
u = (_uc.sample(frac=1.0, random_state=7)
        .groupby("array_id", sort=False).head(CAP)
        .reset_index(drop=True))
u = u.rename(columns={"latitude": "lat", "longitude": "lon"})
u["cell25"] = u["key"]
print("camera coords after cap:", len(u), "arrays", u.array_id.nunique(), flush=True)

MEAN = ee.Reducer.mean()
out = {}
for lab, img, scale in [("swe", ns["SWE"], 500), ("dw", ns["DIST"], 300)]:
    out[lab] = ns["batched"](u, lambda s, i=img, sc=scale: ns["run"](i, s, sc, MEAN, radius=1000),
                             nb=25, label=lab, cache=f"cache_cam_{lab}.json")
    print(lab, "done", len(out[lab]), flush=True)

rows = []
for k in u["key"]:
    rows.append({"key": k,
                 "snowfrac": (out["swe"].get(k) or {}).get("mean"),
                 "dist_water_m": (out["dw"].get(k) or {}).get("mean")})
E = pd.DataFrame(rows).merge(u[["key", "array_id"]], on="key", how="left")
A = E.groupby("array_id")[["snowfrac", "dist_water_m"]].mean().reset_index()
A.to_csv("camera_setb_array.csv", index=False)
print("saved camera_setb_array.csv", len(A),
      A[["snowfrac", "dist_water_m"]].notna().mean().round(3).to_dict(), flush=True)
