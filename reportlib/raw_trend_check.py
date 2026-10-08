"""Compare a fitted camera trend against the raw camera detection rate.

READ THIS BEFORE USING IT AS EVIDENCE. Tested against the 18 established fits,
this screen flags 6 of them -- including white-tailed deer, the control, and
black bear. **A sign disagreement between the naive rate and a fitted trend is
therefore NOT evidence the fit is wrong.** The naive rate ignores effort, the
range mask and covariates, and it is dominated by camera-network growth into a
species' range: moose comes back at +55 %/yr naive and mule deer at +82 %/yr,
which are impossible as population trends.

What the screen is good for, having been calibrated that way:

  * A large SAME-SIGN magnitude gap, where network growth cannot be the
    explanation. Collared peccary is the case that motivated this: fitted
    -54.3 %/yr against a naive -1.7 (p=0.80), a 32-fold gap with both streams
    agreeing the species is roughly flat.
  * The sampling footprint it reports -- `n_detecting_years` above all. That is
    a property of the data rather than a comparison, and it is the column that
    actually separates the suspect Tier B fits (California ground squirrel 6
    detecting years, western gray squirrel 7) from the established ones
    (15-16). A trend estimated over 6 years and printed in the same %/yr column
    as one estimated over 16 is not a peer.

Use `contradicts_sign` as a prompt to look, never as a verdict.
"""
import os

import numpy as np
import pandas as pd
from scipy.stats import linregress

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RAW = os.path.abspath(os.path.join(HERE, "..", "raw_cam_data"))
MIN_DEPLOY_PER_YEAR = 20


def _load():
    dep = pd.read_csv(os.path.join(RAW, "combined_deployments_all.csv"),
                      usecols=["deployment_id", "project_id", "year"],
                      dtype={"project_id": str}
                      ).dropna(subset=["year"]).drop_duplicates("deployment_id")
    seq = pd.read_csv(os.path.join(RAW, "combined_sequences_all.csv"),
                      usecols=["deployment_id", "class", "common_name"], low_memory=True)
    seq = seq[seq["class"] == "Mammalia"].copy()
    seq["sp"] = (seq.common_name.str.strip().str.lower()
                 .str.replace(r"\s+", "_", regex=True).str.replace("'", "", regex=False))
    return dep, seq


def raw_trend(species, dep=None, seq=None):
    """Naive detection-rate trend for one species, with its sampling footprint.

    Restricted to projects that ever recorded the species: the national
    deployment pool spans regions the species cannot occupy, and those zeros
    would manufacture a trend out of network growth alone.
    """
    if dep is None or seq is None:
        dep, seq = _load()
    det = set(seq.loc[seq.sp == species, "deployment_id"])
    if not det:
        return None
    proj = set(dep.loc[dep.deployment_id.isin(det), "project_id"])
    sub = dep[dep.project_id.isin(proj)].copy()
    sub["hit"] = sub.deployment_id.isin(det)
    g = sub.groupby("year").agg(deployments=("hit", "size"), detecting=("hit", "sum"))
    g["rate"] = g.detecting / g.deployments
    g = g[g.deployments >= MIN_DEPLOY_PER_YEAR]
    if len(g) < 4:
        return None
    lr = linregress(g.index.astype(float), np.log(g.rate.clip(lower=1e-4)))
    dd = sub[sub.hit]
    return {"species": species,
            "raw_pct_per_yr": round(100 * (np.exp(lr.slope) - 1), 1),
            "raw_p": round(lr.pvalue, 4),
            "years_used": len(g),
            "detecting_year_span": f"{dd.year.min():.0f}-{dd.year.max():.0f}",
            "n_detecting_years": int(dd.year.nunique()),
            "detecting_deployments": int(len(dd)),
            "projects": int(dd.project_id.nunique()),
            "top_project_share": round(dd.project_id.value_counts().iloc[0] / len(dd), 2)}


def screen(fitted, species_col="species", fitted_col="camera_pct_m"):
    """Flag fitted trends the raw record contradicts.

    `contradicts_sign` is the strong signal: raw and fitted disagree in
    direction. `magnitude_ratio` above ~5 with a flat raw trend is the other.
    """
    dep, seq = _load()
    out = []
    for sp, fit in zip(fitted[species_col], fitted[fitted_col]):
        r = raw_trend(sp, dep, seq)
        if r is None:
            out.append({"species": sp, "fitted_pct_per_yr": fit, "raw_pct_per_yr": np.nan,
                        "verdict": "no raw comparison possible"})
            continue
        r["fitted_pct_per_yr"] = fit
        opp = np.sign(fit) != np.sign(r["raw_pct_per_yr"]) and abs(fit) > 2
        ratio = abs(fit) / max(abs(r["raw_pct_per_yr"]), 0.5)
        r["contradicts_sign"] = bool(opp)
        r["magnitude_ratio"] = round(ratio, 1)
        # Deliberately weak wording: "sign differs" fires on 6 of 18 established
        # fits, so it is a prompt to look, not a finding. The same-sign
        # magnitude gap is the one that isolated peccary.
        r["verdict"] = ("magnitude gap, same sign - INVESTIGATE"
                        if (not opp and ratio > 5 and abs(fit) > 5)
                        else "sign differs - weak, see module docstring" if opp
                        else "consistent")
        out.append(r)
    cols = ["species", "fitted_pct_per_yr", "raw_pct_per_yr", "raw_p", "verdict",
            "contradicts_sign", "magnitude_ratio", "detecting_year_span",
            "n_detecting_years", "detecting_deployments", "projects", "top_project_share"]
    d = pd.DataFrame(out)
    return d[[c for c in cols if c in d.columns]]


def spatial_signal(species_list, cell_m=50_000, min_deploy_per_cell=10):
    """Between-cell variance in camera detection rate, relative to binomial noise.

    A ratio near 1 means the species' camera detections are spatially
    indistinguishable from coin flips at a constant rate -- there is no spatial
    surface for the camera stream to anchor. This matters for theta1, which is
    identified by the correspondence between the two streams' SPATIAL patterns:
    if one stream sees no spatial variation, theta1 is weakly identified and
    integration cannot rescue the species, because the link between the streams
    is the thing being estimated.

    Measured across the 25-species fleet the ratio runs 5.4 to 23.1 (median
    ~10). It is a contributor, not a predictor: American badger is second
    lowest at 5.5 and its chains are bimodal in theta1, but puma is lowest at
    5.4 and converged cleanly, and porcupine is middling at 8.6 and did not.
    Use it to describe a species' camera surface, not to forecast a fit.
    """
    import geopandas as gpd
    dep, seq = _load()
    d2 = pd.read_csv(os.path.join(RAW, "combined_deployments_all.csv"),
                     usecols=["deployment_id", "project_id", "latitude", "longitude"],
                     dtype={"project_id": str}
                     ).dropna(subset=["latitude", "longitude"]).drop_duplicates("deployment_id")
    g = gpd.GeoDataFrame(d2, geometry=gpd.points_from_xy(d2.longitude, d2.latitude),
                         crs=4326).to_crs(5070)
    d2["cell"] = ((g.geometry.y // cell_m).astype(int).astype(str) + "_" +
                  (g.geometry.x // cell_m).astype(int).astype(str)).values
    out = []
    for sp in species_list:
        det = set(seq.loc[seq.sp == sp, "deployment_id"])
        if not det:
            continue
        proj = set(d2.loc[d2.deployment_id.isin(det), "project_id"])
        sub = d2[d2.project_id.isin(proj)].copy()
        sub["hit"] = sub.deployment_id.isin(det)
        c = sub.groupby("cell").agg(n=("hit", "size"), k=("hit", "sum"))
        c = c[c.n >= min_deploy_per_cell]
        if len(c) < 8:
            continue
        p = c.k.sum() / c.n.sum()
        out.append({"species": sp, "cells": len(c), "overall_rate": round(p, 3),
                    "spatial_signal_ratio": round(np.var(c.k / c.n, ddof=1)
                                                  / np.mean(p * (1 - p) / c.n), 1)})
    return pd.DataFrame(out).sort_values("spatial_signal_ratio").reset_index(drop=True)
