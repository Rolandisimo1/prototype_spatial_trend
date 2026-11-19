"""Per-deployment change in fitted occupancy, mapped.

The national-scalar model applies one trend to every cell, so a map of
proportional change in abundance would be one flat colour. Fitted occupancy is
different: cloglog saturates, so the same slope on the log-abundance scale moves
psi by different amounts depending on where a deployment starts. A site already
near psi = 1 can barely rise, while a sparsely used site has room to move. This
figure maps that, which is a statement about what the cameras can detect at each
site as much as about the animals.

Input is `<run>_psi_snapshots.csv`, produced by the extraction step: one row per
camera deployment per camera year of the window, carrying site_id, cell50,
cell100, coordinates and the fitted psi.
"""
import glob
import os

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import figures_win5 as fw
from .figures_maps import conus

PSI_GLOB = "*_psi_snapshots.csv"
REQUIRED = {"site_id", "year", "lon", "lat"}
EPS = 1e-12


def _psi_col(d):
    """The fitted-occupancy column, whatever the extraction named it."""
    for c in ("psi_mean", "psi", "psi_med", "psi_median"):
        if c in d.columns:
            return c
    raise KeyError(f"no fitted-occupancy column in {sorted(d.columns)}")


def _cloglog(p):
    p = np.clip(np.asarray(p, float), EPS, 1 - EPS)
    return np.log(-np.log1p(-p))


def _inv_cloglog(eta):
    return -np.expm1(-np.exp(eta))


def deployment_change(folder=None, products=None, posteriors=None):
    """Change in fitted occupancy at each camera deployment across the window.

    A deployment is fitted in the one year it ran: across the fleet almost no
    site appears in both the first and last year of the window, so a before-and-
    after difference at the same site is unavailable. Instead each deployment is
    projected to both ends of the window, which is exact rather than an
    approximation. The trend is one scalar on the cloglog scale, so

        cloglog(psi_i at year t) = eta_i + beta_yr * standardised(t)

    and subtracting the deployment's own year recovers eta_i, the part of the
    linear predictor carrying its location, habitat and climate. Projecting that
    to the first and last year gives the change the model implies at that site.

    The difference between sites comes entirely from where each one starts. One
    slope on the cloglog scale moves a sparsely used site by many percentage
    points and a site already near certain use by almost none, so this maps the
    saturation, which is a statement about what a camera at that site can detect.
    """
    products = products or os.path.join(fw.PULL, "products")
    posteriors = posteriors or fw.PULL
    rows, skipped = [], []
    for f in sorted(glob.glob(os.path.join(products, PSI_GLOB))):
        d = pd.read_csv(f)
        missing = REQUIRED - set(d.columns)
        if missing:
            raise AssertionError(f"{os.path.basename(f)} lacks {sorted(missing)}")
        run = os.path.basename(f)[: -len("_psi_snapshots.csv")]
        ysf = os.path.join(posteriors, f"year_scaling_{run}.csv")
        pof = glob.glob(os.path.join(posteriors, f"{run.split('_national')[0]}*_posterior_win*.csv"))
        if not (os.path.exists(ysf) and pof):
            # A run whose scaling and posterior live in another pull directory
            # cannot be projected from here. Record it instead of failing, so a
            # five-year call is not blocked by a ten-year run sitting alongside,
            # and return the list so the omission is visible to the caller.
            skipped.append(run)
            continue
        ys = pd.read_csv(ysf)
        cam = ys[ys.coefficient == "year_beta"].iloc[0]
        po = pd.read_csv(pof[0])
        b = po[po.param.astype(str).str.fullmatch("year_beta")]
        if not len(b):
            raise AssertionError(f"{run}: year_beta absent from {os.path.basename(pof[0])}")
        beta = float(b["mean"].iloc[0])
        sp = d.species.iloc[0] if "species" in d.columns else run.split("_v2u")[0]
        pc = _psi_col(d)
        y0, y1 = int(d.year.min()), int(d.year.max())
        z = lambda t: (t - float(cam.year_mean)) / float(cam.year_sd)
        eta = _cloglog(d[pc].to_numpy()) - beta * z(d.year.to_numpy())
        g = pd.DataFrame({
            "species": fw._pretty(sp), "run": run,
            "site_id": d.site_id.to_numpy(), "year_sampled": d.year.to_numpy(),
            "year_first": y0, "year_last": y1,
            "psi_at_first": _inv_cloglog(eta + beta * z(y0)),
            "psi_at_last": _inv_cloglog(eta + beta * z(y1)),
            "psi_as_fitted": d[pc].to_numpy(),
            "lon": d.lon.to_numpy(), "lat": d.lat.to_numpy(),
            "year_beta": beta,
        })
        g["change_pp"] = 100 * (g.psi_at_last - g.psi_at_first)
        rows.append(g.drop_duplicates(subset=["site_id"]))
    if not rows:
        raise AssertionError(
            f"no projectable psi_snapshots in {products}; skipped {skipped}")
    out = pd.concat(rows, ignore_index=True)
    out.attrs["skipped_runs"] = skipped
    return out


def fig_deployment_change(out="fig_deployment_change.png", products=None,
                          posteriors=None, table=None):
    """Map the per-deployment change in fitted occupancy, one panel per species.

    Deployments are drawn as points coloured by change in percentage points on a
    scale centred at zero, so the colour carries sign. Each panel names how many
    deployments were projected, which is every deployment in that species' window.
    """
    t = table if table is not None else deployment_change(
        products=products, posteriors=posteriors)
    if not len(t):
        raise AssertionError(
            "no psi_snapshots files found; this figure needs the extraction "
            f"products in {products or os.path.join(fw.PULL, 'products')}")
    sp = sorted(t.species.unique())
    st = conus()
    lim = float(np.nanpercentile(np.abs(t.change_pp), 98)) or 1.0
    ncol = min(3, len(sp))
    nrow = int(np.ceil(len(sp) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(4.1 * ncol, 2.9 * nrow),
                             squeeze=False)
    for ax in axes.ravel():
        ax.set_axis_off()
    for k, s in enumerate(sp):
        ax = axes[k // ncol, k % ncol]
        ax.set_axis_off()
        st.boundary.plot(ax=ax, color="0.80", linewidth=0.3, zorder=1)
        d = t[t.species == s]
        # conus() is projected; deployment coordinates are degrees. Plotting one
        # onto the other silently stacks every point at the same place, which is
        # what the first render did, so the reprojection is not optional.
        pts = gpd.GeoDataFrame(d, geometry=gpd.points_from_xy(d.lon, d.lat),
                               crs=4326).to_crs(st.crs)
        sc = ax.scatter(pts.geometry.x, pts.geometry.y, c=d.change_pp, cmap="RdBu_r",
                        vmin=-lim, vmax=lim, s=5, linewidths=0, zorder=3)
        ax.set_xlim(*st.total_bounds[[0, 2]]); ax.set_ylim(*st.total_bounds[[1, 3]])
        y0, y1 = int(d.year_first.iloc[0]), int(d.year_last.iloc[0])
        # Not "sampled in both years": almost no deployment ran in both. Each is
        # fitted in its own year and projected to both ends of the window.
        ax.set_title(f"{s}\n{len(d):,} deployments projected to {y0} and {y1}",
                     fontsize=7.6)
    cb = fig.colorbar(sc, ax=axes.ravel().tolist(), fraction=0.022, pad=0.012,
                      shrink=0.6)
    cb.set_label("Change in fitted occupancy\n(percentage points)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    fig.text(0.5, 0.015,
             "Each deployment is fitted in the year it ran, then projected to the first and "
             "last year of the window using the fitted trend. Almost no deployment ran in\n"
             "both years, so this is the change the model implies at each site rather than a "
             "measured before and after. One slope on the cloglog scale moves a sparsely\nused "
             "site by many points and a site near certain use by almost none, so the spread "
             f"within a species is saturation. Scale saturates at {lim:.0f} points.",
             ha="center", va="top", fontsize=6.3, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return out, t
