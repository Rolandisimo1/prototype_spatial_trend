"""Year-by-year abundance surfaces: animation, and a trend that uses every year.

Why this exists. The three-snapshot figure compares 2015 against 2025 and
discards everything between. That is not a cosmetic loss: the fitted surfaces
carry real year-to-year structure, because the iNaturalist covariate array is
indexed by year (xdat_inat[row, covariate, year]) and the covariates do change.
Testing the three snapshots for curvature -- the two equally-spaced 5-year
log-steps must be identical if the surface were a pure exponential ramp -- finds
departures in 1,213 to 2,795 cells per species. Moose is the clearest case: its
mean log-step is +0.003 over 2015-2020 and -0.118 over 2020-2025, so an endpoint
comparison averages a flat period and a steep decline into one misleading number.

So both functions here deliberately use every year:

  per_cell_trend_all_years  -- log-linear slope per cell over the full record,
      which is the standard way to summarise a time series and is far less
      sensitive to what happened to fall on the endpoints.
  animate_surface           -- one frame per year, so year-to-year movement can
      be judged by eye for whether it looks like biology or like noise.

Both require the per-year snapshots. They raise rather than silently falling
back to the 3-year file, because a "yearly" animation built from 3 frames would
look like a smooth ramp whatever the data did.
"""
import glob
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from PIL import Image

from . import style
from .figures_maps import conus, _cells_to_squares

CELL50_M = 50_000
MIN_YEARS = 10

# Three different %/yr quantities exist for each species and they do NOT agree.
# mu is the iNaturalist-OBSERVED intensity, so a trend fitted to it is the iNat
# trend parameter attenuated by theta1 (verified: theta1 x parameter reproduces
# the fitted slope to 0.04 points for bobcat and deer, 0.28 worst case). Showing
# a bare "%/yr" here put coyote at -0.3 in this figure and +2.3 in the trend
# figure with the distinction only in a caption, which is how it should NOT read.
CAMERA_PCT_YR = {"moose": -3.6, "bobcat": 1.0, "coyote": 2.3, "white-tailed_deer": 4.0}
INAT_PARAM_PCT_YR = {"moose": -2.2, "bobcat": -3.4, "coyote": -0.3, "white-tailed_deer": 2.7}


def _yearly(key):
    """Per-year snapshots for one species, masked to its fitted cell list."""
    f = f"{key}_v2u_national_scalar_mu_yearly.csv"
    if not os.path.exists(f):
        raise FileNotFoundError(
            f"{f} not found in {os.getcwd()}. The 3-year snapshot file is not a "
            "substitute: an animation built from 3 frames reads as a smooth ramp "
            "regardless of the data. Request the per-year extraction.")
    mu = pd.read_csv(f)
    fit = pd.read_csv(f"{key}_v2u_national_scalar_fitted_cell50.csv")
    years = sorted(mu.year.unique())
    if len(years) < MIN_YEARS:
        raise AssertionError(f"{key}: only {len(years)} years ({years})")
    extra = set(mu.cell50) - set(fit.cell50)
    if extra:
        raise AssertionError(
            f"{key}: {len(extra)} cells absent from the fitted list -- would paint "
            "area the model never saw.")
    return mu[mu.cell50.isin(set(fit.cell50))].copy(), fit, years


def per_cell_trend_all_years(key):
    """Log-linear slope per cell over every year, as percent per year."""
    mu, fit, years = _yearly(key)
    w = mu.pivot_table(index="cell50", columns="year", values="mu_mean")
    w = w[[y for y in years]].dropna()
    yr = np.asarray(years, dtype=float)
    yr_c = yr - yr.mean()
    slope = (np.log(w.values) @ yr_c) / (yr_c @ yr_c)
    out = pd.DataFrame({"cell50": w.index, "pct_per_year": 100.0 * (np.exp(slope) - 1.0)})
    return out.merge(fit[["cell50", "lon", "lat"]], on="cell50")


def trend_method_comparison(keys):
    """All-years slope against the 2015-2025 endpoint estimate, per species."""
    rows = []
    for key, label in keys:
        allyr = per_cell_trend_all_years(key)
        snap = pd.read_csv(f"{key}_v2u_national_scalar_mu_snapshots_single.csv")
        w = snap.pivot_table(index="cell50", columns="year", values="mu_mean")
        ends = 100.0 * ((w[2025] / w[2015]) ** (1.0 / 10.0) - 1.0)
        j = allyr.set_index("cell50").join(ends.rename("endpoint"), how="inner")
        rows.append({"species": label, "cells": len(j),
                     "all_years_median": round(float(j.pct_per_year.median()), 3),
                     "endpoint_median": round(float(j.endpoint.median()), 3),
                     "median_abs_diff": round(float((j.pct_per_year - j.endpoint).abs().median()), 3),
                     "sign_disagrees": int((np.sign(j.pct_per_year) != np.sign(j.endpoint)).sum())})
    return pd.DataFrame(rows)


def animate_surface(key, label, out=None, fps=2.5):
    """One frame per year. Colour scale fixed across frames so change is real."""
    out = out or f"anim_abundance_{key}.gif"
    mu, fit, years = _yearly(key)
    w = mu.pivot_table(index="cell50", columns="year", values="mu_mean")
    coords = fit.set_index("cell50")[["lon", "lat"]]
    st = conus()
    # One scale for every frame. A per-frame scale would renormalise each year
    # and make a changing surface look static -- the opposite of the question.
    vmax = float(np.nanpercentile(w[years].values, 98))

    frames = []
    for yr in years:
        fig, ax = plt.subplots(figsize=(5.4, 3.5))
        g = _cells_to_squares(coords.join(w[yr].rename("v")).reset_index(), CELL50_M)
        st.plot(ax=ax, facecolor=style.LAND_FILL, edgecolor="0.80", linewidth=0.25, zorder=0)
        g.plot(column="v", ax=ax, cmap="viridis", vmin=0, vmax=vmax, edgecolor="none", zorder=2)
        ax.set_axis_off()
        ax.set_xlim(*st.total_bounds[[0, 2]]); ax.set_ylim(*st.total_bounds[[1, 3]])
        ax.set_title(f"{label} — modelled relative abundance, {yr}", fontsize=9.0)
        # Below the axes, not inside them: placed at the lower right of the map it
        # overlapped the southeastern cells on three of four species.
        fig.text(0.5, 0.045, "yellow high, dark purple low   |   fixed colour scale across years",
                 ha="center", va="top", fontsize=6.0, color="#5f6a73")
        fig.subplots_adjust(bottom=0.11)
        fig.canvas.draw()
        # Canvas.tostring_rgb() was removed in matplotlib 3.10; buffer_rgba is
        # the supported path and needs the alpha channel dropped.
        w_px, h_px = fig.canvas.get_width_height()
        frames.append(Image.frombytes(
            "RGBA", (w_px, h_px), bytes(fig.canvas.buffer_rgba())).convert("RGB"))
        plt.close(fig)

    frames[0].save(out, save_all=True, append_images=frames[1:],
                   duration=int(1000 / fps), loop=0, optimize=True)
    return out, len(frames)


def fig_trajectory(keys, out="fig_trajectory.png"):
    """National mean surface per year against its own straight-line trend.

    The quantitative counterpart to the animation: it shows whether the
    year-to-year departures are structured or scattered. A straight line in log
    space is what the model's single trend parameter would produce on its own,
    so visible departures are the six year-varying covariates acting through
    occ_beta -- NDVI above all, which changes in every sub-cell row every year.
    Those departures are habitat variation, not population dynamics.
    """
    fig, axes = plt.subplots(1, len(keys), figsize=(2.75 * len(keys), 3.0), sharex=True)
    for ax, (key, label) in zip(np.atleast_1d(axes), keys):
        mu, fit, years = _yearly(key)
        nat = mu.groupby("year").mu_mean.mean()
        yr = np.asarray(years, float)
        l = np.log(nat.values)
        X = np.column_stack([np.ones_like(yr), yr - yr.mean()])
        b = np.linalg.lstsq(X, l, rcond=None)[0]
        resid = l - X @ b
        r2 = 1 - resid.var() / l.var()

        ax.plot(yr, np.exp(X @ b), color="0.62", lw=1.3, ls="--", zorder=2,
                label="single-trend fit")
        ax.plot(yr, nat.values, color=style.FOCAL, lw=1.7, marker="o", ms=3.0,
                zorder=3, label="modelled surface")
        ax.set_title(f"{label}\niNat counts {100*(np.exp(b[1])-1):+.2f}%/yr   "
                     f"(R$^2$={r2:.3f})\n"
                     f"cameras {CAMERA_PCT_YR[key]:+.1f}%/yr", fontsize=7.6)
        ax.set_xticks([2008, 2014, 2020, 2025])
        ax.tick_params(labelsize=6.6)
        ax.margins(y=0.16)
    np.atleast_1d(axes)[0].set_ylabel("Mean iNaturalist-observed\nintensity", fontsize=7.4)
    np.atleast_1d(axes)[0].legend(frameon=False, fontsize=6.2, loc="best")
    fig.text(0.5, -0.14,
             "Each panel is the mean of that species' fitted iNaturalist intensity across its own "
             "cells, by year. This is NOT the reported population trend:\nit is the iNaturalist-side "
             "trend attenuated by theta1, the exponent linking iNaturalist counts to intensity "
             "(0.52-0.85 here), so it runs shallower\nthan either stream's coefficient. The "
             "camera-anchored trend we report is printed in each title for comparison, and the two "
             "disagree in sign for\nbobcat and coyote. The dashed line is a single log-linear trend. "
             "Departures from it come from the six year-varying covariates, dominated by\nNDVI, so "
             "they are habitat variation rather than population dynamics.",
             ha="center", va="top", fontsize=6.3, color="#5f6a73")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out
