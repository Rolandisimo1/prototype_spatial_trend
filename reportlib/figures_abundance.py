"""Relative abundance surfaces at three snapshot years, and the change between them.

Supersedes fig_abundance in figures_maps.py, which read the presence-mask (v2b)
snapshots and reconstructed the plotted extent from a mask definition. Two
changes, both correctness rather than style:

  1. Reads the union-mask (v2u) fits -- the reported design.
  2. Masks to <species>_v2u_national_scalar_fitted_cell50.csv, the actual cell
     list the model ran on, rather than rebuilding a mask here. An earlier
     figure clipped to a reconstructed extent and painted area the model never
     saw; reading the fitted list makes that impossible.

Intensity is relative within a species and NOT comparable between species, so
each gets its own scale. The change panel is percent change per cell over the
decade and uses the shared diverging scale (green high, red low, grey midpoint).
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import geopandas as gpd
from matplotlib.colors import TwoSlopeNorm

from . import style
from .figures_maps import conus, _cells_to_squares

SPECIES = [("moose", "Moose"), ("bobcat", "Bobcat"), ("coyote", "Coyote"),
           ("white-tailed_deer", "White-tailed deer")]
YEARS = [2015, 2020, 2025]
CELL50_M = 50_000

# Camera-anchored trend (year_beta, %/yr) -- the REPORTED quantity. mu is the
# iNaturalist-observed intensity, so the change between snapshots is driven by
# the iNaturalist-side trend instead. Both are shown so the figure cannot be
# read as though the change column were the camera result.
CAMERA_PCT_YR = {"moose": -3.6, "bobcat": 1.0, "coyote": 2.3,
                 "white-tailed_deer": 4.0}


def _load(key):
    """Snapshots joined to the fitted cell list. Fails loudly on any mismatch."""
    mu = pd.read_csv(f"{key}_v2u_national_scalar_mu_snapshots_single.csv")
    fit = pd.read_csv(f"{key}_v2u_national_scalar_fitted_cell50.csv")
    missing = set(mu.cell50) - set(fit.cell50)
    if missing:
        raise AssertionError(
            f"{key}: {len(missing)} snapshot cells are absent from the fitted cell "
            "list -- the surface would paint area the model never saw.")
    years = sorted(mu.year.unique())
    if years != YEARS:
        raise AssertionError(f"{key}: snapshot years {years}, expected {YEARS}")
    return mu[mu.cell50.isin(set(fit.cell50))].copy(), fit


def abundance_change_table():
    """Per-species decadal change in mean relative intensity across fitted cells."""
    rows = []
    for key, label in SPECIES:
        mu, fit = _load(key)
        w = mu.pivot_table(index="cell50", columns="year", values="mu_mean")
        pc = 100.0 * (w[YEARS[-1]] / w[YEARS[0]] - 1.0)
        pc = pc.replace([np.inf, -np.inf], np.nan).dropna()
        rows.append({"species": label, "fitted_cells": len(fit),
                     "median_pct_change": round(float(pc.median()), 2),
                     "q25": round(float(pc.quantile(.25)), 2),
                     "q75": round(float(pc.quantile(.75)), 2),
                     "cells_increasing": int((pc > 0).sum()),
                     "cells_decreasing": int((pc < 0).sum())})
    return pd.DataFrame(rows)


def fig_abundance_union(out="fig_abundance_union.png"):
    st = conus()
    fig, axes = plt.subplots(len(SPECIES), 4, figsize=(11.4, 9.4))

    for r, (key, label) in enumerate(SPECIES):
        mu, fit = _load(key)
        w = mu.pivot_table(index="cell50", columns="year", values="mu_mean")
        coords = fit.set_index("cell50")[["lon", "lat"]]

        # one scale per species across its three years: absolute intensity is
        # relative within a species and meaningless between them
        vmax = float(np.nanpercentile(w[YEARS].values, 98))
        for c, yr in enumerate(YEARS):
            ax = axes[r, c]
            g = _cells_to_squares(
                coords.join(w[yr].rename("v")).reset_index(), CELL50_M)
            st.plot(ax=ax, facecolor=style.LAND_FILL, edgecolor="0.80",
                    linewidth=0.25, zorder=0)
            g.plot(column="v", ax=ax, cmap="viridis", vmin=0, vmax=vmax,
                   edgecolor="none", zorder=2)
            ax.set_axis_off()
            ax.set_xlim(*st.total_bounds[[0, 2]]); ax.set_ylim(*st.total_bounds[[1, 3]])
            if r == 0:
                ax.set_title(str(yr), fontsize=9.0)
            if c == 0:
                ax.text(-0.02, 0.5, label, transform=ax.transAxes, rotation=90,
                        ha="right", va="center", fontsize=8.8)

        pc = (100.0 * (w[YEARS[-1]] / w[YEARS[0]] - 1.0)).replace(
            [np.inf, -np.inf], np.nan)
        lim = float(np.nanpercentile(np.abs(pc.dropna()), 98)) or 1.0
        ax = axes[r, 3]
        g = _cells_to_squares(coords.join(pc.rename("v")).reset_index(), CELL50_M)
        st.plot(ax=ax, facecolor=style.LAND_FILL, edgecolor="0.80",
                linewidth=0.25, zorder=0)
        g.plot(column="v", ax=ax, cmap=style.DIVERGING,
               norm=TwoSlopeNorm(vcenter=0, vmin=-lim, vmax=lim),
               edgecolor="none", zorder=2)
        ax.set_axis_off()
        ax.set_xlim(*st.total_bounds[[0, 2]]); ax.set_ylim(*st.total_bounds[[1, 3]])
        if r == 0:
            ax.set_title(f"iNaturalist-side change\n{YEARS[0]} to {YEARS[-1]}", fontsize=8.6)
        med = float(pc.dropna().median())
        ax.text(0.5, -0.04,
                f"median {med:+.1f}%   |   cameras {CAMERA_PCT_YR[key]:+.1f}%/yr",
                transform=ax.transAxes, ha="center", va="top",
                fontsize=7.0, color="#5f6a73")

    fig.suptitle("Where each species occurs, and how that changed over the decade",
                 fontsize=10.4, y=0.995)
    fig.text(0.5, 0.055,
             "Relative abundance from the integrated model, mapped on the 50 km cells each fit "
             "actually used. Yellow is high, dark purple is low. Intensity is relative within a\nspecies and not "
             "comparable between species, so each row has its own scale; the change column is "
             "percent change per cell, green for increase.\n"
             "The change column is the modelled iNaturalist-observed intensity, which the "
             "iNaturalist-side trend drives; the camera-anchored trend we report is printed "
             "beside each median.\nOne national trend applies everywhere, so within-row variation "
             "there is that trend acting on a varying baseline, not regional variation in trend.",
             ha="center", va="top", fontsize=6.6, color="#5f6a73")
    fig.tight_layout(rect=[0.01, 0.075, 1, 0.985])
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=240)
    plt.close(fig)
    return out
