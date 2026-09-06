"""Spatially varying model fields: the baseline surface and the climate responses.

These are habitat and climate associations, not population trends, so they are
estimated from the full record deliberately -- more years of data make them
better, whereas the trend is compromised by the sparse early camera years.

All fields come from the national-model fits. The ecoregion fits carry a
localised spatial-field convergence failure in one contiguous patch of cells;
the national fits do not (largest spatial-field R-hat 1.03 to 1.10 across the
three species), so they are the sound basis for these maps.
"""
import os

import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from . import inputs, style
from .figures_maps import (GRID_CRS, CELL100, conus, union_mask, _cells_to_squares,
                           SPECIES)

CAR_FILE = {
    "bobcat": "bobcat_v2b_national_scalar_spatial_car_fields_single.csv",
    "white-tailed_deer": "white-tailed_deer_v2b_national_scalar_spatial_car_fields_single.csv",
    "moose": "moose_v2b_national_scalar_spatial_car_fields_single.csv",
}


def car_fields(species):
    """Spatial fields for one species, clipped to the lower 48 and its range mask."""
    d = pd.read_csv(os.path.join(inputs.REPO, CAR_FILE[species]))
    if d[["lon", "lat"]].isna().any().any():
        raise AssertionError(f"{species}: CAR file has missing lon/lat")
    g = _cells_to_squares(d, CELL100)
    keep = g.geometry.intersects(union_mask(species)) & g.geometry.intersects(
        conus().union_all())
    return g[keep].copy()


def car_significance_table():
    """How many cells each field resolves away from zero, per species and field."""
    rows = []
    for key, label in SPECIES:
        g = car_fields(key)
        for field, fl in [("link_occ_intercept", "Spatial baseline"),
                          ("MWMT_effect", "Summer warmth response"),
                          ("MCMT_effect", "Winter cold response")]:
            sig = ((g[f"{field}_q025"] > 0) | (g[f"{field}_q975"] < 0)).sum()
            rows.append({"species": label, "field": fl, "cells_in_range": len(g),
                         "cells_excluding_zero": int(sig),
                         "pct_excluding_zero": round(100 * sig / len(g), 1)})
    return pd.DataFrame(rows)


UNRESOLVED_ALPHA = 0.38


def _draw_field(ax, g, field, norm, st):
    """Draw one spatial field, encoding resolution by opacity.

    Cells whose 95% interval excludes zero are drawn at full opacity; the rest
    are faded. Outlining the resolved cells was tried first and fails when they
    are the majority -- at 65% coverage the outline becomes a grid over the whole
    map and stops carrying information. Opacity works at any coverage and is one
    encoding for both the strong baseline field and the weak climate fields.
    """
    st.plot(ax=ax, facecolor=style.LAND_FILL, edgecolor="0.80",
            linewidth=0.3, zorder=0)
    resolved = (g[f"{field}_q025"] > 0) | (g[f"{field}_q975"] < 0)
    for subset, alpha in [(g[~resolved], UNRESOLVED_ALPHA), (g[resolved], 1.0)]:
        if len(subset):
            subset.plot(column=f"{field}_mean", ax=ax, cmap=style.DIVERGING, norm=norm,
                        edgecolor="none", alpha=alpha, zorder=2)
    ax.set_axis_off()
    ax.set_xlim(*st.total_bounds[[0, 2]])
    ax.set_ylim(*st.total_bounds[[1, 3]])
    return g


def fig_spatial_baseline(out="fig_spatial_baseline.png"):
    """Where habitat alone mispredicts occurrence -- the spatial residual field.

    This is the strongest of the three spatial fields -- most cells resolve away
    from zero -- so it is shown on its own rather than alongside the much weaker
    climate responses.
    """
    st = conus()
    gs, vmax = {}, 0.0
    for key, _ in SPECIES:
        g = car_fields(key)
        gs[key] = g
        vmax = max(vmax, float(np.abs(g["link_occ_intercept_mean"]).max()))
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.0))
    counts = {}
    for ax, (key, label) in zip(axes, SPECIES):
        g = _draw_field(ax, gs[key], "link_occ_intercept", norm, st)
        sig = ((g["link_occ_intercept_q025"] > 0) | (g["link_occ_intercept_q975"] < 0)).sum()
        counts[label] = (int(sig), len(g))
        ax.set_title(f"{label}   (clear in {sig} of {len(g)} cells)", fontsize=8.4)

    sm = plt.cm.ScalarMappable(cmap=style.DIVERGING, norm=norm)
    cb = fig.colorbar(sm, ax=axes.tolist(), fraction=0.019, pad=0.010, shrink=0.66,
                      anchor=(0.0, 0.72))
    cb.set_label("Spatial baseline\n(cloglog scale; green = higher)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    fig.suptitle("Places the habitat map alone gets wrong",
                 fontsize=9.4, y=1.02)
    fig.text(0.5, 0.02,
             "The model first predicts where each species should live from habitat alone -- "
             "forest type, farmland, terrain, soil, human population. This map shows\n"
             "what that prediction misses. GREEN means the species turns up more often than its "
             "habitat would suggest; RED means less often; GREY means the habitat\n"
             "prediction was about right. The pattern points to things the habitat layers do not "
             "measure -- harvest pressure, disease, barriers to movement, recent range\n"
             "expansion or contraction. Faded cells are places where the data cannot yet tell "
             "the difference. Lower 48 only, inside each species' range.",
             ha="center", va="top", fontsize=6.4, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, counts


def fig_climate_response(out="fig_climate_response.png"):
    """Spatially varying responses to summer warmth and winter cold.

    Reported with the significance count on every panel because the honest
    headline is that these fields resolve away from zero in only a small
    minority of cells: the spatial pattern is suggestive, the per-cell evidence
    is weak, and a map alone would overstate it.
    """
    st = conus()
    gs = {key: car_fields(key) for key, _ in SPECIES}
    fields = [("MWMT_effect", "Response to summer warmth"),
              ("MCMT_effect", "Response to winter cold")]

    # One scale per climate variable, shared across the three species in that
    # row. A single scale across all six panels was tried first and is set by one
    # extreme cell, which flattens every other panel to near-white. The scale is
    # also robust (98th percentile) for the same reason; the number of cells that
    # saturate is counted and disclosed in the caption rather than hidden.
    norms, clipped = {}, 0
    for field, _ in fields:
        allv = np.concatenate([np.abs(gs[k][f"{field}_mean"].values) for k, _ in SPECIES])
        vmax = float(np.percentile(allv, 98))
        norms[field] = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
        clipped += int((allv > vmax).sum())

    fig, axes = plt.subplots(2, 3, figsize=(11.0, 5.6))
    for i, (field, flabel) in enumerate(fields):
        for j, (key, slabel) in enumerate(SPECIES):
            ax = axes[i, j]
            g = _draw_field(ax, gs[key], field, norms[field], st)
            sig = ((g[f"{field}_q025"] > 0) | (g[f"{field}_q975"] < 0)).sum()
            ax.set_title(f"{slabel}   ({sig} of {len(g)} resolved)", fontsize=7.8)
            if j == 0:
                ax.text(-0.04, 0.5, flabel, transform=ax.transAxes, rotation=90,
                        ha="right", va="center", fontsize=8.4)
        sm = plt.cm.ScalarMappable(cmap=style.DIVERGING, norm=norms[field])
        cb = fig.colorbar(sm, ax=axes[i, :].tolist(), fraction=0.016, pad=0.008,
                          shrink=0.80)
        cb.set_label(f"{flabel.replace('Response to ', '').capitalize()}\n"
                     "(green = higher occurrence)", fontsize=6.6)
        cb.ax.tick_params(labelsize=5.8)
    self_clipped = clipped
    fig.suptitle("Climate responses vary in space, but are weakly resolved at any one place",
                 fontsize=9.4, y=0.99)
    fig.text(0.5, 0.055,
             "Each cell has its own response coefficient, smoothed toward its neighbours. "
             "Opacity marks resolution: full colour where the 95% interval\nexcludes zero, "
             "faded where it does not, with those counts in each title. Because only a small "
             "minority of cells resolve, read these as spatial\npatterns worth investigating "
             "rather than established local effects. Each row has its own colour scale, "
             f"saturating at the 98th percentile\n({self_clipped} of "
             f"{2 * sum(len(gs[k]) for k, _ in SPECIES)} cell values exceed it). "
             "Clipped to the lower 48 and each species' range mask.",
             ha="center", va="top", fontsize=6.4, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, car_significance_table()
