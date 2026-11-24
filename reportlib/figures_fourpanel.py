"""One sheet per species: the two year surfaces, their difference, and the cameras.

Four panels for a species whose trend resolved. The first three come from the
fitted abundance surface; the fourth is the camera arrays, which is a different
measurement of the same change and is the one a field ecologist checks first.

Why four panels rather than the three-panel surface figure plus a separate array
map: the comparison that matters is between the modelled change and what the
cameras saw, and that comparison is only available when the two are on one page.

The surface panels share one colour scale so the level difference between the
two years IS the change and can be read directly. The change panel and the
camera panel each carry their own scale and their own units, because one is a
percent per year in abundance and the other a percentage-point change in
detection probability. They are not interchangeable and are not drawn as if
they were.
"""
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm, TwoSlopeNorm

from .figures_deployment import DARK_BG, DARK_INK, DARK_RULE, DARK_SUB, trend_status
from .figures_products import (PULL, PRODUCTS, _squares, abundance_change, conus,
                               reported_trend, surface_panels)


def _blank(ax):
    ax.set_axis_off()
    ax.set_facecolor(DARK_BG)


def fig_species_four(species, out=None, panels=None, change=None, trends=None,
                     arrays=None, resolved=None, products=PRODUCTS, pull=PULL,
                     posteriors=None, slide=False):
    """Four panels for one species: two year surfaces, the change, the cameras.

    `slide=True` scales the type up and drops the footer, matching the other
    projection variants in this package.
    """
    fs = 1.45 if slide else 1.0
    pn = (panels if panels is not None else surface_panels(products, pull))
    ch = (change if change is not None else abundance_change(products, pull))
    tr = (trends if trends is not None else reported_trend(pull))
    if resolved is None:
        resolved = trend_status(posteriors)
    d = pn[pn.species == species]
    dc = ch[ch.species == species]
    if not len(d):
        raise ValueError(f"no fitted surface for {species!r}")
    st = conus()
    g, gc = _squares(d, crs_to=st.crs), _squares(dc, crs_to=st.crs)
    y0, y1 = int(d.year_first.iloc[0]), int(d.year_last.iloc[0])
    rep = tr.get(species, {})

    fig, axes = plt.subplots(2, 2, figsize=(13.4, 10.2), facecolor=DARK_BG)
    fig.subplots_adjust(hspace=0.16, wspace=0.03)
    for ax in axes.ravel():
        _blank(ax)

    lo = max(min(d.rel_first.min(), d.rel_last.min()), 1e-6)
    norm = LogNorm(vmin=lo, vmax=1.0)
    for j, (col, yr) in enumerate([("rel_first", y0), ("rel_last", y1)]):
        ax = axes[0, j]
        st.boundary.plot(ax=ax, color=DARK_RULE, linewidth=0.5, zorder=1)
        g.plot(column=col, ax=ax, cmap="viridis", norm=norm, linewidth=0, zorder=3)
        ax.set_title(f"Relative abundance, {yr}", fontsize=14 * fs,
                     color=DARK_INK, pad=8)
    cb = fig.colorbar(plt.cm.ScalarMappable(cmap="viridis", norm=norm),
                      ax=axes[0, :].tolist(), fraction=0.030, pad=0.012, shrink=0.78)
    cb.set_label("relative abundance, shared scale across both years",
                 fontsize=10 * fs, color=DARK_SUB)
    cb.ax.tick_params(labelsize=9 * fs, colors=DARK_INK)

    ax = axes[1, 0]
    lim = float(np.abs(gc.pct_per_year).max()) or 1.0
    st.boundary.plot(ax=ax, color=DARK_RULE, linewidth=0.5, zorder=1)
    gc.plot(column="pct_per_year", ax=ax, cmap="RdBu_r",
            norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim), linewidth=0, zorder=3)
    rp = rep.get("pct_per_year")
    ax.set_title("Change in abundance"
                 + (f"   reported {rp:+.1f}%/yr" if rp is not None else ""),
                 fontsize=14 * fs, color=DARK_INK, pad=8)
    cb2 = fig.colorbar(plt.cm.ScalarMappable(
        cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)),
        ax=ax, fraction=0.030, pad=0.012, shrink=0.80, orientation="horizontal")
    cb2.set_label("% per year", fontsize=10 * fs, color=DARK_SUB)
    cb2.ax.tick_params(labelsize=9 * fs, colors=DARK_INK)

    ax = axes[1, 1]
    import geopandas as gpd
    t = arrays[arrays.species == species]
    pts = gpd.GeoDataFrame(t, geometry=gpd.points_from_xy(t.lon, t.lat),
                           crs=4326).to_crs(st.crs)
    alim = float(np.abs(t.change_pp).max()) or 1.0
    st.boundary.plot(ax=ax, color=DARK_RULE, linewidth=0.5, zorder=1)
    size = 16 + 70 * (t.n_deployments / float(t.n_deployments.max())) ** 0.5
    sc = ax.scatter(pts.geometry.x, pts.geometry.y, c=t.change_pp, cmap="RdBu_r",
                    vmin=-alim, vmax=alim, s=size, linewidths=0.3,
                    edgecolors=DARK_BG, zorder=3)
    ax.set_title(f"What the cameras saw   {len(t)} arrays, "
                 f"{int(t.n_deployments.sum()):,} camera-years",
                 fontsize=14 * fs, color=DARK_INK, pad=8)
    cb3 = fig.colorbar(sc, ax=ax, fraction=0.030, pad=0.012, shrink=0.80,
                       orientation="horizontal")
    cb3.set_label("change in detection probability (pp)", fontsize=10 * fs,
                  color=DARK_SUB)
    cb3.ax.tick_params(labelsize=9 * fs, colors=DARK_INK)

    for a in axes.ravel():
        a.set_xlim(*st.total_bounds[[0, 2]])
        a.set_ylim(*st.total_bounds[[1, 3]])

    fig.suptitle(species, fontsize=20 * fs, color=DARK_INK, y=0.965)
    if not slide:
        fig.text(0.5, 0.045,
                 "Top row: the fitted relative-abundance surface at each end of the "
                 "window, on one shared scale, so the difference between the panels is "
                 "the change.\nBottom left: that change as a rate, anchored to the "
                 "reported camera trend. Bottom right: the same window measured at the "
                 "camera arrays, as the change in\nan array's probability of detecting "
                 "the species. The two bottom panels are different quantities in "
                 "different units and share no scale.",
                 ha="center", va="top", fontsize=8.6, color=DARK_SUB)
    out = out or f"fig_four_{species.lower().replace(' ', '_')}.png"
    fig.savefig(out, facecolor=DARK_BG, bbox_inches="tight", pad_inches=0.22, dpi=150)
    plt.close(fig)
    return out


def fig_four_all(species_list, stem="fig_four", **kw):
    """One four-panel sheet per resolved species."""
    return [fig_species_four(s_, out=f"{stem}_{s_.lower().replace(' ', '_')}.png", **kw)
            for s_ in species_list]
