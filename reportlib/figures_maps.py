"""Map figures: abundance surfaces, regional trends, and agency comparison.

Grid CRS was determined empirically, not assumed: joining cell100 to the
CAR-field file (which carries lon/lat for the same cells) and reprojecting gives
a median discrepancy of 0 m against EPSG:5070 over all 893 cells. The 100 km
lattice spacing is exactly 100,000 m in both axes and lon/lat projects onto the
reference point, so those coordinates are cell centres.
"""
import os
import numpy as np
import pandas as pd
import geopandas as gpd
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from matplotlib.patches import Patch

from . import inputs, style

GRID_CRS = 5070
CELL100 = 100_000
CELL50 = 50_000
NON_CONUS = {"HI", "PR", "VI", "GU", "MP", "AS", "AK"}
NO_TREND_GREY = "#e6e6e6"   # interval overlaps zero: no direction claimed.
                            # Lighter than the colour scale's own grey midpoint,
                            # so "no claim" is not confused with "near zero".
RAW_AGENCY_MAP = {"moose": "fig9_moose_agency_trend_map.png",
                  "white-tailed_deer": "fig10_wtd_agency_trend_map.png"}

SPECIES = [("bobcat", "Bobcat"), ("white-tailed_deer", "White-tailed deer"),
           ("moose", "Moose")]
MU_FILE = {"bobcat": "bobcat_v2b_national_scalar_mu_snapshots_single.csv",
           "white-tailed_deer": "white-tailed_deer_v2b_national_scalar_mu_snapshots_single.csv",
           "moose": "moose_v2b_national_scalar_mu_snapshots_single.csv"}
ECO_FILE = {"bobcat": ("pull", "bobcat_v2b_ecoregion_ecoregion_posterior.csv"),
            "white-tailed_deer": ("pull", "white-tailed_deer_v2b_ecoregion_ecoregion_posterior.csv"),
            "moose": ("repo", "moose_v2b_ecoregion_posterior_single.csv")}


def _split_by_significance(m, col="abs_p_negative", lo=0.05, hi=0.95):
    """Split regions into those with a clear direction and those overlapping zero."""
    unclear = (m[col] > lo) & (m[col] < hi)
    return m[~unclear], m[unclear]


def conus():
    """Lower-48 state polygons in the grid CRS."""
    s = gpd.read_file(inputs.states_shapefile()).to_crs(GRID_CRS)
    return s[~s["STUSPS"].isin(NON_CONUS)].copy()


def union_mask(species):
    """Dissolved geometry of the union range mask (IUCN polygon OR iNat presence).

    Reconstructed from the per-cell comparison table and checked against the
    published union counts; all three species reproduce exactly.
    """
    d = pd.read_csv(os.path.join(inputs.REPO, f"{species}_mask_compare_cell100.csv"))
    keep = d[d["in_iucn"] | d["in_presence"]].copy()
    published = pd.read_csv(os.path.join(inputs.REPO, "mask_comparison_final.csv")) \
                  .set_index("species").loc[species, "c100_union"]
    if len(keep) != int(published):
        raise AssertionError(
            f"{species}: union mask rebuilt as {len(keep)} cell100 but the published "
            f"count is {int(published)}. Do not map an unverified mask."
        )
    h = CELL100 / 2
    sq = [gpd.points_from_xy([x], [y])[0].buffer(0) for x, y in zip(keep.x100, keep.y100)]
    from shapely.geometry import box
    g = gpd.GeoDataFrame(keep, geometry=[box(x - h, y - h, x + h, y + h)
                                         for x, y in zip(keep.x100, keep.y100)],
                         crs=GRID_CRS)
    return g.union_all()


def _cells_to_squares(df, size, crs_from=4326):
    from shapely.geometry import box
    g = gpd.GeoDataFrame(df.copy(), geometry=gpd.points_from_xy(df.lon, df.lat),
                         crs=crs_from).to_crs(GRID_CRS)
    h = size / 2
    g["geometry"] = [box(p.x - h, p.y - h, p.x + h, p.y + h) for p in g.geometry]
    return g


def fig_abundance(out="fig_abundance.png"):
    """Relative abundance surfaces at 2015 / 2020 / 2025, plus change.

    The mapped extent is the model's fitted footprint, which for these fits is
    the presence mask -- not the union mask. A union-mask surface cannot be drawn
    until a fit exists under that mask, because the extra cells have no model
    output. Each species has its own colour scale; absolute intensities are not
    comparable between species.
    """
    st = conus()
    years = [2015, 2020, 2025]
    # Dedicated colourbar columns: squeezing a colourbar between panels puts it
    # on top of the neighbouring panel's title.
    fig, axes = plt.subplots(3, 6, figsize=(12.2, 7.2),
                             gridspec_kw={"width_ratios": [1, 1, 1, 0.05, 1, 0.05]})
    for r in range(3):
        for c in (3, 5):
            axes[r, c].set_axis_off()
    for i, (key, label) in enumerate(SPECIES):
        mu = pd.read_csv(os.path.join(inputs.REPO, MU_FILE[key]))
        missing = sorted(set(years) - set(mu["year"].unique()))
        if missing:
            raise AssertionError(f"{key}: snapshot years {missing} absent from {MU_FILE[key]}")
        g = _cells_to_squares(mu[mu.year.isin(years)], CELL50)
        vmin = max(g["mu_mean"][g["mu_mean"] > 0].min(), 1e-6)
        vmax = g["mu_mean"].max()
        norm = LogNorm(vmin=vmin, vmax=vmax)
        for j, yr in enumerate(years):
            ax = axes[i, j]
            st.boundary.plot(ax=ax, color="0.78", linewidth=0.3, zorder=1)
            d = g[g.year == yr]
            d.plot(column="mu_mean", ax=ax, cmap="viridis", norm=norm,
                   linewidth=0, zorder=3)
            ax.set_axis_off()
            ax.set_xlim(*st.total_bounds[[0, 2]]); ax.set_ylim(*st.total_bounds[[1, 3]])
            if i == 0:
                ax.set_title(str(yr), fontsize=8.6)
            if j == 0:
                ax.text(-0.03, 0.5, label, transform=ax.transAxes, rotation=90,
                        ha="right", va="center", fontsize=8.4, fontweight="bold")
        # change panel
        w = mu.pivot(index="cell50", columns="year", values="mu_mean")
        pct = (100 * (w[2025] / w[2015] - 1)).rename("pct").reset_index()
        coords = mu[mu.year == 2025][["cell50", "lon", "lat"]]
        gc = _cells_to_squares(pct.merge(coords, on="cell50"), CELL50)
        lim = float(np.nanpercentile(gc["pct"].abs(), 99)) or 1.0
        ax = axes[i, 4]
        st.boundary.plot(ax=ax, color="0.78", linewidth=0.3, zorder=1)
        gc.plot(column="pct", ax=ax, cmap=style.DIVERGING,
                norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim), linewidth=0, zorder=3)
        ax.set_axis_off()
        ax.set_xlim(*st.total_bounds[[0, 2]]); ax.set_ylim(*st.total_bounds[[1, 3]])
        if i == 0:
            ax.set_title("Change 2015 to 2025", fontsize=8.6)
        med = float(np.nanmedian(gc["pct"]))
        ax.text(0.5, -0.02, f"median {med:+.1f}%", transform=ax.transAxes,
                ha="center", va="top", fontsize=6.8, color="#22282e")

        # A map without a scale is not readable as a quantity. One intensity
        # colourbar per row (the three years share a scale) and one for change.
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        cbi = fig.colorbar(sm, cax=fig.add_axes(_inset(axes[i, 3])))
        cbi.set_label("count intensity", fontsize=6.0)
        cbi.ax.tick_params(labelsize=5.2)
        smc = plt.cm.ScalarMappable(cmap=style.DIVERGING,
                                    norm=TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim))
        cbc = fig.colorbar(smc, cax=fig.add_axes(_inset(axes[i, 5])))
        cbc.set_label("% change", fontsize=6.0)
        cbc.ax.tick_params(labelsize=5.2)

    fig.suptitle("Relative abundance surfaces, 10 and 5 years ago and today", fontsize=9.4, y=0.99)
    fig.text(0.5, 0.035,
             "Yellow = higher modelled iNaturalist count intensity, dark blue = lower (log scale). "
             "Each row has its own scale, so compare years within a\nspecies, not between species. "
             "Extent is each model's fitted footprint (presence mask). Right column: per-cell "
             "percent change, red = decline.",
             ha="center", va="top", fontsize=6.6, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out


def _inset(ax, wfrac=0.42, hfrac=0.62):
    """Position for a slim colourbar centred in a reserved gridspec column."""
    b = ax.get_position()
    return [b.x0 + b.width * (1 - wfrac) / 2, b.y0 + b.height * (1 - hfrac) / 2,
            b.width * wfrac, b.height * hfrac]


def _ecoregions():
    e = gpd.read_file(inputs.ECOREGION_SHP).to_crs(GRID_CRS)
    e["region_key"] = e["NA_L1NAME"].str.upper()
    return e.dissolve(by="region_key").reset_index()[["region_key", "geometry"]]


def fig_regional_trend(window="full", out=None):
    """Regional trend by ecoregion, clipped to the lower 48 and the union range mask.

    Clipping matters: an unclipped ecoregion polygon colours ground the species
    does not occupy. One shared scale across species, since all three are in the
    same log-trend units.
    """
    out = out or f"fig_regional_trend_{window}.png"
    wlabel = inputs.WINDOWS[window][1]
    st = conus(); conus_geom = st.union_all()
    eco = _ecoregions()
    panels = []
    for key, label in SPECIES:
        tag, fname = ECO_FILE[key]
        d0 = inputs.PULL if tag == "pull" else inputs.REPO
        path = (os.path.join(d0, fname) if window == "full"
                else inputs.windowed(fname, window, [d0]))
        d = pd.read_csv(path)
        d["region_key"] = d["region_name"].str.upper()
        m = eco.merge(d, on="region_key", how="inner")
        m["geometry"] = m.geometry.intersection(conus_geom).intersection(union_mask(key))
        m = m[~m.geometry.is_empty]
        panels.append((label, m))

    vmax = max(float(m["abs_trend_mean"].abs().max()) for _, m in panels)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.2))
    for ax, (label, m) in zip(axes, panels):
        st.boundary.plot(ax=ax, color="0.80", linewidth=0.3, zorder=1)
        # A region whose interval overlaps zero gets no colour at all: colouring
        # it and then hatching over the top still reads as a direction.
        clear, unclear = _split_by_significance(m)
        if len(clear):
            clear.plot(column="abs_trend_mean", ax=ax, cmap=style.DIVERGING, norm=norm,
                       edgecolor="white", linewidth=0.4, zorder=3)
        if len(unclear):
            unclear.plot(ax=ax, facecolor=NO_TREND_GREY, edgecolor="white",
                         linewidth=0.4, zorder=3)
        ax.set_title(label, fontsize=8.6)
        ax.set_axis_off()
        ax.set_xlim(*st.total_bounds[[0, 2]]); ax.set_ylim(*st.total_bounds[[1, 3]])
    sm = plt.cm.ScalarMappable(cmap=style.DIVERGING, norm=norm)
    cb = fig.colorbar(sm, ax=axes.tolist(), fraction=0.020, pad=0.010, shrink=0.62,
                      anchor=(0.0, 0.72))
    cb.set_label(f"Regional trend, {wlabel}\n(log scale; green = increasing)", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    fig.text(0.5, -0.02,
             "Clipped to the lower 48 and to each species' union range mask (IUCN polygon or "
             "iNaturalist presence). Grey: the 95% interval overlaps zero, so no direction "
             "is claimed.",
             ha="center", fontsize=6.5, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, panels


def fig_agency_comparison(species, label, out=None):
    """Our regional trend beside the original state-agency map, unmodified.

    The agency panel is the previously produced raw state-by-state figure,
    inserted as-is. No re-aggregation to ecoregions, no derived direction score,
    and no agreement statistic: the two maps are at different spatial
    resolutions and are shown side by side for visual reference only.
    """
    out = out or f"fig_agency_{species}.png"
    raw = os.path.join(inputs.REPO, RAW_AGENCY_MAP[species])
    if not os.path.exists(raw):
        raise FileNotFoundError(f"original agency map missing: {raw}")

    st = conus(); conus_geom = st.union_all()
    eco = _ecoregions()
    tag, fname = ECO_FILE[species]
    d = pd.read_csv(os.path.join(inputs.PULL if tag == "pull" else inputs.REPO, fname))
    d["region_key"] = d["region_name"].str.upper()
    m = eco.merge(d, on="region_key", how="inner")
    m["geometry"] = m.geometry.intersection(conus_geom).intersection(union_mask(species))
    m = m[~m.geometry.is_empty]

    vmax = float(m["abs_trend_mean"].abs().max())
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    panel = f"_panel_model_{species}.png"
    fig, ax = plt.subplots(figsize=(5.4, 3.4))
    st.boundary.plot(ax=ax, color="0.80", linewidth=0.3, zorder=1)
    clear, unclear = _split_by_significance(m)
    if len(clear):
        clear.plot(column="abs_trend_mean", ax=ax, cmap=style.DIVERGING, norm=norm,
                   edgecolor="white", linewidth=0.4, zorder=3)
    if len(unclear):
        unclear.plot(ax=ax, facecolor=NO_TREND_GREY, edgecolor="white",
                     linewidth=0.4, zorder=3)
    ax.set_axis_off()
    ax.set_xlim(*st.total_bounds[[0, 2]]); ax.set_ylim(*st.total_bounds[[1, 3]])
    ax.set_title("Our model: trend by ecoregion", fontsize=8.8, loc="center")
    sm = plt.cm.ScalarMappable(cmap=style.DIVERGING, norm=norm)
    cb = fig.colorbar(sm, ax=ax, fraction=0.026, pad=0.01, shrink=0.78)
    cb.set_label("18-year trend\n(log scale; green = increasing)", fontsize=6.4)
    cb.ax.tick_params(labelsize=5.6)
    fig.text(0.5, 0.02, "Grey: 95% interval overlaps zero. Clipped to the lower 48\n"
                        "and the union range mask.",
             ha="center", va="top", fontsize=6.0, color="#5f6a73")
    fig.savefig(panel, bbox_inches="tight", pad_inches=0.05, dpi=300)
    plt.close(fig)

    # No added title band: the original agency figure carries its own title, and
    # a PIL-drawn caption renders at a fixed pixel size that is illegible next
    # to a 300-dpi panel.
    _compose_side_by_side(panel, raw, out)
    os.remove(panel)
    return out, m


def _compose_side_by_side(left_png, right_png, out, gap=40):
    """Place two rendered PNGs side by side at matched height.

    Used rather than re-plotting because the agency map must appear exactly as
    it was originally produced.
    """
    from PIL import Image
    li, ri = Image.open(left_png).convert("RGB"), Image.open(right_png).convert("RGB")
    h = max(li.height, ri.height)
    def fit(im):
        return im.resize((int(im.width * h / im.height), h), Image.LANCZOS)
    li, ri = fit(li), fit(ri)
    canvas = Image.new("RGB", (li.width + gap + ri.width, h), "white")
    canvas.paste(li, (0, 0))
    canvas.paste(ri, (li.width + gap, 0))
    canvas.save(out)
    return out

# Mask-construction categories, drawn in distinct colours so the figure explains
# how the modelled area was defined rather than arguing for a choice.
MASK_CAT = [("both", "In the IUCN range map\nand has iNaturalist records", "#2c6c9c"),
            ("iucn_only", "IUCN range map only\n(no records: a true absence)", "#a6cee3"),
            ("inat_only", "Added by iNaturalist records\n(outside the IUCN map)", "#e08214")]


def fig_mask_construction(out="fig_mask_construction.png"):
    """How the modelled area was defined, per species.

    Each 100 km cell is shown by which source placed it in range: the IUCN range
    map, the iNaturalist records, or both. This is a methods figure -- it shows
    what area the models cover and where that area came from.
    """
    from shapely.geometry import box
    st = conus()
    h = CELL100 / 2
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.1))
    counts = {}
    for ax, (key, label) in zip(axes, SPECIES):
        d = pd.read_csv(os.path.join(inputs.REPO, f"{key}_mask_compare_cell100.csv"))
        d = d[d["in_iucn"] | d["in_presence"]].copy()
        d["cat"] = np.where(d["in_iucn"] & d["in_presence"], "both",
                    np.where(d["in_iucn"], "iucn_only", "inat_only"))
        g = gpd.GeoDataFrame(d, geometry=[box(x - h, y - h, x + h, y + h)
                                          for x, y in zip(d.x100, d.y100)], crs=GRID_CRS)
        st.boundary.plot(ax=ax, color="0.72", linewidth=0.35, zorder=4)
        for cat, _, col in MASK_CAT:
            sub = g[g.cat == cat]
            if len(sub):
                sub.plot(ax=ax, facecolor=col, edgecolor="none", zorder=2)
        counts[label] = {c: int((d.cat == c).sum()) for c, _, _ in MASK_CAT}
        ax.set_title(f"{label}  ({len(d)} cells)", fontsize=8.6)
        ax.set_axis_off()
        ax.set_xlim(*st.total_bounds[[0, 2]]); ax.set_ylim(*st.total_bounds[[1, 3]])

    handles = [Patch(facecolor=col, edgecolor="none", label=lab) for _, lab, col in MASK_CAT]
    fig.legend(handles=handles, loc="lower center", ncol=3, frameon=False,
               fontsize=6.8, bbox_to_anchor=(0.5, -0.16))
    fig.suptitle("How the modelled area was defined: 100 km cells by source", fontsize=9.4, y=1.02)
    fig.text(0.5, -0.26,
             "A cell enters the model if either source places it in range. Orange cells are real "
             "iNaturalist observations that fall outside the IUCN\nmap; pale blue cells are inside "
             "the map with no records, and contribute genuine absence information.",
             ha="center", va="top", fontsize=6.4, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, pd.DataFrame(counts).T
