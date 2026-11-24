"""Fleet figures built from the extraction products: surfaces, fields, covariates.

A note on what `mu_snapshots` contains, because the name invites a wrong reading
and an older module in this package acts on that wrong reading. `mu` is the
expected number of iNaturalist records in a cell-year, not abundance. Its
temporal slope is theta1 * total_var_beta, confirmed empirically across eleven
species: that reading predicts the observed per-cell change to a mean error of
0.97 percentage points per year, against 6.42 for the reading that treats `mu`
as the abundance surface, and the two species whose covariates barely move match
it to two decimal places.

The shared relative-abundance surface is recovered by dividing log mu by theta1.
The additive constant that goes with it is not identified from the product files
alone, so these surfaces are relative within a species and carry no absolute
scale; the pattern and the change are unaffected, since a constant in log space
cancels from both.
"""
import glob
import os
import re

import geopandas as gpd
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LogNorm, TwoSlopeNorm

from .figures_deployment import DARK_BG, DARK_INK, DARK_RULE, DARK_SUB
from .figures_maps import conus

PULL = ("/Users/rwkays/claude_code/data_integration_arielle/"
        "hazel_pull_20261115_win5/posteriors")
PRODUCTS = os.path.join(PULL, "products")
CELL50 = 50_000
RUN_RE = (r"(.+?)_(v2u|v2un)_national_scalar_(win5|win10)_ctr(_excl)?_"
          r"(mu_snapshots|spatial_car_fields|occbeta|psi_snapshots)\.csv$")


def _pretty(s):
    return s.replace("_", " ").capitalize()


def inventory(products=PRODUCTS, window="win5"):
    """Map species -> {product type -> path} for one window."""
    inv = {}
    for f in sorted(glob.glob(os.path.join(products, "*.csv"))):
        m = re.match(RUN_RE, os.path.basename(f))
        if m and m.group(3) == window:
            inv.setdefault(m.group(1), {})[m.group(5)] = f
    return inv


def _theta1(species, pull=PULL):
    hits = glob.glob(os.path.join(pull, f"{species}*_posterior_win5.csv"))
    if not hits:
        raise AssertionError(f"{species}: no posterior file in {pull}")
    po = pd.read_csv(hits[0])
    r = po[po.param.astype(str).str.fullmatch("theta1")]
    if not len(r):
        raise AssertionError(f"{species}: theta1 absent, cannot recover abundance")
    return float(r["mean"].iloc[0])


def abundance_surface(species, products=PRODUCTS, pull=PULL):
    """Relative abundance per cell-year, recovered from the iNaturalist surface.

    Returns long-form rows with `rel_abundance`, normalised so the species' own
    maximum over the window is 1. Absolute values carry no meaning; only the
    pattern across cells and the change across years do.
    """
    inv = inventory(products)
    mu = pd.read_csv(inv[species]["mu_snapshots"])
    th1 = _theta1(species, pull)
    if abs(th1) < 0.05:
        raise AssertionError(
            f"{species}: theta1 = {th1:.3f} is too near zero to invert; the "
            "iNaturalist surface carries no abundance signal to recover")
    mu = mu[mu.mu_mean > 0].copy()
    mu["log_rel"] = np.log(mu.mu_mean) / th1
    # Anchor every year to the reported camera trend, not only the change panel.
    # Left as fitted, consecutive year surfaces step at the iNaturalist rate, so
    # a 2021-against-2025 pair would contradict the change panel beside it.
    po = pd.read_csv(glob.glob(os.path.join(pull, f"{species}*_posterior_win5.csv"))[0])
    par = lambda p: float(po.loc[po.param.astype(str).str.fullmatch(p), "mean"].iloc[0])
    ys = pd.read_csv(glob.glob(os.path.join(pull, f"year_scaling_{species}*.csv"))[0])
    sdc = float(ys.loc[ys.coefficient == "year_beta", "year_sd"].iloc[0])
    sdi = float(ys.loc[ys.coefficient == "total_var_beta", "year_sd"].iloc[0])
    y0 = mu.year.min()
    step = mu.year - y0
    mu["log_rel_anchored"] = (mu.log_rel
                              - par("total_var_beta") * (step / sdi)
                              + par("year_beta") * (step / sdc))
    for col in ("log_rel", "log_rel_anchored"):
        mu[col] -= mu["log_rel_anchored"].max()
    mu["rel_abundance"] = np.exp(mu.log_rel_anchored)
    mu["species"] = _pretty(species)
    mu["theta1"] = th1
    return mu


def _squares(df, size=CELL50, crs_to=None):
    from shapely.geometry import box
    g = gpd.GeoDataFrame(df.copy(), geometry=gpd.points_from_xy(df.lon, df.lat),
                         crs=4326).to_crs(crs_to)
    h = size / 2
    g["geometry"] = [box(p.x - h, p.y - h, p.x + h, p.y + h) for p in g.geometry]
    return g


def reported_trend(pull=PULL):
    """Per species, the reported camera trend in % per year and whether it resolves."""
    out = {}
    for f in sorted(glob.glob(os.path.join(pull, "*_posterior_win5.csv"))):
        d = pd.read_csv(f)
        b = d[d.param.astype(str).str.fullmatch("year_beta")]
        if not len(b):
            continue
        run = d.run.iloc[0]
        ys = pd.read_csv(os.path.join(pull, f"year_scaling_{run}.csv"))
        sd = float(ys.loc[ys.coefficient == "year_beta", "year_sd"].iloc[0])
        pct = lambda v: 100 * (np.exp(float(v) / sd) - 1)
        lo, hi = pct(b["q2.5"].iloc[0]), pct(b["q97.5"].iloc[0])
        out[_pretty(d.species.iloc[0])] = {
            "pct_yr": pct(b["mean"].iloc[0]), "lo": lo, "hi": hi,
            "resolved": bool(lo > 0 or hi < 0)}
    return out


def abundance_change(products=PRODUCTS, pull=PULL):
    """Per cell, the annual rate of change in recovered relative abundance.

    Two versions, and the difference is the whole point.

    The surface's SPATIAL structure is shared: `link_occ_intercept`, both
    climate fields and all nine `occ_beta` covariate effects appear in the
    camera occupancy block and the iNaturalist intensity block alike, so both
    datasets and the covariates inform where the animals are. Only the year
    coefficient differs between the two blocks -- the camera submodel carries
    `year_beta`, the intensity carries `total_var_beta = year_beta + year_var`.

    So the surface as fitted changes at the iNaturalist-side rate, which is not
    the trend we report. Year enters the log-intensity linearly and identically
    in every cell, so swapping one coefficient for the other is a constant
    shift: it leaves the covariate-driven spatial variation untouched and moves
    the whole surface onto the reported trend. `pct_per_year` is that anchored
    version; `as_fitted_pct_yr` is kept so the two can be compared.
    """
    rows = []
    for sp in sorted(inventory(products)):
        a = abundance_surface(sp, products, pull)
        th1 = float(a.theta1.iloc[0])
        po = pd.read_csv(glob.glob(os.path.join(pull, f"{sp}*_posterior_win5.csv"))[0])
        par = lambda p: float(po.loc[po.param.astype(str).str.fullmatch(p), "mean"].iloc[0])
        ys = pd.read_csv(glob.glob(os.path.join(pull, f"year_scaling_{sp}*.csv"))[0])
        sdc = float(ys.loc[ys.coefficient == "year_beta", "year_sd"].iloc[0])
        sdi = float(ys.loc[ys.coefficient == "total_var_beta", "year_sd"].iloc[0])
        w = a.pivot_table(index="cell50", columns="year", values="log_rel")
        yrs = sorted(w.columns)
        span = yrs[-1] - yrs[0]
        dlog = w[yrs[-1]] - w[yrs[0]]
        dlog_anchored = dlog - par("total_var_beta") * (span / sdi) \
                             + par("year_beta") * (span / sdc)
        xy = a.drop_duplicates("cell50").set_index("cell50")[["lon", "lat"]]
        g = pd.DataFrame({
            "species": _pretty(sp), "cell50": w.index, "theta1": th1,
            "log_rel_first": w[yrs[0]].to_numpy(), "log_rel_last": w[yrs[-1]].to_numpy(),
            "as_fitted_pct_yr": 100 * (np.exp(dlog / span) - 1).to_numpy(),
            "pct_per_year": 100 * (np.exp(dlog_anchored / span) - 1).to_numpy(),
            "year_first": yrs[0], "year_last": yrs[-1],
        }).join(xy, on="cell50")
        rows.append(g.dropna(subset=["lon", "lat"]))
    return pd.concat(rows, ignore_index=True)


PULL10 = "/Users/rwkays/claude_code/data_integration_arielle/ctr_fleet_pull_20261115"
FULL_RECORD = {"Bobcat": (1.02, -0.52, 2.70),
               "White-tailed deer": (4.00, 3.04, 4.97)}


def window_trends(species=("Bobcat", "White-tailed deer"), pull=PULL, pull10=PULL10):
    """Trend per window for the species that have more than one.

    Reads whichever windows exist rather than assuming a set, so the five-year
    row fills itself in when that fit lands instead of needing a code change.
    An unavailable window is returned with NaN estimates and `available=False`
    so the figure can say it is missing rather than silently omitting it.
    """
    out = []
    for sp in species:
        if sp in FULL_RECORD:
            m, lo, hi = FULL_RECORD[sp]
            out.append({"species": sp, "window": "full record\n2008-2025",
                        "trend": m, "lo": lo, "hi": hi, "available": True})
        key = sp.lower().replace(" ", "_")
        f10 = os.path.join(pull10, f"{key}_v2u_national_scalar_posterior_win10.csv")
        for lab, path, sfx in [("10-year\n2016-2025", f10, "win10_ctr"),
                               ("5-year\n2021-2025",
                                os.path.join(pull, f"{key}_v2u_national_scalar_"
                                                   "posterior_win5.csv"), "win5_ctr")]:
            ysf = glob.glob(os.path.join(os.path.dirname(path),
                                         f"year_scaling_{key}*{sfx}.csv"))
            if not (os.path.exists(path) and ysf):
                out.append({"species": sp, "window": lab, "trend": np.nan,
                            "lo": np.nan, "hi": np.nan, "available": False})
                continue
            po = pd.read_csv(path)
            sd = float(pd.read_csv(ysf[0]).query("coefficient == 'year_beta'")
                       ["year_sd"].iloc[0])
            r = po[po.param.astype(str).str.fullmatch("year_beta")].iloc[0]
            pct = lambda v: 100 * (np.exp(float(v) / sd) - 1)
            out.append({"species": sp, "window": lab, "trend": pct(r["mean"]),
                        "lo": pct(r["q2.5"]), "hi": pct(r["q97.5"]), "available": True})
    d = pd.DataFrame(out)
    d["resolved"] = (d.lo > 0) | (d.hi < 0)
    return d


def fig_window_trends(out="fig_window_trends.png", table=None, **kw):
    """Trend by window, for the species fitted on more than one."""
    d = table if table is not None else window_trends(**kw)
    sp = list(dict.fromkeys(d.species))
    fig, axes = plt.subplots(1, len(sp), figsize=(4.6 * len(sp), 3.4),
                             facecolor=DARK_BG, squeeze=False)
    for k, s_ in enumerate(sp):
        ax = axes[0, k]
        ax.set_facecolor(DARK_BG)
        g = d[d.species == s_].reset_index(drop=True)
        y = np.arange(len(g))[::-1]
        ax.axvline(0, color=DARK_RULE, lw=1.0, zorder=1)
        for yi, r in zip(y, g.itertuples()):
            if not r.available:
                ax.text(0, yi, "  not yet fitted", color="#8d949b", fontsize=9.5,
                        va="center", ha="left", style="italic")
                continue
            c = "#4da3e0" if r.resolved else "#9aa4ad"
            ax.plot([r.lo, r.hi], [yi, yi], color=c, lw=3.0, solid_capstyle="round")
            ax.plot(r.trend, yi, "o", color=c, ms=7)
            ax.text(r.hi, yi + 0.26, f"{r.trend:+.2f}", color=c, fontsize=9.5,
                    va="bottom", ha="right")
        ax.set_yticks(y); ax.set_yticklabels(g.window, fontsize=10, color=DARK_INK)
        ax.set_title(s_, fontsize=13.5, color=DARK_INK, pad=8)
        ax.set_xlabel("change in abundance (% per year)", fontsize=10, color=DARK_SUB)
        ax.tick_params(colors=DARK_INK, labelsize=9.5)
        for sp_ in ax.spines.values():
            sp_.set_color(DARK_RULE)
        ax.set_ylim(-0.7, len(g) - 0.3)
    fig.text(0.5, -0.10,
             "The same data over three windows. Filled blue marks an interval excluding zero; "
             "grey marks one that spans it. The five-year window is the reported one.",
             ha="center", va="top", fontsize=10, color=DARK_SUB)
    fig.savefig(out, facecolor=DARK_BG, bbox_inches="tight", pad_inches=0.22, dpi=160)
    plt.close(fig)
    return out, d


def win10_surface(species, products=PRODUCTS, pull10=PULL10):
    """Anchored relative-abundance surface per year for a ten-year fit.

    Same construction as the five-year version, but the ten-year posteriors sit
    in a different pull directory and their files carry a `_win10` suffix, so
    the paths cannot be shared. Only the centred ten-year runs are used; the
    older uncentred ten-year fits in earlier pulls are superseded and must not
    be mixed into a deliverable.
    """
    f = [x for x in glob.glob(os.path.join(products, "*win10*mu_snapshots.csv"))
         if os.path.basename(x).startswith(species + "_")]
    if not f:
        raise AssertionError(f"{species}: no ten-year mu_snapshots in {products}")
    mu = pd.read_csv(f[0])
    po = pd.read_csv(os.path.join(pull10, f"{species}_v2u_national_scalar_posterior_win10.csv"))
    ys = pd.read_csv(os.path.join(
        pull10, f"year_scaling_{species}_v2u_national_scalar_win10_ctr.csv"))
    par = lambda p: float(po.loc[po.param.astype(str).str.fullmatch(p), "mean"].iloc[0])
    sdc = float(ys.loc[ys.coefficient == "year_beta", "year_sd"].iloc[0])
    sdi = float(ys.loc[ys.coefficient == "total_var_beta", "year_sd"].iloc[0])
    th1 = par("theta1")
    mu = mu[mu.mu_mean > 0].copy()
    step = mu.year - mu.year.min()
    mu["log_rel"] = (np.log(mu.mu_mean) / th1
                     - par("total_var_beta") * (step / sdi)
                     + par("year_beta") * (step / sdc))
    mu["log_rel"] -= mu.log_rel.max()
    mu["species"] = _pretty(species)
    mu["trend_pct_yr"] = 100 * (np.exp(par("year_beta") / sdc) - 1)
    lo = 100 * (np.exp(float(po.loc[po.param.astype(str).str.fullmatch("year_beta"),
                                    "q2.5"].iloc[0]) / sdc) - 1)
    hi = 100 * (np.exp(float(po.loc[po.param.astype(str).str.fullmatch("year_beta"),
                                    "q97.5"].iloc[0]) / sdc) - 1)
    mu["lo"], mu["hi"] = lo, hi
    return mu


def fig_win10_windows(out="fig_win10_windows.png", species=("bobcat", "white-tailed_deer"),
                      years=(2016, 2020, 2025), products=PRODUCTS, pull10=PULL10):
    """Three time slices from the ten-year fit, plus the trend, for the two species that have one.

    Shown as an example of what a longer window buys, not as a fleet result:
    only these two species have a centred ten-year fit.
    """
    st = conus()
    fig, axes = plt.subplots(len(species), len(years) + 1,
                             figsize=(5.0 * (len(years) + 1), 5.0 * len(species)),
                             facecolor=DARK_BG, squeeze=False)
    fig.subplots_adjust(hspace=0.30, wspace=0.02)
    for i, sp in enumerate(species):
        a = win10_surface(sp, products, pull10)
        w = a.pivot_table(index="cell50", columns="year", values="log_rel")
        xy = a.drop_duplicates("cell50").set_index("cell50")[["lon", "lat"]]
        miss = [y for y in years if y not in w.columns]
        if miss:
            raise AssertionError(f"{sp}: ten-year fit has no year {miss}")
        base = pd.DataFrame({"cell50": w.index}).join(xy, on="cell50").dropna()
        lo = max(float(np.exp(w[list(years)].to_numpy()).min()), 1e-6)
        norm = LogNorm(vmin=lo, vmax=1.0)
        for j, y in enumerate(years):
            d = base.copy(); d["rel"] = np.exp(w[y]).reindex(d.cell50).to_numpy()
            g = _squares(d, crs_to=st.crs)
            ax = axes[i, j]; ax.set_axis_off(); ax.set_facecolor(DARK_BG)
            st.boundary.plot(ax=ax, color=DARK_RULE, linewidth=0.5, zorder=1)
            g.plot(column="rel", ax=ax, cmap="viridis", norm=norm, linewidth=0, zorder=3)
            ax.set_title((f"{_pretty(sp)}\n{y}" if j == 0 else f"\n{y}"),
                         fontsize=14, color=DARK_INK, pad=8)
        yrs = sorted(w.columns); span = yrs[-1] - yrs[0]
        d = base.copy()
        d["pct"] = (100 * (np.exp((w[yrs[-1]] - w[yrs[0]]) / span) - 1)
                    .reindex(d.cell50).to_numpy())
        g = _squares(d, crs_to=st.crs)
        ax = axes[i, len(years)]; ax.set_axis_off(); ax.set_facecolor(DARK_BG)
        st.boundary.plot(ax=ax, color=DARK_RULE, linewidth=0.5, zorder=1)
        v = float(np.abs(d.pct).max()) or 1.0
        g.plot(column="pct", ax=ax, cmap="RdBu_r",
               norm=TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v), linewidth=0, zorder=3)
        t, l_, h_ = a.trend_pct_yr.iloc[0], a.lo.iloc[0], a.hi.iloc[0]
        resolved = (l_ > 0) or (h_ < 0)
        ax.set_title(f"\ntrend, {yrs[0]} to {yrs[-1]}\n{t:+.2f}%/yr  [{l_:+.2f}, {h_:+.2f}]"
                     + ("" if resolved else "   spans zero"),
                     fontsize=12.5, color=DARK_INK if resolved else "#e0a060", pad=8)
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        cb = fig.colorbar(sm, ax=[axes[i, j] for j in range(len(years))],
                          orientation="horizontal", fraction=0.038, pad=0.02, shrink=0.34)
        cb.set_label("relative abundance, shared scale across all three years",
                     fontsize=10.5, color=DARK_SUB, labelpad=3)
        cb.ax.tick_params(labelsize=9.5, colors=DARK_INK)
        cb.outline.set_edgecolor(DARK_RULE)
        sm2 = plt.cm.ScalarMappable(cmap="RdBu_r",
                                    norm=TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v))
        cb2 = fig.colorbar(sm2, ax=ax, orientation="horizontal", fraction=0.038,
                           pad=0.02, shrink=0.70)
        cb2.set_label("change in abundance (% per year)", fontsize=10.5,
                      color=DARK_SUB, labelpad=3)
        cb2.ax.tick_params(labelsize=9.5, colors=DARK_INK)
        cb2.outline.set_edgecolor(DARK_RULE)
        xmin, ymin, xmax, ymax = g.total_bounds
        side = max(xmax - xmin, ymax - ymin) * 1.12
        cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
        for j in range(len(years) + 1):
            axes[i, j].set_xlim(cx - side / 2, cx + side / 2)
            axes[i, j].set_ylim(cy - side / 2, cy + side / 2)
    fig.suptitle("A longer window: ten years of relative abundance, for the two species "
                 "that have one", fontsize=18, color=DARK_INK, y=0.997)
    fig.text(0.5, 0.030,
             "Only bobcat and white-tailed deer have a centred ten-year fit, so these are "
             "examples of what a longer window shows rather than a fleet result. All panels "
             "share one colour\nscale within a species and are anchored to that species' "
             "reported camera trend, so the difference between the year panels is the trend. "
             "The reported window remains the\nfive-year one; these ten-year estimates are "
             "shown for comparison, and both span zero.",
             ha="center", va="top", fontsize=10.5, color=DARK_SUB, linespacing=1.7)
    fig.savefig(out, facecolor=DARK_BG, bbox_inches="tight", pad_inches=0.26, dpi=150)
    plt.close(fig)
    return out


def surface_panels(products=PRODUCTS, pull=PULL):
    """First-year and last-year anchored surfaces, per cell, for every species."""
    rows = []
    for sp in sorted(inventory(products)):
        a = abundance_surface(sp, products, pull)
        w = a.pivot_table(index="cell50", columns="year", values="log_rel_anchored")
        yrs = sorted(w.columns)
        xy = a.drop_duplicates("cell50").set_index("cell50")[["lon", "lat"]]
        rows.append(pd.DataFrame({
            "species": _pretty(sp), "cell50": w.index,
            "year_first": yrs[0], "year_last": yrs[-1],
            "rel_first": np.exp(w[yrs[0]]).to_numpy(),
            "rel_last": np.exp(w[yrs[-1]]).to_numpy(),
        }).join(xy, on="cell50").dropna(subset=["lon", "lat"]))
    return pd.concat(rows, ignore_index=True)


def fig_abundance_three(out="fig_abundance_three.png", species=None, page=None,
                        n_pages=None, products=PRODUCTS, pull=PULL, change=None,
                        trends=None, panels=None, slide=False):
    """Two year surfaces and the trend map, three panels per species.

    The two year panels share one colour scale, so the level difference between
    them IS the trend and can be read directly rather than inferred. Both are
    anchored to the reported camera trend, so they agree with the third panel.

    `slide=True` renders the projection variant: the explanatory caption is
    dropped and every font scales up. A print caption set at 10.5 pt is
    unreadable projected, and on a slide the speaker says what it says.
    """
    k = 1.75 if slide else 1.0
    pn = panels if panels is not None else surface_panels(products, pull)
    ch = change if change is not None else abundance_change(products, pull)
    tr = trends if trends is not None else reported_trend(pull)
    sp = species or sorted(pn.species.unique())
    pn, ch = pn[pn.species.isin(sp)], ch[ch.species.isin(sp)]
    st = conus()
    fig, axes = plt.subplots(len(sp), 3,
                             figsize=(5.6 * 3, (4.1 if slide else 5.0) * len(sp)),
                             facecolor=DARK_BG, squeeze=False)
    fig.subplots_adjust(hspace=0.34, wspace=0.02)
    for i, s_ in enumerate(sp):
        d, dc = pn[pn.species == s_], ch[ch.species == s_]
        g = _squares(d, crs_to=st.crs)
        gc = _squares(dc, crs_to=st.crs)
        rep = tr.get(s_, {})
        y0, y1 = int(d.year_first.iloc[0]), int(d.year_last.iloc[0])
        lo = max(min(d.rel_first.min(), d.rel_last.min()), 1e-6)
        norm = LogNorm(vmin=lo, vmax=1.0)
        for j, (col, lab) in enumerate([("rel_first", y0), ("rel_last", y1)]):
            ax = axes[i, j]; ax.set_axis_off(); ax.set_facecolor(DARK_BG)
            st.boundary.plot(ax=ax, color=DARK_RULE, linewidth=0.5, zorder=1)
            g.plot(column=col, ax=ax, cmap="viridis", norm=norm, linewidth=0, zorder=3)
            ax.set_title(f"{s_ if j == 0 else ''}\n{lab}" if j == 0 else f"\n{lab}",
                         fontsize=14 * k, color=DARK_INK, pad=8)
        ax = axes[i, 2]; ax.set_axis_off(); ax.set_facecolor(DARK_BG)
        st.boundary.plot(ax=ax, color=DARK_RULE, linewidth=0.5, zorder=1)
        v = float(np.abs(dc.pct_per_year).max()) or 1.0
        gc.plot(column="pct_per_year", ax=ax, cmap="RdBu_r",
                norm=TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v), linewidth=0, zorder=3)
        ax.set_title(f"\ntrend, {y0} to {y1}\nmedian {dc.pct_per_year.median():+.1f}%/yr "
                     f"· reported {rep.get('pct_yr', float('nan')):+.1f}%/yr"
                     + ("" if rep.get("resolved") else ", spans zero"),
                     fontsize=12.5 * k, color=DARK_INK, pad=8)
        sm = plt.cm.ScalarMappable(cmap="viridis", norm=norm)
        cb = fig.colorbar(sm, ax=[axes[i, 0], axes[i, 1]], orientation="horizontal",
                          fraction=0.038, pad=0.02, shrink=0.44)
        cb.set_label("relative abundance, shared scale across both years",
                     fontsize=10.5 * k, color=DARK_SUB, labelpad=3)
        cb.ax.tick_params(labelsize=9.5 * k, colors=DARK_INK)
        cb.outline.set_edgecolor(DARK_RULE)
        sm2 = plt.cm.ScalarMappable(cmap="RdBu_r",
                                    norm=TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v))
        cb2 = fig.colorbar(sm2, ax=axes[i, 2], orientation="horizontal",
                           fraction=0.038, pad=0.02, shrink=0.70)
        cb2.set_label("change in abundance (% per year)", fontsize=10.5 * k,
                      color=DARK_SUB, labelpad=3)
        cb2.ax.tick_params(labelsize=9.5 * k, colors=DARK_INK)
        cb2.outline.set_edgecolor(DARK_RULE)
        xmin, ymin, xmax, ymax = g.total_bounds
        side = max(xmax - xmin, ymax - ymin) * 1.16
        cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
        for j in range(3):
            axes[i, j].set_xlim(cx - side / 2, cx + side / 2)
            axes[i, j].set_ylim(cy - side / 2, cy + side / 2)
    tag = f"  ({page} of {n_pages})" if page and n_pages else ""
    # No suptitle on a slide: the slide carries its own title, and at projection
    # type size this one lands on top of the panel titles.
    if not slide:
        fig.suptitle("Relative abundance at the start and end of the window, and the "
                     "trend" + tag, fontsize=18, color=DARK_INK, y=0.997)
    if not slide:
        fig.text(0.5, 0.030,
                 "The shared surface, which both the camera and the iNaturalist likelihood "
                 "inform along with nine covariates and two climate fields. The two year panels "
                 "share one colour scale,\nso the difference between them is the trend and can be "
                 "read directly. All three panels are anchored to the reported camera trend: the "
                 "surface as fitted steps at the\niNaturalist rate, which is a different quantity. "
                 "Values are relative within a species, so colour is not comparable between "
                 "species and no absolute density is implied.\nThe window begins in 2021, so there "
                 "is no 2020 surface on these fits.",
                 ha="center", va="top", fontsize=10.5, color=DARK_SUB, linespacing=1.7)
    fig.savefig(out, facecolor=DARK_BG, bbox_inches="tight", pad_inches=0.26, dpi=150)
    plt.close(fig)
    return out, pn


def fig_abundance_three_pages(stem="fig_abundance_three", per_page=3,
                              products=PRODUCTS, pull=PULL):
    pn = surface_panels(products, pull)
    ch = abundance_change(products, pull)
    tr = reported_trend(pull)
    sp = sorted(pn.species.unique())
    pages = [sp[i:i + per_page] for i in range(0, len(sp), per_page)]
    return [fig_abundance_three(out=f"{stem}_p{i}.png", species=g, page=i,
                                n_pages=len(pages), panels=pn, change=ch, trends=tr)[0]
            for i, g in enumerate(pages, 1)], pn


def fig_abundance_surface(out="fig_abundance_surface.png", species=None, page=None,
                          n_pages=None, products=PRODUCTS, pull=PULL, change=None,
                          trends=None):
    """Relative abundance and its rate of change, two panels per species.

    One surface rather than one per year: across the window every species keeps
    its shape to a Spearman correlation of at least 0.997, so a five-year
    sequence would be the same picture five times. The change panel carries what
    does differ between years.
    """
    ch = change if change is not None else abundance_change(products, pull)
    tr = trends if trends is not None else reported_trend(pull)
    sp = species or sorted(ch.species.unique())
    ch = ch[ch.species.isin(sp)]
    st = conus()
    fig, axes = plt.subplots(len(sp), 2, figsize=(6.6 * 2, 5.0 * len(sp)),
                             facecolor=DARK_BG, squeeze=False)
    fig.subplots_adjust(hspace=0.34, wspace=0.02)
    for i, s_ in enumerate(sp):
        d = ch[ch.species == s_]
        gsq = _squares(d, crs_to=st.crs)
        rep = tr.get(s_, {})
        for j in (0, 1):
            ax = axes[i, j]
            ax.set_axis_off(); ax.set_facecolor(DARK_BG)
            st.boundary.plot(ax=ax, color=DARK_RULE, linewidth=0.5, zorder=1)
        # left: the surface itself, on a log scale since abundance spans orders
        lo = max(np.exp(d.log_rel_last).min(), 1e-6)
        gsq["rel"] = np.exp(gsq.log_rel_last)
        gsq.plot(column="rel", ax=axes[i, 0], cmap="viridis",
                 norm=LogNorm(vmin=lo, vmax=1.0), linewidth=0, zorder=3)
        axes[i, 0].set_title(f"{s_}\nrelative abundance, {int(d.year_last.iloc[0])} "
                             f"· {len(d):,} grid cells",
                             fontsize=13, color=DARK_INK, pad=8)
        # right: the rate of change, centred on zero
        v = float(np.abs(d.pct_per_year).max()) or 1.0
        gsq.plot(column="pct_per_year", ax=axes[i, 1], cmap="RdBu_r",
                 norm=TwoSlopeNorm(vmin=-v, vcenter=0.0, vmax=v), linewidth=0, zorder=3)
        disagree = (rep and np.sign(d.pct_per_year.median()) != np.sign(rep["pct_yr"]))
        flag = "   differs in sign from the reported trend" if disagree else ""
        axes[i, 1].set_title(
            f"change in abundance{flag}\n"
            f"median {d.pct_per_year.median():+.1f}%/yr · camera trend "
            f"{rep.get('pct_yr', float('nan')):+.1f}%/yr"
            + ("" if rep.get("resolved") else ", spans zero"),
            fontsize=13, color=DARK_INK if not disagree else "#e0a060", pad=8)
        for j, (cmap, vlo, vhi, lab) in enumerate([
                ("viridis", lo, 1.0, "relative abundance (species maximum = 1)"),
                ("RdBu_r", -v, v, "change in abundance (% per year)")]):
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=(
                LogNorm(vmin=vlo, vmax=vhi) if j == 0
                else TwoSlopeNorm(vmin=vlo, vcenter=0.0, vmax=vhi)))
            cb = fig.colorbar(sm, ax=axes[i, j], orientation="horizontal",
                              fraction=0.040, pad=0.02, shrink=0.62)
            cb.set_label(lab, fontsize=10.5, color=DARK_SUB, labelpad=3)
            cb.ax.tick_params(labelsize=9.5, colors=DARK_INK)
            cb.outline.set_edgecolor(DARK_RULE)
        for j in (0, 1):
            xmin, ymin, xmax, ymax = gsq.total_bounds
            side = max(xmax - xmin, ymax - ymin) * 1.16
            cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
            axes[i, j].set_xlim(cx - side / 2, cx + side / 2)
            axes[i, j].set_ylim(cy - side / 2, cy + side / 2)
    tag = f"  ({page} of {n_pages})" if page and n_pages else ""
    fig.suptitle("Where each species is common, and how that is changing" + tag,
                 fontsize=18, color=DARK_INK, y=0.997)
    fig.text(0.5, 0.030,
             "Relative abundance on the model's shared surface, which both the camera and the "
             "iNaturalist likelihood inform, read out through the iNaturalist link by dividing "
             "the fitted log intensity by theta1. Values are "
             "relative within a\nspecies, so colour is not comparable between panels and no "
             "absolute density is implied. One year is shown because every species holds its "
             "spatial shape across the window\nto a rank correlation of at least 0.997, which "
             "makes a year-by-year sequence the same picture repeated.\n"
             "The change panel is anchored to the reported camera trend. The surface's "
             "spatial structure is shared by both likelihoods -- the same spatial intercept "
             "field, climate fields and nine\ncovariate effects appear in the camera and "
             "iNaturalist blocks alike, so both datasets and the covariates inform where the "
             "animals are. Only the year coefficient differs between\nthe two blocks, and "
             "swapping the iNaturalist one for the camera one shifts the whole surface "
             "without touching its spatial variation. Each title gives the reported trend "
             "beside\nthe surface median; they differ by covariate-driven variation between "
             "cells.",
             ha="center", va="top", fontsize=10.5, color=DARK_SUB, linespacing=1.7)
    fig.savefig(out, facecolor=DARK_BG, bbox_inches="tight", pad_inches=0.26, dpi=150)
    plt.close(fig)
    return out, ch


def fig_abundance_pages(stem="fig_abundance_surface", per_page=3, products=PRODUCTS,
                        pull=PULL):
    ch = abundance_change(products, pull)
    tr = reported_trend(pull)
    sp = sorted(ch.species.unique())
    pages = [sp[i:i + per_page] for i in range(0, len(sp), per_page)]
    return [fig_abundance_surface(out=f"{stem}_p{i}.png", species=g, page=i,
                                  n_pages=len(pages), change=ch, trends=tr)[0]
            for i, g in enumerate(pages, 1)], ch


def surface_similarity(products=PRODUCTS, pull=PULL):
    """How much the abundance surface's shape changes across the window.

    The trend is one scalar, so a year-to-year sequence of surfaces risks being
    the same picture rescaled. Annual covariates are the only thing that can
    change the shape. This reports the rank correlation between the first and
    last year, so a figure showing every year can be justified or dropped.
    """
    rows = []
    for sp in sorted(inventory(products)):
        a = abundance_surface(sp, products, pull)
        w = a.pivot_table(index="cell50", columns="year", values="log_rel")
        yrs = sorted(w.columns)
        rows.append({"species": _pretty(sp), "years": len(yrs),
                     "spearman_first_last": w[yrs[0]].corr(w[yrs[-1]], method="spearman"),
                     "pearson_log_first_last": w[yrs[0]].corr(w[yrs[-1]])})
    return pd.DataFrame(rows)
