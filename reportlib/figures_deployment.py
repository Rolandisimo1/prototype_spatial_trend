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
DARK_BG, DARK_INK, DARK_SUB, DARK_RULE = "#1b1f24", "#e8ecef", "#9aa4ad", "#49515a"


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


DEPLOYMENTS = ("/Users/rwkays/claude_code/data_integration_arielle/raw_cam_data/"
               "combined_deployments_all.csv")


def array_change(table=None, deployments=DEPLOYMENTS, **kw):
    """Aggregate per-deployment change to the camera array, i.e. the subproject.

    An array is the unit a field team actually operates: one project's
    subproject, typically about a dozen cameras placed together. The product
    files carry no array field, so deployments are matched to their array on
    exact coordinates, which resolves every row.

    Each array gets the mean change across its deployments and the centroid of
    their positions. Averaging within an array is reasonable because the
    deployments in one array sit close together and so share almost all of the
    linear predictor apart from local habitat.
    """
    t = table if table is not None else deployment_change(**kw)
    dep = pd.read_csv(deployments, low_memory=False)
    dep["array_id"] = (dep.project_id.astype(str) + " | " +
                       dep.subproject_name.fillna("(no subproject)").astype(str))
    for d, a, b in ((dep, "longitude", "latitude"), (t, "lon", "lat")):
        d["lon6"] = d[a].round(6)
        d["lat6"] = d[b].round(6)
    key = dep.drop_duplicates(["lon6", "lat6"]).set_index(["lon6", "lat6"])
    j = t.join(key[["array_id", "project_id", "subproject_name"]], on=["lon6", "lat6"])
    unmatched = int(j.array_id.isna().sum())
    if unmatched:
        raise AssertionError(
            f"{unmatched} of {len(j)} deployments did not match an array on "
            "coordinates; the aggregation would silently drop them")
    g = (j.groupby(["species", "array_id"])
           .agg(change_pp=("change_pp", "mean"),
                change_pp_min=("change_pp", "min"),
                change_pp_max=("change_pp", "max"),
                psi_as_fitted=("psi_as_fitted", "median"),
                n_deployments=("site_id", "size"),
                year_first=("year_first", "first"), year_last=("year_last", "first"),
                lon=("lon", "mean"), lat=("lat", "mean"),
                project_id=("project_id", "first"),
                subproject_name=("subproject_name", "first"))
           .reset_index())
    return g


def trend_status(posteriors=None):
    """Per species: does the camera trend's 95% interval exclude zero?

    An array map paints a point estimate at every array, which reads as a
    confident pattern even when the species' own interval spans zero. House rule
    is that a trend overlapping zero is drawn grey rather than coloured, so this
    supplies the flag.
    """
    posteriors = posteriors or fw.PULL
    out = {}
    for f in sorted(glob.glob(os.path.join(posteriors, "*_posterior_win5.csv"))):
        d = pd.read_csv(f)
        b = d[d.param.astype(str).str.fullmatch("year_beta")]
        if not len(b):
            continue
        lo, hi = float(b["q2.5"].iloc[0]), float(b["q97.5"].iloc[0])
        out[fw._pretty(d.species.iloc[0])] = bool(lo > 0 or hi < 0)
    return out


def fig_array_change(out="fig_array_change.png", table=None, pad=0.08,
                     species=None, page=None, n_pages=None, resolved=None,
                     posteriors=None, slide=False, **kw):
    """Map the change in fitted occupancy per camera array, zoomed to each species.

    Each species is drawn to its own extent and scaled to its own range, because a
    species occupying a few dozen arrays is invisible at continental scale and a
    shared colour scale leaves every small-magnitude species uniformly pale. The
    consequence is that panels are comparable in neither position nor colour, so
    each panel carries its own colour bar in its own units rather than sharing one.

    Four species per figure. Eleven panels on one sheet left each colour bar's
    label overlapping the title of the panel below it, and shrinking the type to
    fit was what made the first version unreadable.

    `slide=True` renders the projection variant: the caption and the legend are
    dropped and the type scales up. A single-species slide call needs neither,
    and at projection size the caption is a grey band nobody can read.
    """
    fs = 1.5 if slide else 1.0
    import geopandas as gpd
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.lines import Line2D
    from mpl_toolkits.axes_grid1.inset_locator import inset_axes
    t = table if table is not None else array_change(**kw)
    cam_max = float(t.n_deployments.max())       # global, so the legend reads true
    if species is not None:
        t = t[t.species.isin(species)]
    st = conus()
    sp = sorted(t.species.unique())
    lims = {s_: float(np.abs(t.loc[t.species == s_, "change_pp"]).max()) or 1.0
            for s_ in sp}
    if resolved is None:
        resolved = trend_status(posteriors)
    ncol = 2 if len(sp) > 1 else 1
    # The +1 reserves a panel for the marker-area legend. A slide drops the
    # legend, so reserving the slot there leaves half the frame empty and
    # shrinks the map to a quarter of the screen.
    nrow = int(np.ceil((len(sp) + (0 if slide else 1)) / ncol))
    fig, axes = plt.subplots(nrow, ncol,
                             figsize=(7.2 * ncol * (1.5 if slide else 1.0),
                                      6.4 * nrow * (1.2 if slide else 1.0)),
                             facecolor=DARK_BG)
    axes = np.atleast_2d(axes)
    fig.subplots_adjust(hspace=0.30, wspace=0.05)
    for ax in axes.ravel():
        ax.set_axis_off()
        ax.set_facecolor(DARK_BG)
    for k, s_ in enumerate(sp):
        ax = axes[k // ncol, k % ncol]
        ax.set_axis_off()
        ax.set_facecolor(DARK_BG)
        d = t[t.species == s_]
        pts = gpd.GeoDataFrame(d, geometry=gpd.points_from_xy(d.lon, d.lat),
                               crs=4326).to_crs(st.crs)
        st.boundary.plot(ax=ax, color=DARK_RULE, linewidth=0.6, zorder=1)
        size = 18 + 90 * (d.n_deployments / cam_max) ** 0.5
        lim = lims[s_]
        # A species whose trend interval spans zero is drawn grey. Colour on such
        # a species would paint a confident spatial pattern from a slope the data
        # do not distinguish from no change.
        ok = bool(resolved.get(s_, False))
        if ok:
            sc = ax.scatter(pts.geometry.x, pts.geometry.y, c=d.change_pp,
                            cmap="RdBu_r", vmin=-lim, vmax=lim, s=size,
                            linewidths=0.35, edgecolors=DARK_BG, zorder=3)
            cvals, cticks = (-lim, lim), [-lim, 0, lim]
        else:
            # Unresolved species are keyed on MAGNITUDE over a light ramp. A
            # signed grey map puts one end in near-black, which disappears on this
            # background -- elk, the only species with a positive slope, vanished.
            # Sign is uniform within a species here, so magnitude loses nothing.
            ramp = LinearSegmentedColormap.from_list("dim", ["#5a6168", "#e9edf0"])
            sc = ax.scatter(pts.geometry.x, pts.geometry.y, c=np.abs(d.change_pp),
                            cmap=ramp, vmin=0, vmax=lim, s=size,
                            linewidths=0.35, edgecolors=DARK_BG, zorder=3)
            cvals, cticks = (0, lim), [0, lim]
        xmin, ymin, xmax, ymax = pts.total_bounds
        side = max(max(xmax - xmin, 1e4), max(ymax - ymin, 1e4)) * (1 + 2 * pad)
        cx, cy = (xmin + xmax) / 2, (ymin + ymax) / 2
        ax.set_xlim(cx - side / 2, cx + side / 2)
        ax.set_ylim(cy - side / 2, cy + side / 2)
        # Own line, not trailing the species name: on a two-column sheet the
        # trailing version runs straight into the next panel's title.
        flag = "" if ok else "\ntrend not distinguishable from zero"
        ax.set_title(f"{s_}{flag}\n{len(d)} arrays · "
                     f"{int(d.n_deployments.sum()):,} camera-years "
                     f"· median {d.change_pp.median():+.1f} pp",
                     fontsize=15 * fs, color=DARK_INK if ok else DARK_SUB, pad=10)
        cax = inset_axes(ax, width="62%", height="2.8%", loc="lower center",
                         borderpad=2.2)
        cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
        cb.set_ticks(cticks)
        cb.set_ticklabels([f"{v:+.1f}".replace("+0.0", "0") for v in cticks]
                          if ok else ["0", f"{lim:.1f}"])
        cb.ax.tick_params(labelsize=11 * fs, colors=DARK_INK, length=3, pad=2)
        cb.outline.set_edgecolor(DARK_RULE)
        cb.set_label("change in detection probability (pp)" if ok
                     else "size of change (pp), direction unresolved",
                     fontsize=10.5 * fs, color=DARK_SUB, labelpad=4)
    if not slide:
        leg = [Line2D([], [], marker="o", color="none", markerfacecolor="0.78",
                      markersize=np.sqrt(18 + 90 * f ** 0.5) * 0.80,
                      label=f"{int(round(f * cam_max))} camera-years")
               for f in (0.05, 0.4, 1.0)]
        # A single-species call leaves no spare slot, so the legend must fall back
        # into the only panel rather than indexing a row that does not exist.
        spare = [axes[r, c] for r in range(axes.shape[0]) for c in range(axes.shape[1])
                 if r * axes.shape[1] + c >= len(sp)]
        host_ax = spare[0] if spare else axes[0, 0]
        host_ax.legend(handles=leg, loc="upper center", frameon=False, fontsize=13,
                       labelcolor=DARK_INK, handletextpad=1.6, borderpad=0.6,
                       labelspacing=2.0, title="Marker area",
                       title_fontproperties={"size": 14})
        host_ax.get_legend().get_title().set_color(DARK_INK)
    tag = f"  ({page} of {n_pages})" if page and n_pages else ""
    if not slide:
        fig.suptitle("Where each camera array's chance of detecting the species changed"
                     + tag, fontsize=18, color=DARK_INK, y=0.995)
    if not slide:
        fig.text(0.5, 0.072,
                 "Change in the probability that a camera array detects the species, from the "
                 "first to the last year of the five-year window.\n"
                 "Each point is one camera array, meaning one project's subproject, averaged over "
                 "its cameras. Every camera is fitted in the year it ran and projected to both "
                 "ends of the window using the\nfitted trend, so this is the change the model "
                 "implies at that array rather than a measured before and after: almost no camera "
                 "ran in both years. The model applies one\ntrend to every grid cell, so arrays "
                 "differ only through where each one starts, and one slope moves a rarely-visited "
                 "array by many points and an array already near\ncertain detection by almost "
                 "none. Panels are zoomed and scaled independently, so neither position nor colour "
                 "is comparable between them; read each against its own bar.\n"
                 "Species drawn in grey are those whose own trend interval spans zero: the "
                 "arrays still carry the point estimate, but the data do not distinguish that "
                 "slope from no change at all.",
                 ha="center", va="top", fontsize=12, color=DARK_SUB, linespacing=1.7)
    fig.savefig(out, facecolor=DARK_BG, bbox_inches="tight", pad_inches=0.3, dpi=160)
    plt.close(fig)
    return out, t


def fig_array_change_pages(stem="fig_array_change", table=None, per_page=4, **kw):
    """Write the array maps as several sheets, four species each."""
    t = table if table is not None else array_change(**kw)
    sp = sorted(t.species.unique())
    pages = [sp[i:i + per_page] for i in range(0, len(sp), per_page)]
    return [fig_array_change(out=f"{stem}_p{i}.png", table=t, species=g,
                             page=i, n_pages=len(pages))[0]
            for i, g in enumerate(pages, 1)]


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
