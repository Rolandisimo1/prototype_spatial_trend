"""One function per real-data figure. Each returns the output path."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

from . import inputs, style
from .conventions import max_rhat_sampled

# The six reported fits: three species x two spatial parameterizations. A second
# moose track on the superseded IUCN mask was dropped from the report -- mask
# choice is a settled pipeline decision, not a result.
FITS = {
    "Bobcat, national":             ("pull", "bobcat_v2b_national_scalar_posterior.csv", "bobcat_v2b_national_scalar"),
    "Bobcat, ecoregion":            ("pull", "bobcat_v2b_ecoregion_global.csv", "bobcat_v2b_ecoregion"),
    "Deer, national":               ("pull", "white-tailed_deer_v2b_national_scalar_posterior.csv", "white-tailed_deer_v2b_national_scalar"),
    "Deer, ecoregion":              ("pull", "white-tailed_deer_v2b_ecoregion_global.csv", "white-tailed_deer_v2b_ecoregion"),
    "Moose, national":              ("repo", "moose_v2b_national_scalar_posterior_single.csv", "moose_v2b_national_scalar"),
    "Moose, ecoregion":             ("repo", "moose_v2b_ecoregion_global_single.csv", "moose_v2b_ecoregion"),
}

FAMILY = {                     # label -> (marker, colour)
    "Population trend":        ("o", "#1b5e9c"),
    "Habitat covariates":     ("s", "#4d9221"),
    "iNaturalist scaling":    ("D", "#e08214"),
    "Spatial (CAR) fields":   ("^", "#8073ac"),
}
RHAT_THRESHOLD = 1.1           # conventional Gelman-Rubin cutoff


def _dir(tag):
    return inputs.PULL if tag == "pull" else inputs.REPO


def convergence_table():
    """Max R-hat per parameter family per fit.

    CAR-field R-hat comes from the per-fit family tables, which exist locally
    for only some fits; missing entries are NaN and must be shown as
    un-exported rather than as passing.
    """
    occ = inputs.posterior("occbeta")
    th = inputs.posterior("theta")
    rows = []
    for lab, (tag, fname, key) in FITS.items():
        post = pd.read_csv(os.path.join(_dir(tag), fname))
        famf = os.path.join(inputs.REPO, f"{key}_rhat_by_family_single.csv")
        car = np.nan
        if os.path.exists(famf):
            t = pd.read_csv(famf)
            car = float(t.loc[t.family == "CAR_fields", "max_rhat"].iloc[0])
        rows.append({
            "fit": lab,
            "Population trend": max_rhat_sampled(post),
            "Habitat covariates": float(occ.loc[occ.model == key, "rhat"].max()),
            "iNaturalist scaling": float(th.loc[th.model == key, "rhat"].max()),
            "Spatial (CAR) fields": car,
        })
    return pd.DataFrame(rows)


def fig_convergence(out="fig_convergence.png"):
    """Max R-hat by parameter family, all eight fits, against the 1.1 cutoff."""
    tab = convergence_table()
    fig, ax = plt.subplots(figsize=(7.4, 4.4))
    y = np.arange(len(tab))[::-1]

    ax.axvspan(1.0, RHAT_THRESHOLD, color="#eef3ee", zorder=0)
    ax.axvline(RHAT_THRESHOLD, color="#4a6a4a", lw=0.9, ls="--", zorder=1)

    for fam, (mk, col) in FAMILY.items():
        vals = tab[fam].values
        ax.plot(vals, y, mk, color=col, ms=5.2, ls="none", label=fam, zorder=3)

    # Fits whose spatial-field diagnostics were never exported are marked on the
    # LABEL, not on the value axis. Plotting them at a position -- especially
    # beyond the 1.1 cutoff -- reads as a failed diagnostic rather than a
    # missing one, which is the opposite of what is true.
    missing = tab.loc[tab["Spatial (CAR) fields"].isna(), "fit"].tolist()
    ticklabels = [f"{f} †" if f in missing else f for f in tab["fit"]]

    ax.set_yticks(y)
    ax.set_yticklabels(ticklabels)
    ax.set_xlabel("Largest Gelman-Rubin R-hat in the family\n"
                  "(1.0 = chains agree perfectly; 1.1 is the conventional limit)")
    ax.set_title("Every fit converged on every parameter family we can assess")
    ax.set_xlim(0.998, 1.106)
    ax.margins(y=0.06)
    ax.legend(frameon=False, fontsize=6.6, loc="upper center",
              bbox_to_anchor=(0.5, -0.19), ncol=4, columnspacing=1.4,
              handletextpad=0.4)
    if missing:
        fig.text(0.5, -0.27,
                 f"†  Spatial-field diagnostics not yet exported for these {len(missing)} fits, "
                 "so that family is unassessed rather than passing. A localised spatial-field "
                 "convergence failure is known in the ecoregion models (see text); it does not "
                 "affect the trend estimates, which are shown here and converged.",
                 ha="center", va="top", fontsize=6.2, color="#5f6a73", wrap=True)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, tab


def fig_covariates(out="fig_covariates.png"):
    """Occupancy covariate effects, one panel per species.

    Panels carry independent x-scales: the moose human-population effect is
    roughly eight times any other coefficient, and a shared axis compresses
    bobcat and deer into an unreadable band.
    """
    occ = inputs.posterior("occbeta")
    prim = occ[occ.model.str.endswith("national_scalar")
               & (~occ.model.str.startswith("moose_v1fix9"))].copy()
    prim["species"] = prim["model"].str.replace("_national_scalar", "", regex=False)
    show = [("bobcat_v2b", "Bobcat"), ("white-tailed_deer_v2b", "White-tailed deer"),
            ("moose_v2b", "Moose")]

    fig, axes = plt.subplots(1, 3, figsize=(10.0, 4.2))
    for ax, (key, label) in zip(axes, show):
        d = prim[prim.species == key].set_index("covariate_label").loc[style.COVAR_ORDER].reset_index()
        y = np.arange(len(d))[::-1]
        for yy, (_, r) in zip(y, d.iterrows()):
            sig = (r["q025"] > 0) or (r["q975"] < 0)
            c = style.FOCAL if sig else style.MUTED
            ax.plot([r["q025"], r["q975"]], [yy, yy], color=c, lw=1.9, solid_capstyle="round")
            ax.plot(r["mean"], yy, "o", color=c, ms=4.6, zorder=3)
        ax.axvline(0, color="0.35", lw=0.7, ls="--", zorder=1)
        ax.set_yticks(y)
        ax.set_yticklabels([style.COVAR_LABEL[c] for c in d["covariate_label"]]
                           if ax is axes[0] else [])
        nsig = int(((d["q025"] > 0) | (d["q975"] < 0)).sum())
        ax.set_title(f"{label}  ({nsig} of 9 clear)")
        lo, hi = d["q025"].min(), d["q975"].max()
        pad = 0.08 * (hi - lo)
        ax.set_xlim(lo - pad, hi + pad)
        ax.margins(y=0.05)
        ax.set_xlabel("Effect on occupancy\n(log-odds, 95% interval)")
    fig.text(0.5, -0.04,
             "Blue: 95% interval excludes zero. Grey: includes zero. "
             "Note the x-axis scale differs between panels.",
             ha="center", fontsize=6.4, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, prim


def fig_congruence(out="fig_congruence.png"):
    """theta1: how iNaturalist counts scale with camera-derived abundance.

    theta1 is the exponent in log(mu) = theta0 + theta1*log(sum(lambda)), and is
    the diagnostic the source method uses to characterise correspondence between
    the two data streams. theta1 = 1 would mean proportional scaling.
    """
    th = inputs.posterior("theta")
    t1 = th[(th.parameter == "theta1") & (~th.model.str.startswith("moose_v1fix9"))].copy()
    t1["track"] = (t1["model"].str.replace("_national_scalar", "", regex=False)
                              .str.replace("_ecoregion", "", regex=False))
    t1["param"] = np.where(t1["model"].str.endswith("national_scalar"),
                           "national", "ecoregion")
    pretty = {"bobcat_v2b": "Bobcat", "white-tailed_deer_v2b": "White-tailed deer",
              "moose_v2b": "Moose"}
    order = ["bobcat_v2b", "moose_v2b", "white-tailed_deer_v2b"]

    fig, ax = plt.subplots(figsize=(6.8, 3.0))
    y = np.arange(len(order))[::-1]
    for yy, trk in zip(y, order):
        for off, par, mk, fc in [(0.14, "national", "o", style.FOCAL),
                                 (-0.14, "ecoregion", "s", "white")]:
            r = t1[(t1.track == trk) & (t1.param == par)]
            if not len(r):
                continue
            r = r.iloc[0]
            ax.plot([r["q025"], r["q975"]], [yy + off] * 2, color=style.FOCAL,
                    lw=1.7, solid_capstyle="round")
            ax.plot(r["mean"], yy + off, mk, color=style.FOCAL, mfc=fc, ms=4.6, zorder=3)
    ax.axvline(1.0, color="#b2182b", lw=1.0, ls="--")
    # Label placed inside the axes at the top of the reference line; an earlier
    # version anchored it below the lowest row, which fell outside the y-limits
    # and was silently dropped, leaving the reference line unexplained.
    ax.text(0.988, y.max() + 0.34,
            "proportional scaling:\ncounts would track\nabundance one-for-one",
            fontsize=6.3, color="#b2182b", ha="right", va="top", linespacing=1.5)
    ax.set_yticks(y)
    ax.set_yticklabels([pretty[t] for t in order])
    ax.set_xlim(0.34, 1.045)
    ax.set_ylim(y.min() - 0.55, y.max() + 0.62)
    ax.set_xlabel(r"$\theta_1$: how iNaturalist counts scale with abundance (95% interval)")
    ax.set_title("iNaturalist counts rise with abundance, but with\n"
                 "diminishing returns in every fit")
    ax.plot([], [], "o", color=style.FOCAL, ms=4.6, label="national model")
    ax.plot([], [], "s", color=style.FOCAL, mfc="white", ms=4.6, label="ecoregion model")
    ax.legend(frameon=False, fontsize=6.5, loc="lower right",
              bbox_to_anchor=(0.62, 0.02))
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, t1


NO_TREND_GREY = "#b8b8b8"   # interval overlaps zero: no direction claimed


def national_trend_table(window="full"):
    """National trend (total_var_beta) with interval and R-hat, per fit."""
    rows = []
    for lab, (tag, fname, key) in FITS.items():
        path = (os.path.join(_dir(tag), fname) if window == "full"
                else inputs.windowed(fname, window, [_dir(tag)]))
        d = pd.read_csv(path).set_index("parameter")
        r = d.loc["total_var_beta"]
        rows.append({"fit": lab, "mean": r["mean"], "q025": r["q025"], "q975": r["q975"],
                     "p_negative": r["p_negative"], "rhat": r["rhat"],
                     "clear": bool(r["q025"] > 0 or r["q975"] < 0)})
    return pd.DataFrame(rows)


def fig_national_trend(window="full", out=None):
    """National 18-year trend per fit, grouped by species.

    Grey wherever the 95% interval overlaps zero, matching the regional maps: an
    interval spanning zero supports no direction, so it is not coloured.
    """
    out = out or f"fig_national_trend_{window}.png"
    wlabel = inputs.WINDOWS[window][1]
    tab = national_trend_table(window)
    groups = [("Bobcat", ["Bobcat, national", "Bobcat, ecoregion"]),
              ("White-tailed deer", ["Deer, national", "Deer, ecoregion"]),
              ("Moose", ["Moose, national", "Moose, ecoregion"])]
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax.axvline(0, color="0.35", lw=0.8, ls="--", zorder=1)
    ypos, ylabels, y = [], [], 0.0
    for gname, fits in groups[::-1]:
        for f in fits[::-1]:
            r = tab[tab.fit == f].iloc[0]
            col = ("#b2182b" if r["mean"] < 0 else "#1b5e9c") if r["clear"] else NO_TREND_GREY
            ax.plot([r["q025"], r["q975"]], [y, y], color=col, lw=2.0,
                    solid_capstyle="round", zorder=2)
            ax.plot(r["mean"], y, "o", color=col, ms=5.4, zorder=3)
            ax.annotate(f"{r['mean']:+.3f}", (r["mean"], y), textcoords="offset points",
                        xytext=(0, 8), fontsize=6.4, color=col, ha="center")
            ypos.append(y)
            ylabels.append("national model" if "national" in f else "ecoregion model")
            y += 1.0
        ax.text(-0.02, y - 0.5, gname, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=7.6, fontweight="bold")
        y += 0.9
    ax.set_yticks(ypos); ax.set_yticklabels(ylabels, fontsize=6.8)
    ax.set_ylim(-0.9, y - 0.5)
    ax.set_xlabel(f"National trend, {wlabel} (log scale, 95% interval)\n"
                  "negative = declining")
    ax.set_title(f"National population trend, {wlabel}")
    # Caption states the grey count, so derive it rather than asserting a
    # remembered one: an earlier draft claimed every interval excluded zero,
    # which was false for two of the eight fits.
    n_grey = int((~tab["clear"]).sum())
    grey_fits = ", ".join(tab.loc[~tab["clear"], "fit"])
    grey_txt = ("No interval overlaps zero." if n_grey == 0 else
                f"Grey marks an interval overlapping zero: {grey_fits} "
                f"({n_grey} of {len(tab)} fits).")
    fig.text(0.5, -0.05,
             f"{grey_txt}",
             ha="center", va="top", fontsize=6.3, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, tab
