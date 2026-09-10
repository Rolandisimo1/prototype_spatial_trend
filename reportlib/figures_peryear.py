"""Trend figures in per-year units, and the simulation scored on the camera trend.

The two submodels' year covariates are standardised on DIFFERENT distributions:
`year_vals` (iNaturalist side) is a generic 18-point sequence, identical for
every species, one standard deviation per 5.34 years. `year_occ` (camera side)
is z-scored on each species' own camera site-year distribution, so one calendar
year is a different distance for each species -- 0.531 for moose (7 distinct
camera years), 0.274 for bobcat, 0.258 for white-tailed deer, against 0.187 for
`year_vals` throughout.

Raw coefficients are therefore not comparable between the two streams, or
between species on the camera side. Everything here is converted to proportional
change per calendar year first: coefficient x step, then exp() - 1.
"""
import numpy as np
import matplotlib.pyplot as plt

YEAR_VALS_STEP = 0.1873172
OCC_STEP = {"Moose": 0.5309095, "Bobcat": 0.2740821, "White-tailed deer": 0.2580241}
CAM_YEARS = {"Moose": 7, "Bobcat": 16, "White-tailed deer": 16}
CAM, INAT = "#2c6c9c", "#e08214"
SPECIES = ["Bobcat", "White-tailed deer", "Moose"]


def per_year_table(fits):
    """Both streams' trends as percent change per calendar year."""
    import pandas as pd
    rows = []
    for (sp, par), f in fits.items():
        t = pd.read_csv(f).set_index("parameter")
        for pnm, step, stream in [("year_beta", OCC_STEP[sp], "Camera surveys"),
                                  ("total_var_beta", YEAR_VALS_STEP, "iNaturalist")]:
            r = t.loc[pnm]
            rows.append(dict(
                species=sp, parameterisation=par, stream=stream, parameter=pnm,
                coefficient=r["mean"], step_per_year=step,
                pct_per_year=100 * (np.exp(r["mean"] * step) - 1),
                pct_lo=100 * (np.exp(r["q025"] * step) - 1),
                pct_hi=100 * (np.exp(r["q975"] * step) - 1),
                excludes_zero=bool(r["q025"] > 0 or r["q975"] < 0),
                years_covered=CAM_YEARS[sp] if pnm == "year_beta" else 18))
    return pd.DataFrame(rows)


def fig_trend_per_year(tab, out="fig_trend_per_year.png"):
    """Camera and iNaturalist trends side by side, in comparable units."""
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    ax.axvline(0, color="0.55", lw=0.9, zorder=1)
    ylab, ytick = [], []
    for i, sp in enumerate(SPECIES):
        base = (len(SPECIES) - 1 - i) * 2.6
        for k, (stream, col) in enumerate([("Camera surveys", CAM), ("iNaturalist", INAT)]):
            for off, par, mk, fc in [(+0.28, "national", "o", col),
                                     (-0.28, "ecoregion", "s", "white")]:
                r = tab[(tab.species == sp) & (tab.stream == stream) &
                        (tab.parameterisation == par)]
                if not len(r):
                    continue
                r = r.iloc[0]
                y = base + (1 - k) * 0.95 + off * 0.42
                sig = bool(r.excludes_zero)
                ax.plot([r.pct_lo, r.pct_hi], [y, y], color=col,
                        lw=2.1 if sig else 1.1, alpha=1.0 if sig else 0.5,
                        solid_capstyle="round", zorder=3)
                ax.plot(r.pct_per_year, y, mk, color=col, ms=5.2, mfc=fc, mew=1.1,
                        alpha=1.0 if sig else 0.5, zorder=4)
            if off:
                ytick.append(base + (1 - k) * 0.95)
                ylab.append(f"{stream}" + (f"\n({CAM_YEARS[sp]} yr)" if k == 0 else "\n(18 yr)"))
        ax.text(-0.13, base + 0.5, sp, transform=ax.get_yaxis_transform(),
                ha="right", va="center", fontsize=8.6)
    ax.set_yticks(ytick); ax.set_yticklabels(ylab, fontsize=6.4)
    ax.set_xlabel("Change in relative abundance per year (%)")
    ax.set_ylim(-0.9, (len(SPECIES) - 1) * 2.6 + 1.9)
    h = [plt.Line2D([], [], marker="o", color="0.35", ls="-", ms=5, label="national-scalar fit"),
         plt.Line2D([], [], marker="s", color="0.35", ls="-", ms=5, mfc="white",
                    label="ecoregion fit")]
    ax.legend(handles=h, frameon=False, fontsize=6.2, loc="lower right")
    ax.set_title("Only white-tailed deer shows a trend both data streams agree on",
                 fontsize=9.4)
    fig.text(0.5, -0.06,
             "Bold intervals exclude zero; faded ones include it. Both streams are converted to "
             "proportional change per calendar year, which their raw\ncoefficients are not -- the "
             "two submodels standardise their year covariate on different distributions. The "
             "camera surveys cover 16 years for\nbobcat and deer but only 7 for moose (cameras "
             "from 2019), so moose's two rows describe different periods.",
             ha="center", va="top", fontsize=6.4, color="#5f6a73")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return out


def fig_sim_camera_trend(sb, out="fig_sim_camera_trend.png"):
    """Simulated recovery of the camera-anchored trend, against its known truth."""
    order = ["bobcat_like", "intermediate", "deer_like"]
    pretty = {"bobcat_like": "bobcat-like", "intermediate": "intermediate",
              "deer_like": "deer-like"}
    arms = {"camera_occ": ("#2c6c9c", "o", "camera occupancy"),
            "camera_rn": ("#3f9e4d", "s", "camera Royle-Nichols"),
            "array_occ": ("#cb4335", "^", "array occupancy"),
            "array_rn": ("#8073ac", "D", "array Royle-Nichols")}
    yb_truth = 0.06450236525
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(8.8, 3.4))
    x = np.arange(len(order))
    for arm, (c, mk, nm) in arms.items():
        dd = sb[sb.estimator == arm].set_index("abundance").reindex(order)
        axL.plot(x, dd.yb_mean.values, marker=mk, color=c, ms=4.8, lw=1.3, label=nm)
        axR.plot(x, dd.yb_bias.values, marker=mk, color=c, ms=4.8, lw=1.3, label=nm)
    axL.axhline(yb_truth, color="0.25", ls="--", lw=1.0)
    axL.annotate("truth  +0.0645", (x[0], yb_truth), textcoords="offset points",
                 xytext=(2, 5), fontsize=6.2, color="0.25")
    axR.axhline(0, color="0.25", ls="--", lw=1.0)
    for ax, lbl, ti in [(axL, "Posterior mean of $year\\_beta$",
                         "Every arm pulls the camera trend toward zero"),
                        (axR, "Bias in $year\\_beta$",
                         "Array-level occupancy is the worst by a wide margin")]:
        ax.set_xticks(x); ax.set_xticklabels([pretty[o] for o in order], fontsize=6.8)
        ax.set_xlabel("Simulated abundance level"); ax.set_ylabel(lbl)
        ax.set_title(ti, fontsize=8.4)
    axR.legend(frameon=False, fontsize=6.0, loc="lower right")
    fig.text(0.5, -0.07,
             "Truth for the camera-anchored trend is +0.0645 and for the iNaturalist-only term "
             "-0.2440, taken from a real converged fit's posterior means. The two\ncomponents "
             "have opposite signs in the truth, and every arm shrinks the camera component "
             "toward zero and compensates in the other -- so the split is\npoorly recovered by "
             "all four, not just by the biased ones. Bias is the posterior mean minus truth, "
             "averaged over 60 replicates per cell.",
             ha="center", va="top", fontsize=6.4, color="#5f6a73")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return out
