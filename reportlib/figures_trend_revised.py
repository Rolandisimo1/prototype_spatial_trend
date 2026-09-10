"""Trend figures under the camera-anchored reading.

Supersedes the earlier national-trend figure, which reported `total_var_beta`.
That parameter is the trend the iNaturalist count stream observes; the camera
occupancy submodel contains only `year_beta`, and `year_var` is their difference,
identified by iNaturalist alone. Following Goldstein et al., the reported
ecological result is the camera submodel's parameter, so `year_beta` leads and
the other two are shown as context rather than as the headline.
"""
import numpy as np
import matplotlib.pyplot as plt

CAM, INAT, TOTAL = "#2c6c9c", "#e08214", "#7a7a7a"
PARAMS = [("year_beta", "Camera-anchored trend\n$year\\_beta$", CAM),
          ("year_var", "iNaturalist-only increment\n$year\\_var$", INAT),
          ("total_var_beta", "iNaturalist-side total\n$total\\_var\\_beta$", TOTAL)]
SPECIES = ["Bobcat", "White-tailed deer", "Moose"]


def fig_trend_decomposition(dec, out="fig_trend_decomposition.png"):
    """Both trend components, all three species, both parameterisations."""
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.5), sharex=True)
    for ax, sp in zip(axes, SPECIES):
        d = dec[dec.species == sp]
        ax.axvline(0, color="0.55", lw=0.8, zorder=1)
        for k, (pnm, _, col) in enumerate(PARAMS):
            yb = len(PARAMS) - 1 - k
            for off, par, mk, fc in [(+0.15, "national", "o", col),
                                     (-0.15, "ecoregion", "s", "white")]:
                r = d[(d.parameter == pnm) & (d.parameterisation == par)]
                if not len(r):
                    continue
                r = r.iloc[0]
                sig = bool(r.excludes_zero)
                ax.plot([r.q025, r.q975], [yb + off] * 2, color=col,
                        lw=2.0 if sig else 1.1, alpha=1.0 if sig else 0.55,
                        solid_capstyle="round", zorder=3)
                ax.plot(r["mean"], yb + off, mk, color=col, ms=5.0,
                        mfc=fc, mew=1.1, alpha=1.0 if sig else 0.55, zorder=4)
        ax.set_yticks(range(len(PARAMS)))
        ax.set_yticklabels([lbl for _, lbl, _ in PARAMS][::-1], fontsize=6.8)
        ax.set_title(sp, fontsize=8.8)
        ax.set_xlabel("Trend coefficient (log scale, per year unit)")
        ax.set_ylim(-0.6, len(PARAMS) - 0.4)
        if ax is not axes[0]:
            ax.tick_params(labelleft=False)
    h = [plt.Line2D([], [], marker="o", color="0.35", ls="-", ms=5, label="national-scalar fit"),
         plt.Line2D([], [], marker="s", color="0.35", ls="-", ms=5, mfc="white",
                    label="ecoregion fit")]
    axes[0].legend(handles=h, frameon=False, fontsize=6.0, loc="lower left")
    fig.suptitle("Only white-tailed deer has a trend the cameras support on their own",
                 fontsize=9.6, y=1.04)
    fig.text(0.5, -0.06,
             "Bold intervals exclude zero; faded ones include it. The camera occupancy submodel "
             "contains only the top row, so that is the trend the camera surveys\nevidence "
             "directly. The middle row is identified by iNaturalist alone and absorbs any "
             "iNaturalist-specific temporal drift; the bottom row is their sum, which is\nwhat "
             "earlier drafts reported. Both parameterisations are shown: the camera-anchored "
             "estimate barely moves between them, because the camera submodel\nis identical in "
             "both.", ha="center", va="top", fontsize=6.4, color="#5f6a73")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return out


def fig_sim_components(g, out="fig_sim_components.png"):
    """Where the simulated estimator difference actually sits.

    The sweep scored bias on `total_var_beta`. Split into components, the four
    arms agree closely on `year_beta` and diverge on `year_var` -- so the
    headline 'occupancy overstates declines 3-6x' is a statement about the
    iNaturalist-only term, not about the camera-anchored trend.
    """
    order = ["bobcat_like", "intermediate", "deer_like"]
    pretty = {"bobcat_like": "bobcat-like", "intermediate": "intermediate",
              "deer_like": "deer-like"}
    arms = {"camera_occ": ("#2c6c9c", "o", "camera occupancy"),
            "camera_rn": ("#3f9e4d", "s", "camera Royle-Nichols"),
            "array_occ": ("#cb4335", "^", "array occupancy"),
            "array_rn": ("#8073ac", "D", "array Royle-Nichols")}
    cols = [("year_beta_mean", "Camera-anchored trend\n$year\\_beta$"),
            ("year_var_mean", "iNaturalist-only increment\n$year\\_var$"),
            ("total_var_beta_mean", "Their sum\n$total\\_var\\_beta$")]
    fig, axes = plt.subplots(1, 3, figsize=(10.4, 3.4))
    x = np.arange(len(order))
    for ax, (col, lbl) in zip(axes, cols):
        for arm, (c, mk, nm) in arms.items():
            dd = g[g.estimator == arm].set_index("abundance").reindex(order)
            ax.plot(x, dd[col].values, marker=mk, color=c, ms=4.6, lw=1.3, label=nm)
        if col == "total_var_beta_mean":
            ax.axhline(g.tvb_true.iloc[0], color="0.25", ls="--", lw=1.0)
            ax.annotate("truth", (x[0], g.tvb_true.iloc[0]), textcoords="offset points",
                        xytext=(2, 5), fontsize=6.2, color="0.25")
        else:
            ax.axhline(0, color="0.7", lw=0.8)
        ax.set_xticks(x); ax.set_xticklabels([pretty[o] for o in order], fontsize=6.8)
        ax.set_title(lbl, fontsize=8.2)
        ax.set_xlabel("Simulated abundance level")
    axes[0].set_ylabel("Posterior mean across replicates")
    axes[1].legend(frameon=False, fontsize=6.0, loc="lower left")
    fig.suptitle("The simulated estimator difference is in the iNaturalist term, "
                 "not the camera-anchored trend", fontsize=9.6, y=1.04)
    fig.text(0.5, -0.07,
             "Truth is only recorded for the SUM (-0.1795); the split between the two components "
             "was set from a real fit's posterior means but never written out,\nso no truth line "
             "can be drawn on the first two panels and neither arm can be called unbiased there "
             "-- only compared with each other. At deer-like\nabundance the two camera arms "
             "differ by 0.002 on the camera-anchored trend and by 0.455 on the iNaturalist term. "
             "Array-level occupancy is the one\narm that is genuinely off on the camera-anchored "
             "trend.", ha="center", va="top", fontsize=6.4, color="#5f6a73")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return out
