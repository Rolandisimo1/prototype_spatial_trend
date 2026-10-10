"""Why a joint fit can show a trend that neither data stream shows alone.

The mechanism behind western gray squirrel's -13.2 %/yr, which sits outside
both the camera-only estimate (-2.7) and the iNaturalist one (+1.8). It has
nothing to do with either stream's trend: it comes from the two streams
disagreeing about SPACE, combined with cameras moving between cells over time.

Numbers in the schematic are a two-cell illustration, not fitted values; the
real decomposition is in table form. The illustration is built from the same
arithmetic the model does -- effort-weighted average of a spatial field -- so
the sign and the logic are exact even though the magnitudes are chosen.
"""
import numpy as np
import matplotlib.pyplot as plt

from . import style

# Two cells, A and B. Camera effort shifts from A to B across the record;
# true occupancy is flat in time and equal in both cells.
EFFORT = {"early": {"A": 0.8, "B": 0.2}, "late": {"A": 0.2, "B": 0.8}}
TRUTH = {"A": -1.0, "B": -1.0}
FIELDS = [("Cameras alone", {"A": -1.0, "B": -1.0},
           "cameras detect equally at A and B,\nso the fitted surface is flat"),
          ("Cameras + iNaturalist", {"A": -1.4, "B": -0.6},
           "iNaturalist records place far more\nanimals at B, so the shared surface tilts")]


def _weighted(field, period):
    return sum(EFFORT[period][c] * field[c] for c in ("A", "B"))


def fig_composition_mechanism(out="fig_composition_mechanism.png"):
    fig, axes = plt.subplots(1, 2, figsize=(8.6, 3.5), sharey=True)
    OBS, EXP = "#22282e", "#2c6c9c"

    for ax, (title, field, sub) in zip(axes, FIELDS):
        x = [0, 1]
        exp = [_weighted(field, "early"), _weighted(field, "late")]
        obs = [_weighted(TRUTH, "early"), _weighted(TRUTH, "late")]
        ax.plot(x, exp, "-o", color=EXP, lw=2.0, ms=6, zorder=3)
        ax.plot(x, obs, "-s", color=OBS, lw=2.0, ms=5.5, zorder=3)
        # The absorbed quantity is the difference in SLOPES, not the level gap at
        # one endpoint. A first version annotated obs[1]-exp[1] = -0.24, which is
        # the late-year level difference and half the real answer; the expected
        # line rises 0.48 while the data are flat, so the year term supplies -0.48.
        d_exp, d_obs = exp[1] - exp[0], obs[1] - obs[0]
        absorbed = d_obs - d_exp
        if abs(absorbed) > 0.01:
            ax.annotate("", xy=(1.0, exp[1]), xytext=(1.0, exp[0]),
                        arrowprops=dict(arrowstyle="-", color=EXP, lw=0.8, ls=":"))
            ax.annotate("", xy=(1.0, obs[1]), xytext=(1.0, exp[1]),
                        arrowprops=dict(arrowstyle="<->", color="#b2182b", lw=1.4))
        ax.text(0.5, 0.03,
                f"model expects {d_exp:+.2f}\ndata show {d_obs:+.2f}\n"
                f"year trend supplies {absorbed:+.2f}",
                transform=ax.transAxes, ha="center", va="bottom", fontsize=7,
                color="#b2182b" if abs(absorbed) > 0.01 else "#5f6a73")
        ax.set_xticks(x)
        ax.set_xticklabels(["early\n(cameras mostly in A)", "late\n(cameras mostly in B)"])
        ax.set_title(title, fontsize=8.4)
        ax.margins(x=0.30, y=0.28)

    axes[0].set_ylabel("Occupancy at the sampled cells\n(log scale)")
    axes[0].text(0.03, 0.96, "expected by the model", transform=axes[0].transAxes,
                 color=EXP, fontsize=7, va="top")
    axes[0].text(0.03, 0.86, "seen in the camera data", transform=axes[0].transAxes,
                 color=OBS, fontsize=7, va="top")

    fig.text(0.5, -0.06,
             "Two cells, A and B. The species is equally common in both and does not change over "
             "time; only WHERE the cameras are changes.\nLeft: fitted on cameras alone, the surface "
             "is flat, so moving the cameras predicts no change and the year trend stays at zero.\n"
             "Right: iNaturalist reports many more animals at B, so the shared surface tilts and the "
             "model now expects detections to RISE as cameras move to B.\nThey don't, and the year "
             "trend is the only parameter that can absorb the difference - so it goes negative. "
             "No trend exists in either data stream.",
             ha="center", va="top", fontsize=6.4, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.08, dpi=300)
    plt.close(fig)
    return out, {t: round(_weighted(f, "late") - _weighted(f, "early"), 3)
                 for t, f, _ in FIELDS}
