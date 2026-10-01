"""Deer camera-likelihood comparison: detection/non-detection against counts.

Deer is the one fleet species where occupancy saturates -- median fitted psi
0.952, 68% of sites above 0.9 -- so the concern was that modelling camera data
as detected/not-detected discards the information that distinguishes a busy site
from a very busy one, and biases the trend. This figure answers that by fitting
the SAME bundle two ways, changing only the camera likelihood.

Both arms come from the presence-mask deer bundle, so the comparison isolates
the observation model. Numbers are per-chain posterior summaries reported by the
Hazel extraction; the conversion to percent change per calendar year uses deer's
own year_occ standard-deviation step, because the two submodels standardise
their year covariate on different distributions and raw coefficients are not
comparable across streams.
"""
import numpy as np
import matplotlib.pyplot as plt

from . import style
from .conventions import verify_text_within_box

# deer's year_occ step: one standard deviation = 3.876 calendar years
DEER_YEAR_STEP = 0.2580241

# posterior mean and 95% interval on year_beta, presence-mask deer bundle
# Occupancy-arm values are read from the extracted posterior table rather than
# retyped: the rounded 0.155 quoted in a status message gives +4.08%/yr, while
# the table's 0.154717 gives +4.07%/yr, and a figure must not disagree with the
# artifact it cites (reporting-standards rule 7).
ARMS = [
    ("Counts\n(negative binomial)", 0.132, 0.107, 0.157, "#2c6c9c"),
    ("Detection /\nnon-detection", 0.154717, 0.115982, 0.192038, "#e08214"),
]

# corrected estimator sweep, deer-like abundance: the two camera observation
# models differed by 0.002 on year_beta. The real-data gap is 0.023.
SIM_GAP, REAL_GAP = 0.002, 0.155 - 0.132


def pct_per_year(coef, step=DEER_YEAR_STEP):
    """Log-scale coefficient -> percent change per calendar year."""
    return 100.0 * (np.exp(coef * step) - 1.0)


def fig_camcount_comparison(out="fig_deer_camcount_comparison.png"):
    fig, ax = plt.subplots(figsize=(7.4, 1.95))
    y = np.arange(len(ARMS))[::-1]

    for yi, (label, mean, lo, hi, colour) in zip(y, ARMS):
        ax.plot([pct_per_year(lo), pct_per_year(hi)], [yi, yi],
                color=colour, lw=3.0, solid_capstyle="round", zorder=3)
        ax.plot(pct_per_year(mean), yi, "o", color=colour, ms=7.5,
                markeredgecolor="white", markeredgewidth=1.1, zorder=4)
        ax.text(pct_per_year(hi) + 0.12, yi,
                f"{pct_per_year(mean):+.2f}%/yr", va="center", ha="left",
                fontsize=8.4, color=colour)

    ax.axvline(0, color="0.72", lw=0.9, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([a[0] for a in ARMS], fontsize=8.4)
        # tight: two rows in a 1.95in panel, so the data fills the box (§3.5)
    ax.set_ylim(-0.48, len(ARMS) - 0.52)
    ax.set_xlim(0, 6.4)
    ax.set_xlabel("Camera-anchored trend in white-tailed deer abundance\n"
                  "(percent change per year, 95% credible interval)")
    ax.set_title("Treating camera data as counts rather than presence/absence\n"
                 "lowers deer's trend by about 15%, and the intervals overlap",
                 fontsize=9.2)

    fig.text(0.5, -0.46,
             "Both arms are the same deer bundle with only the camera likelihood changed, so the "
             "difference isolates the observation model.\nOccupancy saturates for deer (median "
             "fitted occupancy 0.952), which was the reason to check: presence/absence cannot "
             "distinguish a busy\nsite from a very busy one. It does pull the trend up, but it "
             "neither reverses the sign nor inflates the magnitude. The corrected\nsimulation "
             f"predicted a gap of {SIM_GAP:.3f} on the underlying coefficient; the real gap is "
             f"{REAL_GAP:.3f} -- same direction, ten times larger.",
             ha="center", va="top", fontsize=6.4, color="#5f6a73")

    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out
