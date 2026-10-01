"""Camera-anchored population trends, four species, union mask.

The union mask keeps a 100 km cell if EITHER the species' IUCN range polygon or
its iNaturalist records place it there. This is the reported design and the only
one shown: a records-only mask discards true-absence cells -- places inside the
range where a camera ran and the species was not detected -- and discards their
camera sites with them, which is why it was replaced. The comparison against it
is methods history and is deliberately not plotted.

Values are percent change per calendar year, converted from each fit's own
year_occ standard-deviation step, because the step differs between species and
raw coefficients are not on a common scale.
"""
import numpy as np
import matplotlib.pyplot as plt

from . import style

FOCAL, MUTED = "#2c6c9c", "#9aa0a6"

# species -> (pct/yr, lo, hi, P(increase), camera sites, detecting sites)
TRENDS = [
    ("Moose",             -3.6, -10.4,  2.9, 0.15,  3044,   282),
    ("Bobcat",             1.0,  -0.5,  2.7, 0.89, 23471,  3184),
    ("Coyote",             2.3,   1.1,  3.4, 1.00, 24463,  8118),
    ("White-tailed deer",  4.0,   3.0,  5.0, 1.00, 22238, 17764),
]


def fig_union_trend(out="fig_union_trend.png"):
    """Trend per species. Colour marks whether the interval excludes zero."""
    fig, ax = plt.subplots(figsize=(7.4, 3.0))
    y = np.arange(len(TRENDS))[::-1]

    for yi, (label, mean, lo, hi, p, nsite, ndet) in zip(y, TRENDS):
        clear = (lo > 0) or (hi < 0)
        c = FOCAL if clear else MUTED
        ax.plot([lo, hi], [yi]*2, color=c, lw=3.0, solid_capstyle="round", zorder=3)
        ax.plot(mean, yi, "o", color=c, ms=7.0, markeredgecolor="white",
                markeredgewidth=1.1, zorder=4)
        ax.text(hi + 0.40, yi, f"{mean:+.1f}%/yr", va="center", ha="left",
                fontsize=8.0, color=c)

    ax.axvline(0, color="0.45", lw=1.0, zorder=2)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{t[0]}\n{t[6]:,} detecting sites" for t in TRENDS],
                       fontsize=8.0)
    ax.set_ylim(-0.55, len(TRENDS) - 0.45)
    ax.set_xlim(-12.0, 7.6)
    ax.set_xlabel("Population trend the camera surveys support\n"
                  "(percent change per year, 95% credible interval)")
    ax.set_title("Deer and coyote are increasing; bobcat and moose are not\n"
                 "distinguishable from no change", fontsize=9.2)

    fig.text(0.5, -0.155,
             "Blue where the 95% interval excludes zero, grey where it does not. This is the trend "
             "the camera surveys support on their own, which is the\nquantity we report: the "
             "iNaturalist records extend coverage but their survey effort is unmeasured, so they "
             "corroborate rather than carry the\nresult. Moose's interval is wide because only 282 "
             "of its 3,044 camera deployments ever detected a moose. Modelled area is each species' "
             "range\npolygon combined with its iNaturalist records, so cells where a camera ran and "
             "found nothing still inform the fit.",
             ha="center", va="top", fontsize=6.4, color="#5f6a73")

    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out
