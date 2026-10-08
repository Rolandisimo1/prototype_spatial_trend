"""The 18-species trend table: the primary deliverable, and the stream pattern in it.

fig_fleet_trends is figure 1 of the report. fig_stream_compression is the
method finding that falls out of having 18 species rather than 4.

Both read trend_table/trend_table_v2u_national_scalar_18spp.csv, the extracted
camera-anchored trends. Nothing here re-derives a trend.
"""
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from . import style

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TABLE = os.path.join(HERE, "trend_table",
                     "trend_table_v2u_national_scalar_18spp.csv")

UP, DOWN, FLAT, VAGUE = "#2d7a3e", "#b2182b", "#8d969c", "#c8cfd5"
# A camera interval wider than this cannot support "no change detected" either;
# it is simply uninformative. Set from the fleet: only moose exceeds it.
VAGUE_WIDTH_PP = 8.0


def _mid_lo_hi(s):
    v = re.findall(r"[-+]?\d+\.?\d*", str(s))
    return (float(v[0]), float(v[1]), float(v[2])) if len(v) >= 3 else (np.nan,) * 3


def fleet_table():
    """The 18 species with intervals parsed and a four-way verdict per species.

    The extractor's two-way resolved/not-resolvable split puts a tight interval
    straddling zero (skunk, +-2.1) in the same bucket as moose (+-6.6). Those
    are different results: one is 'no change detected', the other is 'cannot
    tell'. Splitting them is the difference between 6 unresolved species and 1.
    """
    if not os.path.exists(TABLE):
        raise FileNotFoundError(
            f"trend table not found at {TABLE} -- it is written by the cluster "
            "extraction step (trend_table_fleet.R), not derived here.")
    d = pd.read_csv(TABLE)
    for c in ("camera_pct", "inat_pct", "gap_pp"):
        d[[f"{c}_m", f"{c}_lo", f"{c}_hi"]] = pd.DataFrame(
            [_mid_lo_hi(v) for v in d[c]], index=d.index)
    d["width"] = d.camera_pct_hi - d.camera_pct_lo
    spans0 = (d.camera_pct_lo < 0) & (d.camera_pct_hi > 0)
    d["verdict"] = np.where(~spans0 & (d.camera_pct_m > 0), "increase",
                   np.where(~spans0, "decline",
                   np.where(d.width > VAGUE_WIDTH_PP, "not resolvable",
                            "no change detected")))
    d["label"] = (d.species.str.replace("_", " ")
                   .str.replace("white-tailed deer", "White-tailed deer")
                   .str.capitalize())
    return d.sort_values("camera_pct_m", ascending=False).reset_index(drop=True)


def fig_fleet_trends(out="fig_fleet_trends.png"):
    """Camera-anchored population trend for each of the 18 species."""
    d = fleet_table()
    col = {"increase": UP, "decline": DOWN,
           "no change detected": FLAT, "not resolvable": VAGUE}
    y = np.arange(len(d))[::-1]

    fig, ax = plt.subplots(figsize=(7.0, 5.4))
    ax.axvline(0, color="0.70", lw=0.9, zorder=1)
    for yi, (_, r) in zip(y, d.iterrows()):
        c = col[r.verdict]
        ax.plot([r.camera_pct_lo, r.camera_pct_hi], [yi] * 2, color=c, lw=2.3,
                solid_capstyle="round", zorder=2)
        ax.plot(r.camera_pct_m, yi, "o", color=c, ms=5.0,
                mec="white", mew=0.7, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(d.label, fontsize=7.6)
    for tick, v in zip(ax.get_yticklabels(), d.verdict):
        tick.set_color("#22282e" if v in ("increase", "decline") else "#6b747b")
    ax.set_xlabel("Population trend, 2008-2025 (% per year)", fontsize=8.2)

    n = d.verdict.value_counts()
    ax.set_title(f"{n.get('increase', 0)} of 18 species increasing, "
                 f"{n.get('decline', 0)} declining",
                 fontsize=9.6, loc="left")
    for lab, c in [("Increasing", UP), ("Declining", DOWN),
                   ("No change detected", FLAT), ("Not resolvable", VAGUE)]:
        ax.plot([], [], "o-", color=c, ms=5, lw=2.3, label=lab)
    ax.legend(frameon=False, fontsize=6.8, loc="lower right", labelspacing=0.75)
    ax.margins(y=0.02)

    fig.text(0.5, -0.015,
             "Camera-anchored trend (year_beta): what the camera-trap surveys support, with "
             "iNaturalist records contributing spatial\ncoverage and a separate term absorbing "
             "iNaturalist-specific temporal signal. Bars are 95% credible intervals. All 18 fits "
             "converged\n(max R-hat 1.010). 'No change detected' means the interval straddles zero "
             f"but is narrower than {VAGUE_WIDTH_PP:.0f} points; only moose is\nwide enough that no "
             "statement is possible. Union range mask; single continuous chains.",
             ha="center", va="top", fontsize=6.4, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, d


def deming(x, y, lam, n_boot=4000, seed=0):
    """Errors-in-variables slope and bootstrap CI.

    Ordinary least squares is not usable here: BOTH trends are posterior
    estimates carrying real uncertainty, and error in x biases an OLS slope
    toward zero -- which is the very effect being tested. lam is the ratio of
    error variances (x over y), taken from the fitted intervals.
    """
    def _slope(a, b):
        sxx, syy = np.var(a, ddof=1), np.var(b, ddof=1)
        sxy = np.cov(a, b, ddof=1)[0, 1]
        return ((syy - lam * sxx) + np.sqrt((syy - lam * sxx) ** 2
                + 4 * lam * sxy ** 2)) / (2 * sxy)
    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(n_boot):
        i = rng.integers(0, len(x), len(x))
        if np.cov(x[i], y[i], ddof=1)[0, 1] != 0:
            bs.append(_slope(x[i], y[i]))
    return _slope(x, y), np.percentile(bs, 2.5), np.percentile(bs, 97.5)


def fig_stream_compression(out="fig_stream_compression.png"):
    """iNaturalist trend against camera trend across the 18 species.

    Only visible with a fleet: at four species there was no slope to estimate.
    """
    d = fleet_table()
    x, y = d.camera_pct_m.values, d.inat_pct_m.values
    cam_se = (d.camera_pct_hi - d.camera_pct_lo) / 3.92
    inat_se = (d.inat_pct_hi - d.inat_pct_lo) / 3.92
    lam = (cam_se ** 2).mean() / (inat_se ** 2).mean()
    b, blo, bhi = deming(x, y, lam)

    fig, ax = plt.subplots(figsize=(5.9, 5.3))
    lim = [min(x.min(), y.min()) - 2, max(x.max(), y.max()) + 2]
    ax.plot(lim, lim, "-", color="0.62", lw=1.0, zorder=1)
    ax.axhline(0, color="0.86", lw=0.7, zorder=0)
    ax.axvline(0, color="0.86", lw=0.7, zorder=0)
    xs = np.array(lim)
    ax.plot(xs, b * xs + (y.mean() - b * x.mean()), "-", color="#2c6c9c", lw=1.7, zorder=2)

    for xi, yi, xe, ye, lab, v in zip(x, y, cam_se, inat_se, d.label, d.verdict):
        c = {"increase": UP, "decline": DOWN}.get(v, FLAT)
        ax.errorbar(xi, yi, xerr=1.96 * xe, yerr=1.96 * ye, fmt="o", ms=4.2,
                    color=c, ecolor=c, elinewidth=0.9, alpha=0.85, zorder=3)
    # Label only where the streams disagree in sign; labelling all 18 is noise.
    # Fixed per-species offsets: raccoon and black bear sit ~1.5 points apart and
    # a uniform offset overprints them. Six labels, so explicit beats a solver.
    OFF = {"Northern raccoon": (0, 9), "American black bear": (7, 4),
           "Coyote": (7, -2), "Bobcat": (7, 3), "Red fox": (7, -8),
           "Mule deer": (7, 2)}
    for xi, yi, lab in zip(x, y, d.label):
        if np.sign(xi) != np.sign(yi) and abs(xi) > 1:
            dx, dy = OFF.get(lab, (6, -3))
            ax.annotate(lab, (xi, yi), textcoords="offset points", xytext=(dx, dy),
                        fontsize=6.4, color="#22282e",
                        ha="center" if dx == 0 else ("right" if dx < 0 else "left"))
    ax.set_xlim(lim); ax.set_ylim(lim); ax.set_aspect("equal")
    ax.set_xlabel("Camera-trap trend (% per year)", fontsize=8.2)
    ax.set_ylabel("iNaturalist trend (% per year)", fontsize=8.2)
    ax.set_title("iNaturalist recovers about a quarter of the\ntrend the cameras measure",
                 fontsize=9.4, loc="left")
    ax.annotate("1:1", (lim[1] - 1.5, lim[1] - 1.5), fontsize=6.6, color="#5f6a73",
                ha="right", va="bottom", rotation=45)
    # Lower-right quadrant is empty (no species has camera-up/iNat-far-down at
    # that magnitude); on the line it overprinted the raccoon/bear cluster.
    ax.text(0.97, 0.06, f"Deming slope {b:.2f}  [{blo:.2f}, {bhi:.2f}]",
            transform=ax.transAxes, fontsize=7.0, color="#2c6c9c", ha="right")

    fig.text(0.5, -0.01,
             "Each point is one species; bars are 95% credible intervals on both axes. The fitted "
             "line is a Deming (errors-in-variables)\nregression, which accounts for both trends "
             "being estimates -- ordinary least squares would bias the slope toward zero, the\n"
             "very effect under test, and gives 0.27 against Deming's "
             f"{b:.2f}. The interval excludes 1, so the compression is not an\nartifact of "
             "measurement error. Named species are those where the two streams disagree in sign.",
             ha="center", va="top", fontsize=6.4, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, (b, blo, bhi)
