"""Fleet-wide urban-rural stream screen: what each stream can and cannot resolve.

The screen fits a trend gradient along human population separately on each data
stream, per species. The point of this figure is NOT the gradients themselves --
it is that the two streams differ enormously in how precisely they measure the
same quantity, and that the camera arm is too imprecise for most species to
support any comparison at all.

Reads stream_screen/stream_screen_<species>.csv, written by the standard screen.
"""
import glob
import os
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from . import style

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCREEN = os.path.join(HERE, "stream_screen")
XLIM = 60.0   # beyond this, intervals are drawn as off-scale arrows


def _parse(s):
    """Mean, lo, hi from a '+2.09 [-6.85, +11.04]' cell; NaN for a verdict string."""
    if not isinstance(s, str):
        return (np.nan,) * 3
    m = re.findall(r"[-+]?\d+\.?\d*", s)
    return (float(m[0]), float(m[1]), float(m[2])) if len(m) >= 3 else (np.nan,) * 3


def screen_table(stratum="NATIONAL"):
    """One row per species for the requested stratum, gradients parsed out."""
    files = [f for f in sorted(glob.glob(os.path.join(SCREEN, "stream_screen_*.csv")))
             if "concentration" not in f]
    if not files:
        raise FileNotFoundError(f"no screen outputs in {SCREEN}")
    d = pd.concat([pd.read_csv(f) for f in files], ignore_index=True)
    d = d[d.stratum == stratum].copy()
    for c in ("cam_gradient", "inat_gradient"):
        d[[f"{c}_m", f"{c}_lo", f"{c}_hi"]] = pd.DataFrame(
            [_parse(v) for v in d[c]], index=d.index)
    # A degenerate fit is flagged in the text field, not by a numeric sentinel.
    d["degenerate"] = d.cam_gradient.astype(str).str.contains("DEGENERATE")
    d["assessed"] = ~d.agreement.astype(str).str.startswith("NOT ASSESSED")
    return d.sort_values("cam_ci_width")


def fig_screen(out="fig_stream_screen.png", stratum="NATIONAL"):
    """Camera and iNaturalist gradients per species, sorted by camera precision."""
    d = screen_table(stratum).dropna(subset=["cam_gradient_m"])
    y = np.arange(len(d))[::-1]
    CAM, INAT, MUTE = "#2c6c9c", "#e08214", "#b9c2c9"

    fig, ax = plt.subplots(figsize=(7.4, 5.0))
    ax.axvline(0, color="0.72", lw=0.8, zorder=1)

    for yi, (_, r) in zip(y, d.iterrows()):
        for off, lo, hi, m, col in [
            (+0.17, r.cam_gradient_lo, r.cam_gradient_hi, r.cam_gradient_m, CAM),
            (-0.17, r.inat_gradient_lo, r.inat_gradient_hi, r.inat_gradient_m, INAT),
        ]:
            if np.isnan(m):
                continue
            c = col if r.assessed else MUTE
            # Clip to the frame and mark the clip, rather than letting one
            # 257-point interval set the scale for all sixteen species.
            l, h = max(lo, -XLIM), min(hi, XLIM)
            ax.plot([l, h], [yi + off] * 2, color=c, lw=1.9,
                    solid_capstyle="butt" if (lo < -XLIM or hi > XLIM) else "round",
                    zorder=2)
            for edge, val, mark in [(-XLIM, lo, "<"), (XLIM, hi, ">")]:
                if (val < -XLIM and mark == "<") or (val > XLIM and mark == ">"):
                    ax.plot(edge, yi + off, mark, color=c, ms=3.6, zorder=3)
            ax.plot(np.clip(m, -XLIM, XLIM), yi + off, "o", color=c, ms=4.0, zorder=4)

    ax.set_yticks(y)
    ax.set_yticklabels([s.replace("_", " ").capitalize() +
                        ("" if r.assessed else "")
                        for s, r in zip(d.species, d.itertuples())], fontsize=7)
    for tick, a in zip(ax.get_yticklabels(), d.assessed):
        tick.set_color("#22282e" if a else "#6b747b")
    ax.set_xlim(-XLIM - 3, XLIM + 3)
    ax.set_xlabel("Trend gradient along human population\n"
                  "(percentage points of %/yr per +1 standardised unit)", fontsize=7.6)
    n_ass = int(d.assessed.sum())
    ratio = d.cam_ci_width / d.inat_ci_width
    ratio_lo, ratio_hi = ratio.min(), ratio.max()
    ax.set_title(f"The camera stream can resolve an urban gradient for only "
                 f"{n_ass} of {len(d)} species", fontsize=9.0, loc="left")

    ax.plot([], [], "o-", color=CAM, ms=4, lw=1.9, label="Camera traps")
    ax.plot([], [], "o-", color=INAT, ms=4, lw=1.9, label="iNaturalist")
    ax.plot([], [], "o-", color=MUTE, ms=4, lw=1.9,
            label="Comparison not assessed\n(camera interval too wide)")
    # Upper right is the only real whitespace: the top-ranked species have the
    # narrowest intervals, so nothing is drawn beyond +20 in those rows.
    ax.legend(frameon=False, fontsize=6.6, loc="upper right", labelspacing=0.8)
    ax.margins(y=0.02)

    fig.text(0.5, -0.02,
             "Species ordered by camera-interval width, narrowest at top. Arrowheads mark intervals "
             "running past the frame; red squirrel's camera\ninterval is 257 points wide. Both "
             "gradients use a mixed model with a random intercept and random year-slope per cell. "
             "Elk is omitted:\nits camera fit was degenerate. A species is assessed only where the "
             "difference between streams is resolvable, which needs the camera\narm to be precise "
             f"enough -- the iNaturalist arm is {ratio_lo:.0f} to {ratio_hi:.0f} times tighter for every species here.",
             ha="center", va="top", fontsize=6.3, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, d
