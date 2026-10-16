"""Does the joint fit stay anchored to the cameras, and what predicts when it doesn't?

Each species with both a joint and a camera-only fit gives one number: how far
the joint camera trend sits from the camera-only trend, in percentage points per
year. Large values mean adding iNaturalist moved the camera-anchored trend.

SUPERSEDED READING, kept so it is not re-derived. On the first 19 species this
module reported range size as the only significant predictor (rho -0.61) and
camera effort as a non-predictor. That comparison was not like-for-like: camera
effort had been correlated against surface agreement while range size was
correlated against the shift, so two different responses were compared and the
wrong predictor won.

On all 26 species, against the same response (|shift|), the ordering reverses:
detecting deployments rho -0.607 (p=0.001), camera cells -0.415 (p=0.035),
range size -0.394 (p=0.046). Restricted to the 14 fits that resolve a direction,
detecting deployments reaches -0.837 (p=0.0002) and range size -0.679.

The mechanism is precision, not geography. |shift| tracks the width of the joint
camera interval at rho +0.642 over all 26 and +0.818 among resolved fits: the
iNaturalist stream moves the camera trend furthest where the camera data
constrain it least. Range size correlates because small-range species tend to
have few detecting deployments, not because of the Scale 2 field resolution.
Porcupine and badger break the range-size story outright -- both have large
ranges (644, 699 cells), both shift heavily (21.8, 9.8 pp/yr), and both are
NOT RESOLVABLE with camera intervals over 138 pp wide.

Practical consequence: the screen is camera information, which is knowable
before fitting from detecting-deployment counts, and a wide joint camera
interval is itself the warning sign.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

FOCAL, MUTED, ALARM = "#2c6c9c", "#9aa6b1", "#b2182b"
PENDING = "#e08214"


def fig_range_size(D, pending=None, out="fig_range_size_anchoring.png"):
    """D: one row per species with range_cells, shift_joint_minus_camonly_pp.

    pending: optional {label: range_cells} for species awaiting a camera-only
    fit, drawn on the x-axis only, since they have no shift to plot yet.
    """
    d = D.dropna(subset=["range_cells", "shift_joint_minus_camonly_pp"]).copy()
    d["abs_shift"] = d.shift_joint_minus_camonly_pp.abs()
    rho, p = spearmanr(d.range_cells.values, d.abs_shift.values)

    fig, ax = plt.subplots(figsize=(7.2, 4.4))
    moved = d.abs_shift >= 5
    ax.scatter(d.loc[~moved, "range_cells"], d.loc[~moved, "abs_shift"],
               s=34, color=FOCAL, zorder=3, label="stays anchored to cameras")
    ax.scatter(d.loc[moved, "range_cells"], d.loc[moved, "abs_shift"],
               s=58, color=ALARM, zorder=4, label="joint fit moved the trend")

    # Pending species are drawn as x-axis-only marks with NO per-point labels.
    # Rotated labels at these five range sizes collide with each other and with
    # the axis title, so the names go in the caption instead.
    if pending:
        ax.scatter(list(pending.values()), [0.35] * len(pending), marker="^",
                   s=44, color=PENDING, zorder=5,
                   label="awaiting a camera-only fit")

    for _, r in d[moved | (d.range_cells < 130)].iterrows():
        ax.annotate(r.species.replace("_", " "), (r.range_cells, r.abs_shift),
                    xytext=(6, 4), textcoords="offset points", fontsize=6.2,
                    color=ALARM if r.abs_shift >= 5 else "0.35")

    ax.set_xscale("log")
    ax.set_xticks([30, 100, 300, 1000])
    ax.set_xticklabels(["30", "100", "300", "1,000"])
    ax.set_xlabel("Size of the species' modelled range (100 km cells, log scale)")
    ax.set_ylabel("How far adding iNaturalist moved\nthe camera trend (percentage points/yr)")
    ax.set_title("Species with small ranges lose their camera anchor;\n"
                 f"camera effort does not predict this (Spearman {rho:+.2f}, p={p:.3f})")
    ax.legend(frameon=False, fontsize=6.6, loc="upper right")
    ax.margins(0.06)
    fig.text(0.5, -0.04,
             "One point per species with both a joint and a camera-only fit (n=%d). "
             "Red marks a shift of 5 percentage points per year or more.\n"
             "Triangles are species awaiting a camera-only fit, placed by range size "
             "only: %s. Range size is known before fitting, so it screens in advance."
             % (len(d), ", ".join(k.replace("_", " ") for k in (pending or {}))),
             ha="center", va="top", fontsize=6.3, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, float(rho), float(p)
