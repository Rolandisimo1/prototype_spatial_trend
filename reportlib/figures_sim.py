"""One function per simulation figure. Each returns the output path."""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from statsmodels.stats.proportion import proportion_confint

from . import inputs, style
from .conventions import normalize_scenario, verify_null_is_null, verify_text_within


def fig_design(out="fig_sim_design.png"):
    """Schematic of the simulation design: what is generated, what is fit, what is measured.

    Estimator colours are threaded from the results figures, so the four arms
    are visually identifiable here before any numbers appear.
    """
    fig, ax = plt.subplots(figsize=(9.8, 5.0))
    ax.set_xlim(0, 103); ax.set_ylim(0, 60); ax.axis("off")

    INK, SUB, RULE = "#22282e", "#5f6a73", "#c8cfd5"
    TOP_Y, TOP_H = 26.0, 32.0

    def step(x, w, n, title, y, h):
        ax.add_patch(Rectangle((x, y), w, h, facecolor="white",
                               edgecolor=RULE, linewidth=0.9, zorder=1))
        ax.add_patch(Rectangle((x, y + h - 0.55), w, 0.55,
                               facecolor=INK, edgecolor="none", zorder=2))
        ax.add_patch(plt.Circle((x + 3.6, y + h - 4.8), 1.95, facecolor=INK,
                                edgecolor="none", zorder=3))
        ax.text(x + 3.6, y + h - 4.8, str(n), ha="center", va="center",
                color="white", fontsize=7.2, fontweight="bold", zorder=4)
        ax.text(x + 7.0, y + h - 4.8, title, ha="left", va="center",
                fontsize=8.3, fontweight="bold", color=INK, zorder=4)
        return x + 3.6, y + h - 9.6

    # --- 1. truth -----------------------------------------------------------
    tx, ty = step(1, 29, 1, "Simulate a known truth", TOP_Y, TOP_H)
    ax.text(tx, ty, "The true trend is set by us, so every\nestimate can be scored against it.",
            fontsize=6.9, color=INK, va="top", linespacing=1.5)
    for k, (lab, val, col) in enumerate([("real decline", f"{style.TVB_TRUE}", "#b2182b"),
                                          ("no trend", "0", "#4a6a4a")]):
        yy = ty - 8.2 - 3.6 * k
        ax.add_patch(Rectangle((tx, yy - 1.3), 1.05, 2.6, facecolor=col,
                               edgecolor="none", zorder=3))
        ax.text(tx + 2.2, yy, f"{lab}:  trend = {val}", fontsize=6.9,
                color=INK, va="center")
    ax.text(tx, ty - 17.0, "crossed with three abundance levels,\nfrom low to high density",
            fontsize=6.7, color=SUB, va="top", linespacing=1.5)

    # --- 2. observations ----------------------------------------------------
    ox, oy = step(33, 26, 2, "Generate observations", TOP_Y, TOP_H)
    ax.text(ox, oy, "Two data streams, as in the\nreal analysis:", fontsize=6.9,
            color=INK, va="top", linespacing=1.5)
    for k, ln in enumerate(["camera detection histories", "iNaturalist counts per cell"]):
        ax.text(ox + 0.5, oy - 8.2 - 3.4 * k, f"—  {ln}", fontsize=6.9,
                color=INK, va="center")
    ax.text(ox, oy - 17.0, "30 replicate datasets for every\nabundance x scenario combination",
            fontsize=6.7, color=SUB, va="top", linespacing=1.5)

    # --- 3. estimators, as an actual 2x2 -----------------------------------
    EST_X, EST_W = 62.0, 40.0          # box geometry, reused by the containment check
    ex, ey = step(EST_X, EST_W, 3, "Fit four estimators", TOP_Y, TOP_H)
    # Row x column labels identify each arm; the cells carry only colour, which
    # is threaded to the results figures. Naming the arm inside the cell as
    # well is redundant and forces the box wider than the layout allows.
    cw, ch = 11.0, 4.6
    col_x = [ex + 12.0 + cw / 2, ex + 12.0 + cw * 1.5 + 1.6]
    col_headers = [
        ax.text(cx, ey + 0.8, clab, fontsize=6.8, color=INK, ha="center",
                va="center", fontweight="bold")
        for cx, clab in zip(col_x, ["occupancy", "Royle-Nichols"])
    ]
    for r, (rlab, arms) in enumerate([("camera level", ["camera_occ", "camera_rn"]),
                                       ("array level", ["array_occ", "array_rn"])]):
        yy = ey - 2.8 - r * (ch + 1.5)
        ax.text(ex + 10.6, yy - ch / 2, rlab, fontsize=6.8, color=INK,
                ha="right", va="center")
        for cx, arm in zip(col_x, arms):
            ax.add_patch(Rectangle((cx - cw / 2, yy - ch), cw, ch,
                                   facecolor=style.EST_COLOR[arm],
                                   edgecolor="none", alpha=0.95, zorder=3))
    ax.text(ex, ey - 17.0,
            "Every estimator sees identical data, so any\ndifference is the method, not the dataset.",
            fontsize=6.7, color=SUB, va="top", linespacing=1.5)

    # --- 4. scoring ---------------------------------------------------------
    sx, sy = step(14, 72, 4, "Score each estimate against the known truth", 2.0, 20.0)
    scores = [("Bias", "is the number right?"),
              ("Interval width", "how precise?"),
              ("Coverage", "are the intervals honest?"),
              ("Power", "does it find a real trend?"),
              ("False-positive", "does it invent one?")]
    for k, (name, gloss) in enumerate(scores):
        xx = sx + (k % 3) * 23.5
        yy = sy - 0.6 - (k // 3) * 5.6
        ax.text(xx, yy, name, fontsize=7.0, color=INK, fontweight="bold", va="center")
        ax.text(xx, yy - 2.7, gloss, fontsize=6.5, color=SUB, va="center")

    for x0, x1 in [(30.2, 32.8), (59.2, 61.8)]:
        ax.add_patch(FancyArrowPatch((x0, TOP_Y + TOP_H / 2), (x1, TOP_Y + TOP_H / 2),
                                     arrowstyle="-|>", mutation_scale=8,
                                     color=SUB, linewidth=0.9, zorder=5))
    ax.add_patch(FancyArrowPatch((80.5, TOP_Y - 0.4), (80.5, 22.4), arrowstyle="-|>",
                                 mutation_scale=8, color=SUB, linewidth=0.9, zorder=5))
    # The estimator-matrix column headers previously overflowed this box's
    # right border. Fail loudly rather than ship it again. ex is the box's left
    # edge and the box is 40 wide; keep a 0.5-unit inset so the border line
    # itself is not touched.
    verify_text_within(ax, col_headers, EST_X + 0.5, EST_X + EST_W - 0.5,
                       what="estimator-matrix column header")

    # pad_inches is load-bearing here: the step boxes' borders sit at the
    # extreme of the drawn content, and a bare tight bbox crops through them.
    fig.savefig(out, bbox_inches="tight", pad_inches=0.10, dpi=300)
    plt.close(fig)
    return out


def fig_estimator_performance(out="fig_sim_estimators.png"):
    """Bias, interval width and coverage against abundance, by estimator.

    Uses the varying-trend scenario: the null scenario has no trend to be
    biased about.
    """
    ab = normalize_scenario(inputs.sim("abund_summary"))
    v = ab[ab.scenario == "varying"]
    x = np.arange(3)

    fig, axes = plt.subplots(1, 3, figsize=(10.2, 3.6))
    panels = [("bias", "Bias  (estimate − truth)", 0.0),
              ("ci_width", "Credible-interval width", None),
              ("coverage", "Coverage of 95% interval", 0.95)]
    for ax, (col, lab, ref) in zip(axes, panels):
        for est in style.EST_LABEL:
            d = v[v.estimator == est].set_index("abundance").loc[style.ABUND_ORDER]
            ax.plot(x, d[col], "-", color=style.EST_COLOR[est],
                    marker=style.EST_MARKER[est], ms=4.6, lw=1.5, label=style.EST_LABEL[est])
        if ref is not None:
            ax.axhline(ref, color="0.35", lw=0.7, ls="--")
            if col == "coverage":
                ax.text(-0.08, ref, "nominal 0.95", fontsize=6.2, color="0.35",
                        va="bottom", ha="left")
        ax.set_xticks(x)
        ax.set_xticklabels([style.ABUND_LABEL[a] for a in style.ABUND_ORDER])
        ax.set_ylabel(lab)
        ax.margins(x=0.10, y=0.14)

    axes[0].set_title("Occupancy bias grows with abundance;\nRoyle-Nichols stays near zero")
    axes[1].set_title("Royle-Nichols intervals are\n3–9× narrower")
    axes[2].set_title("Occupancy coverage collapses as\nabundance rises; RN's does not")
    axes[0].legend(frameon=False, fontsize=6.4, loc="lower left")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return out


def fig_error_rates(out="fig_sim_error_rates.png"):
    """Power against false-positive rate for each estimator.

    Reports the best available estimate per arm: occupancy arms from the
    corrected-null sweep (n=90), Royle-Nichols arms from the higher-replicate
    run (n=270). No pooling across runs -- dataset reuse between them makes a
    cross-arm pooled contrast correlated.
    """
    nf = normalize_scenario(inputs.sim("nullfix_rows"))
    verify_null_is_null(nf)
    bump = normalize_scenario(inputs.sim("rnnull_rows"))
    verify_null_is_null(bump)

    nul = nf[nf.scenario == "null"].groupby("estimator")["tvb_detected"].agg(["sum", "size"])
    bnul = bump[bump.scenario == "null"].groupby("estimator")["tvb_detected"].agg(["sum", "size"])
    best = {e: (bnul.loc[e] if e in bnul.index else nul.loc[e]) for e in style.EST_LABEL}
    power = nf[nf.scenario == "varying"].groupby("estimator")["tvb_detected"].mean() * 100

    fig, ax = plt.subplots(figsize=(6.8, 4.2))
    ax.axvspan(0, 5, color="#eaf1ea", zorder=0)
    ax.axvline(5, color="#4a6a4a", lw=0.9, ls="--", zorder=1)
    # Placed in the marker-free band between the array-occupancy point (y~21)
    # and the Royle-Nichols points (y~36): anywhere higher collides with the
    # camera-occupancy marker and its value label.
    ax.text(0.35, 29.5, "acceptable false-positive\nrate (5% or below)", fontsize=6.4,
            color="#41613f", ha="left", va="center", linespacing=1.5)

    for est, lab in style.EST_LABEL.items():
        s, n = float(best[est]["sum"]), int(best[est]["size"])
        pct = 100 * s / n
        lo, hi = proportion_confint(s, n, method="wilson")
        ax.plot([100 * lo, 100 * hi], [power[est]] * 2, color=style.EST_COLOR[est],
                lw=1.5, alpha=0.85, zorder=2)
        ax.plot(pct, power[est], style.EST_MARKER[est], color=style.EST_COLOR[est],
                ms=7.5, zorder=3, label=lab)
        ax.annotate(f"{pct:.1f}%", (pct, power[est]), textcoords="offset points",
                    xytext=(0, 9), fontsize=6.4, color=style.EST_COLOR[est], ha="center")

    ax.set_xlabel("False-positive rate (%)\nhow often a trend is reported when there is none")
    ax.set_ylabel("Power (%)\nhow often a real decline is detected")
    ax.set_title("All four estimators find real trends far more often\nthan they invent them")
    ax.set_xlim(0, 12); ax.set_ylim(10, 50)
    ax.legend(frameon=False, fontsize=6.6, loc="lower right")
    fig.text(0.5, -0.10,
             "Horizontal bars are Wilson 95% confidence intervals on the false-positive rate. "
             "Royle-Nichols arms n=270, occupancy arms n=90.",
             ha="center", fontsize=6.2, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return out
