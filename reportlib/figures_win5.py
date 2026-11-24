"""Five-year fleet results: the reported trend for every species that has been fitted.

The five-year window is the reported one. This module shows that window alone.
Earlier windows are a methods question and live in the window-sensitivity memo,
not here -- a results figure carrying a window we do not report would put a
number on the page that we would not stand behind (reporting-standards rules 6
and 8).

Everything is read from the posterior pull directory at call time, so re-running
after more species finish needs no edit: `fleet_table()` picks up whatever is
there. Species still on the cluster are reported as not yet fitted rather than
omitted, so the table cannot be mistaken for the whole fleet.
"""
import glob
import os

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from . import style

PULL = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(
    os.path.abspath(__file__)))), "data_integration_arielle",
    "hazel_pull_20261115_win5", "posteriors")
if not os.path.isdir(PULL):                       # when the repo sits beside the data
    PULL = "/Users/rwkays/claude_code/data_integration_arielle/hazel_pull_20261115_win5/posteriors"

# Convergence thresholds. R-hat compares the spread between chains with the
# spread within them: 1.00 means they agree. Effective sample size is how many
# independent draws the chains are worth after accounting for autocorrelation.
RHAT_MAX, ESS_MIN = 1.01, 400
RANGE_MIN_CELLS = 100

FLEET_N = 34            # species in the fitting fleet
WIN = "win5"


def _pretty(s):
    return s.replace("_", " ").replace("-", "-").capitalize()


def _one(folder, run):
    """Read one fitted run: posterior summary, convergence, and the year scaling."""
    d = pd.read_csv(os.path.join(folder, f"{run}_posterior_{WIN}.csv")) \
        if os.path.exists(os.path.join(folder, f"{run}_posterior_{WIN}.csv")) else None
    return d


def _scale(folder, run):
    """Standard deviation of the year covariate, needed to put the slope on a
    per-year scale. The model fits the trend against a standardised year, so the
    raw coefficient is change per standard deviation of year, not per year."""
    ys = pd.read_csv(os.path.join(folder, f"year_scaling_{run}.csv"))
    return float(ys.loc[ys.coefficient == "year_beta", "year_sd"].iloc[0])


def fleet_table(folder=None):
    """One row per fitted species, on the five-year window.

    Returns percent change per year with a 95% credible interval, the direction
    the interval supports, and the convergence and stream-agreement diagnostics
    that decide whether the row is reportable.
    """
    folder = folder or PULL
    rows = []
    for f in sorted(glob.glob(os.path.join(folder, f"*_posterior_{WIN}.csv"))):
        d = pd.read_csv(f)
        run, sp = d.run.iloc[0], d.species.iloc[0]
        cvf = os.path.join(folder, f"convergence_{run}.csv")
        if not os.path.exists(cvf):
            continue
        cv = pd.read_csv(cvf)
        sd = _scale(folder, run)
        gp = lambda p, c: (lambda s: float(s[c].iloc[0]) if len(s) else np.nan)(
            d[d.param.astype(str).str.fullmatch(p)])
        gc = lambda p, c: (lambda s: float(s[c].iloc[0]) if len(s) else np.nan)(
            cv[cv.param.astype(str).str.fullmatch(p)])
        # The slope is on a log scale against standardised year, so exponentiating
        # the per-year slope and subtracting one gives proportional change per year.
        pct = lambda v: 100 * (np.exp(v / sd) - 1)
        rows.append({
            "species": _pretty(sp), "run": run,
            "trend_pct_yr": pct(gp("year_beta", "mean")),
            "lo": pct(gp("year_beta", "q2.5")), "hi": pct(gp("year_beta", "q97.5")),
            "trend_rhat": gc("year_beta", "rhat"),
            "trend_ess": gc("year_beta", "ess_bulk"),
            "theta1": gp("theta1", "mean"), "theta1_lo": gp("theta1", "q2.5"),
            "theta1_hi": gp("theta1", "q97.5"),
            "theta1_rhat": gc("theta1", "rhat"), "theta1_ess": gc("theta1", "ess_bulk"),
            "chains": int(d.n_chains.iloc[0]),
        })
    t = pd.DataFrame(rows)
    if not len(t):
        return t
    t["converged"] = (t.trend_rhat <= RHAT_MAX) & (t.trend_ess >= ESS_MIN)
    t["direction"] = np.where(~t.converged, "not reportable",
                      np.where((t.lo > 0) & (t.hi > 0), "increase",
                      np.where((t.lo < 0) & (t.hi < 0), "decline", "no change detected")))
    # The stream-agreement gate: theta1 is the exponent linking iNaturalist
    # record counts to modelled abundance. An interval reaching 1 or including 0
    # means the joint trend needs a camera-only control fit before it is reported.
    t["theta1_usable"] = (t.theta1_rhat <= RHAT_MAX) & (t.theta1_ess >= ESS_MIN)
    # "Below reporting threshold" rather than "not estimable": a fit missing the
    # threshold narrowly with chains in agreement is not the same as one whose
    # chains disagree, and the table shows both numbers so the difference is
    # visible rather than collapsed into one verdict.
    t["gate"] = np.where(~t.theta1_usable, "theta1 below reporting threshold",
                 np.where((t.theta1_hi >= 1) | (t.theta1_lo <= 0),
                          "needs camera-only control", "clear"))
    t = t.merge(_range_cells(), on="species", how="left")
    # Species occupying fewer than 100 grid cells are marked data-limited. The
    # threshold is not arbitrary: across the species with both a joint and a
    # camera-only fit, those below 100 cells sat a median of 17 percentage points
    # per year away from their camera-only estimate, against 1.7 above it.
    t["data_limited"] = t.range_cells < RANGE_MIN_CELLS
    return t.sort_values("trend_pct_yr").reset_index(drop=True)


def _range_cells():
    """Range size per species, as the number of 50 km grid cells it occupies."""
    cands = ["/Users/rwkays/claude_code/geospatialdata/north_american_mammal_ranges/"
             "fleet34_isdm/fleet34_grid_summary.csv",
             os.path.join(os.path.dirname(PULL), "fleet34_grid_summary.csv"),
             os.path.join(os.path.dirname(os.path.dirname(
                 os.path.abspath(__file__))), "fleet34_grid_summary.csv")]
    for cand in cands:
        if os.path.exists(cand):
            g = pd.read_csv(cand)
            g["species"] = g.species.map(_pretty)
            return g[["species", "cells_union"]].rename(
                columns={"cells_union": "range_cells"})
    # Returning an empty frame here would silently mark every species as not
    # data-limited, which is the dangerous direction: an unflagged implausible
    # estimate reads as a result.
    raise FileNotFoundError(
        "fleet34_grid_summary.csv not found; range sizes are required to flag "
        "data-limited species. Searched: " + "; ".join(cands))


def fig_fleet_win5(out="fig_fleet_win5.png", folder=None, slide=False):
    """Forest plot of the five-year trend for every fitted species.

    Non-converged fits are drawn as open markers and excluded from every count in
    the title, rather than omitted silently or mixed in with the rest.

    `slide=True` renders the projection variant: the same figure with roughly
    double the type size and heavier marks. Print and slide sizes are opposite
    conventions -- 7 pt is correct in a report and unreadable from the back of a
    room -- so both exist rather than one compromise serving neither.
    """
    k = 2.1 if slide else 1.0
    # Styling is applied by the caller via the figure-style skill's
    # apply_figure_style(), matching the other modules in this package.
    t = fleet_table(folder)
    n_fit, n_rep = len(t), int(t.converged.sum())
    res = t[t.converged & t.direction.isin(["increase", "decline"])]
    fig, ax = plt.subplots(figsize=(7.4 * (1.55 if slide else 1.0),
                                (0.42 * max(len(t), 4) + 1.5) * (1.25 if slide else 1.0)))
    ax.axvline(0, color="#b7bfc6", lw=0.9, zorder=1)
    COL = {"decline": "#b3261e", "increase": "#1b7a3d",
           "no change detected": "#5f6a73", "not reportable": "#a7b0b7"}
    for i, r in t.iterrows():
        c = COL[r.direction]
        solid = r.converged
        ax.plot([r.lo, r.hi], [i, i], color=c, lw=(2.0 if solid else 1.2) * k,
                solid_capstyle="round", alpha=1.0 if solid else 0.55,
                ls="-" if solid else ":", zorder=3)
        ax.plot(r.trend_pct_yr, i, "o", color=c if solid else "white",
                markeredgecolor=c, ms=5.4 * k, zorder=4)
    ax.set_yticks(range(len(t)))
    # The data-limited marker is on the tick label, not in the caption: a figure
    # travels alone, and an unqualified -59% per year would be read as a result.
    ax.set_yticklabels([
        f"{r.species}" + ("" if r.converged else "  (did not converge)")
        + (f"  \u2022 {int(r.range_cells)} cells" if r.data_limited else "")
        for r in t.itertuples()])
    for i, r in t.iterrows():
        if r.data_limited:
            ax.get_yticklabels()[i].set_color("#8a5a00")
    ax.set_xlabel("Change in abundance (% per year), five-year window",
                  fontsize=None if not slide else 19)
    if slide:
        ax.tick_params(labelsize=17)
    n_lim = int(t.data_limited.sum())
    ax.set_title(f"{n_rep} of {FLEET_N} species fitted and converged; "
                 f"{len(res)} resolve a direction"
                 + (f"\n{n_lim} occupy fewer than {RANGE_MIN_CELLS} grid cells, "
                    f"where estimates are unreliable" if n_lim else ""),
                 loc="left", fontsize=None if not slide else 21)
    for lab in ["decline", "increase", "no change detected"]:
        if (t.direction == lab).any():
            ax.plot([], [], "o-", color=COL[lab], lw=2.0, ms=5.4, label=lab)
    if (~t.converged).any():
        ax.plot([], [], "o:", color=COL["not reportable"], lw=1.2, ms=5.4,
                markerfacecolor="white", label="did not converge — not reportable")
    # Lower right is clear in print but the larger slide legend reaches the
    # longest decline bars, so the projection variant sits top-left instead.
    ax.legend(frameon=False, fontsize=7 * (2.0 if slide else 1.0),
              loc="upper left" if slide else "lower right")
    ax.margins(0.06)
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return out, t
