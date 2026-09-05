"""Conventions that have caused a real bug in this project, encoded as checks.

Each function here exists because something silently produced a wrong result.
They are cheap to call and they fail loudly. Call them from every figure
function that depends on the convention.
"""
import numpy as np
import pandas as pd


def normalize_scenario(df):
    """The sweep writes scenario='null' as an empty field, which pandas reads as NaN.

    Exists because a groupby on the raw column silently drops every null-scenario
    row, so a false-positive rate computed without this returns the power column
    instead. Returns a copy with 'null' filled in.
    """
    if "scenario" not in df.columns:
        raise AssertionError("no 'scenario' column; wrong sweep file?")
    out = df.copy()
    out["scenario"] = out["scenario"].fillna("null")
    bad = set(out["scenario"].unique()) - {"null", "varying"}
    if bad:
        raise AssertionError(f"unexpected scenario values: {bad}")
    return out


def verify_null_is_null(df, truth_col="tvb_true"):
    """Under scenario='null' the generating trend must be exactly 0.

    Exists because the sweep's original 'null' zeroed only the spatial deviation
    amplitude and left the national trend in place, so what was reported as a
    false-positive rate was actually power under a spatially-flat truth. The
    fix is guarded on the cluster side too; this is the local mirror of it.
    """
    df = normalize_scenario(df)
    nul = df[df["scenario"] == "null"]
    if len(nul) == 0:
        return
    vals = nul[truth_col].unique()
    if not (len(vals) == 1 and float(vals[0]) == 0.0):
        raise AssertionError(
            f"null scenario has {truth_col}={vals}, expected exactly [0.0]. "
            "This is the null-scenario defect; false-positive rates are invalid."
        )


def verify_snapshot_years(df, requested, year_col="year"):
    """Snapshot extraction drops out-of-range years silently.

    Exists because the extractor filters SNAP <- SNAP[SNAP %in% years_all], so
    asking a 5-year windowed fit for 2015/2020/2025 returns one panel and no
    warning. A figure built on that looks like a successful three-year series.
    """
    got = sorted(df[year_col].unique())
    want = sorted(requested)
    if got != want:
        raise AssertionError(
            f"snapshot years {got} != requested {want}. The extractor drops "
            "out-of-range years silently; check the fit's fitted span."
        )


def verify_full_cell_coverage(df, year_col="year", cell_col="cell50"):
    """Every snapshot year must cover the same cell set.

    Exists because a partial extraction yields year panels with different
    footprints, which renders as apparent range contraction rather than as the
    missing data it is.
    """
    per = df.groupby(year_col)[cell_col].nunique()
    if per.nunique() != 1:
        raise AssertionError(f"uneven cell coverage across years: {per.to_dict()}")


# Derived indicators are probabilities pinned near 0/1, where R-hat is not
# interpretable: between-chain variance in the last decimals is compared against
# a near-zero within-chain variance. bobcat_v2b_national_scalar reports R-hat
# 1.1145 on trend_robust_indicator while every sampled trend parameter in the
# same fit is <= 1.0006. Exclude derived quantities from convergence summaries.
DERIVED_PARAMS = ("trend_robust_indicator", "snr_derived")

TREND_PARAMS = (
    "total_var_beta", "year_beta", "year_var", "sigma_region",
) + tuple(f"year_region[{i}]" for i in range(1, 9))


def max_rhat_sampled(df, param_col="parameter", rhat_col="rhat"):
    """Max R-hat over sampled trend parameters only, excluding derived indicators."""
    d = df[df[param_col].isin(TREND_PARAMS)]
    if len(d) == 0:
        raise AssertionError(f"no trend parameters found in {sorted(df[param_col].unique())}")
    return float(d[rhat_col].max())


def verify_text_within(ax, artists, x0, x1, what="text"):
    """Text artists must sit inside the data-x span [x0, x1].

    Exists because the estimator-matrix column header in fig_design() silently
    overflowed its box border: the header ended at px 2308 while the border was
    at px 2287. It was missed on a first check because text extents were
    measured under matplotlib defaults rather than the report's own font
    settings, which underestimated the width. Measure in the live figure, with
    the real rcParams, and fail rather than ship a figure whose label crosses a
    panel border.
    """
    fig = ax.get_figure()
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    bad = []
    for t in artists:
        bb = t.get_window_extent(renderer).transformed(inv)
        if bb.x0 < x0 or bb.x1 > x1:
            bad.append((t.get_text(), round(bb.x0, 2), round(bb.x1, 2)))
    if bad:
        raise AssertionError(
            f"{what} outside allowed span [{x0}, {x1}]: {bad}. "
            "Widen the container or shorten the label."
        )
