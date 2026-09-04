"""Multivariate inference on site-level diel measures: the PI's fixed 7-predictor model.

Each term maps to one distinct mechanism, so a coefficient is a partial association holding the
other six fixed. Artificial light and developed-land share are excluded from the study by PI
decision and appear nowhere in this module.

peak_h is circular. It enters as the signed circular deviation (in hours, range -12..12) from the
species' own precision-weighted circular mean, computed with circ_diff_h / circ_mean_h. Every
difference, mean and interval taken on peak_h anywhere in this module goes through those two
functions; no linear subtraction of two clock hours occurs.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import dielpipe as dp

MECH = {
    'pop_1km': 'human disturbance',
    'nlcd_5k_crop': 'agriculture',
    'tcc_1km': 'vegetation',
    'rug_5km': 'terrain',
    'tmax_hottest_month': 'summer heat',
    'pred_richness': 'predator richness',
    'log_det_rate': 'own detection rate',
}
COLS = list(MECH.keys())
MAP_MECH = {
    'pop_1km': 'human disturbance',
    'nlcd_1k_crop': 'agriculture',
    'tcc_1km': 'vegetation',
    'rug_1km': 'terrain',
    't_warmmonth': 'summer heat',
}
MAP_COLS = list(MAP_MECH.keys())
MEASURES = ['pct_noct', 'pct_crep', 'pct_noon', 'conc', 'peak_h']
CIRC = {'peak_h'}
PRETTY = {
    'pct_noct': 'percent of activity at night',
    'pct_crep': 'percent of activity at dawn and dusk',
    'pct_noon': 'percent of activity around midday',
    'conc': 'concentration of activity into part of the day',
    'peak_h': 'hour of peak activity',
}


def response(g, measure, ref=None):
    """Response vector and inverse-variance weights for one species and one measure.

    For peak_h the response is the circular deviation from `ref` (default: the precision-weighted
    circular mean of this species), so the regression is on a linear quantity that is correct
    across midnight. `ref` is returned so the same centre can be reused for prediction.
    """
    y = g[measure].values.astype(float)
    se = g[measure + '_se'].values.astype(float)
    w = 1.0 / np.maximum(se ** 2, 1e-9)
    if measure in CIRC:
        if ref is None:
            ref = dp.circ_mean_h(y, w)
        y = dp.circ_diff_h(y, ref)
    return y, w, ref


def design(g, cols, logpop=True):
    X = g[cols].copy().astype(float)
    if logpop and 'pop_1km' in X.columns:
        X['pop_1km'] = np.log1p(X['pop_1km'])
    mu, sd = X.mean(), X.std(ddof=1).replace(0, 1)
    return (X - mu) / sd, mu, sd


def _farthest_point_centres(xk, yk, k, seed=0):
    """Space-filling centre selection: deterministic, and it does not cluster centres where the
    sites happen to be dense, so the spatial basis has comparable reach everywhere."""
    n = len(xk)
    k = int(min(k, n))
    P = np.column_stack([xk, yk])
    start = int(np.argmin(np.hypot(xk - xk.mean(), yk - yk.mean())))
    idx = [start]
    d = np.hypot(P[:, 0] - P[start, 0], P[:, 1] - P[start, 1])
    while len(idx) < k:
        j = int(np.argmax(d))
        idx.append(j)
        d = np.minimum(d, np.hypot(P[:, 0] - P[j, 0], P[:, 1] - P[j, 1]))
    return np.array(idx)


def spatial_basis(xk, yk, k=None, fit_xk=None, fit_yk=None):
    """Low-rank radial-basis spatial smooth plus a linear north-south / east-west trend.

    Bandwidth is the median nearest-neighbour distance among centres, floored at 150 km, matching
    the convention already used elsewhere in this project.
    """
    if fit_xk is None:
        fit_xk, fit_yk = xk, yk
    n = len(fit_xk)
    if k is None:
        k = int(np.clip(n // 10, 4, 15))
    idx = _farthest_point_centres(np.asarray(fit_xk, float), np.asarray(fit_yk, float), k)
    cx, cy = np.asarray(fit_xk, float)[idx], np.asarray(fit_yk, float)[idx]
    D = np.hypot(cx[:, None] - cx[None, :], cy[:, None] - cy[None, :])
    np.fill_diagonal(D, np.inf)
    scale = max(float(np.median(D.min(axis=1))), 150.0) if len(cx) > 1 else 150.0
    d2 = (np.asarray(xk, float)[:, None] - cx[None, :]) ** 2 + \
         (np.asarray(yk, float)[:, None] - cy[None, :]) ** 2
    B = np.column_stack([np.asarray(xk, float) / 1000, np.asarray(yk, float) / 1000,
                         np.exp(-d2 / (2 * scale ** 2))])
    return B, dict(k=len(cx), scale_km=scale)


def wls_partials(y, Xs, w, extra=None, names=None):
    """Weighted least squares with an intercept, returning per-term inference and partial R2.

    partial R2 for term j is the share of the residual weighted sum of squares that adding term j
    removes: (SSE_without_j - SSE_full) / SSE_without_j, refitting without the term rather than
    reading it off a t statistic, so it is exact for the weighted fit.
    """
    Xs = np.asarray(Xs, float)
    A = Xs if extra is None else np.column_stack([Xs, extra])
    Z = sm.add_constant(A)
    m = sm.WLS(y, Z, weights=w).fit()
    sse_full = float(np.sum(w * m.resid ** 2))
    ci = m.conf_int()
    p = Xs.shape[1]
    out = []
    for j in range(p):
        keep = [c for c in range(A.shape[1]) if c != j]
        Zr = sm.add_constant(A[:, keep]) if keep else np.ones((len(y), 1))
        mr = sm.WLS(y, Zr, weights=w).fit()
        sse_r = float(np.sum(w * mr.resid ** 2))
        pr2 = (sse_r - sse_full) / sse_r if sse_r > 0 else np.nan
        out.append(dict(predictor=(names[j] if names else j), beta=float(m.params[j + 1]),
                        se=float(m.bse[j + 1]), p=float(m.pvalues[j + 1]),
                        lo=float(ci[j + 1][0]), hi=float(ci[j + 1][1]),
                        partial_r2=float(pr2)))
    return out, dict(r2=float(m.rsquared), r2_adj=float(m.rsquared_adj), n=int(len(y)),
                     df_model=int(m.df_model))
