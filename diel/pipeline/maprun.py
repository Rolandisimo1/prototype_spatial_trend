"""Stage 4: prediction surfaces from the 5-predictor mapping model.

The mapping model is not the inference model. Predator richness and own detection rate are absent
from the prediction grid: detection rate is a property of the survey rather than of the landscape,
so it cannot be mapped at all, and predator richness was never gridded. Supports are 1 km on both
sides, because the grid carries only 1 km crop and ruggedness and because site-versus-grid
disagreement is severe for skewed layers.

peak_h surfaces are built on the circular deviation from the species' precision-weighted circular
mean and mapped back with modulo-24 arithmetic. Neighbour-to-neighbour differences and the spread
across sites use circ_diff_h and circ_sd_h.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import dielpipe as dp
import dielinf as di

MAPC = ['pop_1km', 'nlcd_1k_crop', 'tcc_1km', 'rug_1km', 't_warmmonth']
SCI = {'American Black Bear': 'Ursus americanus', 'Coyote': 'Canis latrans',
       'Eastern Fox Squirrel': 'Sciurus niger', 'Eastern Gray Squirrel': 'Sciurus carolinensis',
       'Northern Raccoon': 'Procyon lotor', 'White-tailed Deer': 'Odocoileus virginianus'}


def fit_map(g, measure, cols=MAPC):
    """Weighted least squares on the five mapping predictors; returns fit plus the standardisation."""
    y = g[measure].values.astype(float)
    se = g[measure + '_se'].values.astype(float)
    w = 1.0 / np.maximum(se ** 2, 1e-9)
    ref = None
    if measure in di.CIRC:
        ref = dp.circ_mean_h(y, w)
        y = dp.circ_diff_h(y, ref)
    X = g[cols].copy().astype(float)
    X['pop_1km'] = np.log1p(X['pop_1km'])
    mu, sd = X.mean(), X.std(ddof=1).replace(0, 1)
    Z = sm.add_constant(((X - mu) / sd).values)
    m = sm.WLS(y, Z, weights=w).fit()
    return dict(model=m, mu=mu, sd=sd, ref=ref, cols=cols, ylo=float(np.nanmin(y)),
                yhi=float(np.nanmax(y)), Xlo=X.min(), Xhi=X.max(), n=len(g))


def predict_grid(fit, grid):
    X = grid[fit['cols']].copy().astype(float)
    X['pop_1km'] = np.log1p(X['pop_1km'])
    Z = sm.add_constant(((X - fit['mu']) / fit['sd']).values)
    m = fit['model']
    pred = Z @ m.params
    V = np.asarray(m.cov_params(), float)
    se = np.sqrt(np.maximum(np.einsum('ij,jk,ik->i', Z, V, Z), 0))
    inside = ((X >= fit['Xlo']) & (X <= fit['Xhi'])).all(axis=1).values
    return pred, se, inside


def neighbour_swing(grid, values, circular=False):
    """Mean absolute difference between adjacent grid cells, on the surface's own row/column index."""
    key = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(grid.r.values, grid.c.values))}
    diffs = []
    for (r, c), i in key.items():
        for dr, dc in ((0, 1), (1, 0)):
            j = key.get((r + dr, c + dc))
            if j is None:
                continue
            a, bv = values[i], values[j]
            if not (np.isfinite(a) and np.isfinite(bv)):
                continue
            diffs.append(abs(dp.circ_diff_h(a, bv)) if circular else abs(a - bv))
    return float(np.mean(diffs)) if diffs else np.nan, len(diffs)
