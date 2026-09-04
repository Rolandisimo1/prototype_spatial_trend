"""Variograms, the three-way variance partition, and the permutation test for regional structure.

Partition logic. The observed between-site variance of a diel metric contains three parts:
  measurement error  - the delta-method sampling variance of each site's own fitted metric;
  local variation    - real site-to-site differences at separations BELOW the 25 km site radius,
                       estimated from the short-lag variogram intercept above the measurement floor;
  regional structure - real variation that is spatially organised at scales ABOVE 25 km.
The regional share is the fraction of real (measurement-corrected) variance that the variogram shows
accruing between the shortest lag and the sill.

Peak time is CIRCULAR: its semivariance uses squared circular differences in hours, never plain
subtraction (rule 7).
"""
import numpy as np
import pandas as pd
import dielpipe as dp

CIRCULAR = {'peak_h', 'mean_h'}
EDGES = np.array([0, 25, 50, 100, 150, 200, 300, 400, 600, 800, 1200, 1800, 2600, 4000], float)


def pair_diffs(x, y, v, circular=False):
    """All pairwise separations (km) and squared differences."""
    n = len(v)
    i, j = np.triu_indices(n, 1)
    h = np.hypot(x[i] - x[j], y[i] - y[j])
    if circular:
        d = dp.circ_diff_h(v[i], v[j])
    else:
        d = v[i] - v[j]
    return h, d ** 2, i, j


def variogram(x, y, v, circular=False, edges=EDGES, min_pairs=30):
    h, d2, _, _ = pair_diffs(x, y, v, circular)
    rows = []
    for a, b in zip(edges[:-1], edges[1:]):
        m = (h >= a) & (h < b)
        if m.sum() < min_pairs:
            continue
        rows.append(dict(lo=a, hi=b, npairs=int(m.sum()), h_mean=float(h[m].mean()),
                         gamma=float(0.5 * d2[m].mean()),
                         gamma_robust=float(0.5 * np.median(d2[m]) / 0.455)))
    return pd.DataFrame(rows)


def partition(v, se, x, y, circular=False, short_km=25.0, sill_frac=0.9):
    """Split observed between-site variance into measurement / local / regional."""
    ok = np.isfinite(v) & np.isfinite(se)
    v, se, x, y = v[ok], se[ok], x[ok], y[ok]
    n = len(v)
    if n < 12:
        return None
    if circular:
        # circular variance expressed in h^2 via the mean resultant length
        a = v / 24 * 2 * np.pi
        R = np.hypot(np.cos(a).mean(), np.sin(a).mean())
        var_obs = (-2 * np.log(max(R, 1e-12))) * (24 / (2 * np.pi)) ** 2
    else:
        var_obs = float(np.var(v, ddof=1))
    var_meas = float(np.mean(se ** 2))
    var_true = max(var_obs - var_meas, 0.0)

    h, d2, _, _ = pair_diffs(x, y, v, circular)
    short = h < short_km
    npairs_short = int(short.sum())
    # semivariance at short lag, minus the measurement floor, is LOCAL real variance
    g_short = float(0.5 * d2[short].mean()) if npairs_short >= 10 else np.nan
    local = max((g_short - var_meas), 0.0) if np.isfinite(g_short) else np.nan
    local = min(local, var_true) if np.isfinite(local) else np.nan
    regional = max(var_true - local, 0.0) if np.isfinite(local) else np.nan
    denom = var_obs if var_obs > 0 else np.nan
    return dict(n=n, var_obs=var_obs, var_meas=var_meas, var_true=var_true,
                var_local=local, var_regional=regional,
                frac_meas=var_meas / denom, frac_local=local / denom,
                frac_regional=regional / denom,
                frac_regional_of_true=regional / var_true if var_true > 0 else np.nan,
                short_gamma=g_short, npairs_short=npairs_short)


def perm_test_regional(v, x, y, circular=False, short_km=25.0, nperm=999, seed=0):
    """Is similarity at short lag greater than under random reassignment of sites to locations?

    Statistic: mean semivariance among pairs closer than short_km. Under the null of no spatial
    structure, nearby sites are no more similar than distant ones, so the statistic is compared to
    its permutation distribution. A LOW observed value means detectable structure.
    """
    rng = np.random.default_rng(seed)
    ok = np.isfinite(v)
    v, x, y = v[ok], x[ok], y[ok]
    if len(v) < 12:
        return dict(p_short_structure=np.nan, nperm=0, npairs_short=0)
    n = len(v)
    i, j = np.triu_indices(n, 1)
    h = np.hypot(x[i] - x[j], y[i] - y[j])
    m = h < short_km
    if m.sum() < 10:
        return dict(p_short_structure=np.nan, nperm=0, npairs_short=int(m.sum()))
    ii, jj = i[m], j[m]

    def stat(vv):
        d = dp.circ_diff_h(vv[ii], vv[jj]) if circular else (vv[ii] - vv[jj])
        return float(0.5 * np.mean(d ** 2))

    obs = stat(v)
    null = np.array([stat(rng.permutation(v)) for _ in range(nperm)])
    p = (1 + np.sum(null <= obs)) / (nperm + 1)
    return dict(p_short_structure=float(p), obs_short_gamma=obs,
                null_median=float(np.median(null)), nperm=nperm, npairs_short=int(m.sum()))


def boot_partition(v, se, x, y, circular=False, nboot=400, seed=0, short_km=25.0):
    """Bootstrap CI on the regional fraction by resampling SITES with replacement."""
    rng = np.random.default_rng(seed)
    ok = np.isfinite(v) & np.isfinite(se)
    v, se, x, y = v[ok], se[ok], x[ok], y[ok]
    n = len(v)
    if n < 12:
        return (np.nan, np.nan, np.nan)
    out = []
    for _ in range(nboot):
        k = rng.integers(0, n, n)
        # jitter duplicated coordinates so resampled pairs are not all at zero separation
        jx = x[k] + rng.normal(0, 0.01, n)
        jy = y[k] + rng.normal(0, 0.01, n)
        r = partition(v[k], se[k], jx, jy, circular, short_km)
        if r and np.isfinite(r['frac_regional']):
            out.append(r['frac_regional'])
    if len(out) < 30:
        return (np.nan, np.nan, np.nan)
    return (float(np.median(out)), float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5)))
