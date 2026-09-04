"""Hypothesis tests on site-level diel metrics, with the simulated-null test required by rule 8.

Why the simulated null exists. Regressing a curve-shape metric on a RATE-LIKE predictor (a species'
own detection rate, or predator richness built from detection rates) is circular in a specific way:
both sides depend on how many detections a site produced. A site with few detections has a noisier
fitted curve, and that noise does not average out symmetrically in bounded or nonlinear metrics, so
an apparent effect can appear with no real change in curve shape at all. The test: hold curve SHAPE
FIXED at the national mean, simulate counts at each site's OWN observed rate and effort, refit, and
measure the effect again. Whatever the null reproduces is not evidence about behaviour.

Spatial control: the same regression with a smooth function of location partialled out of both
sides, so an effect that is really a north-south or east-west gradient does not read as a covariate
effect.
"""
import numpy as np
import pandas as pd
import dielpipe as dp
import harmfit as hf

CIRC = {'peak_h', 'mean_h'}


def zs(v):
    v = np.asarray(v, float)
    s = np.nanstd(v)
    return (v - np.nanmean(v)) / (s if s > 0 else 1.0)


def std_beta(y, x, circular=False, ref=None):
    """Standardised slope of y on x. Circular y is centred on its circular mean first."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 12:
        return np.nan
    xx, yy = x[ok], y[ok]
    if circular:
        c = dp.circ_mean_h(yy) if ref is None else ref
        yy = dp.circ_diff_h(yy, c)
    sy = np.std(yy)
    if sy == 0:
        return np.nan
    return float(np.polyfit(zs(xx), yy / sy, 1)[0])


def boot_beta(y, x, circular=False, nboot=1000, seed=0):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    n = len(x)
    if n < 12:
        return (np.nan, np.nan, np.nan)
    b = std_beta(y, x, circular)
    out = []
    for _ in range(nboot):
        k = rng.integers(0, n, n)
        v = std_beta(y[k], x[k], circular)
        if np.isfinite(v):
            out.append(v)
    if len(out) < 50:
        return (b, np.nan, np.nan)
    return (b, float(np.percentile(out, 2.5)), float(np.percentile(out, 97.5)))


def spatial_control_beta(y, x, xk, yk, circular=False, n_centers=25, lam=10.0, seed=0):
    """Slope of y on x after removing a smooth spatial trend from BOTH sides."""
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(xk) & np.isfinite(yk)
    if ok.sum() < 25:
        return np.nan
    x, y, xk, yk = x[ok], y[ok], xk[ok], yk[ok]
    if circular:
        y = dp.circ_diff_h(y, dp.circ_mean_h(y))
    rng = np.random.default_rng(seed)
    k = min(n_centers, max(5, len(x) // 4))
    idx = rng.choice(len(x), k, replace=False)
    cx, cy = xk[idx], yk[idx]
    # bandwidth = median distance from each centre to its nearest other centre
    Dc = np.hypot(cx[:, None] - cx[None, :], cy[:, None] - cy[None, :])
    np.fill_diagonal(Dc, np.inf)
    scale = max(float(np.median(Dc.min(axis=1))), 150.0)
    d2 = (xk[:, None] - cx[None, :]) ** 2 + (yk[:, None] - cy[None, :]) ** 2
    S = np.column_stack([np.ones(len(x)), xk / 1000, yk / 1000, np.exp(-d2 / (2 * scale ** 2))])

    def resid(v):
        m = v.mean()
        A = S[:, 1:]
        beta = np.linalg.solve(A.T @ A + lam * np.eye(A.shape[1]), A.T @ (v - m))
        return v - (m + A @ beta)

    rx, ry = resid(x), resid(y)
    sy = np.std(ry)
    if sy == 0 or np.std(rx) == 0:
        return np.nan
    return float(np.polyfit(zs(rx), ry / sy, 1)[0])


def simulate_fixed_shape(bysite, national_coef, metrics, nsim=150, seed=0):
    """Simulate counts from a FIXED national curve shape at each site's own rate and effort.

    The simulated curves depend only on the site's total count, its per-bin effort and its
    dispersion, NOT on any predictor, so one set of draws serves every predictor and every metric.
    Returns {metric: array (nsim, nsites)} plus the site order, so each rule-8 null test is then
    just a regression of already-simulated values on the predictor.
    """
    rng = np.random.default_rng(seed)
    shape = hf.curve_from_coef(national_coef)
    shape = shape / shape.sum()
    sites, mus, alphas, effs = [], [], [], []
    for s, g in bysite.items():
        eff = g['effort_h']; tot = g['count'].sum()
        if tot < 25 or not np.all(np.isfinite(eff)) or np.min(eff) <= 0:
            continue
        w = shape * eff
        sites.append(s); mus.append(w / w.sum() * tot)
        alphas.append(max(g['alpha'], 1e-6)); effs.append(eff)
    if not sites:
        return {}, []
    out = {m: np.full((nsim, len(sites)), np.nan) for m in metrics}
    for j, (mu, a, eff) in enumerate(zip(mus, alphas, effs)):
        n_r = np.clip(1 / a, 1e-6, 1e8)
        p_r = np.clip(1 / (1 + a * mu), 1e-9, 1 - 1e-9)
        Y = rng.negative_binomial(n_r, p_r, size=(nsim, len(mu)))
        for i in range(nsim):
            if Y[i].sum() < 5:
                continue
            r = hf.fit_curve(Y[i], eff)
            if r is None:
                continue
            for m in metrics:
                out[m][i, j] = r[m]
    return out, sites


def null_beta_dist(sim_vals, sites, predictor, metric):
    """Distribution of the standardised slope when curve shape is held fixed."""
    x = np.array([predictor.get(s, np.nan) for s in sites], float)
    circular = metric in CIRC
    betas = []
    for i in range(sim_vals.shape[0]):
        v = sim_vals[i]
        ok = np.isfinite(v) & np.isfinite(x)
        if ok.sum() < 12:
            continue
        b = std_beta(v[ok], x[ok], circular)
        if np.isfinite(b):
            betas.append(b)
    if not betas:
        return dict(null_median=np.nan, null_lo=np.nan, null_hi=np.nan,
                    null_abs_median=np.nan, nsim_ok=0)
    betas = np.array(betas)
    return dict(null_median=float(np.median(betas)), null_lo=float(np.percentile(betas, 2.5)),
                null_hi=float(np.percentile(betas, 97.5)),
                null_abs_median=float(np.median(np.abs(betas))), nsim_ok=len(betas))


def simulated_null(binned_sp, sites, predictor, metric, national_coef, nsim=200, seed=0):
    """Rule-8 null: FIXED national curve shape, each site's own rate and effort, refit, re-measure.

    Returns the distribution of the standardised slope obtained when curve shape does not vary at
    all. If the observed slope sits inside this distribution, the apparent effect is a counting
    artefact rather than evidence of a behavioural difference.
    """
    rng = np.random.default_rng(seed)
    circular = metric in CIRC
    shape = hf.curve_from_coef(national_coef)
    shape = shape / shape.sum()
    recs = []
    for s in sites:
        g = binned_sp.get(s)
        if g is None:
            continue
        eff = g['effort_h']; tot = g['count'].sum()
        if tot < 25 or not np.all(np.isfinite(eff)) or eff.min() <= 0:
            continue
        # expected counts: national SHAPE times this site's own total, weighted by its own effort
        w = shape * eff
        mu = w / w.sum() * tot
        recs.append((s, mu, g['alpha']))
    betas = []
    for _ in range(nsim):
        vals, xs = [], []
        for s, mu, alpha in recs:
            a = max(alpha, 1e-6)
            y = rng.negative_binomial(np.clip(1 / a, 1e-6, 1e8),
                                      np.clip(1 / (1 + a * mu), 1e-9, 1 - 1e-9))
            if y.sum() < 5:
                continue
            r = hf.fit_curve(y, binned_sp[s]['effort_h'])
            if r is None or not np.isfinite(r[metric]):
                continue
            vals.append(r[metric]); xs.append(predictor[s])
        if len(vals) >= 12:
            b = std_beta(np.array(vals), np.array(xs), circular)
            if np.isfinite(b):
                betas.append(b)
    if not betas:
        return dict(null_median=np.nan, null_lo=np.nan, null_hi=np.nan, nsim_ok=0)
    return dict(null_median=float(np.median(betas)), null_lo=float(np.percentile(betas, 2.5)),
                null_hi=float(np.percentile(betas, 97.5)), nsim_ok=len(betas),
                null_abs_median=float(np.median(np.abs(betas))))
