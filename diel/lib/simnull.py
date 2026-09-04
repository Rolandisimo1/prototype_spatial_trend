"""Rule-3 simulated-null test for the rate-like predictors, in the multivariate model.

The problem. A site's detection rate and its predator richness are both built from counts. A site
with few detections yields a noisier fitted curve, and that noise does not average out symmetrically
in bounded or nonlinear measures, so an apparent association between a curve-shape measure and a
count-derived predictor can appear with no behavioural difference at all.

The test. Hold curve SHAPE fixed at the species' national pooled shape. At each site draw counts
from a negative binomial with that fixed shape scaled to the site's OWN total detections and its OWN
per-bin effort, using the site's OWN fitted dispersion. Refit the harmonic curve, recompute the
measure, then refit the SAME 7-predictor weighted model and read off the SAME partial coefficient.
Whatever the null reproduces is arithmetic, not behaviour.

Because the draws depend only on the site's total, effort and dispersion, one set of draws per
species serves every measure and every predictor.

peak_h: each simulated replicate is centred on the circular mean of that replicate's own simulated
peak hours, via circ_mean_h and circ_diff_h. No linear subtraction of clock hours occurs.
"""
import numpy as np
import pandas as pd
import statsmodels.api as sm
import harmfit as hf
import dielpipe as dp
import dielinf as di

BC = hf.BC
ANG = hf.ANG
XD = hf.X


def _gap(c):
    d = np.abs(BC - c)
    return np.minimum(d, 24.0 - d)


NIGHT = BC >= 12.0
CREP = (_gap(0.0) <= 1.0) | (_gap(12.0) <= 1.0)
NOON = _gap(6.0) <= 2.0


def measures_from_beta(beta):
    """The five diel measures, verified to reproduce the delivered curve table exactly.

    peak_h is the mode of the fitted rate curve (the bin centre of its maximum), which is the
    definition used in the delivered curves.
    """
    mu = np.exp(XD @ beta)
    p = mu / mu.sum()
    C = float((p * np.cos(ANG)).sum())
    S = float((p * np.sin(ANG)).sum())
    return (100 * float(p[NIGHT].sum()), 100 * float(p[CREP].sum()), 100 * float(p[NOON].sum()),
            float(np.hypot(C, S)), float(BC[int(np.argmax(mu))]))


MEAS_ORDER = ['pct_noct', 'pct_crep', 'pct_noon', 'conc', 'peak_h']


def simulate_species(sp, sites, bg, harm, national_beta, nsim=400, seed=0):
    """Returns array (nsim, nsites, 5) of measures under the fixed-shape null, and the site order."""
    rng = np.random.default_rng(seed)
    shape = np.exp(XD @ national_beta)
    shape = shape / shape.sum()
    keep, mus, alphas, effs = [], [], [], []
    hrow = harm.set_index('final_array')
    for s in sites:
        g = bg.get((sp, s))
        if g is None:
            continue
        eff = g['effort_h'].values.astype(float)
        tot = float(g['count'].sum())
        if tot < 25 or not np.all(np.isfinite(eff)) or eff.min() <= 0:
            continue
        w = shape * eff
        keep.append(s)
        mus.append(w / w.sum() * tot)
        alphas.append(max(float(hrow.loc[s, 'alpha']), 1e-6))
        effs.append(eff)
    if not keep:
        return np.zeros((0, 0, 5)), []
    out = np.full((nsim, len(keep), 5), np.nan)
    for j, (mu, a, eff) in enumerate(zip(mus, alphas, effs)):
        n_r = np.clip(1.0 / a, 1e-6, 1e8)
        p_r = np.clip(1.0 / (1.0 + a * mu), 1e-9, 1 - 1e-9)
        Y = rng.negative_binomial(n_r, p_r, size=(nsim, len(mu)))
        loge = np.log(eff)
        for i in range(nsim):
            if Y[i].sum() < 25:
                continue
            r = hf.fit_one(Y[i], loge)
            if r is None:
                continue
            out[i, j, :] = measures_from_beta(np.array([r[p] for p in hf.PNAMES]))
    return out, keep


def null_partials(sim, sites, gsp, measure, targets, cols=None):
    """Distribution of the target partial coefficients when curve shape is held fixed."""
    cols = cols or di.COLS
    g = gsp.set_index('final_array').loc[sites]
    ok_cov = g[cols + [measure + '_se']].notna().all(axis=1).values
    Xs = di.design(g.loc[ok_cov], cols)[0].values
    Z = sm.add_constant(Xs)
    w = 1.0 / np.maximum(g.loc[ok_cov, measure + '_se'].values ** 2, 1e-9)
    mi = MEAS_ORDER.index(measure)
    circular = measure in di.CIRC
    ti = [cols.index(t) + 1 for t in targets]
    betas = {t: [] for t in targets}
    for i in range(sim.shape[0]):
        v = sim[i, :, mi][ok_cov]
        if np.isfinite(v).sum() < 0.8 * len(v):
            continue
        m_ok = np.isfinite(v)
        vv = v[m_ok]
        if circular:
            vv = dp.circ_diff_h(vv, dp.circ_mean_h(vv, w[m_ok]))
        if np.std(vv) == 0:
            continue
        try:
            fit = sm.WLS(vv, Z[m_ok], weights=w[m_ok]).fit()
        except Exception:
            continue
        for t, k in zip(targets, ti):
            betas[t].append(float(fit.params[k]))
    return {t: np.array(v) for t, v in betas.items()}
