"""Stage 2: per (species x array) 5-parameter harmonic negative-binomial fit with
offset(log_effort) on the 48 sun-time bins, plus delta-method SEs on every derived
diel metric.

Model:  log E[count_b] = log(effort_b) + b0 + s1 sin(t_b) + c1 cos(t_b)
                                             + s2 sin(2 t_b) + c2 cos(2 t_b)
with t_b = 2*pi*bin_centre/24 in DOUBLE-ANCHORED sun time.

Metrics are functionals of the fitted RATE curve mu_b = exp(eta_b) (effort divided
out, i.e. the offset is set to zero so the curve is a rate per hour):
  pct_noct  100 * sum_night p_b            p = mu / sum(mu)
  pct_crep  100 * sum_crep  p_b            crep = within 1.5 sun-h of an anchor
  conc      resultant length R of p on the 24-h circle (1 - circular variance)
  act       Rowcliffe activity level, mean(mu)/max(mu)
  peak_h    circular mean direction, in sun hours -- CYCLIC, never averaged linearly
SEs come from the delta method on the 5x5 coefficient covariance: for a scalar
functional g(beta), se = sqrt(J Sigma J') with J evaluated by central differences.
Peak timing gets a circular SE from the same Jacobian on the (C,S) resultant.
"""
import numpy as np
import statsmodels.api as sm

NBIN = 48
BC = (np.arange(NBIN) + 0.5) * 0.5
REF_SR, REF_SS = 0.0, 12.0          # double-anchored: sunrise at 0, sunset at 12
IS_NIGHT = BC >= 12.0
IS_CREP = (np.minimum(np.abs(BC - 0.0), np.abs(BC - 24.0)) <= 1.5) | (np.abs(BC - 12.0) <= 1.5)
ANG = BC / 24.0 * 2 * np.pi
TR = ANG.copy()
X = np.column_stack([np.ones(NBIN), np.sin(TR), np.cos(TR), np.sin(2 * TR), np.cos(2 * TR)])
PNAMES = ["b0", "s1", "c1", "s2", "c2"]
METRICS = ["pct_noct", "pct_crep", "conc", "act", "peak_h"]


def curve(beta):
    return np.exp(X @ beta)


def metrics_from_beta(beta):
    mu = curve(beta)
    s = mu.sum()
    p = mu / s
    C = float((p * np.cos(ANG)).sum())
    S = float((p * np.sin(ANG)).sum())
    R = float(np.hypot(C, S))
    peak = float((np.arctan2(S, C) % (2 * np.pi)) / (2 * np.pi) * 24.0)
    return dict(pct_noct=100.0 * float(p[IS_NIGHT].sum()),
                pct_crep=100.0 * float(p[IS_CREP].sum()),
                conc=R,
                act=float(mu.mean() / mu.max()),
                peak_h=peak, _C=C, _S=S)


def _jac(beta, keys, h=1e-5):
    J = np.zeros((len(keys), 5))
    for k in range(5):
        bp = beta.copy(); bp[k] += h
        bm = beta.copy(); bm[k] -= h
        mp, mm = metrics_from_beta(bp), metrics_from_beta(bm)
        for i, key in enumerate(keys):
            J[i, k] = (mp[key] - mm[key]) / (2 * h)
    return J


def fit_one(count, log_effort, min_events=25):
    """Returns dict of coefficients, metrics, SEs, dispersion, convergence flag."""
    n = int(np.asarray(count).sum())
    if n < min_events:
        return None
    y = np.asarray(count, float)
    off = np.asarray(log_effort, float)
    # Poisson start, then NB2 with alpha estimated from the Poisson Pearson dispersion.
    try:
        pm = sm.GLM(y, X, family=sm.families.Poisson(), offset=off).fit()
    except Exception:
        return None
    mu0 = np.clip(pm.mu, 1e-9, None)
    # method-of-moments alpha for NB2: Var = mu + alpha mu^2
    num = ((y - mu0) ** 2 - mu0).sum()
    den = (mu0 ** 2).sum()
    alpha = max(float(num / den), 1e-6) if den > 0 else 1e-6
    res, fam = None, None
    try:
        m = sm.GLM(y, X, family=sm.families.NegativeBinomial(alpha=alpha), offset=off)
        res = m.fit(start_params=pm.params, maxiter=200)
        fam = "nb2"
        if not np.all(np.isfinite(res.bse)) or not res.converged:
            raise RuntimeError("nb did not converge")
    except Exception:
        res, fam, alpha = pm, "poisson_fallback", 0.0
    beta = np.asarray(res.params, float)
    Sig = np.asarray(res.cov_params(), float)
    met = metrics_from_beta(beta)
    keys = ["pct_noct", "pct_crep", "conc", "act"]
    J = _jac(beta, keys)
    var = np.einsum('ik,kl,il->i', J, Sig, J)
    out = dict(n_events=n, family=fam, alpha=alpha,
               dispersion=float(((y - np.clip(res.mu, 1e-9, None)) ** 2 /
                                 np.clip(res.mu, 1e-9, None)).sum() / max(NBIN - 5, 1)))
    for i, pn in enumerate(PNAMES):
        out[pn] = float(beta[i])
        out[pn + "_se"] = float(np.sqrt(max(Sig[i, i], 0)))
    for i, k in enumerate(keys):
        out[k] = float(met[k])
        out[k + "_se"] = float(np.sqrt(max(var[i], 0)))
    # circular SE on peak: delta method on (C,S) -> angle, se_angle = |d theta| in hours
    Jcs = _jac(beta, ["_C", "_S"])
    Vcs = Jcs @ Sig @ Jcs.T
    C, S = met["_C"], met["_S"]
    r2 = C * C + S * S
    if r2 > 1e-12:
        gt = np.array([-S / r2, C / r2])          # d theta / d(C,S)
        se_rad = float(np.sqrt(max(gt @ Vcs @ gt, 0)))
    else:
        se_rad = np.nan
    out["peak_h"] = float(met["peak_h"])
    out["peak_h_se"] = se_rad / (2 * np.pi) * 24.0
    out["conc_R"] = float(np.hypot(C, S))
    return out


def raw_metrics(count, effort_h):
    """Metrics read straight off the effort-corrected raw histogram (no smoothing)."""
    c = np.asarray(count, float); e = np.asarray(effort_h, float)
    ok = e > 0
    r = np.zeros(NBIN); r[ok] = c[ok] / e[ok]
    s = r.sum()
    if s <= 0:
        return dict(pct_noct=np.nan, pct_crep=np.nan, conc=np.nan, act=np.nan, peak_h=np.nan)
    p = r / s
    C = float((p * np.cos(ANG)).sum()); S = float((p * np.sin(ANG)).sum())
    return dict(pct_noct=100.0 * float(p[IS_NIGHT].sum()),
                pct_crep=100.0 * float(p[IS_CREP].sum()),
                conc=float(np.hypot(C, S)),
                act=float(r.mean() / r.max()),
                peak_h=float((np.arctan2(S, C) % (2 * np.pi)) / (2 * np.pi) * 24.0))
