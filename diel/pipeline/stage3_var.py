"""Stage 3: variance partition and empirical variograms on array-level diel metrics.

Partition. For a metric y measured at arrays with known measurement SEs s_i:
  var_obs   = Var(y)                          total between-array spread
  var_meas  = mean(s_i^2)                     measurement error
  var_true  = var_obs - var_meas               real between-array variation
The real part splits by SCALE using the semivariogram. The nugget of a variogram
fitted to PAIRS BEYOND the array radius estimates everything unstructured at
distances the design cannot resolve -- measurement error plus genuine variation
below 25 km. Subtracting the known measurement component leaves LOCAL variation:
  var_local    = max(nugget - var_meas, 0)
  var_regional = max(var_true - var_local, 0)
Reported as fractions of var_obs so the three sum to 1.

Estimator: Cressie-Hawkins robust semivariance is used alongside the classical
Matheron estimator; the classical one is reported and the robust one carried as a
check, because a handful of extreme arrays can dominate a distance band.

Variograms are computed on the CENTRED metric (species mean removed) so all
species share a scale, and peak timing is handled as a CIRCULAR quantity: its
semivariance uses the squared circular difference in hours.
"""
import numpy as np


def circ_diff_h(a, b, period=24.0):
    d = (np.asarray(a) - np.asarray(b)) % period
    return np.where(d > period / 2, d - period, d)


def pair_distances(x_km, y_km):
    x = np.asarray(x_km, float); y = np.asarray(y_km, float)
    dx = x[:, None] - x[None, :]
    dy = y[:, None] - y[None, :]
    return np.sqrt(dx * dx + dy * dy)


def empirical_variogram(x_km, y_km, v, edges, circular=False, min_pairs=30):
    """Classical (Matheron) and robust (Cressie-Hawkins) semivariance per band."""
    D = pair_distances(x_km, y_km)
    v = np.asarray(v, float)
    iu = np.triu_indices(len(v), 1)
    d = D[iu]
    if circular:
        diff = circ_diff_h(v[iu[0]], v[iu[1]])
    else:
        diff = v[iu[0]] - v[iu[1]]
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (d >= lo) & (d < hi)
        n = int(m.sum())
        if n < min_pairs:
            out.append((lo, hi, n, np.nan, np.nan, np.nan))
            continue
        dd = diff[m]
        gamma = 0.5 * float(np.mean(dd ** 2))
        # Cressie-Hawkins robust estimator
        ch = (float(np.mean(np.abs(dd) ** 0.5)) ** 4) / (2 * (0.457 + 0.494 / n + 0.045 / n ** 2))
        out.append((lo, hi, n, float(np.mean(dd[~np.isnan(dd)] * 0 + d[m].mean())), gamma, ch))
    return out


def fit_nugget_sill(bands, max_km=2000.0):
    """Weighted least squares exponential model gamma = c0 + c1 (1 - exp(-h/a))."""
    from scipy.optimize import least_squares
    h = np.array([b[3] for b in bands], float)
    g = np.array([b[4] for b in bands], float)
    n = np.array([b[2] for b in bands], float)
    ok = np.isfinite(h) & np.isfinite(g) & (h <= max_km) & (n > 0)
    h, g, n = h[ok], g[ok], n[ok]
    if len(h) < 4:
        return dict(nugget=np.nan, sill=np.nan, range_km=np.nan, ok=False)
    g0 = float(np.nanmin(g)); s0 = max(float(np.nanmax(g) - g0), 1e-9)
    def resid(p):
        c0, c1, a = p
        pred = c0 + c1 * (1 - np.exp(-h / max(a, 1.0)))
        return np.sqrt(n) * (pred - g)
    try:
        r = least_squares(resid, [g0, s0, 400.0],
                          bounds=([0, 0, 10.0], [np.inf, np.inf, 5000.0]), max_nfev=4000)
        c0, c1, a = r.x
        return dict(nugget=float(c0), sill=float(c0 + c1), range_km=float(a), ok=True)
    except Exception:
        return dict(nugget=np.nan, sill=np.nan, range_km=np.nan, ok=False)


def short_band_gamma(bands, max_h=60.0, min_pairs=40):
    """Semivariance at the shortest RESOLVABLE separation, from data rather than
    from extrapolating a model to h=0.

    For a process with measurement error, unresolved local variation, and regionally
    structured variation with correlation rho(h):
        gamma(h) = var_meas + var_local + var_regional (1 - rho(h))
    At the shortest band, rho ~ 1 whenever the regional range exceeds that distance,
    so gamma(h0) estimates var_meas + var_local. Extrapolating a fitted exponential
    to h=0 instead invents a value below the smallest separation the design contains,
    which is what produced implausibly small nuggets for species whose fitted range
    was itself shorter than the first band.
    """
    cand = [(b[3], b[4], b[2]) for b in bands
            if np.isfinite(b[3]) and np.isfinite(b[4]) and b[3] <= max_h and b[2] >= min_pairs]
    if not cand:
        cand = [(b[3], b[4], b[2]) for b in bands
                if np.isfinite(b[3]) and np.isfinite(b[4]) and b[2] >= min_pairs]
        if not cand:
            return np.nan, np.nan
        cand = cand[:1]
    # inverse-variance (pair-count) weighted mean over the resolvable short bands
    h = np.array([c[0] for c in cand]); g = np.array([c[1] for c in cand]); n = np.array([c[2] for c in cand], float)
    return float(np.sum(g * n) / np.sum(n)), float(np.sum(h * n) / np.sum(n))


def short_band_perm_p(x_km, y_km, v, h_cut, nperm=999, seed=11, circular=False):
    """Is semivariance at h < h_cut lower than under random relabelling of arrays?"""
    v = np.asarray(v, float)
    ok = np.isfinite(v)
    v = v[ok]
    D = pair_distances(np.asarray(x_km)[ok], np.asarray(y_km)[ok])
    iu = np.triu_indices(len(v), 1)
    d = D[iu]
    m = (d > 0) & (d < h_cut)
    if m.sum() < 30:
        return np.nan, int(m.sum())
    def stat(vv):
        dd = circ_diff_h(vv[iu[0]], vv[iu[1]]) if circular else vv[iu[0]] - vv[iu[1]]
        return 0.5 * float(np.mean(dd[m] ** 2))
    obs = stat(v)
    rng = np.random.default_rng(seed)
    null = np.array([stat(rng.permutation(v)) for _ in range(nperm)])
    return float((1 + np.sum(null <= obs)) / (nperm + 1)), int(m.sum())


def partition(v, se, x_km, y_km, edges, circular=False):
    v = np.asarray(v, float); se = np.asarray(se, float)
    good = np.isfinite(v) & np.isfinite(se)
    v, se, x_km, y_km = v[good], se[good], np.asarray(x_km)[good], np.asarray(y_km)[good]
    if circular:
        mu = np.angle(np.mean(np.exp(1j * v / 24 * 2 * np.pi)))
        vc = circ_diff_h(v, mu / (2 * np.pi) * 24)
        var_obs = float(np.mean(vc ** 2))
        vv = vc
    else:
        var_obs = float(np.var(v, ddof=1))
        vv = v - v.mean()
    var_meas = float(np.mean(se ** 2))
    # CIRCULAR metrics must keep circular pair differences: centring alone is not
    # enough, because two peak times either side of the wrap become +11 and -11 and
    # a plain subtraction reports 22 h where the true separation is 2 h.
    bands = empirical_variogram(x_km, y_km, vv, edges, circular=circular)
    fit = fit_nugget_sill(bands)
    g0, h0 = short_band_gamma(bands)
    var_true = max(var_obs - var_meas, 0.0)
    # anchor local+measurement on the observed short-distance semivariance
    anchor = g0 if np.isfinite(g0) else var_meas
    var_local = min(max(anchor - var_meas, 0.0), var_true)
    var_reg = max(var_true - var_local, 0.0)
    tot = var_obs if var_obs > 0 else 1.0
    p_short, npair_short = short_band_perm_p(
        x_km, y_km, vv, max(h0 * 1.6, 60.0) if np.isfinite(h0) else 60.0, circular=circular)
    return dict(n=int(len(v)), var_obs=var_obs, var_meas=var_meas, var_true=var_true,
                var_local=var_local, var_regional=var_reg,
                frac_meas=var_meas / tot, frac_local=var_local / tot,
                frac_regional=var_reg / tot,
                short_gamma=g0, short_h_km=h0, p_short_structure=p_short,
                npairs_short=npair_short,
                nugget_modelfit=fit['nugget'], sill=fit['sill'], range_km=fit['range_km'],
                sd_obs=np.sqrt(var_obs), sd_regional=np.sqrt(var_reg)), bands
