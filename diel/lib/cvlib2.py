"""Spatial block cross-validation at 400 km against the null of one nationwide curve per species.

Folds are contiguous 400 km spatial blocks, so held-out sites are geographically separated from the
sites used to fit. This is the relevant test for a map: it asks whether the fitted surface predicts
where you have not sampled, not whether it interpolates among neighbours.

Null model: a single value per species (the effort-weighted national mean of the metric), i.e. no
spatial variation at all. Skill is the proportional reduction in held-out squared error relative to
that null,
    skill = 1 - SSE(model) / SSE(null),
so skill <= 0 means the spatial model predicts held-out blocks no better than one national number.
Circular metrics use squared circular differences in hours (rule 7).

Competing predictors:
  position   - smooth function of location (thin-plate-like radial basis on x_km, y_km)
  covariate  - site covariates (human population, agriculture, tree cover, temperature, ruggedness)
Both are fitted with ridge-penalised least squares, with the penalty chosen by an inner
cross-validation on the training blocks only, so no held-out information reaches model selection.
"""
import numpy as np
import pandas as pd
import dielpipe as dp

BLOCK_KM = 400.0


def blocks(x, y, block_km=BLOCK_KM):
    return (np.floor(x / block_km).astype(int).astype(str) + '_' +
            np.floor(y / block_km).astype(int).astype(str))


def _rbf(x, y, cx, cy, scale):
    d2 = (x[:, None] - cx[None, :]) ** 2 + (y[:, None] - cy[None, :]) ** 2
    return np.exp(-d2 / (2 * scale ** 2))


def _ridge_fit(X, v, lam):
    n, p = X.shape
    A = X.T @ X + lam * np.eye(p)
    return np.linalg.solve(A, X.T @ v)


def _circ_resid(pred, obs, circular):
    return dp.circ_diff_h(pred, obs) if circular else (pred - obs)


def _nn_scale(cx, cy, floor=150.0):
    """Median distance from each centre to its nearest other centre.

    This is the radial-basis bandwidth: it sets how far one centre's influence reaches, so it has to
    be a real inter-centre distance. Computing it by sorting x and y independently and taking
    consecutive differences would pair unrelated points and return roughly half the true value,
    giving a basis that wiggles at a finer scale than the data can support.
    """
    if len(cx) < 2:
        return floor
    D = np.hypot(cx[:, None] - cx[None, :], cy[:, None] - cy[None, :])
    np.fill_diagonal(D, np.inf)
    return max(float(np.median(D.min(axis=1))), floor)


def _design_from(kind, fit_df, apply_df, cov_cols, n_centers=40, seed=0):
    """Build a design matrix for apply_df using parameters estimated ONLY from fit_df.

    Everything the design depends on is set by the training rows: which locations serve as radial
    basis centres, the bandwidth, and the covariate centring and scaling. Held-out rows are then
    projected onto that fixed design. This keeps the fold strictly out of sample, including the
    distribution of its locations and covariates and not merely its response values.
    """
    if kind == 'position':
        rng = np.random.default_rng(seed)
        k = min(n_centers, max(6, len(fit_df) // 3))
        idx = rng.choice(len(fit_df), k, replace=False)
        cx, cy = fit_df.x_km.values[idx], fit_df.y_km.values[idx]
        scale = _nn_scale(cx, cy)
        return np.column_stack([np.ones(len(apply_df)), apply_df.x_km / 1000, apply_df.y_km / 1000,
                                _rbf(apply_df.x_km.values, apply_df.y_km.values, cx, cy, scale)])
    mu = fit_df[cov_cols].mean(); sd = fit_df[cov_cols].std().replace(0, 1)
    Z = ((apply_df[cov_cols] - mu) / sd).fillna(0).values
    return np.column_stack([np.ones(len(apply_df)), Z, Z ** 2])


def _design_full(kind, d, cov_cols, n_centers=40, seed=0):
    """Design over all rows, for FINAL fits where every row is training data (no held-out fold)."""
    return _design_from(kind, d, d, cov_cols, n_centers=n_centers, seed=seed)


def block_cv(df, metric, cov_cols, circular=False, block_km=BLOCK_KM, seed=0,
             lams=(0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)):
    """Leave-one-block-out CV. Returns per-site held-out predictions for each competitor."""
    d = df.dropna(subset=[metric, 'x_km', 'y_km']).copy().reset_index(drop=True)
    if len(d) < 30:
        return None
    d['block'] = blocks(d.x_km.values, d.y_km.values, block_km)
    bl = list(pd.unique(d.block))
    if len(bl) < 4:
        return None
    v = d[metric].values.astype(float)
    bidx = {b: np.where(d.block.values == b)[0] for b in bl}
    kinds = ['position']
    if len(cov_cols) and not d[cov_cols].isna().all().any():
        kinds.append('covariate')

    def ridge_pred(kind, tr_i, te_i, lam):
        """Ridge with an UNPENALISED intercept, on a design built from the training rows only.

        The intercept is fitted as the training mean and the penalty applied only to the slopes, so
        heavy shrinkage collapses the model onto the training mean (the null) rather than onto zero.
        Penalising the intercept would make a strongly regularised model predict near-zero percent
        nocturnal, which is far worse than the null and would make every skill look catastrophic.
        """
        tr_df, te_df = d.iloc[tr_i], d.iloc[te_i]
        Atr = _design_from(kind, tr_df, tr_df, cov_cols, seed=seed)[:, 1:]
        Ate = _design_from(kind, tr_df, te_df, cov_cols, seed=seed)[:, 1:]
        y = v[tr_i]
        m = y.mean()
        beta = np.linalg.solve(Atr.T @ Atr + lam * np.eye(Atr.shape[1]), Atr.T @ (y - m))
        return m + Ate @ beta

    rows = []
    for b in bl:
        te_i = bidx[b]
        tr_i = np.concatenate([bidx[o] for o in bl if o != b])
        if len(tr_i) < 20 or len(te_i) < 1:
            continue
        null_pred = dp.circ_mean_h(v[tr_i]) if circular else float(v[tr_i].mean())
        preds = {'null': np.full(len(te_i), null_pred)}
        inner = [o for o in bl if o != b]
        for kind in kinds:
            best, best_err = lams[0], np.inf
            if len(inner) >= 3:
                for lam in lams:
                    errs = []
                    for ib in inner:
                        i_te = bidx[ib]
                        i_tr = np.concatenate([bidx[o] for o in inner if o != ib])
                        if len(i_tr) < 15 or len(i_te) < 1:
                            continue
                        try:
                            p = ridge_pred(kind, i_tr, i_te, lam)
                            errs.append(np.mean(_circ_resid(p, v[i_te], circular) ** 2))
                        except Exception:
                            pass
                    if errs and np.mean(errs) < best_err:
                        best_err, best = float(np.mean(errs)), lam
            try:
                preds[kind] = ridge_pred(kind, tr_i, te_i, best)
            except Exception:
                continue
        sub = d.iloc[te_i]
        for kind, p in preds.items():
            rows.append(pd.DataFrame(dict(
                species=sub.species.values if 'species' in sub else metric,
                metric=metric, model=kind, block=b, final_array=sub.final_array.values,
                obs=v[te_i], pred=np.asarray(p, float), x_km=sub.x_km.values, y_km=sub.y_km.values)))
    return pd.concat(rows, ignore_index=True) if rows else None


def skill_from_preds(P, circular=False, nboot=1000, seed=0):
    """Skill relative to the null, with a block-bootstrap CI (blocks resampled, not sites)."""
    out = []
    nullP = P[P.model == 'null']
    for model, g in P.groupby('model'):
        if model == 'null':
            continue
        m = g.merge(nullP[['final_array', 'pred']].rename(columns={'pred': 'null_pred'}),
                    on='final_array', how='inner')
        e_m = _circ_resid(m.pred.values, m.obs.values, circular) ** 2
        e_n = _circ_resid(m.null_pred.values, m.obs.values, circular) ** 2
        skill = 1 - e_m.sum() / e_n.sum() if e_n.sum() > 0 else np.nan
        rng = np.random.default_rng(seed)
        ub = m.block.unique()
        # precompute per-block error sums; the block bootstrap then only resamples these
        bi = {b: np.where(m.block.values == b)[0] for b in ub}
        sm_ = np.array([e_m[bi[b]].sum() for b in ub])
        sn_ = np.array([e_n[bi[b]].sum() for b in ub])
        pick = rng.integers(0, len(ub), (nboot, len(ub)))
        num = sm_[pick].sum(1); den = sn_[pick].sum(1)
        bs = 1 - num[den > 0] / den[den > 0]
        lo, hi = (np.percentile(bs, [2.5, 97.5]) if len(bs) > 50 else (np.nan, np.nan))
        # median distance from each held-out site to the nearest training site
        dists = []
        for b in ub:
            te = m[m.block == b]; tr = m[m.block != b]
            if not len(tr):
                continue
            for _, r in te.iterrows():
                dists.append(np.min(np.hypot(tr.x_km - r.x_km, tr.y_km - r.y_km)))
        out.append(dict(model=model, n_sites=len(m), nblocks=len(ub), skill=skill, lo=lo, hi=hi,
                        median_heldout_dist_km=float(np.median(dists)) if dists else np.nan,
                        beats_null=bool(np.isfinite(lo) and lo > 0)))
    return pd.DataFrame(out)
