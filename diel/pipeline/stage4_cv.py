"""Stage 4: spatial block cross-validation of array-level diel metrics.

Competitors, all predicting a HELD-OUT array's metric value:
  (a) nationwide-curve null  -- one number per species: the training-fold mean of the
      metric. This is the honest competitor named in the brief: a single continental
      diel curve carries no spatial information, so its prediction is a constant.
  (b) position-only smooth   -- thin-plate-style smooth of the metric on (x_km, y_km),
      here a ridge-penalised low-order polynomial + radial basis on training arrays.
  (c) covariate model        -- ridge regression on population density (urbanisation),
      agriculture (kept SEPARATE from urbanisation, never composited), tree canopy,
      temperature and terrain, plus effort-weighted day-of-year.

Skill = 1 - SSE_model / SSE_null, computed POOLED over held-out arrays. Uncertainty
is bootstrapped over FOLDS (the brief's rule -- arrays inside a block are not
independent), so the CI reflects how much the answer depends on which blocks you drew.

Errors are inverse-variance weighted by each array's measurement SE where a metric
has one: an array whose metric is known to +/- 15 pp should not dominate the score.
Circular metrics (peak timing) use squared circular differences throughout.
"""
import numpy as np
import pandas as pd


def circ_diff_h(a, b, period=24.0):
    d = (np.asarray(a) - np.asarray(b)) % period
    return np.where(d > period / 2, d - period, d)


def assign_blocks(x_km, y_km, block_km=400.0, seed=0):
    bx = np.floor(np.asarray(x_km) / block_km).astype(int)
    by = np.floor(np.asarray(y_km) / block_km).astype(int)
    keys = pd.Series([f"{a}_{b}" for a, b in zip(bx, by)])
    uniq = keys.unique()
    rng = np.random.default_rng(seed)
    rng.shuffle(uniq)
    return keys.values, {k: i for i, k in enumerate(uniq)}


def rbf_design(xy, centres, scale):
    d2 = ((xy[:, None, :] - centres[None, :, :]) ** 2).sum(-1)
    return np.exp(-d2 / (2 * scale ** 2))


def ridge_fit(X, y, w, lam):
    Xw = X * np.sqrt(w)[:, None]
    yw = y * np.sqrt(w)
    A = Xw.T @ Xw + lam * np.eye(X.shape[1])
    A[0, 0] -= lam                      # do not penalise the intercept
    return np.linalg.solve(A, Xw.T @ yw)


def _standardise(tr, te):
    mu = tr.mean(0); sd = tr.std(0); sd[sd == 0] = 1.0
    return (tr - mu) / sd, (te - mu) / sd


def cv_species(x, covcols, block_km=400.0, seed=0, lam_grid=(1., 10., 100., 1000.),
               n_centres=40, circular=False, metric='pct_noct'):
    """Leave-one-block-out CV for one species x one metric. Returns per-array rows."""
    x = x.reset_index(drop=True)
    y = x[metric].values.astype(float)
    se = x[metric + '_se'].values.astype(float)
    ok = np.isfinite(y) & np.isfinite(se) & np.isfinite(x[covcols].values).all(1)
    x, y, se = x[ok].reset_index(drop=True), y[ok], se[ok]
    if len(x) < 30:
        return None
    w = 1.0 / np.clip(se ** 2, np.percentile(se, 10) ** 2, None)
    w = w / w.mean()
    blk, _ = assign_blocks(x.x_km, x.y_km, block_km, seed)
    xy = x[['x_km', 'y_km']].values
    C = x[covcols].values.astype(float)
    rows = []
    for b in pd.unique(blk):
        te = blk == b
        tr = ~te
        if tr.sum() < 25 or te.sum() < 1:
            continue
        # (a) nationwide-curve null: training mean (circular mean for peak timing)
        if circular:
            ang = y[tr] / 24 * 2 * np.pi
            mu0 = (np.angle(np.average(np.exp(1j * ang), weights=w[tr])) % (2 * np.pi)) / (2 * np.pi) * 24
            e_null = circ_diff_h(y[te], mu0)
        else:
            mu0 = float(np.average(y[tr], weights=w[tr]))
            e_null = y[te] - mu0
        # working response for the regressions: centred (circularly if needed)
        yy = circ_diff_h(y, mu0) if circular else y - mu0
        # (b) position-only smooth
        rng = np.random.default_rng(seed + 1)
        idx = rng.choice(np.where(tr)[0], min(n_centres, tr.sum()), replace=False)
        centres = xy[idx]
        scale = np.median(np.sqrt(((centres[:, None, :] - centres[None, :, :]) ** 2).sum(-1))) / 2
        scale = max(scale, 50.0)
        Ptr = np.column_stack([np.ones(tr.sum()), xy[tr] / 1000.0, rbf_design(xy[tr], centres, scale)])
        Pte = np.column_stack([np.ones(te.sum()), xy[te] / 1000.0, rbf_design(xy[te], centres, scale)])
        # (c) covariate model
        Ctr_s, Cte_s = _standardise(C[tr], C[te])
        Xtr = np.column_stack([np.ones(tr.sum()), Ctr_s])
        Xte = np.column_stack([np.ones(te.sum()), Cte_s])
        preds = {}
        for name, (A, B) in dict(position=(Ptr, Pte), covariate=(Xtr, Xte)).items():
            # inner CV over lam on the TRAINING blocks only
            best, bl = np.inf, lam_grid[0]
            itr_blk = blk[tr]
            for lam in lam_grid:
                sse = 0.0
                for ib in pd.unique(itr_blk):
                    m2 = itr_blk == ib
                    if (~m2).sum() < 15 or m2.sum() < 1:
                        continue
                    beta = ridge_fit(A[~m2], yy[tr][~m2], w[tr][~m2], lam)
                    r = yy[tr][m2] - A[m2] @ beta
                    sse += float(np.sum(w[tr][m2] * r ** 2))
                if sse < best:
                    best, bl = sse, lam
            beta = ridge_fit(A, yy[tr], w[tr], bl)
            preds[name] = B @ beta
            preds[name + '_lam'] = bl
        # held-out distance to nearest TRAINING array
        d_near = np.sqrt(((xy[te][:, None, :] - xy[tr][None, :, :]) ** 2).sum(-1)).min(1)
        for j, i in enumerate(np.where(te)[0]):
            rows.append(dict(final_array=x.final_array.values[i], block=b, metric=metric,
                             y=y[i], se=se[i], w=w[i],
                             e_null=e_null[j],
                             e_position=(yy[i] - preds['position'][j]),
                             e_covariate=(yy[i] - preds['covariate'][j]),
                             lam_position=preds['position_lam'], lam_covariate=preds['covariate_lam'],
                             dist_train_km=d_near[j], n_train=int(tr.sum()), n_test=int(te.sum())))
    return pd.DataFrame(rows)


def skill_with_ci(cv, nboot=2000, seed=3):
    """Pooled skill vs the null, with a bootstrap over FOLDS (blocks)."""
    out = {}
    blocks = cv.block.unique()
    rng = np.random.default_rng(seed)
    def sse(d, col):
        return float(np.sum(d.w * d[col] ** 2))
    for model in ['position', 'covariate']:
        s0, s1 = sse(cv, 'e_null'), sse(cv, 'e_' + model)
        pt = 1 - s1 / s0 if s0 > 0 else np.nan
        bs = []
        by = {b: g for b, g in cv.groupby('block')}
        for _ in range(nboot):
            pick = rng.choice(blocks, len(blocks), replace=True)
            d = pd.concat([by[b] for b in pick], ignore_index=True)
            a, c = sse(d, 'e_null'), sse(d, 'e_' + model)
            if a > 0:
                bs.append(1 - c / a)
        bs = np.array(bs)
        out[model] = dict(skill=pt, lo=float(np.percentile(bs, 2.5)),
                          hi=float(np.percentile(bs, 97.5)), nblocks=len(blocks))
    return out
