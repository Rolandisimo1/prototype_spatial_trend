"""Block cross-validation with correct circular handling for peak_h.

This reuses cvlib2's design construction unchanged, so the fold-leak fix is preserved: every design
parameter (radial-basis centres, bandwidth, covariate centring and scaling) is estimated from
training rows only and held-out rows are projected onto that fixed design.

What is added. cvlib2 fitted the ridge intercept as the arithmetic mean of the training response,
which is wrong for peak_h: the arithmetic mean of clock hours near midnight lands on the far side of
the circle. Here a circular response is first converted to a signed circular deviation from the
training set's circular mean (circ_mean_h, circ_diff_h), the ridge model is fitted on that
deviation, and predictions are mapped back to clock hours modulo 24. Errors are then scored with
circ_diff_h. The null remains the training circular mean.
"""
import numpy as np
import pandas as pd
import dielpipe as dp
import cvlib2 as c2

BLOCK_KM = 400.0
LAMS = (0.1, 1.0, 10.0, 100.0, 1000.0, 10000.0)


def block_cv(df, metric, cov_cols, circular=False, block_km=BLOCK_KM, seed=0, lams=LAMS,
             min_sites=30):
    d = df.dropna(subset=[metric, 'x_km', 'y_km']).copy().reset_index(drop=True)
    if len(d) < min_sites:
        return None
    d['block'] = c2.blocks(d.x_km.values, d.y_km.values, block_km)
    bl = list(pd.unique(d.block))
    if len(bl) < 4:
        return None
    v = d[metric].values.astype(float)
    bidx = {b: np.where(d.block.values == b)[0] for b in bl}
    kinds = ['position']
    if len(cov_cols) and not d[cov_cols].isna().all().any():
        kinds.append('covariate')

    def ridge_pred(kind, tr_i, te_i, lam):
        tr_df, te_df = d.iloc[tr_i], d.iloc[te_i]
        Atr = c2._design_from(kind, tr_df, tr_df, cov_cols, seed=seed)[:, 1:]
        Ate = c2._design_from(kind, tr_df, te_df, cov_cols, seed=seed)[:, 1:]
        y = v[tr_i]
        if circular:
            centre = dp.circ_mean_h(y)
            yc = dp.circ_diff_h(y, centre)
        else:
            centre = float(y.mean())
            yc = y - centre
        beta = np.linalg.solve(Atr.T @ Atr + lam * np.eye(Atr.shape[1]), Atr.T @ yc)
        pred = Ate @ beta
        return (centre + pred) % 24.0 if circular else (centre + pred)

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
                            errs.append(np.mean(c2._circ_resid(p, v[i_te], circular) ** 2))
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
                obs=v[te_i], pred=np.asarray(p, float), x_km=sub.x_km.values, y_km=sub.y_km.values,
                lam=best)))
    return pd.concat(rows, ignore_index=True) if rows else None


skill_from_preds = c2.skill_from_preds
