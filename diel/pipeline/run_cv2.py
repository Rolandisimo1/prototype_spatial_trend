"""Both predictor sets on IDENTICAL rows, so a skill difference is attributable to predictor
membership alone. Rows are restricted to sites complete for the UNION of the two sets."""
import os, sys, numpy as np, pandas as pd
sys.path.append(os.getcwd())
import cvlib3 as c3, dielinf as di

site = pd.read_csv('site_analysis_table.csv')
rel = pd.read_csv('rel.csv')
site['pop_1km'] = np.log1p(site['pop_1km'])
relmap = {(r.species, r.metric): bool(r.passes) for r in rel.itertuples()}

SETS = {'reduced7': di.COLS,
        'superseded10': ['pop_1km', 'ag_5km', 'tcc_1km', 'tmax_hottest_month', 'rug_5km',
                         'nlcd_5k_forest', 'elev', 'nlcd_5k_developed', 'nlcd_5k_crop', 'ntl_1km']}
UNION = sorted(set(SETS['reduced7']) | set(SETS['superseded10']))

allp, alls = [], []
for sp in sorted(site.species.unique()):
    g0 = site[site.species == sp]
    for mt in di.MEASURES:
        if not relmap.get((sp, mt), False):
            continue
        circ = mt in di.CIRC
        g = g0.dropna(subset=UNION + [mt, 'x_km', 'y_km'])
        for name, cols in SETS.items():
            P = c3.block_cv(g, mt, cols, circular=circ, seed=3)
            if P is None or not len(P):
                continue
            P['species'] = sp
            P['predictor_set'] = name
            allp.append(P)
            S = c3.skill_from_preds(P, circular=circ, nboot=1000, seed=3)
            S['species'] = sp
            S['metric'] = mt
            S['predictor_set'] = name
            S['n_predictors'] = len(cols)
            alls.append(S)
        print('done', sp, mt, flush=True)

pd.concat(allp, ignore_index=True).to_csv('stage3_cv_predictions_common.csv', index=False)
pd.concat(alls, ignore_index=True).to_csv('stage3_block_cv_summary_common.csv', index=False)
print('CV2 DONE')
