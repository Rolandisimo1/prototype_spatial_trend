"""Block CV for the 5-predictor MAPPING set, on 1 km supports aggregated from the cameras.

Stage 3 cross-validates the 7-predictor inference model, but the surfaces in Stage 4 come from the
5-predictor mapping model. This run asks the out-of-sample question of the model that is actually
mapped."""
import os, sys, numpy as np, pandas as pd
sys.path.append(os.getcwd())
import cvlib3 as c3, dielinf as di

md = pd.read_csv('mapdat_for_cv.csv')
rel = pd.read_csv('rel.csv')
MAPC = ['pop_1km', 'nlcd_1k_crop', 'tcc_1km', 'rug_1km', 't_warmmonth']
md['pop_1km'] = np.log1p(md['pop_1km'])
relmap = {(r.species, r.metric): bool(r.passes) for r in rel.itertuples()}

alls = []
for sp in sorted(md.species.unique()):
    g0 = md[md.species == sp]
    for mt in di.MEASURES:
        if not relmap.get((sp, mt), False):
            continue
        circ = mt in di.CIRC
        g = g0.dropna(subset=MAPC + [mt, 'x_km', 'y_km'])
        P = c3.block_cv(g, mt, MAPC, circular=circ, seed=3)
        if P is None or not len(P):
            continue
        S = c3.skill_from_preds(P, circular=circ, nboot=1000, seed=3)
        S['species'] = sp
        S['metric'] = mt
        S['predictor_set'] = 'mapping5'
        alls.append(S)
        print('done', sp, mt, flush=True)

pd.concat(alls, ignore_index=True).to_csv('stage4_mapping_cv.csv', index=False)
print('CV3 DONE')
