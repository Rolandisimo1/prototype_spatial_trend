"""Fit every usable species x site diel curve, before and after the exclusions."""
import numpy as np, pandas as pd, sys, os, warnings
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
warnings.filterwarnings('ignore')
import dielpipe as dp, harmfit as hf

THR = 25


def run(path, tag):
    T = pd.read_parquet(path)
    tot = T.groupby(['species', 'final_array'])['count'].sum()
    use = tot[tot >= THR].index
    T = T.set_index(['species', 'final_array']).sort_index()
    rows = []
    for i, (sp, ar) in enumerate(use):
        g = T.loc[(sp, ar)].sort_values('bin')
        r = hf.fit_curve(g['count'].values, g.effort_h.values)
        if r is None:
            continue
        raw = dp.metrics_anchored(g['count'].values / g.effort_h.values)
        r.update(species=sp, final_array=ar, ndep=float(g.ndep.iloc[0]),
                 effort_h=float(g.effort_h.sum()))
        for k, v in raw.items():
            r['raw_' + k] = v
        rows.append(r)
        if (i + 1) % 400 == 0:
            print('  %s %d/%d' % (tag, i + 1, len(use)), flush=True)
    H = pd.DataFrame(rows)
    H.to_csv(f'harmonics_{tag}.csv', index=False)
    print('%s: %d curves, %d converged, %d poisson fallback' %
          (tag, len(H), int(H.converged.sum()), int((H.family == 'poisson_fallback').sum())), flush=True)
    return H


Hc = run('work/binned_clean.parquet', 'clean')
Hp = run('work/binned_preexcl.parquet', 'preexcl')
print('FITS DONE', flush=True)
