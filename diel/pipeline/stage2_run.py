import numpy as np, pandas as pd, sys
sys.path.insert(0, '.')
import harmfit as HF


def run(AG, min_events=25, verbose=True):
    rows = []
    for (sp, arr), x in AG.groupby(['species', 'final_array'], sort=False):
        x = x.sort_values('bin')
        c = np.zeros(HF.NBIN); e = np.zeros(HF.NBIN)
        b = x.bin.values.astype(int)
        c[b] = x['count'].values; e[b] = x.effort_h.values
        if c.sum() < min_events or (e > 0).sum() < HF.NBIN:
            continue
        r = HF.fit_one(c, np.log(e), min_events=min_events)
        if r is None:
            continue
        raw = HF.raw_metrics(c, e)
        r.update({'raw_' + k: v for k, v in raw.items()})
        r['species'] = sp; r['final_array'] = arr
        r['ndep'] = int(x.ndep.max()); r['effort_h'] = float(e.sum())
        rows.append(r)
        if verbose and len(rows) % 300 == 0:
            print(f"  {len(rows)} fits", flush=True)
    return pd.DataFrame(rows)
