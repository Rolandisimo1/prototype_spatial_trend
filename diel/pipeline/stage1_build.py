"""Stage 1: per-deployment sun-time effort (48 bins) over the Aug-Oct active window,
plus clock-time control effort. Emits work/effort_sun.parquet and work/effort_clock.parquet.
"""
import numpy as np, pandas as pd, sys
sys.path.insert(0, '.')
import solarclock as SC
import suneffort as SE

NB = SE.NBIN


def build(dm, out_sun, out_clock):
    """dm: dep_key, sd, ed, latitude, longitude. Aug 1 - Oct 31 clipped."""
    rows_s, rows_c = [], []
    n = len(dm)
    for i, r in enumerate(dm.itertuples(index=False)):
        if i % 2000 == 0:
            print(f"  {i}/{n}", flush=True)
        # clip to Aug-Oct within each calendar year the deployment spans
        segs = []
        for yr in range(r.sd.year, r.ed.year + 1):
            a = max(r.sd, pd.Timestamp(yr, 8, 1))
            b = min(r.ed, pd.Timestamp(yr, 11, 1))
            if b > a:
                segs.append((a, b))
        if not segs:
            continue
        es = np.zeros(NB); ec = np.zeros(NB)
        for a, b in segs:
            days = pd.date_range(a.floor('D'), b.floor('D'), freq='D')
            sr, ss = SC.solar_window(days, np.full(len(days), r.latitude),
                                     np.full(len(days), r.longitude))
            for j, d in enumerate(days):
                lo = max(a, d); hi = min(b, d + pd.Timedelta('1D'))
                dur = (hi - lo).total_seconds() / 3600.0
                if dur <= 0:
                    continue
                if dur >= 24.0 - 1e-9:
                    es += SE.full_day_effort(ss[j] - sr[j])
                else:
                    h0 = (lo - d).total_seconds() / 3600.0
                    es += SE.partial_day_effort(sr[j], ss[j], h0, h0 + dur)
                # clock-time control: uniform per half-hour clock bin
                ck = np.zeros(NB)
                b0 = np.arange(NB) * 0.5
                h0 = (lo - d).total_seconds() / 3600.0
                h1 = h0 + dur
                ck = np.clip(np.minimum(b0 + 0.5, h1) - np.maximum(b0, h0), 0, None)
                ec += ck
        rows_s.append((r.dep_key, es))
        rows_c.append((r.dep_key, ec))

    def flat(rows, col):
        keys = np.repeat([k for k, _ in rows], NB)
        bins = np.tile(np.arange(NB), len(rows))
        vals = np.concatenate([v for _, v in rows])
        return pd.DataFrame({'dep_key': keys, 'bin': bins.astype(np.int16), col: vals})

    flat(rows_s, 'effort_h').to_parquet(out_sun, index=False)
    flat(rows_c, 'effort_h').to_parquet(out_clock, index=False)
    print("done", len(rows_s))
