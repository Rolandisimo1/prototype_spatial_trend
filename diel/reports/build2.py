"""Rebuild the sun-anchored binned count/effort table from source.

Anchoring: sunrise -> sun time 0, sunset -> sun time 12, next sunrise -> 24. Each solar day is
split into 48 half-hour sun-time bins; the day arc's 24 bins each carry day_length/24 CLOCK hours
and the night arc's 24 bins each carry night_length/24 clock hours, so per-bin effort is NOT
uniform in sun-anchored time. Effort is the exact clock measure of the intersection of the
deployment's active interval with each bin's clock-time preimage, so partial first/last days and
the Aug 1 - Oct 31 window truncation are handled exactly.
"""
import numpy as np, pandas as pd, sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import dielpipe as dp
import solarclock as sc

SRC = '/Users/rwkays/claude_code/activity_patterns/activity_data/trend_data_arielle_used/'
OUT = 'work'
os.makedirs(OUT, exist_ok=True)
NB = 48

sp13 = pd.read_csv("/Users/rwkays/.claude-science/orgs/4ac09562-3ef9-4985-8243-8cd6d9ebdcd2/artifacts/proj_47a049eaa5c6/8d3339e0-5f93-42d6-a38b-ddb752bc908e/vb92ad9c8_species_final13.csv")
SPECIES13 = list(sp13.species)

# ---------------------------------------------------------------- deployments
dep = pd.read_csv(SRC + 'combined_deployments_all.csv', dtype={'project_id': str, 'deployment_id': str})
dep['dep_key'] = dep.project_id + '|' + dep.deployment_id
dep = dep.drop_duplicates('dep_key', keep='first').reset_index(drop=True)
dep['t0'] = pd.to_datetime(dep.start_date, format='ISO8601', utc=True).dt.tz_localize(None)
dep['t1'] = pd.to_datetime(dep.end_date, format='ISO8601', utc=True).dt.tz_localize(None)
dep = dep.dropna(subset=['t0', 't1', 'latitude', 'longitude'])
dep = dep[dep.t1 > dep.t0].reset_index(drop=True)

# clip each deployment to the Aug 1 - Oct 31 window of every year it touches
pieces = []
for yr in range(int(dep.t0.dt.year.min()), int(dep.t1.dt.year.max()) + 1):
    a = pd.Timestamp(f'{yr}-08-01'); b = pd.Timestamp(f'{yr}-11-01')
    m = (dep.t0 < b) & (dep.t1 > a)
    if not m.any():
        continue
    s = dep.loc[m, ['dep_key', 't0', 't1', 'latitude', 'longitude']].copy()
    s['a'] = s.t0.clip(lower=a); s['b'] = s.t1.clip(upper=b)
    pieces.append(s[s.b > s.a])
win = pd.concat(pieces, ignore_index=True)
win['nday'] = (win.b.dt.normalize() - win.a.dt.normalize()).dt.days + 1
print('deployment-year ASO windows: %d over %d deployments' % (len(win), win.dep_key.nunique()), flush=True)

deps = pd.Index(win.dep_key.unique())
win['di'] = deps.get_indexer(win.dep_key)
EFF = np.zeros((len(deps), NB))

# ---- expand to solar days and accumulate exact per-bin clock overlap
idx = np.repeat(np.arange(len(win)), win.nday.to_numpy())
offs = np.concatenate([np.arange(n) for n in win.nday.to_numpy()])
day0 = win.a.dt.normalize().to_numpy()[idx] + offs.astype('timedelta64[D]')
lat = win.latitude.to_numpy()[idx]; lon = win.longitude.to_numpy()[idx]
A = win.a.to_numpy()[idx]; Bv = win.b.to_numpy()[idx]; di = win.di.to_numpy()[idx]
print('camera-days to bin: %d' % len(day0), flush=True)

HR = np.timedelta64(3600, 's')
CH = 400_000
for lo in range(0, len(day0), CH):
    sl = slice(lo, min(lo + CH, len(day0)))
    d0 = day0[sl]
    sr, ss = sc.solar_window(pd.DatetimeIndex(d0) + pd.Timedelta(hours=12), lat[sl], lon[sl])
    sr2, _ = sc.solar_window(pd.DatetimeIndex(d0) + pd.Timedelta(hours=36), lat[sl], lon[sl])
    t_sr = d0 + (sr * 3600).astype('timedelta64[s]')
    t_ss = d0 + (ss * 3600).astype('timedelta64[s]')
    t_sr2 = d0 + np.timedelta64(1, 'D') + (sr2 * 3600).astype('timedelta64[s]')
    dl = (t_ss - t_sr) / HR                      # day length, clock hours
    nl = (t_sr2 - t_ss) / HR                     # night length, clock hours
    # 49 bin edges in absolute time, monotone from sunrise to next sunrise
    j = np.arange(25)
    edges = np.empty((len(d0), NB + 1), 'datetime64[s]')
    edges[:, :25] = t_sr[:, None] + ((dl[:, None] / 24 * j[None, :]) * 3600).astype('timedelta64[s]')
    edges[:, 25:] = t_ss[:, None] + ((nl[:, None] / 24 * j[None, 1:]) * 3600).astype('timedelta64[s]')
    a = np.maximum(edges[:, :-1], A[sl][:, None])
    b = np.minimum(edges[:, 1:], Bv[sl][:, None])
    ov = np.clip((b - a) / HR, 0, None)
    np.add.at(EFF, di[sl], ov)
    print('  effort chunk %d/%d' % (sl.stop, len(day0)), flush=True)

tot = EFF.sum(1)
print('effort built. total camera-hours %.1f' % tot.sum(), flush=True)
r = EFF.max(1) / np.where(EFF.min(1) > 0, EFF.min(1), np.nan)
print('within-camera per-bin effort max/min: median %.4f  p99 %.4f  max %.4f'
      % (np.nanmedian(r), np.nanpercentile(r, 99), np.nanmax(r)), flush=True)
np.save(f'{OUT}/EFF.npy', EFF)
pd.DataFrame({'dep_key': deps, 'eff_tot_h': tot}).to_parquet(f'{OUT}/eff_dep.parquet')

# ---------------------------------------------------------------- detections
seq = pd.read_csv(SRC + 'combined_sequences_all.csv',
                  usecols=['project_id', 'deployment_id', 'common_name', 'start_time'],
                  dtype={'project_id': str, 'deployment_id': str, 'common_name': str, 'start_time': str})
seq = seq[seq.common_name.isin(SPECIES13)].copy()
seq['dep_key'] = seq.project_id + '|' + seq.deployment_id
seq['ts'] = pd.to_datetime(seq.start_time, format='ISO8601', utc=True).dt.tz_localize(None)
seq = seq.dropna(subset=['ts'])

# Restrict to the Aug 1 - Oct 31 analysis window and to deployments that carry ASO effort.
# The recorded end_date bounds EFFORT, not detection eligibility: 0.98% of ASO records
# postdate their deployment's recorded end date, and the superseded run retained them.
# Dropping them would remove real detections, so they are kept and counted in their bin.
seq = seq[seq.ts.dt.month.isin(dp.ASO_MONTHS)]
seq = seq[seq.dep_key.isin(set(win.dep_key))]
print('ASO 13-sp sequences on deployments with ASO effort: %d' % len(seq), flush=True)

det = dp.independent_detections(seq[['dep_key', 'common_name', 'ts']])
print('independent detections (30-min running anchor): %d' % len(det), flush=True)

det = det.merge(dep[['dep_key', 'latitude', 'longitude']], on='dep_key', how='left')
d0 = det.ts.dt.normalize()
sr, ss = sc.solar_window(d0 + pd.Timedelta(hours=12), det.latitude.to_numpy(), det.longitude.to_numpy())
srp, ssp = sc.solar_window(d0 - pd.Timedelta(hours=12), det.latitude.to_numpy(), det.longitude.to_numpy())
srn, _ = sc.solar_window(d0 + pd.Timedelta(hours=36), det.latitude.to_numpy(), det.longitude.to_numpy())
h = (det.ts - d0).dt.total_seconds().to_numpy() / 3600.0
dl = (ss - sr) % 24.0
st = np.empty(len(det))
in_day = (h >= sr) & (h < ss)
st[in_day] = (h[in_day] - sr[in_day]) / dl[in_day] * 12.0
post = h >= ss
nl_post = (srn[post] + 24.0 - ss[post])
st[post] = 12.0 + (h[post] - ss[post]) / nl_post * 12.0
pre = h < sr
nl_pre = (sr[pre] + 24.0 - ssp[pre])
st[pre] = 12.0 + (h[pre] + 24.0 - ssp[pre]) / nl_pre * 12.0
det['sun_t'] = st % 24.0
det['bin'] = np.clip((det.sun_t // 0.5).astype(int), 0, NB - 1)
det[['dep_key', 'common_name', 'ts', 'sun_t', 'bin', 'latitude', 'longitude']].to_parquet(f'{OUT}/det.parquet')
print('night share of detections: %.2f%%' % (100 * dp.IS_NIGHT_A[det.bin.to_numpy()].mean()), flush=True)
print('BUILD2 DONE', flush=True)
