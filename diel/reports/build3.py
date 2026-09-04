"""Build the sun-anchored binned count/effort table over the UNION of the combined_* source files
and the swapped published archives.

Several projects (the Carolina datasets and a continental field-station network) were replaced
wholesale by their original published archives earlier in this project, which gave those cameras new
dep_keys. Those records live only in the swapped files, so the full 490-site universe requires the
union. Where a dep_key appears in both, the SWAPPED record wins: the swap was deliberate and the
archive copy is authoritative.

Timestamps in both files are LOCAL CLOCK time. The swapped archive was already converted with a
daylight-saving-aware rule, so no further longitude or flat-hour correction is applied to any record;
solarclock.py resolves the legal time zone and DST per site and date.

Anchoring: sunrise -> sun time 0, sunset -> 12, next sunrise -> 24, in 48 half-hour bins. Per-bin
effort is the exact clock measure of the intersection of the deployment's active interval with each
bin's clock preimage, so it is NOT uniform across sun-time bins (rule 4).
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
DCOLS = ['dep_key', 'project_id', 'deployment_id', 'start_date', 'end_date', 'latitude', 'longitude']


def _prep_dep(df):
    df = df.copy()
    df['dep_key'] = df.project_id.astype(str) + '|' + df.deployment_id.astype(str)
    df['t0'] = pd.to_datetime(df.start_date, format='ISO8601', utc=True).dt.tz_localize(None)
    df['t1'] = pd.to_datetime(df.end_date, format='ISO8601', utc=True).dt.tz_localize(None)
    df = df.dropna(subset=['t0', 't1', 'latitude', 'longitude'])
    return df[df.t1 > df.t0][['dep_key', 't0', 't1', 'latitude', 'longitude']]


# ---------------------------------------------------------------- deployments: swapped wins
dcomb = _prep_dep(pd.read_csv(SRC + 'combined_deployments_all.csv',
                              dtype={'project_id': str, 'deployment_id': str}))
dswap = _prep_dep(pd.read_csv("/Users/rwkays/.claude-science/orgs/4ac09562-3ef9-4985-8243-8cd6d9ebdcd2/artifacts/proj_47a049eaa5c6/0e9cb66e-1ee5-4cfc-8378-cc6bfd67f2f5/v3c5b7b6d_deployments_swapped.csv", dtype={'project_id': str, 'deployment_id': str}))
dcomb['src'] = 'combined'; dswap['src'] = 'swapped'
dep = pd.concat([dswap, dcomb], ignore_index=True).drop_duplicates('dep_key', keep='first')
print('deployments: combined %d, swapped %d, union %d (swapped preferred on %d shared keys)'
      % (dcomb.dep_key.nunique(), dswap.dep_key.nunique(), len(dep),
         len(set(dcomb.dep_key) & set(dswap.dep_key))), flush=True)

# ---------------------------------------------------------------- ASO effort windows
pieces = []
for yr in range(int(dep.t0.dt.year.min()), int(dep.t1.dt.year.max()) + 1):
    a = pd.Timestamp(f'{yr}-08-01'); b = pd.Timestamp(f'{yr}-11-01')
    m = (dep.t0 < b) & (dep.t1 > a)
    if not m.any():
        continue
    s = dep.loc[m].copy()
    s['a'] = s.t0.clip(lower=a); s['b'] = s.t1.clip(upper=b)
    pieces.append(s[s.b > s.a])
win = pd.concat(pieces, ignore_index=True)
win['nday'] = (win.b.dt.normalize() - win.a.dt.normalize()).dt.days + 1
print('ASO deployment-year windows %d over %d deployments' % (len(win), win.dep_key.nunique()), flush=True)

deps = pd.Index(win.dep_key.unique())
win['di'] = deps.get_indexer(win.dep_key)
EFF = np.zeros((len(deps), NB))
idx = np.repeat(np.arange(len(win)), win.nday.to_numpy())
offs = np.concatenate([np.arange(n) for n in win.nday.to_numpy()])
day0 = win.a.dt.normalize().to_numpy()[idx] + offs.astype('timedelta64[D]')
lat = win.latitude.to_numpy()[idx]; lon = win.longitude.to_numpy()[idx]
A = win.a.to_numpy()[idx]; Bv = win.b.to_numpy()[idx]; di = win.di.to_numpy()[idx]
HR = np.timedelta64(3600, 's')
for lo in range(0, len(day0), 400_000):
    sl = slice(lo, min(lo + 400_000, len(day0)))
    d0 = day0[sl]
    sr, ss = sc.solar_window(pd.DatetimeIndex(d0) + pd.Timedelta(hours=12), lat[sl], lon[sl])
    sr2, _ = sc.solar_window(pd.DatetimeIndex(d0) + pd.Timedelta(hours=36), lat[sl], lon[sl])
    t_sr = d0 + (sr * 3600).astype('timedelta64[s]')
    t_ss = d0 + (ss * 3600).astype('timedelta64[s]')
    t_sr2 = d0 + np.timedelta64(1, 'D') + (sr2 * 3600).astype('timedelta64[s]')
    dl = (t_ss - t_sr) / HR; nl = (t_sr2 - t_ss) / HR
    j = np.arange(25)
    edges = np.empty((len(d0), NB + 1), 'datetime64[s]')
    edges[:, :25] = t_sr[:, None] + ((dl[:, None] / 24 * j[None, :]) * 3600).astype('timedelta64[s]')
    edges[:, 25:] = t_ss[:, None] + ((nl[:, None] / 24 * j[None, 1:]) * 3600).astype('timedelta64[s]')
    aa = np.maximum(edges[:, :-1], A[sl][:, None])
    bb = np.minimum(edges[:, 1:], Bv[sl][:, None])
    np.add.at(EFF, di[sl], np.clip((bb - aa) / HR, 0, None))
    print('  effort %d/%d' % (sl.stop, len(day0)), flush=True)

r = EFF.max(1) / np.where(EFF.min(1) > 0, EFF.min(1), np.nan)
print('per-bin effort max/min within camera: median %.4f p99 %.4f max %.4f'
      % (np.nanmedian(r), np.nanpercentile(r, 99), np.nanmax(r)), flush=True)
np.save(f'{OUT}/EFF3.npy', EFF)
pd.DataFrame({'dep_key': deps}).to_parquet(f'{OUT}/eff3_index.parquet')

# ---------------------------------------------------------------- sequences: swapped wins
def _prep_seq(path, cols):
    s = pd.read_csv(path, usecols=cols, dtype={c: str for c in cols})
    s = s[s.common_name.isin(SPECIES13)].copy()
    s['dep_key'] = s.project_id + '|' + s.deployment_id
    s['ts'] = pd.to_datetime(s.start_time, format='ISO8601', utc=True).dt.tz_localize(None)
    return s.dropna(subset=['ts'])[['dep_key', 'common_name', 'ts']]


scomb = _prep_seq(SRC + 'combined_sequences_all.csv',
                  ['project_id', 'deployment_id', 'common_name', 'start_time'])
sswap = _prep_seq("/Users/rwkays/.claude-science/orgs/4ac09562-3ef9-4985-8243-8cd6d9ebdcd2/artifacts/proj_47a049eaa5c6/7a93fc21-69a0-4c6f-8e66-9a3a3727d988/v0c996e8a_sequences_swapped.csv", ['project_id', 'deployment_id', 'common_name', 'start_time'])
# a swapped deployment supersedes its combined counterpart entirely
swap_deps = set(dswap.dep_key)
scomb = scomb[~scomb.dep_key.isin(swap_deps)]
seq = pd.concat([sswap, scomb], ignore_index=True)
print('sequences: swapped %d + combined-not-swapped %d = %d' % (len(sswap), len(scomb), len(seq)), flush=True)

# One sequence can carry several rows (different species, or age/sex/group-size classes).
# Collapse to one row per sequence x species BEFORE the 30-minute anchor or counts inflate.
n0 = len(seq)
seq = seq.drop_duplicates(['dep_key', 'common_name', 'ts'])
print('collapsed %d duplicate sequence x species rows -> %d' % (n0 - len(seq), len(seq)), flush=True)

seq = seq[seq.ts.dt.month.isin(dp.ASO_MONTHS)]
seq = seq[seq.dep_key.isin(set(win.dep_key))]
print('ASO records on deployments with ASO effort: %d' % len(seq), flush=True)

det = dp.independent_detections(seq)
print('independent detections (30-min running anchor): %d' % len(det), flush=True)

# ---------------------------------------------------------------- sun-time bin per detection
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
st[post] = 12.0 + (h[post] - ss[post]) / (srn[post] + 24.0 - ss[post]) * 12.0
pre = h < sr
st[pre] = 12.0 + (h[pre] + 24.0 - ssp[pre]) / (sr[pre] + 24.0 - ssp[pre]) * 12.0
det['sun_t'] = st % 24.0
det['bin'] = np.clip((det.sun_t // 0.5).astype(int), 0, NB - 1)
det.to_parquet(f'{OUT}/det3.parquet')
print('night share of detections %.2f%%' % (100 * dp.IS_NIGHT_A[det.bin.to_numpy()].mean()), flush=True)
print('BUILD3 DONE', flush=True)
