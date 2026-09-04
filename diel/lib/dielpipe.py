"""Sun-anchored binning of camera-trap detections and the exact analytic per-bin effort offset.

Double anchoring: local clock time is mapped monotonically onto a 24 h "sun time" scale so that
sunrise always lands on REF_SR and sunset always on REF_SS. The map is piecewise linear with two
segments (the day arc and the night arc), so an hour of clock time occupies a DIFFERENT amount of
sun time depending on whether it falls in the day or the night arc, and by how much the local day
length differs from the reference. Effort is therefore not uniform across sun-time bins even though
it is uniform in clock time. Effort per bin is computed as the exact clock-time measure of the
preimage of each bin interval, which is what the count models must carry as an offset.
"""
import numpy as np
import pandas as pd
import sys, os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import solarclock as sc

NBIN = 48
BINW = 0.5
REF_SR = 6.2202
REF_SS = 17.9998
BC = (np.arange(1, NBIN + 1) - 0.5) * BINW          # bin centres, sun-time hours
BIN_LO = np.arange(NBIN) * BINW
BIN_HI = BIN_LO + BINW
IS_NIGHT = (BC < REF_SR) | (BC >= REF_SS)           # cvlib.R definition, 24 night / 24 day bins
IS_CREP = (np.abs(BC - REF_SR) <= 1.5) | (np.abs(BC - REF_SS) <= 1.5)
ANG = BC / 24.0 * 2 * np.pi
REF_DAY = REF_SS - REF_SR                           # 11.7796 h of sun time spans the day arc
REF_NIGHT = 24.0 - REF_DAY                          # 12.2204 h of sun time spans the night arc
ASO_MONTHS = (8, 9, 10)

# ---- anchoring actually used by this analysis: sunrise = 0, sunset = 12 (rule 4) ----
# Day arc is sun time [0, 12) = bins 0..23; night arc is [12, 24) = bins 24..47.
SR_A, SS_A = 0.0, 12.0
IS_NIGHT_A = BC >= SS_A
IS_DAY_A = ~IS_NIGHT_A
NOON_A = 6.0                                        # solar noon, midpoint of the day arc


def _circ_gap(x, c):
    d = np.abs(np.asarray(x, float) - c)
    return np.minimum(d, 24.0 - d)


IS_CREP_A = (_circ_gap(BC, SR_A) <= 1.0) | (_circ_gap(BC, SS_A) <= 1.0)   # 2 h bracketing each
IS_NOON_A = _circ_gap(BC, NOON_A) <= 2.0                                  # 4 h around solar noon


def metrics_anchored(v):
    """Diel-shape metrics on a 48-bin rate curve under sunrise=0/sunset=12 anchoring.

    peak_h and mean_h are CIRCULAR quantities on a 24 h sun-time scale and must never be
    averaged, differenced or plotted linearly (rule 7).
    """
    v = np.maximum(np.asarray(v, float), 0)
    s = v.sum()
    if not np.isfinite(s) or s <= 0:
        return dict(pct_noct=np.nan, pct_crep=np.nan, conc=np.nan, peak_h=np.nan,
                    mean_h=np.nan, act=np.nan, pct_noon=np.nan, breadth=np.nan)
    p = v / s
    C = float((p * np.cos(ANG)).sum()); Sn = float((p * np.sin(ANG)).sum())
    R = float(np.hypot(C, Sn))
    mu = (np.arctan2(Sn, C) % (2 * np.pi)) / (2 * np.pi) * 24.0
    mx = float(v.max())
    # activity breadth: exponential of Shannon entropy of the normalized curve, in bins
    q = p[p > 0]
    breadth = float(np.exp(-(q * np.log(q)).sum()))
    return dict(pct_noct=100 * float(p[IS_NIGHT_A].sum()),
                pct_crep=100 * float(p[IS_CREP_A].sum()),
                conc=R,
                peak_h=float(BC[int(np.argmax(v))]),
                mean_h=float(mu),
                act=float(v.mean() / mx) if mx > 0 else np.nan,
                pct_noon=100 * float(p[IS_NOON_A].sum()),
                breadth=breadth)


# ---------------------------------------------------------------- forward map
def clock_to_suntime(h, sr, ss):
    """Local clock hour -> sun-time hour by double anchoring (sunrise->REF_SR, sunset->REF_SS)."""
    h = np.asarray(h, float) % 24.0
    sr = np.asarray(sr, float)
    ss = np.asarray(ss, float)
    day_len = (ss - sr) % 24.0
    night_len = 24.0 - day_len
    in_day = ((h - sr) % 24.0) < day_len
    out = np.empty_like(h)
    # day arc: fraction of the way from sunrise to sunset
    fd = ((h - sr) % 24.0) / np.where(day_len > 0, day_len, np.nan)
    out_day = REF_SR + fd * REF_DAY
    # night arc: fraction of the way from sunset to next sunrise
    fn = ((h - ss) % 24.0) / np.where(night_len > 0, night_len, np.nan)
    out_night = (REF_SS + fn * REF_NIGHT) % 24.0
    out = np.where(in_day, out_day, out_night)
    return out % 24.0


def suntime_bin(t):
    b = np.floor(np.asarray(t, float) / BINW).astype(int)
    return np.clip(b, 0, NBIN - 1)


# ---------------------------------------------------------------- exact effort
def _overlap(a_lo, a_hi, b_lo, b_hi):
    return np.clip(np.minimum(a_hi, b_hi) - np.maximum(a_lo, b_lo), 0, None)


def bin_effort_one_day(day_len):
    """Clock hours falling in each of the 48 sun-time bins, for one camera-day.

    day_len: local day length in clock hours (sunset - sunrise). Returns shape (n, 48).
    The day arc occupies sun time [REF_SR, REF_SS) and carries day_len clock hours;
    the night arc occupies the complement and carries 24 - day_len clock hours.
    """
    day_len = np.atleast_1d(np.asarray(day_len, float))
    night_len = 24.0 - day_len
    kd = (day_len / REF_DAY)[:, None]        # clock hours per sun-hour, day arc
    kn = (night_len / REF_NIGHT)[:, None]    # clock hours per sun-hour, night arc
    lo = BIN_LO[None, :]
    hi = BIN_HI[None, :]
    # day-arc sun-time interval is a single contiguous span [REF_SR, REF_SS)
    ov_day = _overlap(lo, hi, REF_SR, REF_SS)
    # night arc wraps: [REF_SS, 24) plus [0, REF_SR)
    ov_night = _overlap(lo, hi, REF_SS, 24.0) + _overlap(lo, hi, 0.0, REF_SR)
    return ov_day * kd + ov_night * kn


def deployment_daylengths(lat, lon, dates):
    """Local day length (clock h) per (site, date)."""
    sr, ss = sc.solar_window(dates, lat, lon)
    return (ss - sr) % 24.0, sr, ss


# ---------------------------------------------------------------- detections
def independent_detections(df, gap_min=30.0, by=("dep_key", "common_name")):
    """30-minute RUNNING-ANCHOR thinning within camera x species.

    A record is kept when it is more than gap_min after the last KEPT record (the anchor moves to
    each kept record), not merely more than gap_min after the previous raw record.
    """
    d = df.sort_values(list(by) + ["ts"], kind="mergesort").reset_index(drop=True)
    grp = d.groupby(list(by), sort=False).ngroup().to_numpy()
    t = d["ts"].to_numpy("datetime64[s]").astype("int64") / 60.0
    keep = np.zeros(len(d), bool)
    anchor = -np.inf
    prev_g = -1
    for i in range(len(d)):
        if grp[i] != prev_g:
            anchor = -np.inf
            prev_g = grp[i]
        if t[i] - anchor > gap_min:
            keep[i] = True
            anchor = t[i]
    return d.loc[keep].reset_index(drop=True)


def curve_metrics(v):
    """Diel-shape metrics on a 48-bin rate curve. Peak/mean are CIRCULAR."""
    v = np.maximum(np.asarray(v, float), 0)
    s = v.sum()
    if not np.isfinite(s) or s <= 0:
        return dict(pct_noct=np.nan, pct_crep=np.nan, conc=np.nan,
                    peak_h=np.nan, mean_h=np.nan, act=np.nan, pct_noon=np.nan)
    p = v / s
    C = float((p * np.cos(ANG)).sum())
    Sn = float((p * np.sin(ANG)).sum())
    R = float(np.hypot(C, Sn))
    mu = (np.arctan2(Sn, C) % (2 * np.pi)) / (2 * np.pi) * 24.0
    mx = v.max()
    # share of activity within 2 h of solar noon on the sun-time scale
    noon = (REF_SR + REF_SS) / 2.0
    dnoon = np.abs(BC - noon)
    dnoon = np.minimum(dnoon, 24 - dnoon)
    return dict(pct_noct=100 * float(p[IS_NIGHT].sum()),
                pct_crep=100 * float(p[IS_CREP].sum()),
                conc=R,
                peak_h=float(BC[int(np.argmax(v))]),
                mean_h=float(mu),
                act=float(v.mean() / mx) if mx > 0 else np.nan,
                pct_noon=100 * float(p[dnoon <= 2.0].sum()))


def circ_diff_h(a, b):
    d = (np.asarray(a, float) - np.asarray(b, float)) % 24.0
    return np.where(d > 12, d - 24, d)


def circ_mean_h(x, w=None):
    x = np.asarray(x, float)
    m = np.isfinite(x)
    if w is None:
        w = np.ones_like(x)
    w = np.asarray(w, float)
    if m.sum() == 0:
        return np.nan
    a = x[m] / 24.0 * 2 * np.pi
    ww = w[m]
    C = (ww * np.cos(a)).sum()
    Sn = (ww * np.sin(a)).sum()
    return (np.arctan2(Sn, C) % (2 * np.pi)) / (2 * np.pi) * 24.0


def circ_sd_h(x):
    """Circular SD in hours (Fisher), for spread of peak times."""
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return np.nan
    a = x / 24.0 * 2 * np.pi
    R = np.hypot(np.cos(a).mean(), np.sin(a).mean())
    R = min(max(R, 1e-12), 1 - 1e-12)
    return np.sqrt(-2 * np.log(R)) / (2 * np.pi) * 24.0
