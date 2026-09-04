"""Exact analytic per-bin effort in DOUBLE-ANCHORED sun time.

Sun time: sunrise -> 0, sunset -> 12, next sunrise -> 24. Day is stretched/compressed
onto [0,12), night onto [12,24). 48 half-hour sun-time bins.

Per-bin effort is NOT uniform in sun time even though cameras run continuously in
clock time: 12 sun-hours of day occupy (sunset-sunrise) clock hours, so each day bin
is worth (daylen)/24 clock hours while each night bin is worth (24-daylen)/24.

Full calendar days are handled with that closed form; the (at most two) partial days
at a deployment's boundaries are handled by exact interval intersection against each
bin's clock-time footprint.
"""
import numpy as np
import pandas as pd

NBIN = 48
BINW = 24.0 / NBIN          # 0.5 sun-hours
BC = (np.arange(NBIN) + 0.5) * BINW   # bin centres in sun time
IS_NIGHT_BIN = BC >= 12.0


def suntime(h, sr, ss):
    """Clock hour -> double-anchored sun time. Arrays."""
    h = np.asarray(h, float); sr = np.asarray(sr, float); ss = np.asarray(ss, float)
    dl = ss - sr
    nl = 24.0 - dl
    st = np.empty_like(h)
    day = (h >= sr) & (h < ss)
    st[day] = 12.0 * (h[day] - sr[day]) / dl[day]
    ev = h >= ss
    st[ev] = 12.0 + 12.0 * (h[ev] - ss[ev]) / nl[ev]
    mo = h < sr
    st[mo] = 12.0 + 12.0 * (h[mo] + 24.0 - ss[mo]) / nl[mo]
    return st


def bin_clock_edges(sr, ss):
    """(48,) clock-time start and end of each sun bin for ONE day. Night bins may wrap."""
    dl = ss - sr
    nl = 24.0 - dl
    b = np.arange(NBIN + 1) * BINW           # sun-time edges 0 .. 24
    ck = np.where(b <= 12.0,
                  sr + b * dl / 12.0,
                  ss + (b - 12.0) * nl / 12.0)
    return ck[:-1] % 24.0, ck[1:] % 24.0


def _overlap(a0, a1, b0, b1):
    """Overlap length of [a0,a1) with [b0,b1); both may wrap past 24."""
    segs_a = [(a0, a1)] if a1 > a0 else [(a0, 24.0), (0.0, a1)]
    segs_b = [(b0, b1)] if b1 > b0 else [(b0, 24.0), (0.0, b1)]
    tot = 0.0
    for x0, x1 in segs_a:
        for y0, y1 in segs_b:
            tot += max(0.0, min(x1, y1) - max(x0, y0))
    return tot


def partial_day_effort(sr, ss, a, b):
    """(48,) clock hours of effort for a day active only over clock interval [a,b)."""
    e0, e1 = bin_clock_edges(sr, ss)
    return np.array([_overlap(a, b, e0[i], e1[i]) for i in range(NBIN)])


def full_day_effort(dl):
    """(48,) clock hours for one fully-active calendar day with day length dl."""
    out = np.empty(NBIN)
    out[~IS_NIGHT_BIN] = dl / 24.0
    out[IS_NIGHT_BIN] = (24.0 - dl) / 24.0
    return out


def independent_events(ts, gap_min=30.0):
    """Running-anchor filter: keep a detection only if >= gap_min after the last KEPT one.

    ts must be sorted. Returns a boolean mask.
    """
    t = np.asarray(ts, dtype='datetime64[s]').astype('int64') / 60.0
    keep = np.zeros(len(t), bool)
    anchor = -np.inf
    for i in range(len(t)):
        if t[i] - anchor >= gap_min:
            keep[i] = True
            anchor = t[i]
    return keep
