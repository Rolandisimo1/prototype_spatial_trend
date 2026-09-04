"""Sunrise/sunset in LOCAL CLOCK hours — the units camera timestamps actually use.

Why this exists: computing SOLAR sunrise and comparing it to a CLOCK timestamp is
wrong by the site's longitude offset within its time zone (up to an hour) plus
daylight saving (another hour). That bias is systematic, it concentrates at the
western edge of each zone, and it makes normal crepuscular activity look
nocturnal — which reads as a broken camera clock.

Time-zone boundaries follow legal borders, not meridians, and Arizona does not
observe daylight saving. Only the IANA database gets these right, so the offset is
looked up per site and per date.
"""
import numpy as np
import pandas as pd
from zoneinfo import ZoneInfo
from timezonefinder import TimezoneFinder

_TF = TimezoneFinder()
_ZCACHE = {}


def site_timezone(lat, lon):
    """IANA zone name for one site, cached on rounded coordinates."""
    key = (round(float(lat), 2), round(float(lon), 2))
    if key not in _ZCACHE:
        _ZCACHE[key] = _TF.timezone_at(lat=key[0], lng=key[1])
    return _ZCACHE[key]


def utc_offsets(timestamps, lat, lon):
    """UTC offset in hours per record, honouring legal zones and DST."""
    ts = pd.DatetimeIndex(pd.to_datetime(timestamps))
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    zones = np.array([site_timezone(a, o) for a, o in zip(lat, lon)])
    out = np.full(len(ts), np.nan)
    for z in pd.unique(zones):
        m = zones == z
        if z is None:
            out[m] = np.round(lon[m] / 15.0)
            continue
        tzi = ZoneInfo(z)
        # The hour repeated at the DST fall-back is genuinely ambiguous; resolve to
        # standard time (ambiguous=False) and skip forward over the spring-gap hour.
        # A one-hour error on those rare records is far smaller than the zone-wide
        # bias this function exists to remove.
        loc = pd.DatetimeIndex(ts[m]).tz_localize(
            tzi, ambiguous=False, nonexistent="shift_forward")
        out[m] = np.array([o.total_seconds() / 3600.0 for o in loc.map(lambda x: x.utcoffset())])
    return out


def solar_window(timestamps, lat, lon):
    """(sunrise, sunset) in local clock hours. NOAA declination/equation of time."""
    ts = pd.DatetimeIndex(pd.to_datetime(timestamps))
    lat = np.asarray(lat, float)
    lon = np.asarray(lon, float)
    doy = ts.dayofyear.values.astype(float)
    ga = 2 * np.pi / 365.0 * (doy - 1)
    decl = (0.006918 - 0.399912 * np.cos(ga) + 0.070257 * np.sin(ga)
            - 0.006758 * np.cos(2 * ga) + 0.000907 * np.sin(2 * ga)
            - 0.002697 * np.cos(3 * ga) + 0.00148 * np.sin(3 * ga))
    eqt = 229.18 * (0.000075 + 0.001868 * np.cos(ga) - 0.032077 * np.sin(ga)
                    - 0.014615 * np.cos(2 * ga) - 0.040849 * np.sin(2 * ga))
    latr = np.radians(lat)
    cosw = np.clip((np.cos(np.radians(90.833)) - np.sin(latr) * np.sin(decl))
                   / (np.cos(latr) * np.cos(decl)), -1, 1)
    w = np.degrees(np.arccos(cosw))
    noon = 12.0 - lon / 15.0 - eqt / 60.0 + utc_offsets(ts, lat, lon)
    return (noon - w / 15.0) % 24, (noon + w / 15.0) % 24


def is_night(timestamps, lat, lon):
    """Boolean: is each timestamp between sunset and sunrise at its own site?"""
    ts = pd.DatetimeIndex(pd.to_datetime(timestamps))
    h = ts.hour.values + ts.minute.values / 60.0 + ts.second.values / 3600.0
    sr, ss = solar_window(ts, lat, lon)
    return (h < sr) | (h >= ss)
