"""Albers equal-area conic used by this project's maps, plus the inverse of the screen
transform stored in viewer_states.json.

viewer_states.json holds state boundary paths in SCREEN space, not lon/lat. It also holds the
transform that produced them, so the projection does not need reverse-engineering:
    screen_x = ox + s * (X - x0)
    screen_y = oy + s * (y1 - Y)
where X, Y come from albers() below. Use draw_states(ax, states_json) for outlines and
albers(lon, lat) for the data, and both land on the same axes.
"""
import numpy as np


def albers(lon, lat):
    """Albers equal-area conic, CONUS parameters. Returns x, y in arbitrary units."""
    lon = np.radians(np.asarray(lon, float)); lat = np.radians(np.asarray(lat, float))
    lon0, lat0 = np.radians(-96.0), np.radians(37.5)
    p1, p2 = np.radians(29.5), np.radians(45.5)
    n = 0.5 * (np.sin(p1) + np.sin(p2))
    C = np.cos(p1) ** 2 + 2 * n * np.sin(p1)
    rho = np.sqrt(C - 2 * n * np.sin(lat)) / n
    rho0 = np.sqrt(C - 2 * n * np.sin(lat0)) / n
    theta = n * (lon - lon0)
    return rho * np.sin(theta), rho0 - rho * np.cos(theta)


def screen_to_albers(sx, sy, transform):
    """Invert the stored screen transform. `transform` is the dict in viewer_states.json."""
    ox, oy = transform["ox"], transform["oy"]
    s, x0, y1 = transform["s"], transform["x0"], transform["y1"]
    return (np.asarray(sx, float) - ox) / s + x0, y1 - (np.asarray(sy, float) - oy) / s


def parse_svg_paths(path_strings, transform):
    """Parse 'M x,y x,y ...' path strings into (x, y) arrays in Albers space."""
    out = []
    for ps in path_strings:
        for sub in ps.split("M")[1:]:
            pts = []
            for tok in sub.replace("Z", "").strip().split():
                if "," not in tok:
                    continue
                a, b = tok.split(",")[:2]
                try:
                    pts.append((float(a), float(b)))
                except ValueError:
                    continue
            if len(pts) > 2:
                arr = np.array(pts)
                ax_, ay_ = screen_to_albers(arr[:, 0], arr[:, 1], transform)
                out.append((ax_, ay_))
    return out


def draw_states(ax, states_json, lw=0.5, color="#9aa0a6", zorder=1):
    """Draw lower-48 state outlines in Albers space onto `ax`. Returns the ring count."""
    for ax_, ay_ in parse_svg_paths(states_json["state_paths"], states_json["transform"]):
        ax.plot(ax_, ay_, lw=lw, color=color, zorder=zorder, solid_joinstyle="round")
    return len(states_json["state_paths"])
