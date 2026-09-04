#!/usr/bin/env python
"""Verify the built viewer the way a BROWSER would, not the way Python would.

The bug this exists to catch: Python's json.loads accepts bare NaN/Infinity;
JavaScript's JSON.parse does not. A payload that parses fine here can throw in
the browser and blank the entire page. Every check below is chosen to fail on
something a browser would reject.
"""
import json, re, sys, math

def check(path="diel_activity_viewer.html"):
    h = open(path).read()
    fails, warns = [], []

    m = re.search(r'<script id="payload" type="application/json">(.*?)</script>', h, re.S)
    if not m:
        return ["no payload block found"], []
    raw = m.group(1)

    # 1. STRICT JSON — the check that was missing. parse_constant fires on
    #    NaN / Infinity / -Infinity, which JSON.parse rejects outright.
    def boom(c):
        raise ValueError(f"bare {c} in payload — JSON.parse would throw and blank the page")
    try:
        D = json.loads(raw, parse_constant=boom)
    except ValueError as e:
        return [f"STRICT JSON: {e}"], []

    # 2. no bare non-finite tokens anywhere in the serialized text
    for tok in ("NaN", "Infinity"):
        n = len(re.findall(r'(?<![A-Za-z_"])' + tok + r'(?![A-Za-z_"])', raw))
        if n:
            fails.append(f"{n} bare {tok} token(s) in payload text")

    # 3. payload cannot terminate its own <script> element
    if "</script" in raw.lower():
        fails.append("payload contains a literal </script>")

    # 4. every surface cell has geometry AND a clickable curve
    for sp in D["surfaces"]:
        cells = set(D["surfaces"][sp][D["metrics"][0]["key"]])
        missing_geom = [c for c in cells if c not in D["paths"]]
        if missing_geom:
            fails.append(f"{sp}: {len(missing_geom)} cells lack geometry")
        cc = D.get("cell_curves", {}).get(sp, {})
        missing_curve = [c for c in cells if c not in cc]
        if missing_curve:
            fails.append(f"{sp}: {len(missing_curve)} cells lack a click curve")
        for met in D["surfaces"][sp]:
            for cid, rec in D["surfaces"][sp][met].items():
                if any(v is None or (isinstance(v, float) and not math.isfinite(v)) for v in rec):
                    fails.append(f"{sp}/{met}/{cid}: non-finite value in surface record")
                    break

    # 5. no path outside the declared viewBox
    W, H = D["vb"]
    for cid, d in D["paths"].items():
        xs = [float(a) for a, _ in re.findall(r'(-?\d+\.?\d*),(-?\d+\.?\d*)', d)]
        ys = [float(b) for _, b in re.findall(r'(-?\d+\.?\d*),(-?\d+\.?\d*)', d)]
        if xs and (min(xs) < -1 or max(xs) > W + 1 or min(ys) < -1 or max(ys) > H + 1):
            fails.append(f"path {cid} outside viewBox")
            break

    # 6. every DOM id the JS reaches for exists in the markup
    js = h.split("<script>")[-1]
    ids_js = set(re.findall(r"getElementById\('([^']+)'\)", js))
    ids_html = set(re.findall(r'id="([^"]+)"', h))
    for i in sorted(ids_js - ids_html):
        fails.append(f"JS references missing DOM id: {i}")

    # 7. self-contained: no external resources
    ext = re.findall(r'(?:src|href)="(https?://[^"]+)"', h)
    if ext:
        fails.append(f"external references present: {ext[:3]}")

    # 8. tables wrapped so they cannot overflow their card
    ntab = len(re.findall(r'<table', h))
    nwrap = len(re.findall(r'class="tw"', h))
    if nwrap < ntab:
        warns.append(f"{ntab} tables but only {nwrap} scroll wrappers")

    # 9. placeholder banner consistency
    banner = "synthetic surfaces for layout review" in h.lower()
    if D.get("placeholder") and not banner:
        fails.append("placeholder=True but no banner shown")
    if banner and not D.get("placeholder"):
        fails.append("banner shown but placeholder=False")

    return fails, warns

if __name__ == "__main__":
    path = sys.argv[1] if len(sys.argv) > 1 else "diel_activity_viewer.html"
    fails, warns = check(path)
    for w in warns:
        print("WARN:", w)
    for f in fails:
        print("FAIL:", f)
    print(f"\n{len(fails)} failures, {len(warns)} warnings")
    sys.exit(1 if fails else 0)
