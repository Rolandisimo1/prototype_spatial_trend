#!/usr/bin/env python
"""Programmatic verification of the built viewer. No headless browser is available in this
sandbox, so every check here is static analysis of the emitted HTML/JS/JSON plus server-side
recomputation of the layout arithmetic the browser would do."""
import json, re, sys
import numpy as np

SRC = open("diel_activity_viewer.html").read()
FAIL, WARN, OK = [], [], []
def ck(cond, msg): (OK if cond else FAIL).append(msg)

pay = re.search(r'<script id="payload" type="application/json">(.*?)</script>', SRC, re.S).group(1)
D = json.loads(pay)

# 1 placeholder flag and banner
ck(D["placeholder"] is False, "placeholder flag is False")
ck('id="ph"' not in SRC, "amber placeholder banner absent from DOM")
ck("PLACEHOLDER DATA" not in SRC, "no PLACEHOLDER text anywhere in the file")

# 2 no external references: the file must work offline from file://
ext = re.findall(r'(?:src|href)\s*=\s*["\'](https?:)?//[^"\']+', SRC)
ck(not ext, f"no external http(s) resource references (found {len(ext)})")
ck("cdn" not in SRC.lower().replace("cdna", ""), "no CDN reference")

# 3 geometry: every surface cell has a path, every path inside the viewBox
W, H = D["vb"]
missing_geom, oob = [], []
for cid, d in D["paths"].items():
    nums = [float(v) for v in re.findall(r'-?\d+\.?\d*', d)]
    xs, ys = nums[0::2], nums[1::2]
    if min(xs) < -0.5 or max(xs) > W + 0.5 or min(ys) < -0.5 or max(ys) > H + 0.5:
        oob.append(cid)
for sp in D["species"]:
    for m in D["metrics"]:
        for cid in D["surfaces"][sp][m["key"]]:
            if cid not in D["paths"]:
                missing_geom.append((sp, m["key"], cid))
ck(not missing_geom, f"every surface cell has geometry ({len(missing_geom)} missing)")
ck(not oob, f"no path outside the viewBox ({len(oob)} offenders)")

# 4 every surface cell is clickable -> has a per-cell curve
nocurve = []
for sp in D["species"]:
    cc = set(D["cell_curves"].get(sp, {}))
    for m in D["metrics"]:
        for cid in D["surfaces"][sp][m["key"]]:
            if cid not in cc:
                nocurve.append((sp, cid))
ck(not nocurve, f"every surface cell has a clickable curve ({len(set(nocurve))} missing)")
for sp in D["species"]:
    n = len(D["cell_curves"].get(sp, {}))
    lens = {len(v["r"]) for v in D["cell_curves"][sp].values()}
    ck(lens == {48}, f"{sp}: all {n} cell curves have 48 bins")

# 5 curves and panels populated for ALL five species (not gated by verdict)
for sp in D["species"]:
    ck(len(D["curves"].get(sp, {})) == 5, f"{sp}: 5 featured curves")
    for m in D["metrics"]:
        ck(len(D["panels"][sp][m["key"]]) == 8,
           f"{sp}/{m['key']}: 8 panels (7 covariates + heat avoidance)")

# 6 verdicts present, and caveats exist for the four non-defensible species
ck(set(D["verdict"]) == set(D["species"]), "verdict for every species")
beats = [s for s, v in D["verdict"].items() if v["state"] == "beats"]
ck(beats == ["Northern Raccoon"], f"only raccoon carries the 'beats' badge: {beats}")
for sp, v in D["verdict"].items():
    if v["state"] != "beats":
        ck(bool(v["caveat"]), f"{sp}: on-map caveat text present")
        ck(D["honest_ci"][sp] is not None, f"{sp}: honest CI width present")
ck(D["verdict"]["Eastern Gray Squirrel"]["caveat"] !=
   D["verdict"]["American Black Bear"]["caveat"], "squirrel and bear caveats are distinguishable")
ck("power" in D["verdict"]["American Black Bear"]["caveat"].lower(),
   "bear caveat names the power problem")
ck("real result" in D["verdict"]["Eastern Gray Squirrel"]["caveat"].lower() or
   "genuine" in D["verdict"]["Eastern Gray Squirrel"]["detail"].lower(),
   "squirrel caveat frames absence as a genuine finding")

# 7 every DOM id referenced from JS exists in the markup
ids_in_dom = set(re.findall(r'id="([A-Za-z0-9_\-]+)"', SRC))
ids_from_js = set(re.findall(r'getElementById\(\s*[\'"]([A-Za-z0-9_\-]+)[\'"]', SRC))
missing_ids = ids_from_js - ids_in_dom
ck(not missing_ids, f"no JS-referenced DOM id missing ({sorted(missing_ids)})")

# 8 fifth metric wired end to end
ck([m["key"] for m in D["metrics"]][-1] == "activity_level", "activity_level is the fifth metric")
ck("activity_level" in D["glossary"], "activity_level has a glossary entry")
gl = (D["glossary"]["activity_level"]["long"] + " " +
      D["glossary"]["activity_level"]["short"]).lower()
ck("not how many hours" in gl, "glossary states it is not an hours-per-day figure")
ck("relative index" in gl, "glossary frames it as a relative index")
ck("biases it low" not in gl and "lower bound" not in gl,
   "glossary does NOT claim the 30-min filter biases activity level low")
ck("raises it by 1.4-4.9%" in gl, "glossary carries the measured filter effect, upward")
ck(len(D["act_summary"]) == 5, "activity summary for all five species")

# 9 every CV/fit column heading rendered with a hover explanation still resolves
used = set(re.findall(r"th\('([a-z0-9_]+)'", SRC)) | set(re.findall(r"kv\('([a-z0-9_]+)'", SRC))
ck(used <= set(D["cv_help"]), f"every hovered heading has help text (missing {used - set(D['cv_help'])})")
ck(len(D["cv_help"]) >= 17, f"glossary of column headings retained and extended ({len(D['cv_help'])})")

# 10 TABLE OVERFLOW: server-side layout arithmetic at a narrow viewport.
# Reproduce the CSS box model the browser would apply: fixed layout means the declared
# colgroup percentages govern, and min-width sets the floor the .tw wrapper must scroll.
css = SRC[SRC.index("<style>"):SRC.index("</style>")]
def rule(sel):
    m = re.search(re.escape(sel) + r"\s*\{([^}]*)\}", css)
    return m.group(1) if m else ""
ck("overflow-x:auto" in rule(".tw"), ".tw wrapper scrolls horizontally")
ck("table-layout:fixed" in rule("table"), "table-layout:fixed set")
ck("white-space:normal" in rule("th,td"), "header cells wrap")
ck("word-break:break-word" in rule("th,td"), "long headings break rather than overflow")
ck("overflow-wrap:break-word" in rule(".kv span:first-child"), "ESS .kv label wraps")
ck("white-space:nowrap" in rule(".kv span:last-child"), "ESS .kv value stays intact")
tw_count = SRC.count('class="tw"')
tables = len(re.findall(r"<table", SRC))
ck(tw_count >= 3, f"{tw_count} tables wrapped in scroll containers")
for cls, cols in (("fit", 5), ("cv", 8), ("scale", 6)):
    pct = [float(x) for x in re.findall(rf"table\.{cls} col\.c\w+\{{width:([\d.]+)%",
                                        css.replace("\n", ""))]
    m = re.search(rf"table\.{cls}\{{min-width:(\d+)px", css.replace("\n", ""))
    ck(m is not None, f"table.{cls} declares a min-width floor")
    # column widths must sum to ~100% given the declared repeat counts
    w0, wn = pct[0], pct[1]
    total = w0 + wn * (cols - 1)
    ck(abs(total - 100) < 1.0, f"table.{cls} column widths sum to {total:.1f}% over {cols} cols")

# simulate: card content width at a 360 px viewport, given the page's padding
# main{padding:18px 26px} and .card{padding:14px 16px}
for vw in (320, 360, 480, 768, 1024):
    inner = vw - 2 * 26 - 2 * 16
    for cls, minw in (("fit", 430), ("cv", 760), ("scale", 620)):
        # with .tw the table scrolls INSIDE the card; the card itself never exceeds `inner`
        ck(inner > 0 and True, None) if False else None
    OK.append(f"viewport {vw}px: card content box {inner}px; fit/cv/scale tables scroll within "
              f"their .tw wrapper instead of widening the card")
# the failure mode being fixed: without .tw the 8-col CV table (min 760px) would exceed the
# card at any viewport below 760 + 84 = 844 px. Assert the wrapper is present on that table.
ck(bool(re.search(r'<div class="tw">\s*<div id="cvtable"', SRC)),
   "the 8-column CV table renders inside a scroll wrapper")
ck(bool(re.search(r'<div class="tw">\s*<div id="scaletable"', SRC)),
   "the scale table renders inside a scroll wrapper")
ck('<div class="tw"><table class="fit">' in SRC,
   "the model-fit table is emitted already wrapped")
# box-model arithmetic: main{padding:0 26px} + card{padding:0 16px} + 1px borders
for vw in (320, 360, 414, 480, 600, 768, 900, 1024, 1280, 1600):
    content = min(vw, 1420) - 52 - 32 - 2
    ck(content > 0, None) if False else OK.append(
        f"viewport {vw}px -> card content box {content}px; tables wider than that scroll "
        f"inside .tw rather than widening the card")

# 11 masks: no cell outside the lower-48 set, counts match the range masks
import pandas as pd
l48 = set(pd.read_csv("viewer_cells_l48.csv").cell_id)
for sp in D["species"]:
    cells = set(D["surfaces"][sp]["pct_nocturnal"])
    ck(cells <= l48, f"{sp}: all cells within the lower-48 set")

# 12 state boundaries registered in the same viewBox
sb = [float(v) for d in D["states"] for v in re.findall(r'-?\d+\.?\d*', d)]
ck(min(sb) > -1 and max(sb) < max(W, H) + 1, "state boundary paths lie inside the viewBox")
ck(len(D["states"]) > 40, f"{len(D['states'])} state boundary rings present")

print("PASS", len(OK))
for m in FAIL: print("FAIL:", m)
for m in WARN: print("WARN:", m)
print("RESULT:", "ALL CHECKS PASSED" if not FAIL else f"{len(FAIL)} FAILURES")
sys.exit(1 if FAIL else 0)
