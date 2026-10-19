#!/usr/bin/env python3
"""score_agency_concordance_windows.py -- CANONICAL builder of the model-vs-agency concordance tables (2026-10-19).

The concordance table cited in the paper (table_model_vs_agency_concordance.csv, first written 2026-10-13) had no
builder on disk, locally, on Hazel or in session history. This script replaces it: it regenerates that 5-year table
from its two sources and writes the 5-year vs 10-year comparison alongside. No new fits.

usage: python3 score_agency_concordance_windows.py [trend_table.csv]
  default trend table: trend_table/trend_table_v2u_national_scalar_33spp.csv (18-yr fits; fisher not in it -> "no fit")

Inputs
  agency_trends_combined.csv   one row per jurisdiction x species; trend_5yr_direction / trend_10yr_direction as the
                               agency states them (increasing / decreasing / stable / mixed / unknown)
  trend table                  model_direction = cam_direction (decline / increase / NOT RESOLVABLE),
                               model_camera_pct = camera_pct (camera-anchored year_beta, %/yr)

Rule
  rows      region_level == "state" only (US states). Provinces, sub-state units and regions are not scored.
  scored    direction != unknown
  majority  plurality of increasing / decreasing / (stable + mixed); ties -> stable/mixed; none if 0 scored
  concordance  model NOT RESOLVABLE -> model not resolvable; species absent from the trend table -> no model fit;
               agency stable/mixed or none -> agency no direction; same sign -> agree; opposite -> differ
  pop_estimate_records  state rows with metric_category "population/survey estimate" (informational, not scored)

Reconciliation with the 2026-10-13 published table (kept as table_model_vs_agency_concordance.published_20261013.csv)
  The states-only rule reproduces agency_n_scored, agency_majority and concordance for all 34 species. This resolves
  the pronghorn row first reconstructed with states + provinces: the 12th scored row was Saskatchewan (province,
  5-yr stable), which the published table did not score; table_agency_state_assessment_counts.csv likewise counts
  pronghorn as 11 STATES with a trend.
  One informational difference remains and could not be resolved: published pop_estimate_records for american_badger
  is 2, which are its BC and ON PROVINCE rows (both 5-yr unknown, so nothing scored changes), whereas the published
  counts for pronghorn and grey wolf exclude their province rows. No single level rule reproduces both, so the
  regenerated value (states only: badger 0) is canonical going forward.

Notes column (cases that must not be read as corroboration either way)
  black-tailed_jackrabbit  The only scored agency record at either window is Nebraska's July Rural Mail Carrier
                           Survey, whose index POOLS black-tailed and white-tailed jackrabbits. Neither its 5-yr +28%
                           nor its 10-yr -26% speaks to black-tailed jackrabbit; the 5-yr "differ" and the 10-yr
                           "agree" are both set aside.

Window
  The 5-year agency window is primary (METHODS_agency_trend_search.md): the camera record is mostly 2019-2025, and
  the windows reverse for pronghorn (Montana 2019-2021). The 10-year column is reported as a sensitivity check; it does
  not remove disagreement (5-yr: 1 differ; 10-yr: 3 differ).
"""
import csv, collections, os, sys

AGENCY = "agency_trends_combined.csv"
TREND = sys.argv[1] if len(sys.argv) > 1 else "trend_table/trend_table_v2u_national_scalar_33spp.csv"
PUBLISHED = "table_model_vs_agency_concordance.published_20261013.csv"
NOTES = {
    "black-tailed_jackrabbit": "agency record set aside at both windows: Nebraska's roadside index pools black-tailed "
                               "and white-tailed jackrabbits (5-yr +28%, 10-yr -26%); not evidence for or against "
                               "the model",
}

R = list(csv.DictReader(open(AGENCY)))
M = {x["species"]: x for x in csv.DictReader(open(TREND))}
key = lambda x: x["sp_key"] or x["species"].replace(" ", "_")
species = sorted({key(x) for x in R} | set(M))


def majority(c):
    i, d, s = c["increasing"], c["decreasing"], c["stable"] + c["mixed"]
    if i + d + s == 0:
        return "none", (0, 0, 0)
    m = "increase" if i > d and i > s else "decline" if d > i and d > s else "stable/mixed"
    return m, (i, s, d)


def conc(model, agency):
    if model == "NOT RESOLVABLE": return "model not resolvable"
    if model == "no fit": return "no model fit"
    if agency in ("stable/mixed", "none"): return "agency no direction"
    return "agree" if model == agency else "differ"


five, both = [], []
for sp in species:
    rows = [x for x in R if key(x) == sp and x["region_level"] == "state"]
    m = M.get(sp)
    md = m["cam_direction"] if m else "no fit"
    mp = m["camera_pct"] if m else ""
    agw = {}
    for w in ("5yr", "10yr"):
        c = collections.Counter(x["trend_%s_direction" % w] for x in rows if x["trend_%s_direction" % w] != "unknown")
        agw[w] = majority(c)
    pop = sum(1 for x in rows if x["metric_category"] == "population/survey estimate")
    (m5, (i5, s5, d5)), (m10, (i10, s10, d10)) = agw["5yr"], agw["10yr"]
    five.append({"species": sp, "model_direction": md, "model_camera_pct": mp, "agency_n_scored": i5 + s5 + d5,
                 "agency_majority": m5, "pop_estimate_records": pop, "concordance": conc(md, m5),
                 "note": NOTES.get(sp, "")})
    both.append({"species": sp, "model_direction": md, "model_camera_pct": mp,
                 "agency_5yr_n_scored": i5 + s5 + d5, "agency_5yr_inc_stab_dec": "%d/%d/%d" % (i5, s5, d5),
                 "agency_5yr_majority": m5, "concordance_5yr": conc(md, m5),
                 "agency_10yr_n_scored": i10 + s10 + d10, "agency_10yr_inc_stab_dec": "%d/%d/%d" % (i10, s10, d10),
                 "agency_10yr_majority": m10, "concordance_10yr": conc(md, m10), "note": NOTES.get(sp, "")})

# reproduce check against the published 2026-10-13 table
if os.path.exists(PUBLISHED):
    pub = {x["species"]: x for x in csv.DictReader(open(PUBLISHED))}
    cols = ("model_direction", "model_camera_pct", "agency_n_scored", "agency_majority", "pop_estimate_records",
            "concordance")
    diffs = [(r["species"], c, pub[r["species"]][c], str(r[c])) for r in five for c in cols
             if r["species"] in pub and pub[r["species"]][c] != str(r[c])]
    assert set(pub) == {r["species"] for r in five}, "species set differs from published table"
    print("vs published 2026-10-13 table: %d species, %d differing cells" % (len(pub), len(diffs)))
    for d in diffs:
        print("  %s %s: published %s -> regenerated %s" % d)
    assert all(d[:2] == ("american_badger", "pop_estimate_records") for d in diffs), "unexpected difference"

for fn, rows in (("table_model_vs_agency_concordance.csv", five),
                 ("table_model_vs_agency_concordance_by_window.csv", both)):
    with open(fn, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), lineterminator="\n"); w.writeheader(); w.writerows(rows)
    print("wrote", fn)

for w in ("5yr", "10yr"):
    print(w, dict(collections.Counter(r["concordance_%s" % w] for r in both)))
print("species whose concordance changes with the window:")
for r in both:
    if r["concordance_5yr"] != r["concordance_10yr"]:
        print("  %-26s model %-8s %-24s | 5yr %-6s %-20s | 10yr %-6s %-20s %s" % (
            r["species"], r["model_direction"], r["model_camera_pct"], r["agency_5yr_inc_stab_dec"],
            r["concordance_5yr"], r["agency_10yr_inc_stab_dec"], r["concordance_10yr"], "[note]" if r["note"] else ""))
