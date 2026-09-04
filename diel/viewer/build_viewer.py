#!/usr/bin/env python
"""
Build the interactive diel-activity viewer as ONE self-contained HTML file.

No external JS/CSS/CDN and no tile server: the map is inline SVG drawn from the
100 km grid polygons, so the page works offline, from a file:// path, or on any
static host (GitHub Pages, lab web space).

Data contract (viewer_contract.json declares it; see SCHEMAS below):
  surfaces : species, metric, cell_id, value, ci_lo, ci_hi, ci_width, in_mask
  curves   : species, location, bin, sun_hour, rate, lo, hi
  panels   : species, metric, panel_id, title, x, y, lo, hi, xlabel, ylabel, caption
  fit_model      : per species x model fit numbers
  fit_validation : per species x cv_scheme, incl. null baseline
  fit_ess        : per species effective sample size + CI inflation

Usage:
  python build_viewer.py --placeholder      # build with PLACEHOLDER_* inputs
  python build_viewer.py --real             # build with the real Phase 1-3 outputs
"""
import argparse, json, os
import numpy as np
import pandas as pd

OUT = "diel_activity_viewer.html"

# ---------------------------------------------------------------- projection
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


# ---------------------------------------------------------------- colour ramps
def _ramp(stops):
    def f(t):
        t = 0.0 if not np.isfinite(t) else min(max(t, 0.0), 1.0)
        k = t * (len(stops) - 1)
        i = int(np.floor(k)); i = min(i, len(stops) - 2)
        u = k - i
        a, b = stops[i], stops[i + 1]
        return "#%02x%02x%02x" % tuple(int(round(a[j] + (b[j] - a[j]) * u)) for j in range(3))
    return f

# viridis-like, colour-blind safe, for metric values
RAMP_VALUE = _ramp([(68,1,84),(72,40,120),(62,74,137),(49,104,142),(38,130,142),
                    (31,158,137),(53,183,121),(109,205,89),(180,222,44),(253,231,37)])
# single-hue for uncertainty: light = precise, dark = uncertain
RAMP_CI = _ramp([(247,247,247),(217,217,217),(189,189,189),(150,150,150),
                 (115,115,115),(82,82,82),(37,37,37)])
# CYCLIC ramp for PEAK TIMING, which is an hour on a 24 h circle, not a magnitude.
# This is a correctness fix, not decoration: coyote peaks near midnight, so 72% of its cells
# fall in 0-3 h or 21-24 h. On a linear ramp a cell at 23.8 h and one at 0.2 h -- 24 minutes
# apart in reality -- land at OPPOSITE ends of the colour scale and read as a maximal
# difference. A cyclic ramp closes the circle so they render nearly identically. The first
# and last stops are the same colour by construction.
RAMP_CYC = _ramp([(60,45,120),(35,95,165),(30,150,160),(90,185,110),(200,195,60),
                  (225,140,55),(200,70,80),(140,55,120),(60,45,120)])
CYCLIC_METRICS = {"peak_hour"}


# ---------------------------------------------------------------- plain-language help
# Covariate names as they should appear to a reader who is not on this project.
GLOSSARY = {
 "ntl": {"label": "Artificial night light",
   "short": "Brightness of artificial lighting at night (satellite radiance, 1 km).",
   "long": "Measured by the VIIRS satellite sensor as average night-time radiance in a 1 km "
           "neighbourhood of each camera. This is our measure of urbanisation, chosen over "
           "impervious surface because light is a direct mechanism that can change when an animal "
           "is active, not merely a correlate of development."},
 "crop": {"label": "Cropland",
   "short": "Percent of land within 1 km that is cultivated crops.",
   "long": "From the National Land Cover Database (cultivated crops class). Held separate from "
           "pasture and from urbanisation throughout, because row-crop agriculture, grazing land "
           "and cities act on animals through different mechanisms."},
 "pasture": {"label": "Pasture and hay",
   "short": "Percent of land within 1 km that is pasture or hayfield.",
   "long": "From the National Land Cover Database (pasture/hay class) - semi-natural grazing and "
           "forage land, distinct from row-crop agriculture."},
 "canopy": {"label": "Tree canopy cover",
   "short": "Percent tree canopy within 1 km.",
   "long": "From the National Land Cover Database tree-canopy-cover product; a continuous measure "
           "of how wooded the surroundings are."},
 "tmax": {"label": "Daytime high temperature during the survey",
   "short": "Average daily high temperature (degrees C) on the days each camera was actually running.",
   "long": "Not a climate average. For every camera we took the daily maximum temperature from "
           "Daymet for each day that camera was actually recording, then averaged over that "
           "camera's own survey period. So a camera running through a hot September gets a hot "
           "value, and a nearby camera running in a cool October gets a cooler one. This is what "
           "lets us ask whether animals avoid the hottest hours when conditions are genuinely hot, "
           "rather than only whether hot places differ from cool places."},
 "elev": {"label": "Elevation",
   "short": "Height above sea level (m) at the camera.",
   "long": "From the Copernicus 30 m digital elevation model."},
 "pop": {"label": "Human population density",
   "short": "People per pixel within 1 km (GHSL).",
   "long": "Kept separate from night light because dispersed rural housing and dense city cores "
           "affect wildlife differently. Keeping them apart matters: within a species the two push "
           "in OPPOSITE directions (for coyote, night light -6.9 percentage points on % nocturnal "
           "across its range while population is +7.2), so a combined 'human footprint' index would "
           "cancel them out and report nothing."},
 "activity_level": {"label": "Activity level (proportion of day active)",
   "short": "Where animals are active for more of the day - NOT how many hours they are active.",
   "long": "An analogue of the Rowcliffe et al. (2014) activity-level statistic: the mean of the "
           "fitted daily activity curve divided by its maximum. A flat curve (active at a similar "
           "rate around the clock) scores near 1; a sharply peaked curve scores low. It is the one "
           "metric here that is NOT purely about timing - the other four are computed on curves "
           "normalised to sum to 1, so an animal active four hours a day and one active twenty "
           "hours a day score identically on them if the timing matches. "
           "<b>Read it as a relative index across space within one species, not as an absolute "
           "figure.</b> Three reasons: Rowcliffe's estimator assumes an animal is detectable "
           "whenever it is active, which camera placement (on-trail vs off-trail), detection-zone "
           "geometry and reduced infrared range at night all violate; array-level differences in "
           "detectability mean absolute levels are not comparable between projects; and we compute "
           "it from binned fitted rates rather than Rowcliffe's circular-kernel fit to event times, "
           "so it is an analogue of the published statistic rather than the identical estimator. "
           "Absolute values are therefore not directly comparable to published Rowcliffe figures. "
           "The 30-minute independence filter applied upstream raises it by 1.4-4.9% depending on "
           "species (largest for bear and deer), which is small next to both the spatial variation "
           "mapped here and the credible intervals."},
}

# Every cross-validation / precision column, explained for a non-statistician.
CV_HELP = {
 "cv_scheme": ("How the test data were held back",
   "Spatial block: we hide whole geographic regions (300-500 km squares), fit the model on what is "
   "left, and try to predict the hidden regions. Held-out cameras end up 440-530 km from the nearest "
   "camera used for fitting, so this asks whether the map works in country we never sampled."),
 "regions_tested": ("Regions tested",
   "How many held-out geographic blocks cleared the minimum-detection threshold and could be scored. "
   "For black bear only 17 of 39 blocks were scorable, which is why its result is underpowered "
   "rather than negative."),
 "skill_unit": ("Skill vs one nationwide rhythm",
   "How much the map reduces prediction error in withheld REGIONS compared with using a single "
   "nationwide daily rhythm for the species. Above 0 means the map helps; 0 means it adds nothing. "
   "Every held-out region counts equally here, so one very large camera array cannot dominate."),
 "skill_ci": ("95% confidence interval",
   "Range of skill values consistent with the data, from 4000 bootstrap resamples of the held-out "
   "regions. If this interval contains 0, the map has not been shown to beat one nationwide rhythm. "
   "Every interval on this page contains 0."),
 "folds_pos": ("Folds positive",
   "How many of the five held-out region sets the map predicted better than one nationwide rhythm. "
   "Four out of five sounds encouraging but is well within what chance produces when the true skill "
   "is zero."),
 "folds": ("Number of held-back subsets",
   "The data were split into this many parts; each part took a turn as the hidden test set."),
 "shape_correlation": ("Predicted vs observed daily rhythm (r)",
   "How closely the predicted 24-hour activity shape matched the observed one at held-back sites. "
   "1.0 is perfect, 0 is no better than a flat line. This scores the SHAPE of the rhythm, not how "
   "many animals were photographed."),
 "null_shape_correlation": ("Same score using one nationwide rhythm instead of a map",
   "The score you get by ignoring location entirely: predict a single average daily rhythm for the "
   "species and apply it everywhere. The map is only worth having if it beats this number."),
 "metric_mae": ("Typical error in the metric",
   "Average absolute error in the mapped metric at held-back sites, in the metric's own units "
   "(e.g. percentage points for % nocturnal). Smaller is better."),
 "null_metric_mae": ("Typical error using one nationwide rhythm instead of a map",
   "The same error, for the single nationwide average rhythm. If this is smaller than the model's "
   "error, the map is not adding information."),
 "beats_null": ("Is the map better than one nationwide rhythm?",
   "Yes means the map predicted withheld regions better than assuming the species behaves the same "
   "way everywhere. No means you would do just as well ignoring location. This is the single most "
   "important number on the page: it is the test of whether a map is worth having at all."),
 "nominal_n_deployments": ("Cameras",
   "The raw number of camera deployments for this species."),
 "n_arrays": ("Camera arrays",
   "Cameras are deployed in clusters (median spacing about 50 m). Cameras in one cluster see the "
   "same trail, habitat and often the same individual animals, so they are not independent "
   "samples of behaviour."),
 "effective_n": ("Effective sample size",
   "How many genuinely independent samples the clustered design is worth. It is far closer to the "
   "number of arrays than to the number of cameras - which is why we do not treat every camera as "
   "an independent data point."),
 "ci_width_naive": ("Error bar if cameras were treated as independent",
   "Width of the uncertainty interval you would report if you wrongly assumed every camera was an "
   "independent sample. Misleadingly narrow."),
 "ci_width_hierarchical": ("Honest error bar",
   "Width of the uncertainty interval from the model that accounts for camera clustering. This is "
   "the one we report."),
 "ci_inflation_factor": ("How much wider the honest error bar is",
   "The ratio of the two. A value of 2 means naive analysis would have claimed twice the precision "
   "the data actually support."),
 "dev_explained": ("Variation explained",
   "Share of the variation in detection counts the model accounts for."),
 "edf_diel": ("Flexibility used for the daily rhythm",
   "Effective degrees of freedom in the time-of-day curve - higher means a more wiggly rhythm was "
   "needed to fit the data."),
 "edf_spatial": ("Flexibility used for the map",
   "Effective degrees of freedom in the spatial surface - higher means activity varied over "
   "finer geographic detail."),
 "dispersion": ("Spread relative to expectation",
   "About 1 means the counts vary as the statistical model assumes; well above 1 means extra "
   "variability, well below 1 means less."),
 "skill_noct": ("Block-CV skill (% nocturnal)",
   "How much the map reduces prediction error in withheld REGIONS compared with using one "
   "nationwide rhythm. Above 0 means the map helps; at or below 0 means it does not. This is the "
   "number the verdict badge is based on."),
 "honest_ci_width_pp": ("Honest interval width",
   "Width, in percentage points, of the credible interval for % nocturnal at a NEW array - the "
   "prediction situation that matters if you want to use the map somewhere we did not sample. It "
   "accounts for camera clustering and for array-to-array variation."),
 "n_units": ("Scorable blocks or arrays",
   "How many held-out units actually cleared the minimum-detection threshold and could be scored. "
   "For black bear only 17 of 39 geographic blocks were scorable, which is why its result is "
   "underpowered rather than negative."),
 "n_events": ("Detections",
   "Number of independent detection events (30-minute independence filter) contributing."),
 "med_dist_to_train_km": ("Distance from training cameras",
   "Median distance from a held-out camera to the NEAREST camera used to fit the model. In the "
   "block scheme this is 440-530 km, so the test genuinely asks about unsampled country."),
 "mi_0_50": ("Similarity at 0-50 km",
   "Moran's I of observed % nocturnal between camera arrays less than 50 km apart. Near 0 means "
   "nearby arrays are no more alike than distant ones; higher means real local structure."),
 "mi_50_200": ("Similarity at 50-200 km",
   "The same measure for arrays 50-200 km apart."),
 "mi_200_800": ("Similarity at 200-800 km",
   "The same measure for arrays 200-800 km apart. Values near zero mean that at the scale a "
   "continental map can resolve, there is nothing left to map."),
 "dev_space_k10": ("Space alone, coarse map",
   "Deviance explained by location alone with a deliberately coarse spatial spline."),
 "dev_space_k120": ("Space alone, fine map",
   "Deviance explained by location alone with a much more flexible spatial spline. That this keeps "
   "rising means the spatial signal lives at finer scales than a continental map resolves."),
 "activity_level": ("Activity level",
   "Mean of the fitted daily activity curve divided by its maximum - relative across space within "
   "a species, not an absolute hours-per-day figure."),
}

# Metric key in the viewer  ->  metric key in skill_bootstrap_cis.csv. Activity level was not among
# the four metrics cross-validated, and the page says so rather than borrowing another metric's number.
BOOT_METRIC = {"pct_nocturnal": "noct", "crepuscular": "crep",
               "concentration": "conc", "peak_hour": "peak", "activity_level": None}



# ---------------------------------------------------------------- covariate take-home text
# The PI's instruction: every caption OPENS with the relationship in plain language, before any
# statistic, and only panels whose credible band actually excludes zero are shown.
SPP_PLURAL = {"Northern Raccoon": "Raccoons", "White-tailed Deer": "White-tailed deer",
              "Coyote": "Coyotes", "Eastern Gray Squirrel": "Eastern gray squirrels",
              "American Black Bear": "Black bears"}
COV_PHRASE = {
 "cov_ntl":     ("where nights are brighter with artificial light", "the dimmest tenth of sites", "the brightest tenth"),
 "cov_crop":    ("where more of the surrounding land is row crops", "the least cropped sites", "the most cropped"),
 "cov_pasture": ("where more of the land is pasture or hayfield", "the least pastured sites", "the most pastured"),
 "cov_pop":     ("where more people live nearby", "the least populated sites", "the most populated"),
 "cov_canopy":  ("where tree cover is denser", "the most open sites", "the most wooded"),
 "cov_tmax":    ("where daytime highs during the survey were hotter", "the coolest surveys", "the hottest"),
 "cov_elev":    ("at higher elevations", "the lowest sites", "the highest"),
}
MET_PHRASE = {
 "pct_nocturnal": ("more nocturnal", "less nocturnal", "% nocturnal"),
 "crepuscular":   ("more concentrated around dawn and dusk", "less concentrated around dawn and dusk", "crepuscular index"),
 "concentration": ("more tightly peaked in their daily rhythm", "more spread through the day", "activity concentration"),
 "peak_hour":     ("peaking later relative to sunrise", "peaking earlier relative to sunrise", "peak timing"),
 "activity_level":("active over more of the day", "active over less of the day", "activity level"),
}


def panel_significant(g):
    """Does this panel's credible band exclude zero effect ANYWHERE across the covariate range?

    For a fitted curve, 'no effect' means a flat line, so the test is whether the band at some
    point sits entirely above the band at another point (max of the lower bounds above min of the
    upper bounds). For the three-estimator heat panel the relevant estimator is the within-array
    one, which is the identified contrast; its own interval must exclude zero.
    """
    kind = str(g.panel_kind.iloc[0]) if "panel_kind" in g.columns else "curve"
    lo = g.lo.values.astype(float); hi = g.hi.values.astype(float)
    if kind == "estimators":
        return bool(lo[2] > 0 or hi[2] < 0)
    return bool(np.nanmax(lo) > np.nanmin(hi))


def panel_takehome(sp, metric, pid, g):
    """One plain-language sentence describing the relationship, to OPEN the caption."""
    y = g.y.values.astype(float); lo = g.lo.values.astype(float); hi = g.hi.values.astype(float)
    who = SPP_PLURAL.get(sp, sp)
    if pid == "heat_avoid":
        w, wl, wh = y[2], lo[2], hi[2]
        if wl > 0:
            return (f"{who} shift activity INTO the three hottest hours of the day when the survey "
                    f"period is hotter (+{w:.2f} percentage points, 90% CI {wl:.2f} to {wh:.2f}), "
                    f"comparing the same cameras across their own cool and hot days.")
        if wh < 0:
            return (f"{who} pull activity OUT of the three hottest hours when the survey period is "
                    f"hotter ({w:.2f} percentage points, 90% CI {wl:.2f} to {wh:.2f}), comparing the "
                    f"same cameras across their own cool and hot days.")
        return (f"{who} show no clear heat avoidance once you compare the same cameras across their "
                f"own cool and hot days ({w:+.2f} percentage points, 90% CI {wl:.2f} to {wh:.2f}). "
                f"The naive estimate ({y[0]:+.2f}) mostly reflects hot PLACES differing from cool "
                f"places rather than animals responding to heat.")
    phrase, lowlab, hilab = COV_PHRASE.get(pid, ("across this predictor", "the low end", "the high end"))
    up, down, mlab = MET_PHRASE.get(metric, ("higher", "lower", metric))
    y0, y1 = float(y[0]), float(y[-1])
    d = y1 - y0
    rng = float(y.max() - y.min())
    imax, imin = int(y.argmax()), int(y.argmin())
    interior_peak = 0 < imax < len(y) - 1
    interior_dip = 0 < imin < len(y) - 1
    # A hump or a dip is NOT a direction. Asserting one from the endpoints would misdescribe the
    # relationship, so those panels get a shape sentence instead of a direction sentence.
    if abs(d) < 0.6 * rng and (interior_peak or interior_dip):
        noun = str(g.title.iloc[0])
        noun = noun[0].lower() + noun[1:]
        if interior_peak:
            return (f"For {who.lower()}, {mlab} peaks at intermediate levels of {noun} rather than "
                    f"rising or falling steadily: {y0:.2f} at {lowlab}, up to {y.max():.2f} part-way "
                    f"through the range, back to {y1:.2f} at {hilab}.")
        return (f"For {who.lower()}, {mlab} dips in the middle of the {noun} range rather than "
                f"moving steadily in one direction: {y0:.2f} at {lowlab}, down to {y.min():.2f} "
                f"part-way through, back up to {y1:.2f} at {hilab}.")
    direction = up if d > 0 else down
    return (f"{who} are {direction} {phrase}: {mlab} runs from {y0:.2f} at {lowlab} to {y1:.2f} at "
            f"{hilab}, a change of {d:+.2f} across the range.")



def _json_strict(obj):
    """Serialize to STRICTLY VALID JSON for the browser.

    json.dumps emits bare NaN/Infinity by default. JavaScript's JSON.parse REJECTS
    those tokens: one anywhere in the payload throws at parse time and the ENTIRE
    page renders blank -- tabs empty, map missing -- with no visible error. Python's
    json.loads accepts them, so a Python-side parse check does NOT catch this.
    This shipped once (six NaNs in the scale table). Non-finite floats become null.
    """
    import math

    def clean(o):
        if isinstance(o, float):
            return None if (math.isnan(o) or math.isinf(o)) else o
        if isinstance(o, dict):
            return {k: clean(v) for k, v in o.items()}
        if isinstance(o, (list, tuple)):
            return [clean(v) for v in o]
        return o

    s = json.dumps(clean(obj), allow_nan=False)   # raises rather than emitting NaN
    return s.replace("</", "<\\/")               # never close the script element early


# ---------------------------------------------------------------- support tiers for the MAP TAB
# The map tab must not present five equivalent maps. Ordering and marking come from
# covariate_model_cv_fixed.csv, columns beats_fix / skill_fix (the corrected constrained model;
# the *_raw columns are superseded and are not read anywhere in this build).
#
# Tier rule, stated on the page as well as here:
#   strong  : beats a single nationwide curve on 4 or 5 of its 5 measures  -> squirrel (4, +0.255)
#   partial : 1-3 measures, on a sample large enough to rank              -> coyote (2), raccoon (1),
#                                                                           deer (1, mean skill < 0)
#   weak    : any count, but resting on too few held-out arrays to rank   -> bear (24 arrays, 4 folds)
# Bear clears 2 of 5 but on 24 arrays over 4 folds, so it is marked weakly supported rather than
# being ranked alongside squirrel. Within a tier, order is by number of measures beaten, then by
# mean skill.
SUPPORT_TIER_LABEL = {
    "strong":  "Defensible national map",
    "partial": "Partly supported",
    "weak":    "Too little held-out data to judge",
}
SUPPORT_GLYPH = {"strong": "\u25c9", "partial": "\u25d4", "weak": "\u25cb"}
WEAK_MIN_ARRAYS = 40          # below this the species is placed in the 'weak' tier regardless of hits

SUPPORT_NOTE = {
 "Eastern Gray Squirrel":
   "The one species with a defensible national map. It beats a single nationwide curve on 4 of its "
   "5 measures, with a mean skill of +0.255 \u2014 by far the largest margin here. The reason is not "
   "that squirrels vary regionally: essentially none of their between-array variation is regional "
   "(0.5%). They respond to conditions we actually measured, so a map built from those conditions "
   "transfers to country we never sampled.",
 "Coyote":
   "Beats a nationwide curve on 2 of 5 measures (% nocturnal and activity level), with a small mean "
   "margin of +0.035 across 253 camera arrays. Coyotes do have regional structure (42% of their "
   "between-array variation), but it is not captured by the conditions in this model, so most of "
   "their measures do not transfer.",
 "Northern Raccoon":
   "Beats a nationwide curve on 1 of 5 measures (activity concentration), mean margin +0.035 across "
   "386 arrays. None of the raccoon's between-array variation is regional \u2014 it is all local, "
   "below 25 km \u2014 so there is little for a national surface to carry.",
 "White-tailed Deer":
   "The largest sample here (532 arrays) and the weakest result: 1 of 5 measures beaten and a "
   "mean skill of \u22120.044, i.e. on average worse than assuming deer keep the same hours "
   "everywhere. Deer have the most regional structure of any species (59% of between-array "
   "variation), so this is not an absence of pattern \u2014 it is a pattern driven by something this "
   "model does not measure.",
 "American Black Bear":
   "Beats a nationwide curve on 2 of 5 measures, but on only 24 held-out arrays across 4 folds \u2014 "
   "a fifth to a tenth of every other species. Treat this as too little data to judge, not as "
   "evidence comparable to the gray squirrel's.",
}

# One sentence per species for the deviation readout, keyed off the recomputed decomposition.
DEV_NOTE = {
 "White-tailed Deer":
   "Almost everything you see in a deer curve is the nationwide deer rhythm. Very little of it is "
   "about the place you clicked.",
 "Coyote":
   "Most of a coyote curve is the nationwide coyote rhythm; about a fifth is specific to the place.",
 "Eastern Gray Squirrel":
   "About a fifth of a squirrel curve is specific to the place \u2014 and it is the one species where "
   "that remainder is predictable from local conditions.",
 "Northern Raccoon":
   "About a third of a raccoon curve is specific to the place, but that remainder is local rather "
   "than regional, so a national map cannot deliver it.",
 "American Black Bear":
   "Roughly half of a bear curve is specific to the place \u2014 the largest local share of the five "
   "\u2014 but bear data are too thin to test whether a map can predict it.",
}


def curve_decomposition(cell_curves_long, species):
    """Split every predicted cell curve into a species-typical part and a cell-specific part.

    Each curve is normalised to sum to 1 over the 48 half-hour bins, so this is a decomposition of
    SHAPE (when the animal is active), not of detection volume. Sum of squares is taken relative to
    a FLAT 24-hour day (1/48 in every bin), which is the "no diel structure at all" reference:

        SS_total   = sum over cells and bins of (curve - flat)^2
        SS_species = n_cells * sum over bins of (species mean curve - flat)^2
        SS_cell    = sum over cells and bins of (curve - species mean curve)^2

    The cross term vanishes exactly because the species mean is the arithmetic mean of the curves,
    so SS_total = SS_species + SS_cell with no residual (asserted below). The reported percentages
    are those two parts over the total.
    """
    out = {}
    means = {}
    for sp in species:
        g = cell_curves_long[cell_curves_long.species == sp]
        piv = g.pivot_table(index="cell25", columns="bin", values="rate_permille",
                            observed=True).dropna()
        M = piv.values.astype(float)
        M = M / M.sum(axis=1, keepdims=True)
        nb = M.shape[1]
        flat = 1.0 / nb
        m = M.mean(axis=0)
        ss_tot = float(((M - flat) ** 2).sum())
        ss_sp = float(M.shape[0] * ((m - flat) ** 2).sum())
        ss_cell = float(((M - m) ** 2).sum())
        assert abs(ss_tot - ss_sp - ss_cell) < 1e-9 * max(ss_tot, 1.0), \
            f"{sp}: decomposition does not close ({ss_tot} vs {ss_sp}+{ss_cell})"
        out[sp] = {"n_cells": int(M.shape[0]),
                   "pct_species": 100.0 * ss_sp / ss_tot,
                   "pct_cell": 100.0 * ss_cell / ss_tot}
        means[sp] = m
    return out, means



# ---------------------------------------------------------------- 25 km covariate payload
TIER_ORDER = ["Interpolation", "Near", "Moderate extrapolation", "Severe extrapolation"]
MET_SHORT = {"pct_nocturnal": "noct", "crepuscular": "crep", "concentration": "conc",
             "peak_hour": "peak", "activity_level": "act"}
WQ = {"pct_nocturnal": 1, "crepuscular": 1, "concentration": 100, "peak_hour": 1,
      "activity_level": 100}
COEFQ = 200

NBIN = 48
BCEN = np.array([(i + 0.5) * 0.5 for i in range(NBIN)])
REF_SR, REF_SS = 6.2202, 17.9998


def _reconstruct(s1, c1, s2, c2):
    """Normalised 48-bin activity share from the four non-intercept harmonic coefficients.

    Identical arithmetic to the browser's covCurve(): every metric on this page is a property of
    curve SHAPE, so the intercept carries no information and is not shipped.
    """
    t = BCEN / 24.0 * 2 * np.pi
    v = (np.outer(np.atleast_1d(s1), np.sin(t)) + np.outer(np.atleast_1d(c1), np.cos(t))
         + np.outer(np.atleast_1d(s2), np.sin(2 * t)) + np.outer(np.atleast_1d(c2), np.cos(2 * t)))
    v = v - v.max(axis=1, keepdims=True)
    e = np.exp(v)
    return e / e.sum(axis=1, keepdims=True)


def _metrics_from_curves(Pm):
    """The five diel metrics from normalised 48-bin curves; matches the browser's covMetrics()."""
    isn = ((BCEN < REF_SR) | (BCEN >= REF_SS)).astype(float)
    isc = ((np.abs(BCEN - REF_SR) <= 1.5) | (np.abs(BCEN - REF_SS) <= 1.5)).astype(float)
    ang = BCEN / 24.0 * 2 * np.pi
    cc, ss = Pm @ np.cos(ang), Pm @ np.sin(ang)
    return {"pct_nocturnal": 100 * (Pm @ isn), "crepuscular": 100 * (Pm @ isc),
            "concentration": np.sqrt(cc ** 2 + ss ** 2),
            "peak_hour": BCEN[Pm.argmax(axis=1)],
            "activity_level": Pm.mean(axis=1) / Pm.max(axis=1)}


def enforce_range_quantised(q, obs_ranges, max_steps=60):
    """Re-impose the observed-range constraint IN QUANTISED SPACE, where the browser reads it.

    This exists because of a bug that shipped once. Upstream constrains each cell's coefficient
    vector to the joint region real arrays occupy, but the payload then rounds those coefficients to
    integers (x COEFQ). Rounding moves a vector by up to half a quantisation step, which is enough
    to push a cell that was exactly on the boundary back outside -- so a constraint satisfied in
    continuous space is NOT automatically satisfied in the space the page actually draws from.
    The fix is the same one upstream uses, applied after rounding: shrink the offending vector
    toward the species' mean vector by the SMALLEST amount that brings every metric back inside,
    and round again, iterating until the rounded result itself passes.

    Peak timing is tested with a one-bin tolerance on the 24-hour circle. The fitted curve lives on
    a 0.5 h grid, so a predicted peak within one bin of an observed peak (measured the short way
    around midnight) is inside the observed support, not outside it: 23.75 h against observed peaks
    at 23.25 h and 0.25 h sits in the gap between two occupied bins.

    Returns the corrected integer array and the number of cells that needed correcting.
    """
    q = q.astype(float).copy()
    m = q.mean(axis=0)

    def outside(qi):
        MM = _metrics_from_curves(_reconstruct(*(qi / COEFQ).T))
        bad = np.zeros(len(qi), bool)
        for mk, (lo, hi, obs) in obs_ranges.items():
            v = MM[mk]
            if mk == "peak_hour":
                d = np.abs(v[:, None] - obs[None, :])
                d = np.minimum(d, 24.0 - d)
                bad |= d.min(axis=1) > 0.51        # one 0.5 h bin, circularly
            else:
                bad |= (v < lo - 1e-9) | (v > hi + 1e-9)
        return bad

    qi = np.round(q)
    bad = outside(qi)
    n_fixed = int(bad.sum())
    lam = np.zeros(len(q))
    step = 0
    while bad.any() and step < max_steps:
        lam[bad] += 1.0 / max_steps
        qi = np.round(m + (1.0 - lam)[:, None] * (q - m))
        bad = outside(qi)
        step += 1
    assert not bad.any(), \
        f"{int(bad.sum())} cells still outside the observed range after {max_steps} shrink steps"
    return qi.astype(int), n_fixed


def build_cov_payload(coefs, gridcov, conf, preds, cv, tierval, species, metrics,
                      to_view, mets_keys, harmonics=None):
    """Assemble the 25 km 'local conditions' payload: geometry, coefficients, confidence, CV.

    Shipping FOUR harmonic coefficients per species x cell (quantised by COEFQ) rather than five
    pre-computed metric surfaces is what lets the full 25 km grid fit the size budget; the browser
    derives every metric and the click curve from them.
    """
    cells = sorted(set(coefs.cell25), key=lambda c: (int(c.split("_")[0][1:]), int(c.split("_")[1])))
    cidx = {c: i for i, c in enumerate(cells)}
    g = gridcov.set_index("cell25")
    missing = [c for c in cells if c not in g.index]
    assert not missing, f"{len(missing)} coefficient cells absent from the grid covariate table"
    gg = g.loc[cells]

    # ---- geometry: project each 25 km cell's four corners into the SAME viewbox transform as the
    # 100 km position map, then ship [x0, y0, dx1, dy1, dx2, dy2, dx3, dy3] at 0.1 px resolution.
    import pyproj
    tr = pyproj.Transformer.from_crs("EPSG:5070", "EPSG:4326", always_xy=True)
    X = gg.X_5070_km.values * 1000.0
    Y = gg.Y_5070_km.values * 1000.0
    half = 12500.0
    corner = []
    for dx, dy in ((-half, half), (half, half), (half, -half), (-half, -half)):
        lo, la = tr.transform(X + dx, Y + dy)
        corner.append(to_view(lo, la))
    # Round the ABSOLUTE corner positions, then take differences. Rounding each corner's OFFSET
    # instead (the earlier encoding) leaves a hairline gap along every shared edge: at this grid the
    # true cell width is 4.94 view units while the centre spacing is 4.95, so an independently
    # rounded width of 4.9 falls 0.05 units short of the neighbour's edge and the shortfall
    # accumulates into visible vertical seams every few columns. Two adjacent cells share a corner
    # position, so rounding the position itself makes both agree exactly and the tiling closes.
    R = [np.round(np.asarray(c) * 10).astype(np.int64) for c in corner]
    x0i, y0i = R[0][0], R[0][1]
    paths = np.column_stack([x0i, y0i] + [R[k][j] - (x0i if j == 0 else y0i)
                                          for k in (1, 2, 3) for j in (0, 1)]).tolist()

    # ---- coefficients, quantised. An all-zero vector marks a cell outside this species' range
    # mask; the browser uses that as an exact presence test, so assert no in-range cell is all-zero.
    CF, CONF, SCORE, CIW, SEAM = {}, {}, {}, {}, {}
    tier_counts, tier_pct, fixed = {}, {}, {}
    curves_mean = {}
    for sp in species:
        arr = np.zeros((len(cells), 4), int)
        sub = coefs[coefs.species == sp]
        ii = np.array([cidx[c] for c in sub.cell25])
        q = sub[["s1", "c1", "s2", "c2"]].values * COEFQ
        if harmonics is not None:
            OBSC = {"pct_nocturnal": "obs_noct", "crepuscular": "obs_crep",
                    "concentration": "obs_conc", "peak_hour": "obs_peak",
                    "activity_level": "obs_act"}
            ah = harmonics[harmonics.species == sp]
            ranges = {}
            for mk in mets_keys:
                o = ah[OBSC[mk]].dropna().values.astype(float)
                ranges[mk] = (float(o.min()), float(o.max()), o)
            q, n_fix = enforce_range_quantised(q, ranges)
            fixed[sp] = n_fix
        else:
            q = np.round(q).astype(int)
        zero = (q == 0).all(axis=1)
        assert not zero.any(), \
            f"{sp}: {int(zero.sum())} in-range cells quantise to an all-zero coefficient vector"
        arr[ii] = q
        CF[sp] = arr.tolist()
        Pm = _reconstruct(*(arr[ii] / COEFQ).T)
        curves_mean[sp] = Pm.mean(axis=0)

        cs = conf[conf.species == sp].set_index("cell25")
        cs = cs.loc[[c for c in sub.cell25]]
        tv = np.zeros(len(cells), int)
        sv = np.zeros(len(cells), int)
        tv[ii] = [TIER_ORDER.index(t) + 1 for t in cs.tier]
        sv[ii] = np.round(cs.score.values).astype(int)
        CONF[sp], SCORE[sp] = tv.tolist(), sv.tolist()
        tc = {t: int((tv[ii] == k + 1).sum()) for k, t in enumerate(TIER_ORDER)}
        tier_counts[sp] = tc
        tier_pct[sp] = {t: round(100.0 * v / max(len(ii), 1), 1) for t, v in tc.items()}
        SEAM[sp] = [int(bool(x)) for x in
                    gg.seam_interpolated.values[np.isin(np.arange(len(cells)), ii)]] \
            if "seam_interpolated" in gg else []

        CIW[sp] = {}
        for mk in mets_keys:
            p = preds[(preds.species == sp) & (preds.metric == MET_SHORT[mk])].set_index("cell25")
            w = np.zeros(len(cells), int)
            common = [c for c in sub.cell25 if c in p.index]
            jj = np.array([cidx[c] for c in common])
            pw = (p.loc[common].ci_hi.values - p.loc[common].ci_lo.values) * WQ[mk]
            w[jj] = np.round(np.nan_to_num(pw, nan=0.0)).astype(int)
            CIW[sp][mk] = w.tolist()

    # ---- cross-validation, read from the beats_fix / skill_fix columns ONLY
    CV = {sp: {} for sp in species}
    n_beat = 0
    for r in cv.itertuples():
        mk = next((k for k, v in MET_SHORT.items() if v == r.metric), None)
        if mk is None:
            continue
        CV[str(r.species)][mk] = {
            "skill": round(float(r.skill_fix), 4), "lo": round(float(r.fix_lo), 4),
            "hi": round(float(r.fix_hi), 4), "folds": f"{int(r.folds_pos_fix)}/{int(r.n_folds)}",
            "beats": bool(r.beats_fix), "inside": round(float(r.skill_fix_inside), 4),
            "n_te": int(r.n_arrays), "n_folds": int(r.n_folds)}
        n_beat += int(bool(r.beats_fix))
    n_tot = sum(len(v) for v in CV.values())

    # Non-finite interval bounds DO occur here (a band with too few folds to resample). An explicit
    # JSON null is fine -- JSON.parse accepts it and the table renders an em dash -- but a bare NaN
    # token would throw at parse time and blank the whole page, so convert rather than pass through.
    import math as _math
    def _fin(v):
        if isinstance(v, float):
            return round(v, 3) if _math.isfinite(v) else None
        return v
    tvv = [{k: _fin(v) for k, v in r.items()} for r in tierval.to_dict("records")]
    pool = [r for r in tvv if str(r["species"]).startswith("ALL SPECIES")]
    tierpool = [next((r["nmae_cov"] for r in pool if r["tier"] == t), None) for t in TIER_ORDER]

    return {"cells": cells, "paths": paths, "coefs": CF, "conf": CONF, "score": SCORE,
            "ciw": CIW, "cv": CV, "coefq": COEFQ, "wq": WQ, "seam": SEAM,
            "tiers": TIER_ORDER, "n_beat": n_beat, "n_tot": n_tot,
            "tier_counts": tier_counts, "tier_pct": tier_pct,
            "tierval": tvv, "tierpool": tierpool, "n_range_fixed": fixed}, curves_mean



# ---------------------------------------------------------------- support ranking for the map tab
def build_support(cv, species):
    """Per-species map-support tier and ordering, from beats_fix / skill_fix only."""
    sup = {}
    for sp in species:
        s = cv[cv.species == sp]
        beats = int(s.beats_fix.sum())
        n = int(len(s))
        mean_skill = float(s.skill_fix.mean())
        n_arrays = int(s.n_arrays.max())
        n_folds = int(s.n_folds.max())
        tier = "weak" if n_arrays < WEAK_MIN_ARRAYS else ("strong" if beats >= 4 else "partial")
        sup[sp] = {"beats": beats, "n": n, "mean_skill": round(mean_skill, 3),
                   "n_arrays": n_arrays, "n_folds": n_folds, "tier": tier,
                   "tier_label": SUPPORT_TIER_LABEL[tier], "glyph": SUPPORT_GLYPH[tier],
                   "note": SUPPORT_NOTE.get(sp, "")}
    rank = {"strong": 0, "partial": 1, "weak": 2}
    order = sorted(species, key=lambda s: (rank[sup[s]["tier"]], -sup[s]["beats"],
                                           -sup[s]["mean_skill"]))
    return sup, order


# ---------------------------------------------------------------- the results section (static HTML)
def results_html(dec, vp, cvsup, hyp, nulltest, sqx, msel, rangetab, png_b64):
    """Plain-language results carrying every number, with the scale figure embedded as base64.

    User-facing rule followed throughout: no 'null', 'baseline', 'MAE' or 'beta' in any label a
    reader sees. Where a statistic is unavoidable it is spelled out ("in standard deviations of the
    response, per standard deviation of the predictor").
    """
    SHORT = {"White-tailed Deer": "White-tailed deer", "Northern Raccoon": "Raccoons",
             "Eastern Gray Squirrel": "Gray squirrels", "Coyote": "Coyotes",
             "American Black Bear": "Black bears"}
    order = sorted(dec, key=lambda s: dec[s]["pct_cell"])
    rows = "".join(
        f"<tr><td>{s}</td><td>{dec[s]['pct_species']:.1f}%</td><td>{dec[s]['pct_cell']:.1f}%</td>"
        f"<td>{dec[s]['n_cells']:,}</td></tr>" for s in order)

    vpd = {r["species"]: r for r in vp.to_dict("records")}
    vprows = ""
    for s in ["White-tailed Deer", "Coyote", "Northern Raccoon", "Eastern Gray Squirrel"]:
        r = vpd.get(s)
        if not r or not np.isfinite(float(r.get("pct_regional_mappable", float("nan")))):
            continue
        vprows += (f"<tr><td>{s}</td><td>{float(r['pct_measurement_error']):.1f}%</td>"
                   f"<td>{float(r['pct_local_under25km']):.1f}%</td>"
                   f"<td>{float(r['pct_regional_mappable']):.1f}%</td>"
                   f"<td>{int(r['n_arrays'])}</td></tr>")

    # density / predator effects: bootstrap-significant, in standardized units
    def eff(hypo, sp, met, term):
        q = hyp[(hyp.hypothesis == hypo) & (hyp.species == sp) & (hyp.metric == met)
                & (hyp.term == term)]
        return None if not len(q) else q.iloc[0]

    dens = [("Eastern Gray Squirrel", "crep", "more of their activity into the hours around dawn "
             "and dusk"),
            ("White-tailed Deer", "crep", "more of their activity into the hours around dawn and "
             "dusk"),
            ("White-tailed Deer", "noct", "less of their activity into the night"),
            ("Coyote", "noct", "more of their activity into the night")]
    drows = ""
    for sp, met, phrase in dens:
        r = eff("H1_density", sp, met, "log_det")
        if r is None:
            continue
        drows += (f"<tr><td>{SHORT[sp]}</td><td style=\"text-align:left\">{phrase}</td>"
                  f"<td>{float(r.beta_std):+.2f}</td>"
                  f"<td>{float(r.beta_std_spatial):+.2f}</td></tr>")
    puma = eff("H2_carnivores", "White-tailed Deer", "noct", "cam_pres_Puma")
    coyn = eff("H1_density", "Coyote", "noct", "log_det")

    # null test: which measures survive it
    ntg = nulltest.groupby("metric").artifact_ratio.mean().to_dict()
    NTLAB = {"act": "activity level", "peak": "peak timing", "conc": "activity concentration",
             "noct": "% nocturnal", "crep": "activity around dawn and dusk"}
    ntrows = ""
    for m in ["noct", "crep", "conc", "peak", "act"]:
        keep = int(nulltest[nulltest.metric == m].interpretable.sum())
        ntrows += (f"<tr><td>{NTLAB[m]}</td><td>{100 * ntg[m]:.0f}%</td>"
                   f"<td>{'usable' if keep >= 4 else 'discarded'}</td></tr>")

    sqd = {r["metric"]: r for r in sqx.to_dict("records")}
    sqrows = ""
    for m in ["% nocturnal", "crepuscular", "concentration", "peak hour", "activity level"]:
        r = sqd.get(m)
        if r is None:
            continue
        same = str(r["differs"]).strip().lower() == "no"
        sqrows += (f"<tr><td>{m}</td><td>{float(r['east_mean']):.2f}</td>"
                   f"<td>{float(r['west_mean']):.2f}</td>"
                   f"<td>{float(r['diff_in_east_SDs']):+.2f}</td>"
                   f"<td>{float(r['perm_p']):.3f}</td>"
                   f"<td class=\"{'ok' if same else 'no'}\">"
                   f"{'indistinguishable' if same else 'differs'}</td></tr>")
    n_same = sum(1 for r in sqx.to_dict("records") if str(r["differs"]).strip().lower() == "no")

    mrows = ""
    for r in msel.to_dict("records"):
        mrows += (f"<tr><td>{r['metric']}</td><td>{float(r['reliability']):.2f}</td>"
                  f"<td>{int(r['beats_null_of_5'])} of 5</td>"
                  f"<td>{float(r['mean_cv_skill']):+.3f}</td>"
                  f"<td>{r['reconstruction_r']}</td></tr>")

    rtrows = ""
    for r in rangetab:
        flag = "" if r["n_outside"] == 0 else " class=\"no\""
        rtrows += (f"<tr><td>{r['species']}</td><td style=\"text-align:left\">{r['label']}</td>"
                   f"<td>{r['obs_lo']:.2f} to {r['obs_hi']:.2f}</td>"
                   f"<td>{r['pred_lo']:.2f} to {r['pred_hi']:.2f}</td>"
                   f"<td{flag}>{r['n_outside']}</td></tr>")

    return f"""
  <div class="card resultscard" id="results" style="margin-top:22px">
    <h2>What we found</h2>
    <p class="lead">Four results. The first two say what a national map can and cannot be; the
    third and fourth say what does move an animal's daily rhythm.</p>

    <h3 class="rh">1. A daily rhythm is mostly the species, not the place</h3>
    <p>Every curve this page can draw was split into two parts: the rhythm this species keeps
    everywhere in the country, and the part specific to the individual place. The split is on
    curve shape, measured against a completely flat day.</p>
    <div class="tw"><table class="res"><colgroup><col class="c0"><col class="cn"><col class="cn">
      <col class="cn"></colgroup>
      <tr><th>Animal</th><th>Same everywhere</th><th>Specific to the place</th>
      <th>Places compared</th></tr>{rows}</table></div>
    <p>Read the extremes. For white-tailed deer, <b>{dec['White-tailed Deer']['pct_species']:.0f}%
    of the curve at any place in the country is simply the nationwide deer rhythm</b> &mdash; click
    two cells a thousand kilometres apart and you are looking at almost the same picture. For black
    bears it is about half and half. This is why the curve panel now draws the nationwide rhythm
    behind the place you clicked and shades the gap between them: the gap is the only part that is
    about your click.</p>

    <h3 class="rh">2. The scale the variation lives at differs by species &mdash; and neither
    pattern gives a predictive map</h3>
    <p>The place-specific remainder is real. Splitting the differences between camera arrays into
    measurement noise, variation among nearby arrays (under 25&nbsp;km) and variation between
    regions shows two completely different regimes.</p>
    <div class="tw"><table class="res"><colgroup><col class="c0"><col class="cn"><col class="cn">
      <col class="cn"><col class="cn"></colgroup>
      <tr><th>Animal</th><th>Measurement noise</th><th>Differs between neighbours (under
      25&nbsp;km)</th><th>Differs between regions</th><th>Camera arrays</th></tr>
      {vprows}</table></div>
    <p><b>Deer and coyotes are regional animals</b> &mdash; 59% and 42% of their between-array
    differences are structure a map could in principle draw. <b>Raccoons and gray squirrels are
    not</b>: 0.0% and 0.5% regional, essentially all of their variation sits between neighbouring
    arrays.</p>
    <p><b>And yet the ranking of map quality is upside down.</b> Gray squirrel, with no regional
    pattern at all, is the best-predicted species (beats a single nationwide curve on
    {cvsup['Eastern Gray Squirrel']['beats']} of 5 measures, average margin
    {cvsup['Eastern Gray Squirrel']['mean_skill']:+.3f}). White-tailed deer, with the most regional
    pattern and the largest sample, is the worst ({cvsup['White-tailed Deer']['beats']} of 5,
    average margin {cvsup['White-tailed Deer']['mean_skill']:+.3f} &mdash; on average worse than
    assuming deer keep the same hours everywhere). The two regimes fail for opposite reasons: deer
    have regional structure driven by something we never measured (state-level hunting rules, deer
    management and predator communities are all plausible and none is in this model), while
    squirrels have no regional pattern but do respond to conditions we did measure, which is
    exactly what transfers to unsampled country. <b>It would be wrong to summarise this as "the
    variation is local"</b> &mdash; that is true of raccoons and squirrels and false of deer and
    coyotes.</p>
    <figure class="scalefig">
      <img alt="Three panels. Left: differences between camera arrays grow with distance for deer
      and coyote but not for raccoon and squirrel. Middle: share of between-array variance that is
      regional, local or measurement noise, by species. Right: regional structure plotted against
      out-of-sample skill, showing gray squirrel high-skill with no regional structure and
      white-tailed deer below zero with the most regional structure."
      src="data:image/png;base64,{png_b64}">
      <figcaption>How alike two camera arrays are as they get further apart (left); where each
      species' variation sits (middle); and why more regional structure did not buy better
      prediction (right).</figcaption>
    </figure>

    <h3 class="rh">3. More animals and more predators do shift the clock &mdash; but do not help
    predict it</h3>
    <p>Effects below survived three hurdles in sequence: a scrambled-data check that the effect is
    not an artefact of how the measure is computed, a resampling interval that excludes no effect,
    and refitting with a spatial term so a regional gradient cannot masquerade as a local one.
    Sizes are given the same way for every row: <b>how many standard deviations the rhythm moves
    per one standard deviation of the predictor</b>.</p>
    <div class="tw"><table class="res"><colgroup><col class="c0"><col class="c1"><col class="cn">
      <col class="cn"></colgroup>
      <tr><th>Animal</th><th>Where the animal is more often photographed, it shifts&hellip;</th>
      <th>Size</th><th>With a spatial term added</th></tr>{drows}</table></div>
    <p>The largest single effect is a predator one: where a puma was photographed by the same
    cameras, deer put <b>{abs(float(puma.beta_std)):.2f} standard deviations less</b> of their
    activity into the night ({float(puma.beta):.1f} percentage points,
    {float(puma.boot_lo):.1f} to {float(puma.boot_hi):.1f}). <b>This is not independent of the
    density result above:</b> deer detection rate and puma presence push % nocturnal the same way
    and the two are correlated across sites, so they should be read as one entangled signal, not
    two confirmations. Coyotes shift {float(coyn.beta_std):+.2f} standard deviations more nocturnal
    where they are more often photographed &mdash; a real but small effect.</p>
    <p>Three of the five measures had to be discarded for this analysis. Rebuilding each measure
    from scrambled detection times reproduces a large share of the observed effect &mdash; the
    measure responds to how many photographs there are, not only to when they were taken:</p>
    <div class="tw"><table class="res"><colgroup><col class="c0"><col class="cn"><col class="cn">
      </colgroup><tr><th>Measure</th><th>Share of the effect that scrambled data reproduces</th>
      <th>Verdict</th></tr>{ntrows}</table></div>
    <p><b>And none of it helps predict.</b> Adding detection rate or any carnivore combination to
    the model changes out-of-sample accuracy by a median of
    {hyp[hyp.hypothesis == 'H2_carnivores'].cv_skill_delta.median():+.3f} across
    {int((hyp.hypothesis == 'H2_carnivores').sum())} carnivore combinations tested. These are
    explanations of what the rhythm responds to, not ingredients that make the map work.</p>

    <h3 class="rh">4. A natural experiment: squirrels moved 3,000&nbsp;km keep their hours</h3>
    <p>Western gray squirrels were introduced to the Pacific Northwest, a different climate,
    different forest and a different predator community from the eastern range. If a daily rhythm
    were mostly a response to local conditions, they should have changed. On
    {n_same} of 5 measures they are statistically indistinguishable from
    {int(sqd['% nocturnal']['n_east'])} native eastern arrays, including the headline one: their
    {float(sqd['% nocturnal']['west_mean']):.1f}% nocturnal sits at the
    {float(sqd['% nocturnal']['west_mean_pctile_of_east']):.0f}th percentile of the native
    distribution.</p>
    <div class="tw"><table class="res"><colgroup><col class="c0"><col class="cn"><col class="cn">
      <col class="cn"><col class="cn"><col class="cn"></colgroup>
      <tr><th>Measure</th><th>Native east</th><th>Introduced west</th>
      <th>Difference (in eastern spreads)</th><th>Chance of this by luck</th><th>Verdict</th></tr>
      {sqrows}</table></div>
    <p><b>The caveat is load-bearing.</b> All {int(sqd['% nocturnal']['n_west'])} western arrays are
    suburban, so part of what looks like constancy could be a suburban setting reproducing suburban
    behaviour. Matched against population-comparable eastern arrays the % nocturnal difference is
    {float(sqd['% nocturnal']['diff_vs_population_matched_east']):+.2f} percentage points, similar
    to the unmatched figure. Two of five measures did change: activity around dawn and dusk fell
    {abs(float(sqd['crepuscular']['diff_in_east_SDs'])):.2f} eastern spreads
    (chance {float(sqd['crepuscular']['perm_p']):.3f}) and the rhythm is more tightly packed
    ({float(sqd['concentration']['diff_in_east_SDs']):+.2f},
    chance {float(sqd['concentration']['perm_p']):.3f}). Six arrays is a small experiment.</p>

    <details style="margin-top:14px"><summary>How the five measures compare, and the range check
    on every displayed value</summary>
      <p class="note" style="margin-top:8px">Repeatability is how much of the difference between
      arrays is real rather than measurement noise. "Beats one nationwide curve" counts species.
      Fidelity is how well a measure survives being rebuilt from the smoothed curve.</p>
      <div class="tw"><table class="res"><colgroup><col class="c0"><col class="cn"><col class="cn">
        <col class="cn"><col class="cn"></colgroup>
        <tr><th>Measure</th><th>Repeatability</th><th>Beats one nationwide curve</th>
        <th>Average margin</th><th>Fidelity</th></tr>{mrows}</table></div>
      <p class="note" style="margin-top:12px">Every value drawn on the local-conditions map is
      constrained to the range real camera arrays actually produced for that animal and measure.
      The table below is regenerated at build time from the shipped payload, so it describes what
      the page draws rather than what an upstream file claimed.</p>
      <div class="tw"><table class="res"><colgroup><col class="c0"><col class="c1"><col class="cn">
        <col class="cn"><col class="cn"></colgroup>
        <tr><th>Animal</th><th>Measure</th><th>Real arrays span</th><th>Map draws</th>
        <th>Outside</th></tr>{rtrows}</table></div>
    </details>
  </div>
"""



def range_check(covpay, harmonics, species, mets):
    """Is every value the map draws inside the range real camera arrays produced?

    Recomputed from the SHIPPED, QUANTISED payload rather than from an upstream file, because
    quantisation is applied after the constraint and could in principle push a value out.

    Peak timing is tested CIRCULARLY. It is an hour on a 24-hour clock, so a linear min/max test
    calls 23.75 h "outside" a range whose maximum is 23.25 h even though the two are 30 minutes
    apart across midnight. Both tests are reported: the linear one to stay comparable with the
    upstream table, the circular one because it is the correct test for an angle.
    """
    OBS = {"pct_nocturnal": "obs_noct", "crepuscular": "obs_crep", "concentration": "obs_conc",
           "peak_hour": "obs_peak", "activity_level": "obs_act"}
    rows = []
    for sp in species:
        B = np.array(covpay["coefs"][sp], float) / covpay["coefq"]
        keep = ~(B == 0).all(axis=1)
        MM = _metrics_from_curves(_reconstruct(*B[keep].T))
        a = harmonics[harmonics.species == sp]
        for m in mets:
            mk = m["key"]
            obs = a[OBS[mk]].dropna().values.astype(float)
            v = MM[mk]
            lo, hi = float(obs.min()), float(obs.max())
            outside = (v < lo - 1e-9) | (v > hi + 1e-9)
            n_lin = int(outside.sum())
            if mk == "peak_hour" and n_lin:
                # ONLY the values the linear test flagged get re-tested. On a 24-hour circle a
                # predicted peak of 23.75 h against an observed maximum of 23.25 h is 30 minutes
                # away across midnight, not outside the range, so the circular distance to the
                # nearest observed peak is the correct test for exactly those cases. Values the
                # linear test did not flag already sit inside the observed span.
                w = v[outside]
                d = np.abs(w[:, None] - obs[None, :])
                d = np.minimum(d, 24.0 - d)                  # circular distance, hours
                # Tolerance is ONE bin (0.5 h) plus rounding slack. Predicted and observed peaks both
                # live on the same 0.5 h bin-centre grid, so a predicted peak 0.5 h from an observed
                # one is at the adjacent grid point -- interpolated between two observed values, not
                # extrapolated past them. Bear's 23.75 h cells sit exactly between observed peaks at
                # 23.25 h and 0.25 h. Same tolerance as the enforcement step, so the check and the
                # constraint cannot disagree.
                n_circ = int((d.min(axis=1) > 0.51).sum())
            else:
                n_circ = n_lin
            rows.append({"species": sp, "metric": mk, "label": m["label"],
                         "n_cells": int(len(v)), "obs_lo": lo, "obs_hi": hi,
                         "pred_lo": float(v.min()), "pred_hi": float(v.max()),
                         "n_outside": n_circ, "n_outside_linear": n_lin})
    return rows


def build(surfaces, curves, panels, fit_model, fit_validation, fit_ess,
          grid_geojson, contract, methods_href="methods.html", placeholder=False,
          l48_cells="viewer_cells_l48.csv", states_json="viewer_states.json",
          cell_curves=None, act_summary=None,
          coefs=None, gridcov=None, conf=None, preds=None, cv=None, tierval=None,
          cellcurves_long=None, harmonics=None, varpart=None, hyp=None, nulltest=None,
          sqx=None, msel=None, scale_png=None, out=OUT):
    grid = json.load(open(grid_geojson))

    # ---- restrict to grid cells overlapping lower-48 LAND
    keep_ids = None
    if l48_cells and os.path.exists(l48_cells):
        keep_ids = set(pd.read_csv(l48_cells)["cell_id"])
        grid["features"] = [f for f in grid["features"]
                            if f["properties"]["cell_id"] in keep_ids]
        surfaces = surfaces[surfaces.cell_id.isin(keep_ids)].copy()

    # ---- project every grid cell once, to a shared SVG viewbox
    polys = {}
    for f in grid["features"]:
        ring = np.array(f["geometry"]["coordinates"][0], float)
        x, y = albers(ring[:, 0], ring[:, 1])
        polys[f["properties"]["cell_id"]] = np.c_[x, y]
    allpts = np.vstack(list(polys.values()))
    x0, x1 = allpts[:, 0].min(), allpts[:, 0].max()
    y0, y1 = allpts[:, 1].min(), allpts[:, 1].max()
    W, H = 980.0, 620.0
    pad = 12.0
    s = min((W - 2 * pad) / (x1 - x0), (H - 2 * pad) / (y1 - y0))
    ox = pad + ((W - 2 * pad) - (x1 - x0) * s) / 2.0
    oy = pad + ((H - 2 * pad) - (y1 - y0) * s) / 2.0

    def to_view(lon, lat):
        ax_, ay_ = albers(lon, lat)
        return ox + (ax_ - x0) * s, oy + (y1 - ay_) * s

    def path_of(cid):
        p = polys[cid]
        sx = ox + (p[:, 0] - x0) * s
        sy = oy + (y1 - p[:, 1]) * s          # flip y for screen coords
        return "M" + " ".join(f"{a:.1f},{b:.1f}" for a, b in zip(sx, sy)) + "Z"

    cell_paths = {cid: path_of(cid) for cid in polys}

    species_all = [sp["key"] for sp in contract["species"]]
    sci = {sp["key"]: sp["sci"] for sp in contract["species"]}
    metrics = contract["metrics"]
    mets_keys = [m["key"] for m in metrics]

    # ---- MAP SUPPORT: ordering and tier come from beats_fix / skill_fix
    SUP, sp_order = build_support(cv, species_all)
    species = sp_order              # species tabs are emitted in support order, best first

    # ---- PART 1: the decomposition, recomputed from covariate_cellcurves_fixed.csv
    DEC, cc_mean_25 = curve_decomposition(cellcurves_long, species_all)

    # ---- pack the 100 km position surfaces
    S = {}
    for sp in species:
        S[sp] = {}
        for m in metrics:
            sub = surfaces[(surfaces.species == sp) & (surfaces.metric == m["key"])]
            S[sp][m["key"]] = {
                r.cell_id: [round(float(r.value), 3), round(float(r.ci_width), 3),
                            round(float(r.ci_lo), 3), round(float(r.ci_hi), 3)]
                for r in sub.itertuples()}

    dom = {}
    for m in metrics:
        v = surfaces[surfaces.metric == m["key"]]["value"].astype(float)
        w = surfaces[surfaces.metric == m["key"]]["ci_width"].astype(float)
        dom[m["key"]] = {"vlo": float(np.nanpercentile(v, 2)), "vhi": float(np.nanpercentile(v, 98)),
                         "wlo": float(np.nanpercentile(w, 2)), "whi": float(np.nanpercentile(w, 98))}
    domsp = {}
    for sp in species:
        domsp[sp] = {}
        for m in metrics:
            sub = surfaces[(surfaces.species == sp) & (surfaces.metric == m["key"])]
            v = sub["value"].astype(float); w = sub["ci_width"].astype(float)
            domsp[sp][m["key"]] = {
                "vlo": float(np.nanpercentile(v, 2)), "vhi": float(np.nanpercentile(v, 98)),
                "wlo": float(np.nanpercentile(w, 2)), "whi": float(np.nanpercentile(w, 98))}

    # ---- covariate panels, unchanged filtering rule
    P = {sp: {} for sp in species}
    PCOUNT = {sp: {} for sp in species}
    for sp in species:
        for m in metrics:
            sub = panels[(panels.species == sp) & (panels.metric == m["key"])]
            arr, n_total, hidden = [], 0, []
            for pid, g in sub.groupby("panel_id", sort=False):
                g = g.sort_values("x")
                n_total += 1
                kind = str(g.panel_kind.iloc[0]) if "panel_kind" in g.columns else "curve"
                if not panel_significant(g):
                    hidden.append(str(g.title.iloc[0]))
                    continue
                take = panel_takehome(sp, m["key"], pid, g)
                arr.append({"id": pid, "title": str(g.title.iloc[0]), "kind": kind,
                            "xlabel": str(g.xlabel.iloc[0]), "ylabel": str(g.ylabel.iloc[0]),
                            "take": take, "caption": take + "  " + str(g.caption.iloc[0]),
                            "x": [round(float(v), 4) for v in g.x],
                            "y": [round(float(v), 4) for v in g.y],
                            "lo": [round(float(v), 4) for v in g.lo],
                            "hi": [round(float(v), 4) for v in g.hi]})
            P[sp][m["key"]] = arr
            PCOUNT[sp][m["key"]] = {"kept": len(arr), "total": n_total, "hidden": hidden}

    FM = {sp: fit_model[fit_model.species == sp].to_dict("records") for sp in species}
    FV = {sp: fit_validation[fit_validation.species == sp].to_dict("records") for sp in species}
    FE = {sp: fit_ess[fit_ess.species == sp].to_dict("records") for sp in species}

    state_paths = []
    if states_json and os.path.exists(states_json):
        sj = json.load(open(states_json))
        if abs(sj.get("vb", [0, 0])[0] - W) < 1 and abs(sj.get("vb", [0, 0])[1] - H) < 1:
            state_paths = sj["state_paths"]

    # ---- per-cell click curves on the 100 km position map, plus the SPECIES MEAN in the same
    # units, so the deviation display has a reference curve that is exactly comparable.
    CC, CCMEAN = {}, {}
    if cell_curves is not None and len(cell_curves):
        sub = cell_curves[cell_curves.cell_id.isin(keep_ids)] if keep_ids is not None else cell_curves
        for (sp, cid), g in sub.groupby(["species", "cell_id"], sort=False):
            g = g.sort_values("bin")
            CC.setdefault(sp, {})[cid] = {
                "r": [int(round(float(v) * 1000)) for v in g.rate],
                "w": int(round(float(g.ci_rel.iloc[0]) * 1000)) if "ci_rel" in g else 180}
        for sp, d in CC.items():
            M = np.array([d[c]["r"] for c in d], float)
            CCMEAN[sp] = [int(round(v)) for v in M.mean(axis=0)]

    HCI = {}
    for sp in species:
        r = fit_ess[fit_ess.species == sp]
        HCI[sp] = float(r.honest_ci_width_pp.iloc[0]) if len(r) and "honest_ci_width_pp" in r else None

    # ---- 25 km covariate payload; COVMEAN is the species mean curve in NORMALISED share units,
    # the reference the deviation view draws behind a clicked cell.
    COV, covmean = build_cov_payload(coefs, gridcov, conf, preds, cv, tierval,
                                     species_all, metrics, to_view, mets_keys,
                                     harmonics=harmonics)
    COVMEAN = {sp: [round(float(v), 6) for v in covmean[sp]] for sp in covmean}

    # cross-check the two independent means of the same 25 km curves
    for sp in species_all:
        d = float(np.abs(np.array(covmean[sp]) - cc_mean_25[sp]).max())
        # Boundary cells are shrunk toward the species mean when quantisation pushes them outside
        # the observed range, so the two means are close but not identical by construction.
        assert d < 5e-3, f"{sp}: payload mean curve disagrees with source curves by {d:.4f}"

    import math as _m
    def _allfin(o, where):
        if o is None:
            return
        if isinstance(o, float):
            assert _m.isfinite(o), f"non-finite value in covariate payload at {where}"
        elif isinstance(o, dict):
            for k, v in o.items(): _allfin(v, f"{where}.{k}")
        elif isinstance(o, (list, tuple)):
            for i, v in enumerate(o): _allfin(v, f"{where}[{i}]")
    _allfin(COV, "cov")

    RANGE = range_check(COV, harmonics, species_all, metrics)

    DECPAY = {sp: {"pct_species": round(DEC[sp]["pct_species"], 1),
                   "pct_cell": round(DEC[sp]["pct_cell"], 1),
                   "n_cells": DEC[sp]["n_cells"], "note": DEV_NOTE.get(sp, "")}
              for sp in DEC}

    payload = {"species": species, "sci": sci, "metrics": metrics, "dom": dom,
               "domsp": domsp, "cov": COV, "paths": cell_paths, "surfaces": S,
               "cell_curves": CC, "cc_mean": CCMEAN, "cov_mean": COVMEAN,
               "panels": P, "pcount": PCOUNT, "states": state_paths,
               "fit_model": FM, "fit_validation": FV, "fit_ess": FE,
               "honest_ci": HCI, "act_summary": act_summary or {},
               "dec": DECPAY, "support": SUP, "sp_order": species,
               "tier_label": SUPPORT_TIER_LABEL,
               "vb": [W, H], "placeholder": bool(placeholder),
               "glossary": GLOSSARY, "cv_help": CV_HELP,
               "ramp_value": [RAMP_VALUE(i / 40) for i in range(41)],
               "ramp_ci": [RAMP_CI(i / 40) for i in range(41)],
               "ramp_cyc": [RAMP_CYC(i / 40) for i in range(41)],
               "cyclic": sorted(CYCLIC_METRICS)}

    png_b64 = ""
    if scale_png and os.path.exists(scale_png):
        import base64
        png_b64 = base64.b64encode(open(scale_png, "rb").read()).decode("ascii")

    RESULTS = results_html(DEC, varpart, SUP, hyp, nulltest, sqx, msel, RANGE, png_b64)

    banner = ("""<div id="ph">PLACEHOLDER DATA — synthetic surfaces for layout review only.
      No result on this page is a scientific finding until rebuilt with the fitted models.</div>"""
              if placeholder else "")


    html = ("""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Mammal diel activity across the lower 48</title>
<style>
:root{--ink:#1a1a1a;--mut:#5b6470;--line:#e2e2e2;--bg:#fff;--accent:#2b6cb0;--warn:#8a4a08}
*{box-sizing:border-box}
body{margin:0;font:15px/1.55 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,Helvetica,Arial,sans-serif;
color:var(--ink);background:var(--bg)}
header{padding:22px 26px 14px;border-bottom:1px solid var(--line)}
h1{margin:0 0 4px;font-size:23px;font-weight:650;letter-spacing:-.2px}
.sub{color:var(--mut);font-size:13.5px}
#ph{background:#fffbeb;border:1px solid #b45309;color:var(--warn);padding:9px 14px;margin:12px 26px 0;
border-radius:6px;font-size:13px;font-weight:600}
main{padding:18px 26px 60px;max-width:1420px}
.tabs{display:flex;flex-wrap:wrap;gap:6px;margin:14px 0 6px}
.tab{padding:7px 13px;border:1px solid #b9c2cc;border-radius:999px;background:#fafafa;cursor:pointer;
font-size:13.5px;white-space:nowrap;color:#28313b;text-align:left}
.tab:hover{background:#f0f4f8;border-color:#8595a5}
.tab[aria-selected=true]{background:var(--accent);color:#fff;border-color:#1d4e80;font-weight:600}
.tab .glyph{font-weight:800;margin-right:5px}
.tab .tsub{display:block;font-size:10.5px;letter-spacing:.02em;font-weight:600;opacity:.95}
/* unsupported measure tabs are struck through AND labelled, never colour alone */
.tab.unsup .tlab{text-decoration:line-through;text-decoration-thickness:1.5px}
.tablab{font-size:11.5px;text-transform:uppercase;letter-spacing:.08em;color:var(--mut);
margin:12px 0 2px;font-weight:700}
.tabhint{font-size:12px;color:var(--mut);margin:0 0 4px}
.layout{display:grid;grid-template-columns:392px minmax(0,1fr);gap:22px;align-items:start;margin-top:18px}
@media(max-width:1120px){.layout{grid-template-columns:1fr}}
.card{border:1px solid var(--line);border-radius:9px;padding:14px 16px;background:#fff}
.card h2{margin:0 0 3px;font-size:16px;font-weight:640}
.card .note{color:var(--mut);font-size:12.5px;margin:0 0 10px}
svg{display:block;max-width:100%;height:auto}
.cell{stroke:#fff;stroke-width:.35;cursor:pointer}
.cell:hover{stroke:#111;stroke-width:1.4}
.cell.sel{stroke:#111;stroke-width:2.1}
.cell.flash{stroke:#b91c1c;stroke-width:3.2}
.statebd{fill:none;stroke:#3a3a3a;stroke-width:.5;opacity:.34;pointer-events:none;stroke-linejoin:round}
.hdr{cursor:help;border-bottom:1px dotted #767f8a}
.hdr sup{color:var(--accent);font-weight:700;margin-left:1px}
details{border-bottom:1px solid var(--line);padding:7px 0}
details summary{cursor:pointer;font-size:13.5px}
.legend{display:flex;align-items:center;gap:9px;margin-top:9px;font-size:12px;color:var(--mut);flex-wrap:wrap}
.swatches{display:flex;height:11px;width:210px;border:1px solid #949ca6;border-radius:2px;overflow:hidden}
.swatches i{flex:1}
.grid2{display:grid;grid-template-columns:repeat(auto-fill,minmax(215px,1fr));gap:12px;margin-top:8px}
.mini{border:1px solid var(--line);border-radius:7px;padding:7px;cursor:zoom-in;background:#fff}
.mini:hover{border-color:var(--accent);box-shadow:0 1px 7px rgba(43,108,176,.14)}
.mini b{display:block;font-size:11.5px;font-weight:620;margin-bottom:3px}
.mini .take{display:block;font-size:11.5px;color:#374151;line-height:1.4;margin-top:4px}
.nopanel{border:1px dashed #a9adb4;border-radius:8px;padding:12px 14px;background:#fafafa;
color:#4b5057;font-size:13px}
/* TABLE OVERFLOW: .tw scrolls on its own, table-layout is fixed and every column has an explicit
   width, so the deliberately wordy headings can never widen a card. */
.tw{overflow-x:auto;-webkit-overflow-scrolling:touch;max-width:100%}
table{border-collapse:collapse;width:100%;font-size:12.5px;margin-top:4px;table-layout:fixed}
th,td{text-align:right;padding:5px 7px;border-bottom:1px solid var(--line);
font-variant-numeric:tabular-nums;white-space:normal;overflow-wrap:break-word;word-break:break-word}
th:first-child,td:first-child{text-align:left}
th{font-weight:640;color:#4a5360;font-size:11.5px;text-transform:uppercase;letter-spacing:.04em;
vertical-align:bottom}
table.fit{min-width:430px}
table.fit col.c0{width:34%}table.fit col.cn{width:16.5%}
table.cv{min-width:640px}
table.cv col.c0{width:20%}table.cv col.c1{width:20%}table.cv col.cn{width:15%}
table.res{min-width:560px}
table.res col.c0{width:23%}table.res col.c1{width:31%}table.res col.cn{width:15.3%}
table.tier{min-width:560px}
table.tier col.c0{width:28%}table.tier col.cn{width:18%}
.kv{display:flex;justify-content:space-between;gap:10px;padding:4px 0;
border-bottom:1px dotted var(--line);font-size:13px;align-items:baseline}
.kv span:first-child{flex:1 1 auto;min-width:0;overflow-wrap:break-word;word-break:break-word}
.kv span:last-child{flex:0 0 auto;white-space:nowrap;font-variant-numeric:tabular-nums;font-weight:600}
.ok{color:#14532d;font-weight:640}.no{color:#8f1414;font-weight:640}
#tip{position:fixed;pointer-events:none;background:#111;color:#fff;padding:6px 9px;border-radius:5px;
font-size:12px;opacity:0;transition:opacity .1s;z-index:60;max-width:250px}
#modal{position:fixed;inset:0;background:rgba(0,0,0,.62);display:none;align-items:center;
justify-content:center;z-index:80;padding:26px}
#modal.on{display:flex}
#modalbox{background:#fff;border-radius:11px;padding:20px;max-width:860px;width:100%}
#modalbox h3{margin:0 0 3px;font-size:17px}
#modalcap{color:#4a5360;font-size:13px;margin:0 0 12px}
.x{float:right;cursor:pointer;color:#4a5360;font-size:20px;line-height:1;border:0;background:0}
a{color:#1c5c9e}
.foot{margin-top:26px;padding-top:14px;border-top:1px solid var(--line);color:var(--mut);font-size:12.5px}
.seg{display:inline-flex;border:1px solid #b9c2cc;border-radius:7px;overflow:hidden;margin-top:6px}
.seg button{border:0;background:#fafafa;padding:6px 12px;cursor:pointer;font-size:13px;color:#28313b}
.seg button[aria-pressed=true]{background:var(--accent);color:#fff;font-weight:600}
.seg button[disabled]{color:#5b6470;cursor:not-allowed}
.mapwrap{position:relative}
.mapfr{border:1px solid var(--line);border-radius:8px;overflow:hidden;background:#fff}
.mapcav{position:absolute;left:0;right:0;bottom:0;background:rgba(255,251,235,.97);
border-top:2px solid #b45309;color:#6d2a0f;font-size:12.5px;font-weight:600;padding:8px 11px;
line-height:1.4;border-radius:0 0 7px 7px}
.mapcav .ciw{display:block;font-weight:700;margin-top:3px;color:#7c3708}
/* ---- SUPPORT MARKING. A surface that does not beat one nationwide curve must not look like one
   that does. Three channels, none of them colour alone: a stamp across the map, a hatch over the
   whole surface, and a worded badge. */
.nsstamp{position:absolute;top:12px;left:50%;transform:translateX(-50%);z-index:4;
background:#fff;border:2px solid #8f1414;color:#6d1010;border-radius:7px;padding:7px 13px;
font-size:13px;font-weight:800;text-align:center;line-height:1.35;max-width:86%;
box-shadow:0 1px 5px rgba(0,0,0,.14)}
.nsstamp .sub2{display:block;font-weight:600;font-size:11.5px;color:#3f2323;margin-top:2px}
.badge{display:inline-block;padding:5px 11px;border-radius:999px;font-size:12.5px;font-weight:700;
margin:2px 0 8px;line-height:1.35;background:#e8eef4;color:#1f2937;border:1px solid #7d8b99}
.badge .bglyph{font-weight:800;margin-right:5px}
.badge.supported{background:#e6f2ef;color:#0c3229;border:1px solid #3f7365}
.badge.unsupported{background:#f6ecea;color:#5c1f14;border:1px solid #9c5b4e}
.badge.weak{background:#f3eee2;color:#4a3708;border:1px solid #8f7a45}
.vnums{font-size:12.5px;color:#333b45;margin:0 0 7px;line-height:1.55}
.vnums b{font-variant-numeric:tabular-nums}
/* ---- headline */
.headline{border:1px solid #b6c1cd;border-left:5px solid #8a4a08;border-radius:9px;
background:#fffdf7;padding:15px 18px;margin:16px 0 4px}
.headline h2{margin:0 0 7px;font-size:19px;font-weight:680;letter-spacing:-.2px;color:#6d2a0f}
.headline p{margin:0 0 8px;font-size:14px;line-height:1.6}
.headline p:last-child{margin-bottom:0}
.hl-why{color:#2c3540;border-top:1px dotted #cfcbc7;padding-top:8px}
/* ---- deviation curve panel: the departure is the object, not the absolute shape */
.selcurve{border:2px solid var(--accent);border-radius:9px;background:#f7fbff;padding:10px 11px;
margin-bottom:12px}
.selcurve .lab{font-size:13px;font-weight:700;color:#17457f;margin-bottom:1px}
.selcurve .hint{font-size:11.5px;color:#3c4653;margin-bottom:4px}
.devkey{display:flex;flex-wrap:wrap;gap:6px 14px;font-size:11.5px;color:#333b45;margin:4px 0 2px}
.devkey span{display:inline-flex;align-items:center;gap:5px;white-space:nowrap}
.devkey i{width:17px;height:0;border-top-width:3px;border-top-style:solid;display:inline-block}
.devkey i.band{height:9px;border:1px solid #7d8b99;border-top:0;background:#c9b8e0}
.decbox{border:1px solid #b6c1cd;border-radius:8px;background:#fbfaf7;padding:9px 11px;
margin:0 0 11px;font-size:12.5px;line-height:1.5;color:#2c3540}
.decbox .big{font-size:15px;font-weight:750;color:#4a2c06}
.decbar{display:flex;height:14px;border:1px solid #7d8b99;border-radius:3px;overflow:hidden;
margin:6px 0 5px}
.decbar i{display:block}
.decbar i.sp{background:#b8c4d2}
.decbar i.cl{background:#7b52ab}
.cmphead{font-size:11.5px;text-transform:uppercase;letter-spacing:.07em;color:var(--mut);
font-weight:700;margin:14px 0 2px}
.cmp{border:1px solid var(--line);border-radius:7px;padding:7px 8px;margin-bottom:8px;cursor:pointer;
background:#fff}
.cmp:hover{border-color:var(--accent);background:#f7fbff}
.cmp .lab{font-size:12px;font-weight:620}
.cmp .go{font-size:11px;color:#1c5c9e}
/* ---- support ranking table on the map tab */
.suprow{display:grid;grid-template-columns:1.4fr 2.6fr;gap:8px 14px;align-items:start;
border-top:1px solid var(--line);padding:9px 0;font-size:12.5px}
.suprow:first-of-type{border-top:0}
.suprow .who{font-weight:700}
.suprow .who .glyph{font-weight:800;margin-right:5px}
.suprow .who .num{display:block;font-weight:600;color:#4a5360;font-size:11.5px;
font-variant-numeric:tabular-nums}
.suprow.here{background:#f2f7fc;border-radius:6px;padding-left:7px;padding-right:7px}
/* ---- confidence key for the local-conditions map */
.conftitle{font-size:12px;font-weight:650;color:#2c3540;margin:9px 0 4px}
.confrow{display:flex;flex-wrap:wrap;gap:8px 15px;align-items:center}
.confk{display:inline-flex;align-items:center;gap:5px;font-size:11.5px;color:#2c3540;white-space:nowrap}
.confk i{width:15px;height:11px;border-radius:2px;display:inline-block;border:1px solid #949ca6}
.confk b{font-variant-numeric:tabular-nums;color:#0f172a}
.covcell{stroke:none}
.covcell:hover{stroke:#111827;stroke-width:.8}
/* ---- results section */
.resultscard{border-color:#b7d7ea;background:#fcfeff}
.resultscard h2{font-size:19px}
.resultscard .lead{font-size:15.5px;line-height:1.6;margin:2px 0 10px;color:#0f172a}
.resultscard p{font-size:14px;line-height:1.62;margin:0 0 10px}
.resultscard .rh{font-size:15.5px;font-weight:700;margin:20px 0 6px;color:#0c3f5e;
border-top:1px solid #b7d7ea;padding-top:12px}
.scalefig{margin:12px 0 14px;border:1px solid var(--line);border-radius:8px;background:#fff;
padding:8px;overflow-x:auto}
.scalefig img{display:block;width:100%;min-width:680px;height:auto}
.scalefig figcaption{font-size:12px;color:#4a5360;margin-top:7px;line-height:1.5}
</style></head><body>
<header>
  <h1>Mammal diel activity across the lower 48</h1>
  <div class="sub">Predicted 24-hour activity from camera-trap detections &middot; sun-anchored time
  &middot; how much of a local rhythm is the species, and how much is the place</div>
</header>
""" + banner + """
<main>
  <div class="headline">
    <h2>A daily rhythm is mostly the species. The part that belongs to the place is real, small,
    and predictable for only one of these five animals</h2>
    <p><b>Click a cell anywhere on the map and most of what you see would look the same anywhere
    else in the country.</b> For white-tailed deer, <b>88% of the curve</b> at any location is
    simply the nationwide deer rhythm; for black bear it is about half. So the curve panel draws
    the nationwide rhythm behind your clicked cell and shades the gap between them &mdash; the gap
    is the part that is about the place you clicked.</p>
    <p><b>Only the gray squirrel has earned a national map.</b> Holding out whole regions and
    predicting them, the squirrel surface beats "assume this animal keeps the same hours
    everywhere" on 4 of its 5 measures. Coyote and raccoon manage 2 and 1; black bear's 2 rest on
    24 held-out arrays; and white-tailed deer &mdash; the largest sample here at 532 arrays &mdash;
    manages 1 and is on average <i>worse</i> than assuming deer are the same everywhere. Species
    and measures below are ordered and marked accordingly.</p>
    <p class="hl-why"><b>The reason is not the one you would guess.</b> Deer and coyotes are the
    two species whose behaviour genuinely varies between regions (59% and 42% of their
    array-to-array differences), and they predict worst. Raccoons and squirrels have essentially no
    regional pattern (0.0% and 0.5%) and squirrel predicts best. Deer's regional structure is
    driven by something this model never measured; squirrels have no regional structure but do
    respond to the conditions we did measure. <a href="#results">The full result is below</a>.</p>
  </div>

  <div class="tablab">Species &mdash; ordered by whether a national map is supported</div>
  <p class="tabhint" id="sphint"></p>
  <div class="tabs" id="sptabs" role="tablist"></div>
  <div class="tablab">Activity measure</div>
  <p class="tabhint" id="mhint"></p>
  <div class="tabs" id="mtabs" role="tablist"></div>

  <div class="layout">
    <div>
      <div class="card">
        <h2>Daily rhythm: the place against the species</h2>
        <p class="note" id="curvenote"></p>
        <div class="seg" role="group" aria-label="Curve display">
          <button id="btndev" aria-pressed="true">Difference from the species rhythm</button>
          <button id="btnabs" aria-pressed="false">Absolute curve</button>
        </div>
        <div id="decreadout"></div>
        <div id="curves"></div>
      </div>
    </div>

    <div>
      <div class="card">
        <h2 id="maptitle">&nbsp;</h2>
        <div id="verdictbadge"></div>
        <p class="vnums" id="verdictdetail"></p>
        <p class="note" id="mapnote">&nbsp;</p>
        <div class="seg" role="group" aria-label="Map type">
          <button id="btncov" aria-pressed="true">Local conditions (the tested map)</button>
          <button id="btnpos" aria-pressed="false">Position only (fails everywhere)</button>
        </div>
        <p class="note" id="mtnote" style="margin:7px 0 0"></p>
        <div class="seg" role="group" aria-label="Map layer">
          <button id="btnval" aria-pressed="true">Measure value</button>
          <button id="btnci" aria-pressed="false">Uncertainty (interval width)</button>
        </div>
        <div class="seg" style="margin-left:8px" role="group" aria-label="Colour scale">
          <button id="btnpool" aria-pressed="true">Scale: all species</button>
          <button id="btnsp" aria-pressed="false">Stretch to this species</button>
        </div>
        <p class="note" id="scalenote" style="margin:7px 0 0"></p>
        <div class="mapwrap mapfr">
        <svg id="map" viewBox="0 0 980 620" role="img"
             aria-label="Map of the selected activity measure across the lower 48"></svg>
        <div id="mapstamp"></div>
        <div id="mapcaveat"></div>
        </div>
        <div class="legend">
          <span id="leglo"></span><div class="swatches" id="legsw"></div><span id="leghi"></span>
          <span id="legunit"></span>
        </div>
        <div id="conflegend"></div>
        <p class="note" style="margin-top:8px">Masked to each species' union range (IUCN/MDD polygon
        plus extending 100&nbsp;km iNaturalist cells). Cells outside the mask are not drawn.</p>
        <p class="note" id="seamnote" style="margin-top:4px"></p>
        <p class="note" style="margin-top:4px"><b>Range constraint.</b> Predicted curves are
        constrained to the joint region real camera arrays occupy, so no cell is shown a value
        outside the observed across-array range for that species and measure. The check is
        regenerated from the shipped page data at build time and
        <a href="#results">tabulated in the results section</a>.</p>
      </div>
    </div>
  </div>

  <div class="card" id="supcard" style="margin-top:20px">
    <h2>Which species have earned a national map</h2>
    <p class="note">Whole regions were hidden, the model refitted without them, and each map asked
    to predict the hidden regions. A map is worth having only if it beats assuming the animal keeps
    the same hours everywhere. Counted over all five measures per species.</p>
    <div id="suptable"></div>
    <p class="note" style="margin-top:10px" id="supfoot"></p>
  </div>

  <div class="card" style="margin-top:20px">
    <h2>Does the map predict regions we never sampled?</h2>
    <p class="note">The selected measure, for every species. Positive means the map beat assuming
    one nationwide rhythm; the interval is the range consistent with the held-out regions.
    <span style="color:#4a5360">Hover any heading marked <sup>?</sup> for a plain-language
    explanation.</span></p>
    <div class="tw"><div id="cvtable"></div></div>
  </div>

  <div class="card" style="margin-top:20px">
    <h2>What activity changes with: the local-condition relationships</h2>
    <p class="note" id="panelnote"></p>
    <div class="grid2" id="panels"></div>
    <p class="note" id="panelhidden" style="margin-top:10px"></p>
  </div>
""" + RESULTS + """
  <div class="card" id="tiercard" style="margin-top:20px">
    <h2>Does the fading actually track worse predictions?</h2>
    <p class="note">A confidence shading is only worth having if error really does get worse as it
    fades. Each held-out camera array was assigned a band using only the arrays it was <i>not</i>
    fitted with. Error is scaled so 1.00 is what you get by assuming the animal keeps the same
    hours everywhere &mdash; <b>below 1.00 the map helps, above 1.00 it hurts</b>. Empty cells mean
    no held-out array fell in that band.</p>
    <div class="tw"><div id="tiertable"></div></div>
    <p class="note" style="margin-top:10px" id="tiernote"></p>
  </div>

  <div class="card" style="margin-top:20px">
    <h2>Model fit</h2>
    <p class="note" id="fitnote">Numbers for the selected species.</p>
    <div id="fitcards"></div>
  </div>

  <div class="card" style="margin-top:20px">
    <h2>What the predictors mean</h2>
    <p class="note">Click any entry for the full definition and how it was measured.</p>
    <div id="glossary"></div>
  </div>

  <div class="foot">
    <strong>Methods:</strong> <a id="methodslink" href="#">full methods document</a> &mdash;
    modelling framework, sun-time transformation and its effort consequence, spatial autocorrelation
    and effective sample size, range-mask construction, and known biases.
  </div>
</main>

<div id="tip"></div>
<div id="modal"><div id="modalbox">
  <button class="x" id="modalx" aria-label="Close">&times;</button>
  <h3 id="modaltitle"></h3><p id="modalcap"></p><div id="modalbody"></div>
</div></div>
""")


    # JS is assembled as ONE script element (JS_1 + JS_2) so every function shares one scope and the
    # strict verifier's DOM-id scan sees the whole program.
    JS_1 = """
const D = JSON.parse(document.getElementById('payload').textContent);
document.getElementById('methodslink').href = __METHODS_HREF__;
let SP = D.sp_order[0], MET = D.metrics[0].key, LAYER = 'value', PICKED = null, SCALE = 'pooled';
let MAPTYPE = 'cov';          // the tested map leads; 'position' is the one that fails everywhere
let CURVEMODE = 'deviation';  // DEFAULT: show the departure from the species rhythm, not the curve

/* ---------------- local-conditions map: rebuild curves from harmonic coefficients ----------------
   The payload carries four non-intercept harmonic coefficients per species x cell (quantised by
   D.cov.coefq). Every measure here is a property of curve SHAPE, so those four numbers determine
   all five and the 48-bin click curve. Same arithmetic as the fitting pipeline. */
const NBIN = 48;
const BCEN = Array.from({length:NBIN},(_,i)=>(i+0.5)*0.5);
const REF_SR = 6.2202, REF_SS = 17.9998;
const XS1=[],XC1=[],XS2=[],XC2=[],ANG=[],ISN=[],ISC=[];
for(let i=0;i<NBIN;i++){
  const t=BCEN[i]/24*2*Math.PI;
  XS1.push(Math.sin(t)); XC1.push(Math.cos(t));
  XS2.push(Math.sin(2*t)); XC2.push(Math.cos(2*t));
  ANG.push(t);
  ISN.push((BCEN[i]<REF_SR||BCEN[i]>=REF_SS)?1:0);
  ISC.push((Math.abs(BCEN[i]-REF_SR)<=1.5||Math.abs(BCEN[i]-REF_SS)<=1.5)?1:0);
}
function covCurve(sp, ci){
  const C=D.cov, b=C.coefs[sp][ci], q=C.coefq;
  const s1=b[0]/q, c1=b[1]/q, s2=b[2]/q, c2=b[3]/q;
  const e=new Array(NBIN); let mx=-Infinity;
  for(let i=0;i<NBIN;i++){ const v=s1*XS1[i]+c1*XC1[i]+s2*XS2[i]+c2*XC2[i]; e[i]=v; if(v>mx)mx=v; }
  let tot=0; for(let i=0;i<NBIN;i++){ e[i]=Math.exp(e[i]-mx); tot+=e[i]; }
  for(let i=0;i<NBIN;i++) e[i]/=(tot||1);
  return e;
}
function covMetrics(p){
  let noct=0,crep=0,Cc=0,Ss=0,mean=0,mx=-Infinity,mi=0;
  for(let i=0;i<NBIN;i++){
    noct+=ISN[i]*p[i]; crep+=ISC[i]*p[i];
    Cc+=p[i]*Math.cos(ANG[i]); Ss+=p[i]*Math.sin(ANG[i]);
    mean+=p[i]; if(p[i]>mx){mx=p[i]; mi=i;}
  }
  return {pct_nocturnal:100*noct, crepuscular:100*crep,
          concentration:Math.sqrt(Cc*Cc+Ss*Ss), peak_hour:BCEN[mi],
          activity_level:(mean/NBIN)/(mx||1)};
}
const _covCache={};
function covSurface(sp, met){
  const key=sp+'|'+met; if(_covCache[key]) return _covCache[key];
  const C=D.cov, n=C.cells.length, vals=new Float64Array(n), idx=[];
  for(let i=0;i<n;i++){
    // an all-zero coefficient vector marks a cell outside this species' range mask; asserted at
    // build time to be an exact presence test, which saves shipping a separate flag array
    const b=C.coefs[sp][i];
    if(b[0]===0&&b[1]===0&&b[2]===0&&b[3]===0) continue;
    vals[i]=covMetrics(covCurve(sp,i))[met]; idx.push(i);
  }
  return (_covCache[key]={vals:vals, idx:idx});
}
function covPath(i){
  const a=D.cov.paths[i], x0=a[0]/10, y0=a[1]/10;
  return `M${x0.toFixed(1)},${y0.toFixed(1)} ${((a[0]+a[2])/10).toFixed(1)},${((a[1]+a[3])/10).toFixed(1)}`
       + ` ${((a[0]+a[4])/10).toFixed(1)},${((a[1]+a[5])/10).toFixed(1)}`
       + ` ${((a[0]+a[6])/10).toFixed(1)},${((a[1]+a[7])/10).toFixed(1)}Z`;
}
const TIER_ALPHA=[1.0, 0.72, 0.42, 0.20];
const TIER_DESAT=[0.0, 0.30, 0.62, 0.88];
function desat(hex, f){
  const r=parseInt(hex.slice(1,3),16), g=parseInt(hex.slice(3,5),16), b=parseInt(hex.slice(5,7),16);
  const y=0.299*r+0.587*g+0.114*b;
  const m=v=>Math.round(v+(y-v)*f);
  return '#'+[m(r),m(g),m(b)].map(v=>Math.max(0,Math.min(255,v)).toString(16).padStart(2,'0')).join('');
}
const tip = document.getElementById('tip');
const fmt = (v,d=1)=> (v==null||!isFinite(v)) ? '\\u2014' : (+v).toFixed(d);
const sgn = (v,d=3)=> (v==null||!isFinite(v)) ? '\\u2014' : ((v>=0?'+':'')+(+v).toFixed(d));
function isCyclic(met){ return (D.cyclic||[]).indexOf(met)>=0; }
function rampCol(t, ci, met){
  // peak timing is an hour on a circle: a cyclic ramp so 23.8 h and 0.2 h render alike
  const cyc = (!ci && met && isCyclic(met) && D.ramp_cyc);
  const r = ci ? D.ramp_ci : (cyc ? D.ramp_cyc : D.ramp_value);
  let u = Math.max(0,Math.min(1,t));
  return r[Math.round(u*(r.length-1))]; }

/* ---------------- SUPPORT STATE. One function decides it, everywhere. ----------------
   'cov'      : per species x measure, from beats_fix in covariate_model_cv_fixed.csv, downgraded
                to 'weak' when the species' whole result rests on too few held-out arrays.
   'position' : never supported. Position alone beat one nationwide rhythm for no species and no
                measure, so every position surface is stamped. */
function supportOf(sp, met){
  const S=D.support[sp]||{};
  if(MAPTYPE==='position') return {state:'unsupported', tier:S.tier, rec:null};
  const rec=((D.cov.cv[sp]||{})[met])||null;
  if(!rec) return {state:'untested', tier:S.tier, rec:null};
  if(S.tier==='weak') return {state:(rec.beats?'weak':'unsupported'), tier:S.tier, rec:rec};
  return {state:(rec.beats?'supported':'unsupported'), tier:S.tier, rec:rec};
}
function tabs(el, items, get, set){
  el.innerHTML='';
  items.forEach(it=>{
    const b=document.createElement('button');
    b.className='tab'+(it.cls?(' '+it.cls):''); b.setAttribute('role','tab');
    b.innerHTML=it.html; b.setAttribute('aria-selected', get()===it.key);
    if(it.title) b.title=it.title;
    b.onclick=()=>{ set(it.key); render(); };
    el.appendChild(b);
  });
}

/* ---------------- THE DEVIATION CHART: the point of this page ----------------
   Top panel: the species' nationwide rhythm as a reference line, the clicked place over it, and
   the gap between them shaded. Bottom panel: the same gap on its own axis, centred on zero, which
   is what makes a 12%-local deer curve look different from a 50%-local bear curve. */
function devChart(hs, cell, mean, w, h, small, xlabel, ylabel){
  const pl = small?42:60, pr=8, pt=8, gap=small?18:24;
  const hA = Math.round((h-pt-gap-(small?26:44))*0.62), hB = (h-pt-gap-(small?26:44))-hA;
  const yA0=pt, yB0=pt+hA+gap;
  const xmin=hs[0], xmax=hs[hs.length-1];
  const X=v=>pl+(v-xmin)/((xmax-xmin)||1)*(w-pl-pr);
  let aMin=Infinity,aMax=-Infinity,dMax=0;
  const dif=[];
  for(let i=0;i<cell.length;i++){
    aMin=Math.min(aMin,cell[i],mean[i]); aMax=Math.max(aMax,cell[i],mean[i]);
    const d=cell[i]-mean[i]; dif.push(d); dMax=Math.max(dMax,Math.abs(d));
  }
  const padA=(aMax-aMin)*0.08||0.1; aMin-=padA; aMax+=padA;
  dMax=Math.max(dMax*1.12, 1e-6);
  const YA=v=>yA0+(1-(v-aMin)/((aMax-aMin)||1))*hA;
  const YB=v=>yB0+(1-(v+dMax)/(2*dMax))*hB;
  const px=hs.map(X);
  let s='';
  // -- top panel: shaded departure between the two curves, drawn as one closed band
  const up=px.map((x,i)=>`${x.toFixed(1)},${YA(cell[i]).toFixed(1)}`);
  const dn=px.map((x,i)=>`${x.toFixed(1)},${YA(mean[i]).toFixed(1)}`).reverse();
  s+=`<polygon points="${up.concat(dn).join(' ')}" fill="#7b52ab" opacity=".30"></polygon>`;
  s+=`<line x1="${pl}" y1="${yA0}" x2="${pl}" y2="${yA0+hA}" stroke="#949ca6"></line>`
   + `<line x1="${pl}" y1="${yA0+hA}" x2="${w-pr}" y2="${yA0+hA}" stroke="#949ca6"></line>`;
  s+=`<polyline points="${px.map((x,i)=>x.toFixed(1)+','+YA(mean[i]).toFixed(1)).join(' ')}"`
   + ` fill="none" stroke="#4a5360" stroke-width="${small?1.5:2}" stroke-dasharray="7,4"></polyline>`;
  s+=`<polyline points="${px.map((x,i)=>x.toFixed(1)+','+YA(cell[i]).toFixed(1)).join(' ')}"`
   + ` fill="none" stroke="#1d4e80" stroke-width="${small?1.8:2.4}"></polyline>`;
  const fs=small?9:11;
  s+=`<text x="${pl-4}" y="${yA0+7}" font-size="${fs}" text-anchor="end" fill="#4a5360">`
   + `${fmt(aMax,2)}</text>`
   + `<text x="${pl-4}" y="${yA0+hA}" font-size="${fs}" text-anchor="end" fill="#4a5360">`
   + `${fmt(aMin,2)}</text>`;
  // -- bottom panel: the difference alone, zero-centred, so the shape of the departure is readable
  const zp=px.map((x,i)=>`${x.toFixed(1)},${YB(dif[i]).toFixed(1)}`);
  const base=`${px[px.length-1].toFixed(1)},${YB(0).toFixed(1)} ${px[0].toFixed(1)},${YB(0).toFixed(1)}`;
  s+=`<polygon points="${zp.join(' ')} ${base}" fill="#7b52ab" opacity=".42"></polygon>`;
  s+=`<line x1="${pl}" y1="${YB(0)}" x2="${w-pr}" y2="${YB(0)}" stroke="#2c3540"></line>`
   + `<line x1="${pl}" y1="${yB0}" x2="${pl}" y2="${yB0+hB}" stroke="#949ca6"></line>`;
  s+=`<polyline points="${zp.join(' ')}" fill="none" stroke="#5b3b86" stroke-width="${small?1.5:2}">`
   + `</polyline>`;
  s+=`<text x="${pl-4}" y="${yB0+7}" font-size="${fs}" text-anchor="end" fill="#4a5360">`
   + `+${fmt(dMax,2)}</text>`
   + `<text x="${pl-4}" y="${yB0+hB}" font-size="${fs}" text-anchor="end" fill="#4a5360">`
   + `\\u2212${fmt(dMax,2)}</text>`;
  s+=`<text x="${pl}" y="${yB0+hB+fs+3}" font-size="${fs}" fill="#4a5360">${fmt(xmin,0)}</text>`
   + `<text x="${w-pr}" y="${yB0+hB+fs+3}" font-size="${fs}" text-anchor="end" fill="#4a5360">`
   + `${fmt(xmax,0)}</text>`;
  const lfs=small?9.5:11.5;
  s+=`<text x="${(w+pl)/2}" y="${h-3}" font-size="${lfs}" text-anchor="middle" fill="#2c3540">`
   + `${xlabel}</text>`;
  if(!small){
    s+=`<text x="13" y="${yA0+hA/2}" font-size="10.5" text-anchor="middle" fill="#2c3540"`
     + ` transform="rotate(-90 13 ${yA0+hA/2})">${ylabel}</text>`
     + `<text x="13" y="${yB0+hB/2}" font-size="10.5" text-anchor="middle" fill="#2c3540"`
     + ` transform="rotate(-90 13 ${yB0+hB/2})">this place minus species</text>`;
  }
  return `<svg viewBox="0 0 ${w} ${h}" role="img" aria-label="Upper panel: the nationwide rhythm`
   + ` for this species as a dashed line with this place's rhythm over it and the gap shaded.`
   + ` Lower panel: that gap alone, centred on zero.">${s}</svg>`;
}
function devKey(){
  return `<div class="devkey">`
   + `<span><i style="border-top-color:#1d4e80"></i>this place</span>`
   + `<span><i style="border-top-color:#4a5360;border-top-style:dashed"></i>`
   + `nationwide rhythm for this species</span>`
   + `<span><i class="band"></i>the difference \\u2014 the only part about this place</span></div>`;
}
function lineChart(o, w, h, small){
  const xs=o.x||o.h, ys=o.y||o.r, lo=o.lo, hi=o.hi;
  const xmin=Math.min(...xs), xmax=Math.max(...xs);
  const ymin=Math.min(...(lo||ys)), ymax=Math.max(...(hi||ys));
  const pl = small?38:60, pb = small?30:48, pt=6, pr=6;
  const X=v=>pl+(v-xmin)/((xmax-xmin)||1)*(w-pl-pr);
  const Y=v=>pt+(1-(v-ymin)/((ymax-ymin)||1))*(h-pt-pb);
  let band='';
  if(lo&&hi){ const up=xs.map((x,i)=>`${X(x).toFixed(1)},${Y(hi[i]).toFixed(1)}`);
    const dn=xs.map((x,i)=>`${X(x).toFixed(1)},${Y(lo[i]).toFixed(1)}`).reverse();
    band=`<polygon points="${up.concat(dn).join(' ')}" fill="#2b6cb0" opacity=".16"></polygon>`; }
  const line=xs.map((x,i)=>`${X(x).toFixed(1)},${Y(ys[i]).toFixed(1)}`).join(' ');
  const ax=`<line x1="${pl}" y1="${h-pb}" x2="${w-pr}" y2="${h-pb}" stroke="#949ca6"></line>`
          +`<line x1="${pl}" y1="${pt}" x2="${pl}" y2="${h-pb}" stroke="#949ca6"></line>`;
  const fs = small?9:11;
  const labs=`<text x="${pl-4}" y="${pt+7}" font-size="${fs}" text-anchor="end" fill="#4a5360">${fmt(ymax,1)}</text>`
   +`<text x="${pl-4}" y="${h-pb}" font-size="${fs}" text-anchor="end" fill="#4a5360">${fmt(ymin,1)}</text>`
   +`<text x="${pl}" y="${h-pb+fs+3}" font-size="${fs}" fill="#4a5360">${fmt(xmin,1)}</text>`
   +`<text x="${w-pr}" y="${h-pb+fs+3}" font-size="${fs}" text-anchor="end" fill="#4a5360">${fmt(xmax,1)}</text>`;
  const lfs = small?9.5:11.5;
  const xl = o.xlabel || 'x', yl = o.ylabel || 'y';
  const axt = `<text x="${(w+pl)/2}" y="${h-3}" font-size="${lfs}" text-anchor="middle" `
            + `fill="#2c3540">${xl}</text>`
            + `<text x="${small?10:12}" y="${pt+(h-pt-pb)/2}" font-size="${lfs}" `
            + `text-anchor="middle" fill="#2c3540" `
            + `transform="rotate(-90 ${small?10:12} ${pt+(h-pt-pb)/2})">${yl}</text>`;
  return `<svg viewBox="0 0 ${w} ${h}" role="img" aria-label="${yl} against ${xl}">${ax}${band}`
    +`<polyline points="${line}" fill="none" stroke="#2b6cb0" stroke-width="${small?1.6:2.1}"></polyline>`
    +`${labs}${axt}</svg>`;
}
function estChart(o, w, h, small){
  const labs=['Naive','Between-array','Within-array'];
  const ys=o.y, lo=o.lo, hi=o.hi;
  const ymin=Math.min(0,...lo), ymax=Math.max(0,...hi);
  const pl= small?38:64, pb= small?32:52, pt=8, pr=10;
  const X=i=>pl+(ys.length<2?0.5:i/(ys.length-1))*(w-pl-pr);
  const Y=v=>pt+(1-(v-ymin)/((ymax-ymin)||1))*(h-pt-pb);
  let s=`<line x1="${pl}" y1="${Y(0)}" x2="${w-pr}" y2="${Y(0)}" stroke="#949ca6" `
      + `stroke-dasharray="3,3"></line>`;
  ys.forEach((v,i)=>{
    const sig=(lo[i]>0||hi[i]<0), col=sig?'#8a4a08':'#5b6470';
    s += `<line x1="${X(i)}" y1="${Y(lo[i])}" x2="${X(i)}" y2="${Y(hi[i])}" stroke="${col}" `
       + `stroke-width="${small?1.6:2.4}"></line>`
       + `<circle cx="${X(i)}" cy="${Y(v)}" r="${small?3:5}" fill="${col}"></circle>`;
    if(!small) s += `<text x="${X(i)}" y="${h-pb+16}" font-size="11" text-anchor="middle" `
       + `fill="#4a5360">${labs[i]||''}</text>`
       + `<text x="${X(i)}" y="${Y(v)-10}" font-size="11" text-anchor="middle" fill="#2c3540">`
       + `${fmt(v,2)}</text>`;
  });
  const fs= small?9:11;
  s += `<text x="${pl-4}" y="${pt+7}" font-size="${fs}" text-anchor="end" fill="#4a5360">${fmt(ymax,1)}</text>`
     + `<text x="${pl-4}" y="${h-pb}" font-size="${fs}" text-anchor="end" fill="#4a5360">${fmt(ymin,1)}</text>`;
  if(small) labs.forEach((L,i)=>{ s += `<text x="${X(i)}" y="${h-pb+12}" font-size="8.5" `
     + `text-anchor="middle" fill="#4a5360">${L.slice(0,6)}</text>`; });
  const lfs2 = small?9.5:11.5;
  const xl2 = o.xlabel || 'Estimator (how the effect is computed)';
  const yl2 = o.ylabel || 'Effect size (change in metric units)';
  s += `<text x="${(w+pl)/2}" y="${h-3}" font-size="${lfs2}" text-anchor="middle" `
     + `fill="#2c3540">${xl2}</text>`
     + `<text x="${small?10:12}" y="${pt+(h-pt-pb)/2}" font-size="${lfs2}" text-anchor="middle" `
     + `fill="#2c3540" transform="rotate(-90 ${small?10:12} ${pt+(h-pt-pb)/2})">${yl2}</text>`;
  return `<svg viewBox="0 0 ${w} ${h}" role="img" aria-label="${yl2} by ${xl2}">${s}</svg>`;
}
function selectCell(cid, flash){
  const svg=document.getElementById('map');
  svg.querySelectorAll('.cell.sel').forEach(q=>q.classList.remove('sel'));
  const el = (typeof cid==='string' && cid.startsWith('COV:'))
    ? svg.querySelector(`.covcell[data-ci="${cid.slice(4)}"]`)
    : svg.querySelector(`.cell[data-cid="${cid}"]`);
  if(el){ el.classList.add('sel');
    if(flash){ el.classList.add('flash'); setTimeout(()=>el.classList.remove('flash'), 1400); } }
  return !!el;
}
"""


    JS_2 = """
/* ---------------- the support stamp: an unsupported surface must not look supported ------------ */
function drawStamp(){
  const el=document.getElementById('mapstamp');
  const met=D.metrics.find(m=>m.key===MET), S=supportOf(SP,MET), rec=S.rec;
  if(S.state==='supported'){ el.className=''; el.innerHTML=''; return; }
  el.className='nsstamp';
  if(MAPTYPE==='position'){
    el.innerHTML = `NOT SUPPORTED \\u2014 position alone predicts nothing`
      + `<span class="sub2">No species and no measure beat assuming the animal keeps the same`
      + ` hours everywhere. Use the local-conditions map.</span>`;
    return;
  }
  if(S.state==='untested'){
    el.innerHTML = `NOT TESTED for this measure<span class="sub2">No held-out result exists for`
      + ` this combination.</span>`;
    return;
  }
  if(S.state==='weak'){
    el.innerHTML = `WEAKLY SUPPORTED \\u2014 too little held-out data`
      + `<span class="sub2">${SP} beats one nationwide rhythm on ${met.label.toLowerCase()}`
      + ` (${sgn(rec.skill)}), but on only ${rec.n_te} held-out arrays across ${rec.n_folds}`
      + ` folds. Not comparable to the gray squirrel's evidence.</span>`;
    return;
  }
  el.innerHTML = `NOT SUPPORTED \\u2014 this surface does not beat assuming ${SP.toLowerCase()}`
    + ` behave the same everywhere`
    + `<span class="sub2">${met.label}: ${sgn(rec.skill)} (range ${sgn(rec.lo)} to`
    + ` ${sgn(rec.hi)}, ${rec.folds} held-out regions improved). Read the pattern below as`
    + ` decoration, not as a finding.</span>`;
}
function drawVerdict(){
  const bb=document.getElementById('verdictbadge'), vd=document.getElementById('verdictdetail');
  const mc=document.getElementById('mapcaveat');
  const met=D.metrics.find(m=>m.key===MET), S=supportOf(SP,MET), rec=S.rec, sup=D.support[SP]||{};
  const cls={supported:'supported', unsupported:'unsupported', weak:'weak', untested:''}[S.state];
  const gly={supported:'\\u25b2', unsupported:'\\u25cb', weak:'\\u25d4', untested:'\\u2014'}[S.state];
  const word={supported:'Supported', unsupported:'Not supported', weak:'Weakly supported',
              untested:'Not tested'}[S.state];
  bb.innerHTML=`<span class="badge ${cls}"><span class="bglyph" aria-hidden="true">${gly}</span>`
    + `${word} \\u2014 ${met.label.toLowerCase()}, ${MAPTYPE==='position'?'position only'
    :'local conditions'}</span>`;
  let t='';
  if(MAPTYPE==='position'){
    t = `<div class="vnums">This surface uses <b>geographic position alone</b>. Tested by hiding`
      + ` whole regions, it beat a single nationwide rhythm for <b>no species and no measure</b>,`
      + ` so it is shown only for comparison. The local-conditions map is the one with results.`
      + `</div>`;
  } else if(!rec){
    t = `<div class="vnums">No held-out result exists for this species and measure.</div>`;
  } else {
    t = `<div class="vnums">On <b>${met.label}</b>, hiding whole regions and predicting them from`
      + ` local conditions: prediction error is <b>${fmt(Math.abs(rec.skill)*100,0)}%`
      + ` ${rec.skill>=0?'smaller':'larger'}</b> than assuming this animal keeps the same hours`
      + ` everywhere (${sgn(rec.skill)}, range ${sgn(rec.lo)} to ${sgn(rec.hi)},`
      + ` <b>${rec.folds}</b> held-out regions improved, ${rec.n_te} arrays).`
      + (rec.beats ? ` The range stays above zero.` : ` The range includes zero, so this surface`
        + ` is <b>not</b> shown to help and should not be read as a finding.`)
      + `</div>`;
  }
  t += `<div class="vnums" style="color:#3c4653"><b>${SP}:</b> ${sup.glyph} ${sup.tier_label}`
    + ` \\u2014 beats one nationwide rhythm on <b>${sup.beats} of ${sup.n}</b> measures, average`
    + ` margin <b>${sgn(sup.mean_skill)}</b>, ${sup.n_arrays} held-out arrays over`
    + ` ${sup.n_folds} folds. ${sup.note}</div>`;
  vd.innerHTML=t;
  mc.className='mapcav';
  const ciw=D.honest_ci[SP];
  mc.innerHTML = (MAPTYPE==='position')
    ? `Position-only surface, kept for comparison \\u2014 coordinates alone did not predict any`
      + ` measure for any species.`
      + `<span class="ciw">Honest interval at a new array: ${fmt(ciw,1)} percentage points wide`
      + `</span>`
    : `Predictive, not causal \\u2014 these are the rhythms we expect where conditions take these`
      + ` values, not evidence those conditions cause the rhythm.`
      + `<span class="ciw">${word} for ${met.label.toLowerCase()} &middot; fades where conditions`
      + ` were never sampled</span>`;
}
function drawCovMap(){
  const svg=document.getElementById('map'), met=D.metrics.find(m=>m.key===MET);
  const snote=document.getElementById('seamnote');
  if(snote){ const ns=(D.cov.seam&&D.cov.seam[SP])?D.cov.seam[SP].reduce((a,x)=>a+x,0):0;
    snote.innerHTML = ns
      ? `<b>Seam repair.</b> The 25\\u00a0km grid column near 96\\u00b0\\u00a0W (${ns} cells for`
        + ` this species) was missing from an earlier build and has been restored, with its`
        + ` conditions measured directly from the source rasters at the same spatial support as`
        + ` every other cell.`
      : ''; }
  const C=D.cov, S=covSurface(SP,MET), ci=(LAYER==='ci');
  const wq=C.wq[MET]||1, wArr=(C.ciw[SP]||{})[MET]||null;
  let vals=[];
  if(SCALE==='species'){ S.idx.forEach(i=>vals.push(ci?(wArr?wArr[i]/wq:0):S.vals[i])); }
  else { D.species.forEach(sp=>{ const T=covSurface(sp,MET), w=(C.ciw[sp]||{})[MET];
      T.idx.forEach(i=>vals.push(ci?(w?w[i]/wq:0):T.vals[i])); }); }
  vals.sort((a,b)=>a-b);
  const qt=p=>vals.length?vals[Math.min(vals.length-1,Math.max(0,Math.round(p*(vals.length-1))))]:0;
  // a cyclic measure spans its FULL period: percentile-stretching an angle would make a
  // 24-minute difference across midnight read as the largest contrast on the map
  const cyc=(!ci && isCyclic(MET));
  const lo = cyc?0:qt(0.02), hi = cyc?24:qt(0.98);
  const unsup = supportOf(SP,MET).state!=='supported';
  let out='', hatch=[];
  for(const i of S.idx){
    const tier=(C.conf[SP][i]||1)-1;
    const v = ci ? (wArr?wArr[i]/wq:0) : S.vals[i];
    const t = (v-lo)/((hi-lo)||1);
    const col = desat(rampCol(t,ci,MET), TIER_DESAT[tier]);
    out += `<path class="cell covcell" d="${covPath(i)}" fill="${col}"`
        +  ` fill-opacity="${TIER_ALPHA[tier]}" data-ci="${i}" data-tier="${tier}"></path>`;
    if(tier===3) hatch.push(i);
  }
  out += hatch.map(i=>`<path d="${covPath(i)}" fill="url(#hatchsev)" stroke="none"></path>`).join('');
  // an unsupported surface is hatched WHOLE, so it cannot be mistaken for a supported one even in
  // a screenshot, in greyscale, or by a colour-blind reader
  if(unsup) out += `<rect x="0" y="0" width="980" height="620" fill="url(#hatchunsup)"`
                + ` pointer-events="none"></rect>`;
  out += D.states.map(d=>`<path d="${d}" class="statebd"></path>`).join('');
  out = `<defs><pattern id="hatchsev" width="5" height="5" patternUnits="userSpaceOnUse"
          patternTransform="rotate(45)"><line x1="0" y1="0" x2="0" y2="5" stroke="#5b6470"
          stroke-width="1.1"></line></pattern>
         <pattern id="hatchunsup" width="9" height="9" patternUnits="userSpaceOnUse"
          patternTransform="rotate(-45)"><line x1="0" y1="0" x2="0" y2="9" stroke="#8f1414"
          stroke-width="1.5" opacity=".34"></line></pattern></defs>` + out;
  svg.innerHTML=out;
  svg.querySelectorAll('.covcell').forEach(pEl=>{
    const i=+pEl.dataset.ci, tier=+pEl.dataset.tier;
    pEl.onmousemove=e=>{ tip.style.opacity=1; tip.style.left=(e.clientX+13)+'px';
      tip.style.top=(e.clientY+13)+'px';
      const w=wArr?wArr[i]/wq:null;
      tip.innerHTML=`<b>${met.label}</b>: ${fmt(S.vals[i],2)} ${met.unit}`
        + (w!=null?`<br>95% interval width ${fmt(w,2)}`:'')
        + `<br>confidence: <b>${C.tiers[tier]}</b> (score ${C.score[SP][i]}/100)`
        + (unsup?`<br><b>this surface is not supported</b>`:'')
        + `<br><i>click for this place\\u2019s rhythm</i>`; };
    pEl.onmouseleave=()=>tip.style.opacity=0;
    pEl.onclick=()=>{ PICKED='COV:'+i; selectCell(PICKED,false); drawCurves(); };
  });
  document.getElementById('maptitle').textContent =
    `${SP} \\u2014 ${met.label} \\u2014 local conditions`;
  drawVerdict(); drawStamp();
  document.getElementById('mapnote').innerHTML =
    `<i>${D.sci[SP]}</i> &middot; ${met.desc} &middot; ${S.idx.length} cells at 25 km`;
  document.getElementById('legsw').innerHTML =
    Array.from({length:41},(_,i)=>`<i style="background:${rampCol(i/40,ci,MET)}"></i>`).join('');
  document.getElementById('leglo').textContent=fmt(lo,2);
  document.getElementById('leghi').textContent=fmt(hi,2);
  document.getElementById('legunit').textContent = ci ? 'interval width ('+met.unit+')'
    : (cyc ? met.unit+' (colour wraps at 24 h \\u2014 the scale is a circle)' : met.unit);
  const pctv=C.tier_pct[SP]||{};
  document.getElementById('conflegend').innerHTML =
    `<div class="conftitle">How far each place sits outside the conditions we sampled</div>`
    + `<div class="confrow">` + C.tiers.map((tn,k)=>{
        const sw=desat('#2f6b9a',TIER_DESAT[k]);
        return `<span class="confk"><i style="background:${sw};opacity:${TIER_ALPHA[k]}`
             + (k===3?';background-image:repeating-linear-gradient(45deg,#5b6470 0 1px,transparent 1px 5px)':'')
             + `"></i>${tn}<b>${fmt(pctv[tn],0)}%</b></span>`; }).join('')
    + `</div><p class="note" style="margin:5px 0 0">Percentages are this species' share of its own`
    + ` range. Paler and hatched places are drawn but should not be read as findings.</p>`;
  const sn=document.getElementById('scalenote');
  if(sn) sn.innerHTML = (SCALE==='species')
    ? `<b>Colours are stretched to this species' own range</b>, so contrast here does not indicate`
      + ` a large effect.`
    : `Colours use one scale across all five species, so maps are comparable between animals.`;
}
function drawMap(){
  if(MAPTYPE==='cov' && D.cov){ drawCovMap(); return; }
  const sn0=document.getElementById('seamnote'); if(sn0) sn0.innerHTML='';
  const svg=document.getElementById('map'), met=D.metrics.find(m=>m.key===MET);
  const surf=D.surfaces[SP][MET]||{};
  const dm = (SCALE==='species' && D.domsp[SP]) ? D.domsp[SP][MET] : D.dom[MET];
  const ci = LAYER==='ci';
  const lo = ci?dm.wlo:dm.vlo, hi = ci?dm.whi:dm.vhi;
  let out='';
  for(const cid in surf){
    const rec=surf[cid], v = ci?rec[1]:rec[0], t=(v-lo)/((hi-lo)||1);
    out += `<path class="cell" d="${D.paths[cid]}" fill="${rampCol(t,ci,MET)}" data-cid="${cid}"`
        +  ` data-v="${rec[0]}" data-w="${rec[1]}" data-lo="${rec[2]}" data-hi="${rec[3]}"></path>`;
  }
  // the position map is never supported, so it is always hatched whole
  out += `<defs><pattern id="hatchunsup" width="9" height="9" patternUnits="userSpaceOnUse"
          patternTransform="rotate(-45)"><line x1="0" y1="0" x2="0" y2="9" stroke="#8f1414"
          stroke-width="1.5" opacity=".34"></line></pattern></defs>`
       + `<rect x="0" y="0" width="980" height="620" fill="url(#hatchunsup)" pointer-events="none">`
       + `</rect>`;
  out += D.states.map(d=>`<path d="${d}" class="statebd"></path>`).join('');
  svg.innerHTML=out;
  const hasCC = D.cell_curves[SP] && Object.keys(D.cell_curves[SP]).length;
  svg.querySelectorAll('.cell').forEach(p=>{
    p.onmousemove=e=>{ tip.style.opacity=1; tip.style.left=(e.clientX+13)+'px';
      tip.style.top=(e.clientY+13)+'px';
      tip.innerHTML=`<b>${met.label}</b>: ${fmt(p.dataset.v,2)} ${met.unit}`
        +`<br>95% interval ${fmt(p.dataset.lo,2)} \\u2013 ${fmt(p.dataset.hi,2)}`
        +`<br><b>position-only surface \\u2014 not supported</b>`
        +(hasCC?'<br><i>click for this cell\\u2019s rhythm</i>':''); };
    p.onmouseleave=()=>tip.style.opacity=0;
    p.onclick=()=>{ if(!hasCC) return; PICKED=p.dataset.cid; selectCell(PICKED,false); drawCurves(); };
  });
  document.getElementById('maptitle').textContent=`${SP} \\u2014 ${met.label} \\u2014 position only`;
  drawVerdict(); drawStamp();
  const cl=document.getElementById('conflegend'); if(cl) cl.innerHTML='';
  document.getElementById('mapnote').innerHTML =
    `<i>${D.sci[SP]}</i> &middot; ${met.desc} &middot; ${Object.keys(surf).length} cells at 100 km`;
  document.getElementById('legsw').innerHTML =
    Array.from({length:41},(_,i)=>`<i style="background:${rampCol(i/40,ci,MET)}"></i>`).join('');
  document.getElementById('leglo').textContent=fmt(lo,2);
  document.getElementById('leghi').textContent=fmt(hi,2);
  document.getElementById('legunit').textContent = ci ? 'interval width ('+met.unit+')' : met.unit;
  const sn=document.getElementById('scalenote');
  if(sn) sn.innerHTML = (SCALE==='species')
    ? `<b>Colours are stretched to this species' own range</b>, so contrast here does not indicate`
      + ` a large effect.`
    : `Colours use one scale across all five species, so maps are comparable between animals.`;
}
/* ---------------- curves: DEVIATION IS THE DEFAULT ---------------- */
function decReadout(){
  const el=document.getElementById('decreadout'), d=D.dec[SP];
  if(!el) return;
  if(!d){ el.innerHTML=''; return; }
  el.innerHTML = `<div class="decbox">`
    + `<span class="big">${fmt(d.pct_species,0)}% of a ${SP.toLowerCase()} curve is the same`
    + ` everywhere in the country.</span>`
    + `<div class="decbar" role="img" aria-label="${fmt(d.pct_species,0)} percent species rhythm,`
    + ` ${fmt(d.pct_cell,0)} percent specific to the place">`
    + `<i class="sp" style="width:${d.pct_species}%"></i>`
    + `<i class="cl" style="width:${d.pct_cell}%"></i></div>`
    + `Only <b>${fmt(d.pct_cell,0)}%</b> is specific to the place you click`
    + ` (${d.n_cells.toLocaleString()} places compared). ${d.note}</div>`;
}
function drawCurves(){
  const box=document.getElementById('curves'); let out='';
  const cc=D.cell_curves[SP]||{}, dev=(CURVEMODE==='deviation');
  const XL='Sun-anchored hour of day (h, 0\\u201324)';
  const YL='Share of daily activity (% per half-hour bin)';
  document.getElementById('curvenote').innerHTML = dev
    ? `The dashed line is this species' nationwide rhythm; the solid line is the place you clicked;`
      + ` the shaded gap between them is the only part that is about that place. The lower panel`
      + ` shows that gap on its own.`
    : `The clicked place's rhythm on its own. Most of this shape is the species, not the place`
      + ` \\u2014 switch back to the difference view to see which part is which.`;
  decReadout();
  if(PICKED && typeof PICKED==='string' && PICKED.startsWith('COV:') && D.cov){
    const i=+PICKED.slice(4), C=D.cov, tier=(C.conf[SP][i]||1)-1;
    const p=covCurve(SP,i).map(v=>v*100);
    const mean=(D.cov_mean[SP]||[]).map(v=>v*100);
    const warn = tier>=2
      ? `<div class="hint" style="color:#6d2a0f"><b>${C.tiers[tier]}</b> \\u2014 conditions here`
        + ` sit outside what the model learned from, so treat this shape as an illustration.</div>`
      : `<div class="hint">Confidence: <b>${C.tiers[tier]}</b> (score ${C.score[SP][i]}/100).</div>`;
    let chart, key='';
    if(dev && mean.length===p.length){
      chart=devChart(BCEN, p, mean, 0, 356, 300, false, XL, YL); key=devKey();
    } else {
      const w=0.18;
      chart=lineChart({h:BCEN.slice(), r:p, lo:p.map(v=>Math.max(v*(1-w),0)),
                       hi:p.map(v=>v*(1+w)), xlabel:XL, ylabel:YL},356,200,false);
    }
    out += `<div class="selcurve"><div class="lab">Selected place \\u2014 ${C.cells[i]}`
        + ` (local conditions)</div>`
        + `<div class="hint">Rebuilt from this place's conditions. Click any other cell to change`
        + ` it.</div>` + warn + key + chart
        + `<div style="font-size:11.5px;margin-top:2px"><a href="#" id="clearpick">Clear`
        + ` selection</a></div></div>`;
  }
  else if(PICKED&&cc[PICKED]){
    // cell_curves ship per-mille scaled so the 48 bins MEAN 1000; /480 gives percent share of day
    const q=cc[PICKED], r=q.r.map(v=>v/480), w=(q.w||180)/1000;
    const mean=(D.cc_mean[SP]||[]).map(v=>v/480);
    const hs=Array.from({length:r.length},(_,i)=>i*0.5+0.25);
    let chart, key='';
    if(dev && mean.length===r.length){
      chart=devChart(hs, r, mean, 0, 356, 300, false, XL, YL); key=devKey();
    } else {
      chart=lineChart({h:hs, r:r, lo:r.map(v=>Math.max(v*(1-w),0)), hi:r.map(v=>v*(1+w)),
                       xlabel:XL, ylabel:YL},356,200,false);
    }
    out += `<div class="selcurve"><div class="lab">Selected cell \\u2014 ${PICKED}</div>`
        + `<div class="hint">This is the cell outlined on the map.</div>` + key + chart
        + `<div style="font-size:11.5px;margin-top:2px"><a href="#" id="clearpick">Clear`
        + ` selection</a></div></div>`;
  } else {
    const mean=(MAPTYPE==='cov'?(D.cov_mean[SP]||[]).map(v=>v*100)
                              :(D.cc_mean[SP]||[]).map(v=>v/480));
    let chart='';
    if(mean.length) chart=lineChart({h:BCEN.slice(0,mean.length), r:mean, xlabel:XL, ylabel:YL},
                                    356,190,false);
    out += `<div class="selcurve"><div class="lab">The nationwide ${SP.toLowerCase()} rhythm</div>`
        + `<div class="hint">This is the reference curve \\u2014 the rhythm this species keeps`
        + ` across the whole country. <b>Click any cell on the map</b> to draw that place against`
        + ` it and shade the difference.</div>` + chart + `</div>`;
  }
  box.innerHTML=out;
  const cp=document.getElementById('clearpick');
  if(cp) cp.onclick=e=>{ e.preventDefault(); PICKED=null;
    document.querySelectorAll('.cell.sel').forEach(q=>q.classList.remove('sel')); drawCurves(); };
}
function drawSupport(){
  const el=document.getElementById('suptable'); if(!el) return;
  let out='';
  D.sp_order.forEach(sp=>{
    const s=D.support[sp];
    out += `<div class="suprow${sp===SP?' here':''}"><div class="who">`
      + `<span class="glyph" aria-hidden="true">${s.glyph}</span>${sp}`
      + `<span class="num">${s.tier_label} &middot; ${s.beats} of ${s.n} measures &middot;`
      + ` average margin ${sgn(s.mean_skill)} &middot; ${s.n_arrays} arrays,`
      + ` ${s.n_folds} folds</span></div><div>${s.note}</div></div>`;
  });
  el.innerHTML=out;
  const f=document.getElementById('supfoot');
  if(f) f.innerHTML = `Species tabs above follow this order. A measure that does not beat one`
    + ` nationwide rhythm is struck through in the measure tabs, stamped on the map and hatched`
    + ` across the whole surface \\u2014 supported and unsupported maps are never drawn alike.`
    + ` The position-only map is unsupported for every species and measure.`;
}
function th(key, fallback){
  const g=D.cv_help[key];
  if(!g) return `<th>${fallback}</th>`;
  return `<th><span class="hdr" data-k="${key}">${g[0]}<sup>?</sup></span></th>`;
}
function drawCV(){
  const met=D.metrics.find(m=>m.key===MET);
  let out='<table class="cv"><colgroup><col class="c0"><col class="c1"><col class="cn">'
        + '<col class="cn"><col class="cn"></colgroup>'
        + '<tr><th>Species</th><th>Verdict for this measure</th>'
        + th('skill_unit','Margin over one nationwide rhythm') + th('skill_ci','Range')
        + th('folds_pos','Regions improved') + '</tr>';
  D.sp_order.forEach(sp=>{
    const r=(D.cov.cv[sp]||{})[MET], s=D.support[sp];
    const st = r ? (s.tier==='weak' ? (r.beats?'weak':'unsupported')
                                    : (r.beats?'supported':'unsupported')) : 'untested';
    const word={supported:'\\u25b2 supported', unsupported:'\\u25cb not supported',
                weak:'\\u25d4 too little data', untested:'\\u2014 not tested'}[st];
    const cls={supported:'ok', unsupported:'no', weak:'', untested:''}[st];
    if(!r){ out+=`<tr${sp===SP?' style="background:#f2f7fc"':''}><td>${sp}</td>`
      + `<td colspan="4" style="text-align:left;color:#4a5360">not tested for this measure</td></tr>`;
      return; }
    out+=`<tr${sp===SP?' style="background:#f2f7fc"':''}><td>${sp}</td>`
      +`<td style="text-align:left" class="${cls}">${word}</td>`
      +`<td>${sgn(r.skill)}</td><td>${sgn(r.lo)} to ${sgn(r.hi)}</td><td>${r.folds}</td></tr>`;
  });
  out+='</table>';
  out+=`<p class="note" style="margin-top:8px">Showing <b>${met.label}</b>. Across all`
     + ` ${D.cov.n_tot} species-and-measure combinations, <b>${D.cov.n_beat} beat one nationwide`
     + ` rhythm and ${D.cov.n_tot-D.cov.n_beat} do not</b>. Positive margins mean lower prediction`
     + ` error in held-out regions.</p>`;
  document.getElementById('cvtable').innerHTML=out;
}
function drawPanels(){
  const box=document.getElementById('panels');
  const arr=(D.panels[SP]||{})[MET]||[];
  const pc=((D.pcount[SP]||{})[MET])||{kept:0,total:0,hidden:[]};
  const met=D.metrics.find(m=>m.key===MET);
  const chart=(p,w,h,sm)=> (p.kind==='estimators') ? estChart(p,w,h,sm) : lineChart(p,w,h,sm);
  document.getElementById('panelnote').innerHTML =
    `How <b>${SP}</b> activity changes with local conditions \\u2014 ${met.label}. Showing the`
    + ` <b>${pc.kept} of ${pc.total}</b> predictors whose interval excludes no effect somewhere`
    + ` across its range. Each caption starts with what the relationship means; click to enlarge.`;
  if(!arr.length){
    box.innerHTML = `<div class="nopanel"><b>No predictor shows a clear effect on ${met.label} for`
      + ` ${SP}.</b> All ${pc.total} fitted relationships have intervals wide enough to include a`
      + ` flat line, so none is shown. Try another measure.</div>`;
  } else {
    box.innerHTML = arr.map((p,i)=>
      `<div class="mini" data-i="${i}"><b>${p.title}</b>${chart(p,205,110,true)}`
      + `<span class="take">${p.take}</span></div>`).join('');
  }
  document.getElementById('panelhidden').innerHTML = pc.hidden.length
    ? `<b>${pc.hidden.length} panel${pc.hidden.length>1?'s':''} hidden</b> because the interval`
      + ` never excludes a flat line: ${pc.hidden.join(', ')}.`
    : `No panels hidden for this selection.`;
  box.querySelectorAll('.mini').forEach(el=>{
    el.onclick=()=>{ const p=arr[+el.dataset.i];
      document.getElementById('modaltitle').textContent=`${SP} \\u2014 ${p.title}`;
      document.getElementById('modalcap').textContent=p.caption;
      document.getElementById('modalbody').innerHTML=chart(p,800,400,false);
      document.getElementById('modal').classList.add('on'); };
  });
}
function drawFit(){
  const fm=D.fit_model[SP]||[], fe=D.fit_ess[SP]||[];
  let out='<div class="tw"><table class="fit">'
        + '<colgroup><col class="c0"><col class="cn"><col class="cn"><col class="cn">'
        + '<col class="cn"></colgroup>'
        + '<tr><th>Model</th>' + th('dev_explained','Dev. expl.')
        + th('edf_diel','EDF diel') + th('edf_spatial','EDF spatial')
        + th('dispersion','Disp.') + '</tr>';
  fm.forEach(r=>{ out+=`<tr><td>${r.model}</td><td>${fmt(r.dev_explained*100,1)}%</td>`
    +`<td>${fmt(r.edf_diel,1)}</td><td>${fmt(r.edf_spatial,1)}</td>`
    +`<td>${fmt(r.dispersion,2)}</td></tr>`; });
  out+='</table></div>';
  const al=D.act_summary[SP];
  if(al){ out += `<div style="margin-top:11px;padding:8px 10px;border:1px solid var(--line);`
    + `border-radius:7px;background:#fafafa;font-size:12.5px">`
    + `<b>Activity level (mean across mapped cells): ${fmt(al.suntime,3)}</b>`
    + ` &middot; range ${fmt(al.suntime_min,3)}\\u2013${fmt(al.suntime_max,3)}`
    + `<div style="color:#4a5360;margin-top:3px">A relative index for comparing places within this`
    + ` species \\u2014 where animals are active for more of the day, not how many hours. The`
    + ` 30-minute independence filter raises it by about ${fmt(al.filter_pct,1)}% here.</div></div>`; }
  if(fe.length){ const e=fe[0];
    const kv=(key,val)=>{ const g=D.cv_help[key];
      return `<div class="kv"><span class="hdr" data-k="${key}">${g?g[0]:key}<sup>?</sup></span>`
           + `<span>${val}</span></div>`; };
    out+='<div style="margin-top:12px">'
      + kv('nominal_n_deployments', e.nominal_n_deployments)
      + kv('n_arrays', e.n_arrays)
      + kv('effective_n', e.effective_n)
      + kv('ci_width_naive', fmt(e.ci_width_naive,2))
      + kv('ci_width_hierarchical', fmt(e.ci_width_hierarchical,2))
      + kv('ci_inflation_factor', fmt(e.ci_inflation_factor,2)+'\\u00d7')
      +'</div>'; }
  document.getElementById('fitcards').innerHTML=out;
}
function drawTiers(){
  const el=document.getElementById('tiertable'); if(!el||!D.cov||!D.cov.tierval) return;
  const TVv=D.cov.tierval, tiers=D.cov.tiers;
  let out='<table class="tier"><colgroup><col class="c0"><col class="cn"><col class="cn">'
        + '<col class="cn"><col class="cn"></colgroup><tr><th>Animal</th>'
        + tiers.map(t=>`<th>${t}</th>`).join('') + '</tr>';
  const rows={}; TVv.forEach(r=>{ (rows[r.species]=rows[r.species]||{})[r.tier]=r; });
  Object.keys(rows).forEach(sp=>{
    out+=`<tr${sp===SP?' style="background:#f2f7fc"':''}><td>${sp}</td>`;
    tiers.forEach(t=>{ const r=rows[sp][t];
      out += (!r||r.n_arrays===0||r.nmae_cov==null) ? `<td>\\u2014</td>`
        : `<td>${fmt(r.nmae_cov,2)}<span style="color:#4a5360;font-size:11px">`
          + ` (${r.n_arrays})</span></td>`; });
    out+='</tr>';
  });
  el.innerHTML=out+'</table>';
  const tn=document.getElementById('tiernote');
  if(tn) tn.innerHTML = `Each cell shows relative error with the number of held-out arrays in`
    + ` brackets. Pooling all five animals, error rises across the four bands`
    + ` (${(D.cov.tierpool||[]).map(v=>fmt(v,2)).join(' \\u2192 ')}), so the shading reflects`
    + ` measured degradation rather than an asserted threshold.`;
}
function drawGlossary(){
  document.getElementById('glossary').innerHTML = Object.keys(D.glossary).map(k=>{
    const g=D.glossary[k];
    return `<details><summary><b>${g.label}</b> &mdash; ${g.short}</summary>`
         + `<p style="margin:7px 0 2px;color:#333b45;font-size:13px">${g.long}</p></details>`;
  }).join('');
}
function wireHelp(){
  document.querySelectorAll('.hdr').forEach(el=>{
    el.onmouseenter=e=>{ const g=D.cv_help[el.dataset.k]; if(!g) return;
      tip.style.opacity=1; tip.innerHTML=`<b>${g[0]}</b><br>${g[1]}`;
      const r=el.getBoundingClientRect();
      tip.style.left=Math.min(r.left, window.innerWidth-300)+'px';
      tip.style.top=(r.bottom+8)+'px'; };
    el.onmouseleave=()=>tip.style.opacity=0;
  });
}
function render(){
  // species tabs in support order, each carrying its tier word and glyph
  tabs(document.getElementById('sptabs'), D.sp_order.map(s=>{
    const u=D.support[s];
    return {key:s, html:`<span class="glyph" aria-hidden="true">${u.glyph}</span>`
      + `<span class="tlab">${s}</span><span class="tsub">${u.tier_label} \\u00b7 ${u.beats}/${u.n}`
      + ` measures</span>`, title:u.note};
  }), ()=>SP, v=>{SP=v; PICKED=null;});
  // measure tabs: struck through and worded when this species x measure does not beat one
  // nationwide rhythm, so an unsupported selection is visible BEFORE it is chosen
  tabs(document.getElementById('mtabs'), D.metrics.map(m=>{
    const st=supportOf(SP,m.key).state;
    const word={supported:'supported', unsupported:'not supported', weak:'too little data',
                untested:'not tested'}[st];
    return {key:m.key, cls:(st==='supported'?'':'unsup'),
            html:`<span class="tlab">${m.label}</span><span class="tsub">${word}</span>`,
            title:`${m.label}: ${word} for ${SP}`};
  }), ()=>MET, v=>{MET=v;});
  document.getElementById('sphint').innerHTML =
    `\\u25c9 defensible national map \\u00b7 \\u25d4 too little held-out data to judge \\u00b7`
    + ` \\u25cb partly supported. Only the gray squirrel clears 4 of 5 measures.`;
  document.getElementById('mhint').innerHTML = (MAPTYPE==='position')
    ? `Every measure is unsupported on the position-only map.`
    : `Struck-through measures do not beat assuming ${SP.toLowerCase()} behave the same everywhere.`;
  document.getElementById('btnpos').setAttribute('aria-pressed', MAPTYPE==='position');
  document.getElementById('btncov').setAttribute('aria-pressed', MAPTYPE==='cov');
  document.getElementById('btndev').setAttribute('aria-pressed', CURVEMODE==='deviation');
  document.getElementById('btnabs').setAttribute('aria-pressed', CURVEMODE==='absolute');
  const bc=document.getElementById('btncov');
  if(!D.cov){ bc.disabled=true; bc.title='local-conditions payload not built'; }
  document.getElementById('mtnote').innerHTML = (MAPTYPE==='cov')
    ? `Predicted from the <b>conditions</b> at each place \\u2014 woodland, people, cropland,`
      + ` pasture, ruggedness, summer and winter temperature \\u2014 at 25&nbsp;km. This is the map`
      + ` that was tested and partly works.`
    : `Predicted from <b>position</b> alone at 100&nbsp;km. Shown for comparison only: it beat one`
      + ` nationwide rhythm for no species and no measure.`;
  document.getElementById('btnval').setAttribute('aria-pressed', LAYER==='value');
  document.getElementById('btnci').setAttribute('aria-pressed', LAYER==='ci');
  document.getElementById('btnpool').setAttribute('aria-pressed', SCALE==='pooled');
  document.getElementById('btnsp').setAttribute('aria-pressed', SCALE==='species');
  drawMap(); drawCurves(); drawSupport(); drawCV(); drawPanels(); drawFit(); drawTiers();
  drawGlossary(); wireHelp();
}
document.getElementById('btnpos').onclick=()=>{MAPTYPE='position';PICKED=null;render();};
document.getElementById('btncov').onclick=()=>{if(!D.cov)return;MAPTYPE='cov';PICKED=null;render();};
document.getElementById('btndev').onclick=()=>{CURVEMODE='deviation';drawCurves();
  document.getElementById('btndev').setAttribute('aria-pressed',true);
  document.getElementById('btnabs').setAttribute('aria-pressed',false);};
document.getElementById('btnabs').onclick=()=>{CURVEMODE='absolute';drawCurves();
  document.getElementById('btndev').setAttribute('aria-pressed',false);
  document.getElementById('btnabs').setAttribute('aria-pressed',true);};
document.getElementById('btnval').onclick=()=>{LAYER='value';render();};
document.getElementById('btnci').onclick=()=>{LAYER='ci';render();};
document.getElementById('btnpool').onclick=()=>{SCALE='pooled';render();};
document.getElementById('btnsp').onclick=()=>{SCALE='species';render();};
document.getElementById('modalx').onclick=()=>document.getElementById('modal').classList.remove('on');
document.getElementById('modal').onclick=e=>{ if(e.target.id==='modal')
  document.getElementById('modal').classList.remove('on'); };
document.addEventListener('keydown',e=>{ if(e.key==='Escape')
  document.getElementById('modal').classList.remove('on'); });
render();
"""

    html += ('<script id="payload" type="application/json">' + _json_strict(payload) + "</script>\n"
             + "<script>\n"
             + (JS_1 + JS_2).replace("__METHODS_HREF__", json.dumps(methods_href))
             + "\n</script></body></html>")

    with open(out, "w") as fh:
        fh.write(html)
    return out, len(html), DEC, RANGE, SUP


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--indir", default="in")
    ap.add_argument("--out", default=OUT)
    ap.add_argument("--methods", default="methods.html")
    a = ap.parse_args()
    d = lambda f: os.path.join(a.indir, f)
    contract = json.load(open(d("viewer_contract.json")))
    res = build(
        pd.read_csv(d("viewer_surfaces.csv")),
        pd.read_csv(d("viewer_curves.csv")),
        pd.read_csv(d("viewer_panels.csv")),
        pd.read_csv(d("viewer_fit_model.csv")),
        pd.read_csv(d("viewer_fit_validation.csv")),
        pd.read_csv(d("viewer_fit_ess.csv")),
        d("viewer_grid.geojson"), contract, methods_href=a.methods,
        l48_cells=d("viewer_cells_l48.csv"), states_json=d("viewer_states.json"),
        cell_curves=pd.read_csv(d("viewer_cellcurves.csv")),
        act_summary=json.load(open(d("viewer_activity_summary.json"))),
        coefs=pd.read_csv(d("covariate_coefs_fixed.csv")),
        gridcov=pd.read_csv(d("grid_covariates_25km_seamfixed.csv")),
        conf=pd.read_csv(d("national_confidence_fixed.csv")),
        preds=pd.read_csv(d("covariate_predictions_fixed.csv")),
        cv=pd.read_csv(d("covariate_model_cv_fixed.csv")),
        tierval=pd.read_csv(d("tier_validation_fixed.csv")),
        cellcurves_long=pd.read_csv(d("covariate_cellcurves_fixed.csv")),
        harmonics=pd.read_csv(d("array_harmonics.csv")),
        varpart=pd.read_csv(d("variance_partition.csv")),
        hyp=pd.read_csv(d("hypothesis_results_standardized.csv")),
        nulltest=pd.read_csv(d("density_null_test.csv")),
        sqx=pd.read_csv(d("squirrel_natural_experiment.csv")),
        msel=pd.read_csv(d("metric_selection.csv")),
        scale_png=d("spatial_scale_decomposition.png"), out=a.out)
    outp, n, DEC, RANGE, SUP = res
    dec = pd.DataFrame([{"species": s, "pct_species_typical": round(v["pct_species"], 2),
                         "pct_cell_specific": round(v["pct_cell"], 2), "n_cells": v["n_cells"]}
                        for s, v in DEC.items()]).sort_values("pct_cell_specific")
    dec.to_csv("curve_decomposition.csv", index=False)
    pd.DataFrame(RANGE).to_csv("range_check_shipped.csv", index=False)
    n_out = int(sum(r["n_outside"] for r in RANGE))
    print(f"{outp}  {n/1e6:.2f} MB  range_outside={n_out}")
    print(dec.to_string(index=False))


if __name__ == "__main__":
    main()
