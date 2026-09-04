
import html as _h, pandas as pd, numpy as np

def build_report2(figs, H, cn, EX, strong, cv, cmp_vp, guilds, host, outfile="analysis_report.html"):
    IMG = {k: _img(v[0], host) for k, v in figs.items()}
    cvb = cv[cv.beats_null == True].sort_values("skill", ascending=False)

    def fig(key, title, caption, plain=None):
        if key not in IMG: return ""
        s = (f'<figure><img src="{IMG[key]}" alt="{_h.escape(title)}">'
             f'<figcaption><b>{_h.escape(title)}.</b> {caption}</figcaption></figure>')
        if plain: s += f'<div class="plain"><b>Reading this figure:</b> {plain}</div>'
        return s

    S = []
    S.append(f"""
<section id="summary">
<h2><span class="num">1</span>What we find</h2>
<p class="lead">Thirteen North American mammal and bird species, 686 camera sites, 1.22 million
independent detections from August to October. For each species at each site we measure the shape of
the daily activity cycle.</p>

<div class="grid">
  <div class="card"><div class="big">13</div><div class="lbl">species</div></div>
  <div class="card"><div class="big">686</div><div class="lbl">camera sites</div></div>
  <div class="card"><div class="big">1.22M</div><div class="lbl">detections</div></div>
  <div class="card"><div class="big">1,935</div><div class="lbl">species-by-site activity curves</div></div>
</div>

<p><b>Result 1.</b> Each species keeps a recognisable daily schedule across the continent. A raccoon
in Oregon and a raccoon in Georgia are active at close to the same times relative to sunrise and
sunset. Across all 13 species, the spread among sites is small compared with the differences between
species.</p>

<p><b>Result 2.</b> Site-to-site differences within a species are real and connect to measurable
conditions. Fifty-six species-by-condition relationships hold with and without a spatial control.
Black bears are more nocturnal where summers are hot and less nocturnal in rugged terrain; red foxes
are less nocturnal where people are dense; opossums concentrate around dawn and dusk in farmland.</p>

<p><b>Result 3.</b> How much a species varies between sites is set by where its schedule already sits,
not by its body size. A species averaging 50% night activity has room to move in both directions and
varies by 13 to 17 percentage points between sites. A species averaging 3% or 97% has little room and
varies by 2 to 5 points. The rank correlation between site-to-site spread and distance from the
all-day or all-night extremes is +0.92; with body mass it is +0.41 (p = 0.73).</p>
</section>""")

    dl, il, nl = guilds["diurnal"], guilds["inter"], guilds["nocturnal"]
    S.append(f"""
<section id="curves">
<h2><span class="num">2</span>The daily schedules</h2>
<p class="lead">Every species-by-site activity curve in the study, grouped by species.</p>

{fig("fig_curves_13species", "Daily activity of 13 species, and of every site within each species",
 "One panel per species, ordered from most daytime to most night-time. The heavy line is that species' "
 "schedule pooled across all sites; each faint line is one site with at least 50 detections. The "
 "shaded band is night. Time runs from sunrise through midday, sunset and midnight back to sunrise.",
 "The horizontal axis is stretched separately for each site and date so that sunrise and sunset always "
 "fall in the same position. A 6 a.m. detection in Maine in October and a 6 a.m. detection in Arizona "
 "in August sit at different points on this axis, because the sun is in a different place. Without "
 "this adjustment, latitude and season would masquerade as behaviour. Species shapes are distinct and "
 "the faint site lines cluster tightly around each heavy line.")}

<p>At the 20% and 80% night-activity boundaries, the species divide into three groups:</p>
<ul>
<li><b>Daytime ({len(dl)} species):</b> {", ".join(dl)}</li>
<li><b>Both day and night ({len(il)}):</b> {", ".join(il)}</li>
<li><b>Night ({len(nl)}):</b> {", ".join(nl)}</li>
</ul>
</section>""")

    # ---- hypotheses
    h1 = H[(H.covariate=="pop_1km") & (H.metric=="pct_noct")].sort_values("body_kg")
    h3 = H[(H.covariate=="tmax_hottest_month") & (H.metric=="pct_noct")].sort_values("body_kg")
    from scipy.stats import spearmanr
    r1 = spearmanr(h1.body_kg, h1.beta_std); r3 = spearmanr(h3.body_kg, h3.beta_std)
    ag = H[H.covariate=="ag_5km"].pivot(index="species", columns="metric", values="beta_std")
    agp = H[H.covariate=="ag_5km"].pivot(index="species", columns="metric", values="p")
    n_conc = int((agp["conc"]<0.05).sum()); n_noct = int((agp["pct_noct"]<0.05).sum())
    n_bigger = int((ag["conc"].abs() > ag["pct_noct"].abs()).sum())

    tH = h1.copy()
    tH["Species"] = tH.species; tH["Body mass (kg)"] = tH.body_kg
    tH["Sites"] = tH.n; tH["Effect"] = tH.beta_std; tH["p"] = tH.p
    tH["With spatial control"] = tH.beta_std_spatial; tH["p (spatial)"] = tH.p_spatial
    tbl_h1 = _table(tH, ["Species","Body mass (kg)","Sites","Effect","p","With spatial control","p (spatial)"],
                    highlight=lambda r: r["p"] < 0.05)

    S.append(f"""
<section id="predictions">
<h2><span class="num">3</span>Predictions, and how they fared</h2>
<p class="lead">Three predictions were set before the tests were run. Each names a mechanism and a
species ordering, so each can fail.</p>

<h3>Prediction 1: urban areas concentrate activity more than farmland does</h3>
<p>Mechanism: a food source that is available on a predictable schedule should pack an animal's
activity into a narrow window. Urban food is available year-round at consistent times, from bins,
compost and feeders. Crops are seasonal, with a harvest that removes the resource entirely.</p>
<p><b>Outcome: supported in direction, not in strength.</b> Urban density concentrates activity more
than farmland does in {int((out.urban > out.agriculture).sum())} of 13 species. The paired difference
across species is not significant (p = {pw:.2f}), so the ordering holds as a tendency rather than a
established difference. Six species show a significant urban effect and five a significant farmland
effect. Raccoon fits the mechanism most cleanly, with urban concentrating activity (+0.20, p &lt; 0.001)
and farmland loosening it (&minus;0.15, p = 0.001), opposite signs on the two axes within one species.</p>
<p>A second test of the same mechanism: if urban food is more consistent, urban sites should also be
more consistent in behaviour. Among the 8 species with enough sites in both settings, site-to-site
spread is smaller in high-urban than in high-farmland settings for
{int(CV2.urban_more_consistent.sum())} of 8 (mean spread {CV2.sd_urban.mean():.1f} against
{CV2.sd_ag.mean():.1f} points, paired p = {pw2:.2f}).</p>

<h3>Prediction 2: heat suppresses midday activity, more so in large species</h3>
<p>Mechanism: a large body sheds heat slowly, because heat production scales with mass while heat loss
scales with surface area. Midday heat should therefore cost a bear more than a chipmunk, even though
the bear tolerates cold better.</p>
<p><b>Outcome: the mechanism holds for individual species; the size ordering does not.</b> Measuring
the share of activity in the three hours around solar noon, {int((MD.beta_pp<0).sum())} of 13 species
reduce midday activity where summers are hotter. Two survive a spatial control: black bear
(&minus;2.1 points per unit of heat, p &lt; 0.001) and gray squirrel (&minus;2.4 points, p &lt; 0.001).
Those two sit at opposite ends of the size range, which is the difficulty for the prediction: the rank
correlation between body mass and midday loss is {r_heat_prop[0]:+.2f} (p = {r_heat_prop[1]:.2f}) among species with
midday activity to lose.</p>
<p>Black bear also shows the largest heat effect on overall night activity in the study (+0.61 standard
deviations, p &lt; 0.001, unchanged by a spatial control), consistent with the mechanism for that
species.</p>

<h3>Prediction 3: activity shifts into the night where people are dense, most strongly in large species</h3>
<p>Mechanism: large mammals are hunted and are more conspicuous, so avoiding daylight near people
should pay off more for them.</p>
<p><b>Outcome: the species ordering is not supported.</b> Seven of 13 species show a relationship with
human population density, but the direction is inconsistent and does not track body mass (rank
correlation {r_pop[0]:+.2f}, p = {r_pop[1]:.2f}). Chipmunks, opossums and coyotes become more nocturnal
where people are dense; cottontails, red foxes and mule deer become less so. Red fox shows the largest
effect of any species, at 0.63 standard deviations, in the direction opposite to the prediction.</p>
{tbl_h1}
<p class="note">Effect is the change in night activity per one standard deviation of human population
density, in standard deviations of that species' night activity. Highlighted rows are significant at
p &lt; 0.05. The spatial control adds a linear north-south and east-west trend.</p>

{fig("fig_hypotheses", "Testing the three predictions",
 "Panel a pairs each species' urban and farmland effect on how tightly activity is packed into a short "
 "window. Panel b plots the change in midday activity share under heat against body mass, for species "
 "with at least 5% of activity at midday. Panel c plots site-to-site spread in night activity against "
 "how far that species' average sits from fully diurnal or fully nocturnal.",
 "In panel a, a species whose blue point sits right of its orange point concentrates activity more in "
 "urban settings than in farmland. In panel b, points below the line reduce midday activity where "
 "summers are hot; a size ordering would show as a trend across the horizontal axis. Panel c shows what "
 "does order site-to-site variation: a species averaging 50% night activity can shift in either "
 "direction, while one at 3% or 97% has little room.")}

<h3>What orders the variation instead</h3>
<p>Body mass ordered neither the heat response nor the human-density response. Site-to-site spread
tracks how far a species' average sits from fully diurnal or fully nocturnal, with a rank correlation
of +0.92. Fitting spread against both body mass and distance from those extremes leaves body mass at
p = 0.73 and distance from the extremes at p = 0.003. Black bear averages 46% night activity, closest
to an even split of any species here, and has both the largest site-to-site spread and the largest
heat effect.</p>
</section>""")

    # ---- what conditions matter
    ts = strong.head(18).copy()
    ts["Species"] = ts.species; ts["Sites"] = ts.n
    ts["Effect"] = ts.beta_std; ts["With spatial control"] = ts.beta_std_spatial
    tbl_s = _table(ts, ["Species","Condition","Measurement","Sites","Effect","With spatial control"])
    te = EX.copy()
    te["Species"] = te.species; te["Sites"] = te.n
    te["Share that is counting noise (%)"] = te.pct_is_measurement.round(0).astype(int)
    te["Share of the rest explained by habitat (%)"] = te.pct_real_explained.round(0).astype(int)
    tbl_e = _table(te, ["Species","Sites","Share that is counting noise (%)",
                        "Share of the rest explained by habitat (%)"])

    S.append(f"""
<section id="conditions">
<h2><span class="num">4</span>Which conditions move a species' schedule</h2>
<p class="lead">Fifty-six of 195 species-by-condition-by-measurement relationships are significant
both on their own and after adding a spatial trend, with the same sign in both. The largest are below.</p>
{tbl_s}
<p class="note">Effect is in standard deviations of the measurement per standard deviation of the
condition. Conditions are human population density, farmland cover, summer maximum temperature, tree
cover and terrain ruggedness.</p>

<h3>How much of the site-to-site variation these conditions account for</h3>
<p>Two questions can be separated. First, how much of the apparent variation between sites is simply
uncertainty from counting a limited number of photographs. Second, of the variation that remains, how
much do the six habitat measurements explain.</p>
{tbl_e}
<p class="note">The first column is computed from the standard error on each site's own activity
curve. The second is the reduction in residual variance when all six habitat measurements are fitted.</p>

<p>Counting noise accounts for 7% to 47% of the apparent variation depending on species, and is
largest where detections are sparse: wild turkey 44%, chipmunk 47%, red squirrel 32%. Of the variation
that is not counting noise, the six habitat measurements explain a median of 6.9% across species,
ranging from 34% for black bear down to none for wild turkey and red squirrel.</p>

<p>So the site-to-site differences are largely not explained by these habitat layers, all of which are
measured at 1 to 5 km resolution. Whether finer-resolution conditions would explain more is untested
here; the variation that remains is unaccounted for rather than shown to be unexplainable.</p>

{fig("fig_variance_simple", "Site-to-site differences, in three separate questions",
 "Panel a: how much sites differ within each species, with the grey whisker showing how much of that "
 "apparent difference is uncertainty from limited photograph counts. Panel b: the strength of evidence "
 "that nearby sites resemble each other. Panel c: the share of the real site-to-site variation "
 "accounted for by the six habitat measurements.",
 "Reading across the three panels gives each species' situation. Black bear differs most between sites "
 "(17 points), that difference follows geography, and habitat measurements account for a third of it. "
 "Wild turkey differs least, and what difference there is follows geography but is not captured by any "
 "habitat measurement we have.")}
</section>""")

    S.append(f"""
<section id="geography">
<h2><span class="num">5</span>Geographic pattern, species by species</h2>
<p class="lead">If a species' schedule were set by broad geography, two nearby sites would resemble
each other more than two distant sites. Seven of 13 species show this; six do not.</p>

{fig("fig_variograms_grouped", "Whether similarity between sites fades with distance",
 "Each point compares all pairs of sites separated by that distance, showing how different their night "
 "activity is. Top row: species where the curve rises with distance. Bottom row: species where it is "
 "flat. Red dashed line is the total spread among all sites for that species; grey dotted line is the "
 "spread expected from limited photograph counts alone.",
 "A curve that climbs and then flattens at the red line means nearby sites are similar and distant "
 "sites differ, which is the pattern geography produces. A flat curve means two sites 20 km apart "
 "differ about as much as two 2,000 km apart. The seven species in the top row pass a permutation "
 "test for short-distance similarity at p &lt; 0.05; the six in the bottom row do not.")}

<p>Species with geographic structure in night activity: wild turkey, red squirrel, black bear,
white-tailed deer, red fox, eastern cottontail, northern raccoon. Mule deer is borderline (p = 0.14).
Without: eastern chipmunk, gray squirrel, fox squirrel, coyote, opossum.</p>


</section>""")

    # ---- data & methods
    S.append(f"""
<section id="data">
<h2><span class="num">6</span>Data and measurement</h2>

{fig("fig_sampling", "Sampling coverage",
 "Panel a maps the cameras (grey) and the sites they group into (bubbles sized by camera count, "
 "coloured by how far apart the cameras are). Panel b shows cameras per site. Panel c shows survey "
 "effort by month with the analysis window in red. Panels d and e show detections and usable sites per "
 "species. Panel f places each species on the day-night axis.",
 "Coverage is dense in the eastern United States and sparse in the arid West and northern Plains. "
 "Survey effort is concentrated in September and October, so the analysis covers August to October; "
 "outside those months there is too little data to separate geography from season.")}

<h3>How a site's activity curve is built</h3>
<p>A camera records the time of each detection. Those times are converted to a scale where sunrise is
0 and sunset is 12, so that day length and latitude do not distort comparisons. Detections are counted
in 48 half-hour intervals on that scale.</p>

<p>The counts are then described by a curve built from two cycles: one that repeats once per day and
one that repeats twice. A once-per-day cycle captures a single activity peak; adding a twice-per-day
cycle captures two peaks, which is what dawn-and-dusk species show. Five numbers describe the curve:
an overall level, and two pairs describing the position and size of each cycle. Fitting a curve rather
than using the raw counts directly means each site's schedule is described by five numbers instead of
48, which makes sites comparable when some have thousands of detections and others have 25.</p>

<p>From each curve we read off: the share of activity falling at night, the share within 1.5 hours of
sunrise or sunset, and how tightly activity is packed into a short window rather than spread evenly.</p>

{fig("fig_s2_fitted_vs_raw", "Comparing the fitted curves against the raw detection counts",
 "For each species and measurement, the value read from the fitted curve plotted against the value "
 "computed directly from the raw detection counts with no curve fitting.",
 "A fitted curve can smooth away detail. This compares each number against the same number computed "
 "straight from the photographs. Points on the diagonal mean the curve is faithful. Three measurements "
 "sit on the diagonal: night share, dawn-dusk share, and concentration. A fourth, overall activity "
 "level, does not, and it is excluded from the rest of the analysis. The second axis on each panel is "
 "the highest correlation the data could produce given how few detections some sites have; a point "
 "below that line is limited by the fitting, not by the data.")}

{fig("fig_s1_effort_map_timeline", "Survey effort per site, and the effect of the sun-based time scale",
 "Effort per site on the map, deployment timelines, and the distribution of within-camera effort "
 "imbalance across time intervals on the sun-based scale compared with ordinary clock time.",
 "On the sun-based scale a camera does not record for equal lengths of real time in each interval, "
 "because day and night are stretched by different amounts. The imbalance reaches 1.6-fold at some "
 "cameras. Each count is adjusted by the exact recording time in its own interval. On the ordinary "
 "clock the equivalent measure is flat, as expected.")}

<h3>Time zones and camera clocks</h3>
<p>Sunrise and sunset are computed in local clock time using each site's time zone from the standard
time zone database, including daylight saving. Zone boundaries follow state lines rather than
longitude, and Arizona does not observe daylight saving, so the zone is looked up by location rather
than derived from longitude. The calculation agrees with published sunrise and sunset times to within
14 minutes at nine cities across four time zones.</p>

<p>Photograph times were checked project by project using species with fixed schedules: squirrels,
chipmunks and turkeys are active in daylight, and raccoons, opossums and skunks at night. Across the
full dataset the measured values are chipmunk 4.9%, gray squirrel 5.7% and raccoon 91.2% night
activity, against expectations of roughly 5%, 5% and 90%. One project's times were stored in universal
time and were shifted by five hours; eleven projects were taken from their original published archive.
Twenty-two individual cameras with faulty clocks were excluded, identified by strictly daytime species
recorded at night together with a second independent signal, and two of them confirmed from the
photographs.</p>

{fig("camera_trap_diagnostics", "Routine data checks",
 "Camera locations, deployment timelines by project, monthly effort, and the interval between "
 "consecutive deployments at the same location.",
 "The timeline shows which projects ran in which years and reveals deployments spanning unusually long "
 "periods. The gap distribution has two peaks: gaps of a few days, which are one survey split across "
 "records and are merged, and gaps of about a year, which are separate annual surveys at the same "
 "location and are kept as separate deployments within one site.")}
</section>""")

    # ---- prediction / map
    tcv = cvb.copy()
    tcv["Species"] = tcv.species
    tcv["Measurement"] = tcv.metric.map({"pct_noct":"% active at night","pct_crep":"% at dawn/dusk",
                                         "conc":"how concentrated","peak_h":"time of peak"}).fillna(tcv.metric)
    tcv["Predicted from"] = tcv.model.map({"covariate":"habitat at the site",
                                           "position":"map position"}).fillna(tcv.model)
    tcv["Accuracy gain"] = tcv.skill
    tcv["Plausible range"] = [f"{a:.2f} to {b:.2f}" for a,b in zip(tcv.lo, tcv.hi)]
    tbl_cv = _table(tcv, ["Species","Measurement","Predicted from","Accuracy gain","Plausible range"])

    S.append(f"""
<section id="mapping">
<h2><span class="num">7</span>Predicting an unsampled place</h2>
<p class="lead">A map of a species' activity is a prediction for places with no cameras. We tested
that directly.</p>

<p>The test hides every camera site inside a 400 km region, fits the model on the sites that remain,
and predicts the hidden ones. The median hidden site sits
{cv.median_heldout_dist_km.median():.0f} km from the nearest site used for fitting. The comparison is
against the simplest possible alternative: using that species' average schedule everywhere, ignoring
location and habitat entirely. An accuracy gain of zero means the model matches that alternative;
above zero means the model does better; below zero means worse.</p>

{fig("fig_s4_block_cv", "Prediction accuracy for unsampled regions",
 "Accuracy gain for every species and measurement, with plausible ranges. Zero is the accuracy of "
 "using one average schedule for the whole species. Also shown: prediction error against distance "
 "from the nearest fitting site, and the regional blocks used.",
 "Most points sit at or below zero, meaning a species' average schedule predicts a new region as well "
 "as location or habitat does. Eleven of 104 combinations have a plausible range entirely above zero.")}

<p>Those 11 are where a place-specific prediction adds accuracy:</p>
{tbl_cv}

<p>Black bear is the strongest case, and its accuracy comes from habitat at the site rather than map
position: night activity gain 0.43, time of peak 0.56. This matches the heat and terrain effects in
section 3.</p>

{fig("fig_prediction_raster", "Predicted night activity where local conditions carry predictive accuracy",
 "Predicted share of activity at night for black bear and coyote, at 25 km resolution. Open circles are "
 "the camera sites the models were fitted on. Cells are shown only where the local conditions fall "
 "inside the range the model was fitted on and inside the species' known range.",
 "Black bear is predicted from summer heat and terrain ruggedness alone. Fitting all six conditions to "
 "72 sites raised block-tested accuracy from 0.21 to 0.26 but produced a surface where neighbouring "
 "cells differed by 9.6 points on average, against a total spread between sites of 18.1 points; the "
 "two-condition version differs by 2.3 points between neighbours. Coyote uses all six conditions and "
 "its predicted range is narrow, from 57% to 86% night activity.")}

{fig("fig_s5_precision", "Prediction uncertainty compared with the variation being mapped",
 "For each species and measurement, the uncertainty on a prediction for a new site alongside the "
 "range of variation across sampled sites.",
 "For every species and measurement, the uncertainty on a prediction for a new site is between 1.00 "
 "and 1.60 times the full range of variation across sampled sites. A map of predicted values would "
 "show colour differences smaller than the uncertainty on each value.")}

<div class="callout"><h4>Where a map is supported</h4>
Black bear, on night activity and time of peak, predicted from habitat. Both gains have plausible
ranges entirely above zero (0.18 to 0.62, and 0.31 to 0.68), and black bear also has the largest
site-to-site spread of any species in the study. Raccoon concentration and coyote night activity and
concentration meet the same bar with smaller gains.</div>
</section>""")

    S.append("""
<section id="limits">
<h2><span class="num">8</span>Limitations</h2>
<ul>
<li>The analysis covers August to October. Schedules outside those months are untested.</li>
<li>The split between local and regional variation has a wide plausible range for 11 of 13 species, and
shifts by 0.21 on average between two definitions of a site, so it is reported as a range.</li>
<li>Habitat measurements are at 1 to 5 km resolution and explain a median of 6.9% of the site-to-site
variation that is not counting noise.</li>
<li>Variation between sites closer than 25 km cannot be resolved, because sites are defined by grouping
cameras within 25 km.</li>
<li>Four species have curves that do not flatten within the continent, so their correlation distance is
not estimated.</li>
<li>Overall activity level is excluded: the fitted curves do not reproduce raw values for that
measurement in 7 species.</li>
<li>Coyote time of peak has more measurement uncertainty than variation between sites.</li>
</ul>
</section>""")

    body = "".join(S)
    nav = ('<nav><a href="#summary">What we find</a><a href="#curves">Daily schedules</a>'
           '<a href="#predictions">Predictions</a><a href="#conditions">Conditions</a>'
           '<a href="#geography">Geography</a><a href="#data">Data</a>'
           '<a href="#mapping">Prediction</a><a href="#limits">Limitations</a></nav>')
    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Daily activity schedules of 13 North American species</title>
<style>{CSS}</style></head><body>
<header><h1>Daily activity schedules of 13 North American species</h1>
<p>686 camera sites, 1.22 million detections, August to October. Click any figure to enlarge.</p></header>
{nav}<main>{body}</main>
<footer>All figures are embedded in this file; nothing is loaded from the internet.</footer>
<div id="lb"><img id="lbimg" alt="enlarged figure"></div>
<script>{JS}</script></body></html>"""
    with open(outfile,"w") as f: f.write(doc)
    return outfile, len(doc)
