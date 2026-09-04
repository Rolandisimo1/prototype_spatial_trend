
import html as _h, pandas as pd, numpy as np

def build_report(figs, cmp_vp, cv, rel, host, outfile="analysis_report.html"):
    IMG = {k: _img(v[0], host) for k, v in figs.items()}
    cvb = cv[cv.beats_null == True].sort_values("skill", ascending=False)

    def fig(key, title, caption, plain):
        if key not in IMG:
            return ""
        return (f'<figure><img src="{IMG[key]}" alt="{_h.escape(title)}">'
                f'<figcaption><b>{_h.escape(title)}.</b> {caption}</figcaption></figure>'
                f'<div class="plain"><b>In plain terms:</b> {plain}</div>')

    S = []

    # ---------- 1. WHAT WE SET OUT TO DO
    S.append(f"""
<section id="goal">
<h2><span class="num">1</span>What we set out to do, and what we found</h2>
<p class="lead">We asked whether the daily activity schedule of a mammal species &mdash; when it is
awake and moving through the 24-hour cycle &mdash; changes from place to place across the United
States, and whether we can predict that change well enough to draw it on a map.</p>

<p>The short answer is that a species' daily schedule is mostly a property of <em>the species</em>,
not of the place. A raccoon behaves like a raccoon almost everywhere. There is real local variation
on top of that, but it is organised at a scale finer than a continental camera network can resolve,
so it does not turn into a usable national map for most species.</p>

<div class="grid">
  <div class="card"><div class="big">13</div><div class="lbl">species analysed</div></div>
  <div class="card"><div class="big">686</div><div class="lbl">camera arrays</div></div>
  <div class="card"><div class="big">1.22M</div><div class="lbl">independent detections</div></div>
  <div class="card"><div class="big">11 / 104</div><div class="lbl">tests where a map beat a single
     nationwide schedule</div></div>
</div>

<div class="callout good"><h4>The main result, stated without jargon</h4>
If you want to guess what time of day a deer will be active at a place you have never sampled, you
do about as well by using the national average deer schedule as by using anything we can build from
that place's location or its habitat. That holds for 93 of the 104 combinations of species and
measurement we tested. It is a real finding, not a failure to find one: it says daily activity is a
conserved trait of the species, and that local flexibility exists but operates below the scale we
can currently see.</div>

<h3>Why this is worth publishing as a positive result</h3>
<p>Two things make it more than an absence of signal. First, the same pattern holds across all 13
species &mdash; including 8 we added specifically to test whether the first 5 were unusual. They were
not. Second, we can say <em>how</em> coarse the limit is and <em>which</em> species are exceptions:
7 of 13 species do show genuine large-scale geographic structure in their night-time activity, and
black bear and raccoon show it strongly enough to predict.</p>
</section>""")

    # ---------- 2. THE DATA
    S.append(f"""
<section id="data">
<h2><span class="num">2</span>The data, and how we checked it</h2>
<p class="lead">Everything below rests on camera-trap photographs contributed by dozens of separate
projects. Combining them is where most of the risk lies, so this section documents what we have and
what was wrong with it.</p>

{fig("atlas1_sampling", "Where the cameras are, when they ran, and which species they caught",
 "Panel a maps every camera (grey) and the arrays they group into (coloured bubbles; size = number of "
 "cameras, colour = how far apart they are). Panel b shows array size. Panel c shows survey effort by "
 "month, with our analysis window in red. Panels d and e show detections and usable arrays per species. "
 "Panel f places each species on the day&ndash;night axis.",
 "Cameras are not spread evenly. There is a dense cluster in the Carolinas and thin coverage across the "
 "arid West, which limits how far any prediction can travel. Almost all the photographs come from "
 "August to October, which is why we restricted the analysis to those three months &mdash; otherwise "
 "'Montana looks different from Georgia' could just mean 'Montana was sampled in a different season'.")}

{fig("camera_trap_diagnostics", "Routine data-quality plots we now make before any modelling",
 "Deployment locations, deployment timelines by project, monthly effort, and the gap between "
 "consecutive deployments at the same location.",
 "These four plots exist because looking at the data caught errors that summary tables hid. The "
 "timeline in particular revealed projects whose records span years when they were supposed to cover "
 "a single season.")}

{fig("array_rebuild_map", "Rebuilding the unit of analysis",
 "Before and after the array definition was rebuilt, with the size distribution alongside.",
 "An 'array' is a cluster of nearby cameras treated as one sampling unit. The original grouping had a "
 "serious flaw: when it could not find a spatial grouping it fell back on lumping an entire project "
 "together &mdash; in one case cameras spread over 4,000 km were treated as a single site. Since "
 "almost every statistic in the analysis depends on what counts as an independent site, this had to be "
 "fixed before anything else. The rebuilt version groups cameras within 25 km of each other and "
 "nothing wider.")}

<h3>Problems we found and fixed in the source data</h3>
<ol>
<li><b>Timestamps in the wrong time zone.</b> Some projects stored photograph times in universal time
rather than local time, which makes a squirrel look nocturnal. We caught this by using strictly
daytime species as clocks: if squirrels appear to be active at midnight, the clock is wrong, not the
squirrel. Rather than patch the affected data we replaced 11 projects with the original published
archive, which had correct times.</li>
<li><b>Individual broken camera clocks.</b> 22 cameras were dropped. We deliberately did <em>not</em>
shift their times to fit expectations &mdash; doing so would erase exactly the unusual local behaviour
the study is trying to measure. Two of the 22 were confirmed by eye from the photographs.</li>
<li><b>A bug in our own sunrise calculation.</b> We were comparing solar sunrise times against clock
timestamps, which is wrong by up to an hour. That made ordinary dusk activity at sites on the western
edge of a time zone look like night-time activity, and generated false 'broken camera' reports. Fixed
and checked against published sunrise tables at nine cities.</li>
<li><b>Repeated surveys counted as independent sites.</b> Annual surveys revisit the same cameras year
after year. Treating each year as a new site inflates the apparent sample size. The rebuilt arrays
pool them &mdash; correct here, though it would be the wrong choice for a study asking how behaviour
changed between years.</li>
</ol>
</section>""")

    # ---------- 3. SEASON
    S.append(f"""
<section id="season">
<h2><span class="num">3</span>Why the analysis covers only August to October</h2>

{fig("november_assessment", "Testing whether to extend the window into November",
 "Panel a shows monthly effort with November highlighted. Panel b shows what happens to strictly "
 "daytime species around the clock change. Panel c shows how many extra arrays November would add.",
 "November has real data, so it was worth asking. But the clocks change in early November, and about "
 "a third of cameras do not adjust automatically. Panel b is the evidence: strictly daytime species "
 "&mdash; squirrels, chipmunks, turkeys &mdash; suddenly appear to become night-active, but only "
 "<em>after</em> the clock change, not before it. We checked the obvious alternative explanation, that "
 "shorter days simply push animals into darkness, and it accounts for only about 3% of the jump. "
 "Meanwhile November would add only one genuinely new location. Not worth the contamination.")}
</section>""")

    # ---------- 4. MEASURING
    S.append(f"""
<section id="measure">
<h2><span class="num">4</span>How we measure a daily schedule</h2>
<p class="lead">For each species at each array we fit a smooth curve describing activity through the
24-hour cycle, then summarise that curve with a few simple numbers.</p>

<p>Two details matter. First, we measure time relative to <b>sunrise and sunset</b> rather than by the
clock, because 7 a.m. means something different in Maine in October than in Arizona in August.
Second, stretching the day this way has a consequence that is easy to miss: a camera does not collect
equal amounts of data in each stretched time slot. We compute that unevenness exactly and correct for
it. Ignoring it shifts a site's night-time activity estimate by up to 10 percentage points, and the
error changes direction across the equinox &mdash; so it would have looked like a real
north&ndash;south pattern.</p>

{fig("fig_s1_effort_map_timeline", "Survey effort, and the unevenness the sun-based clock introduces",
 "Effort per array on the map, the deployment timeline, and the distribution of within-camera effort "
 "imbalance in sun-based time versus ordinary clock time.",
 "In ordinary clock time a camera records equally in every time slot, so the flat control line is "
 "exactly flat. In sun-based time it does not, and the imbalance reaches 1.6-fold at some cameras. "
 "This is the correction described above.")}

{fig("fig_s2_curves_13species", "The fitted daily schedule of each species, and of each array",
 "One panel per species. The heavy line is the species-wide schedule; each faint line is one array.",
 "This is the raw material of the whole study. The species signatures are unmistakable and match "
 "what is known from field biology: chipmunks and turkeys active in the middle of the day, opossums "
 "almost entirely at night, deer concentrated at dawn and dusk. What matters for our question is the "
 "faint lines: they scatter around the heavy line, but not by much, and not in an obviously "
 "geographic way.")}

{fig("fig_s2_fitted_vs_raw", "Checking that the smooth curves are faithful to the raw photographs",
 "Fitted values plotted against values computed directly from the raw detection histogram.",
 "A smooth curve can lie. This check compares each summary number against the same number computed "
 "straight from the photograph counts with no smoothing. Three of our four measures track the raw "
 "data closely. The fourth, 'activity level', does not: for 7 species it performs worse than its own "
 "noise floor allows, so we excluded it from all geographic analysis rather than report it.")}
</section>""")

    # ---------- 5. HOW MUCH VARIES, AND AT WHAT SCALE
    tvp = cmp_vp.copy()
    tvp["Species"] = tvp.species
    tvp["Measurement noise"] = tvp.frac_meas
    tvp["Local (under 25 km)"] = tvp.frac_local
    tvp["Regional"] = tvp.frac_regional
    tvp["Regional, plausible range"] = [f"{lo:.2f} to {hi:.2f}" for lo, hi in zip(tvp.frac_reg_lo, tvp.frac_reg_hi)]
    tvp["If years not pooled"] = tvp.frac_reg_oneyear
    tbl1 = _table(tvp, ["Species","Measurement noise","Local (under 25 km)","Regional",
                        "Regional, plausible range","If years not pooled"],
                  highlight=lambda r: r["species"] in ("Northern Raccoon","American Black Bear"))

    S.append(f"""
<section id="scale">
<h2><span class="num">5</span>How much does activity vary from place to place, and at what scale?</h2>
<p class="lead">We split the differences between arrays into three parts: noise from having limited
photographs, genuine differences between nearby sites, and genuine differences organised across
large regions. Only the third kind can be mapped.</p>

{fig("fig_s3_variance_partition", "Splitting the variation three ways",
 "Bars show the share of between-array variation attributable to measurement noise (red), local "
 "differences under 25 km (orange), and regionally structured differences (blue). Black bars are the "
 "plausible range on the regional share. The right panel compares the estimate under two defensible "
 "definitions of a site.",
 "The orange band is the important one: for most species, most of the real variation is between "
 "<em>nearby</em> sites, not between regions. Two neighbouring forests can differ more than two "
 "distant states. That is why a national map does not work &mdash; not because animals are identical "
 "everywhere, but because their differences are organised at the wrong scale for a map.")}

<div class="callout bad"><h4>An honest limitation: we cannot put a single number on this split</h4>
We would like to say something like "72% of the variation is local". We cannot. The plausible range
on the regional share covers most of the possible span for 11 of the 13 species, and the estimate
moves by 0.21 on average simply by changing which of two reasonable definitions of a "site" we use
&mdash; black bear reads 0.74 one way and 0.07 the other. Any single percentage would be an artefact
of an arbitrary choice. We report the ranges instead, and we lean on a different test (below) that
<em>is</em> precise enough to act on.</div>

{tbl1}
<p class="note">Highlighted rows are the two species whose regional share is distinguishable from zero.
"If years not pooled" shows the same quantity when repeated annual surveys are treated separately;
the gap between the two columns is the instability described above.</p>

{fig("fig_s3_variograms", "How quickly similarity fades with distance",
 "For each species, how different two arrays' schedules are as a function of how far apart they are.",
 "If geography mattered strongly, these lines would rise steeply and then flatten: nearby sites "
 "similar, distant sites different. For some species they do. For others the line is essentially flat, "
 "meaning two cameras 20 km apart differ about as much as two cameras 2,000 km apart. A flat line is a "
 "direct statement that there is no large-scale geographic pattern to map.")}

{fig("fig_s3_metric_maps", "What the variation looks like drawn on a map",
 "Each array's night-time activity share, coloured relative to that species' own average.",
 "These maps are included precisely because they look unconvincing, and that is the honest result. "
 "For most species the colours are scattered rather than forming coherent regions. Raccoon and black "
 "bear are the exceptions where visible regional patterning appears.")}

<h3>Which species genuinely do have large-scale geographic structure</h3>
<p>Because the percentage split is imprecise, we tested the question directly instead: is the
similarity between nearby arrays greater than you would get by shuffling the arrays at random? That
test is precise enough to give a clean answer. <b>7 of 13 species pass</b>: wild turkey, red squirrel,
black bear, white-tailed deer, red fox, eastern cottontail and northern raccoon. Mule deer is
borderline. Five do not: eastern chipmunk, gray squirrel, fox squirrel, coyote and opossum.</p>
</section>""")

    # ---------- 6. CAN WE PREDICT
    tcv = cvb.copy()
    tcv["Species"] = tcv.species
    tcv["Measurement"] = tcv.metric.map({"pct_noct":"% active at night","pct_crep":"% at dawn/dusk",
                                         "conc":"how concentrated","peak_h":"time of peak"}).fillna(tcv.metric)
    tcv["Predicted from"] = tcv.model.map({"covariate":"local habitat","position":"location alone"}).fillna(tcv.model)
    tcv["Skill"] = tcv.skill
    tcv["Plausible range"] = [f"{a:.2f} to {b:.2f}" for a, b in zip(tcv.lo, tcv.hi)]
    tbl2 = _table(tcv, ["Species","Measurement","Predicted from","Skill","Plausible range"])

    S.append(f"""
<section id="predict">
<h2><span class="num">6</span>Can we predict a place we have never sampled?</h2>
<p class="lead">This is the test that decides whether a map is publishable. We hide entire 400 km
regions, fit the model on what remains, and ask whether it predicts the hidden region better than
simply using that species' nationwide average schedule.</p>

<p>Hiding whole regions rather than individual sites matters. Camera arrays sit close together, so
hiding one array leaves its neighbours in the training data and the test is far too easy. The median
hidden array in our test sits <b>{cv.median_heldout_dist_km.median():.0f} km</b> from the nearest
training data.</p>

{fig("fig_s4_block_cv", "Does a map beat a single nationwide schedule?",
 "Predictive skill for every species and measurement, with plausible ranges. Zero means 'no better "
 "than one nationwide schedule for that species'. Also shown: error against distance from training "
 "data, and the regional blocks used.",
 "Almost everything sits at or below zero. Only 11 of 104 combinations clearly beat the nationwide "
 "average &mdash; and those 11 are the honest contribution: they identify exactly where a "
 "place-specific prediction earns its keep. Black bear is the standout, and its activity is "
 "predictable from local habitat rather than from location on the map.")}

{tbl2}
<p class="note">These are the 11 combinations that beat a single nationwide schedule. "Local habitat"
means predicted from conditions at the site (human population, farmland, tree cover, temperature,
terrain); "location alone" means predicted from map position only.</p>

{fig("fig_s5_precision", "How precisely could a map ever be drawn?",
 "The honest uncertainty on a prediction for a new location, compared against how much variation "
 "the map would actually be showing.",
 "This is the quiet killer. For every species and every measurement, the uncertainty on a prediction "
 "for a new place is <em>larger</em> than the whole range of variation the map would display "
 "(ratio 1.00 to 1.60). In other words, even where a pattern is real, the error bars are wider than "
 "the pattern. A map drawn without those error bars would look far more confident than the data "
 "supports.")}
</section>""")

    # ---------- 7. CHANGES
    S.append(f"""
<section id="changes">
<h2><span class="num">7</span>What changed when we rebuilt the analysis, and why</h2>

{fig("fig_s6_old_vs_new", "Before and after the rebuild",
 "Headline numbers from the previous version of the analysis against the current one.",
 "The analysis was rerun from scratch after the data problems in section 2 were fixed and 8 species "
 "were added. This panel shows what moved.")}

<div class="callout good"><h4>The main conclusion survived a much harder test</h4>
The previous version found 2 of 25 combinations beating the nationwide average; this one finds 11 of
104. The <em>proportion</em> barely moved (8% to 11%) on four times as many tests, with better data
and eight more species. That is genuine confirmation rather than a lucky repeat.</div>

<div class="callout"><h4>A previously reported result was withdrawn</h4>
An earlier version reported that the local-versus-regional split had reversed between analyses. We
now know that quantity is not precisely estimable at all, so the "reversal" was movement in a number
that was never pinned down. It has been withdrawn rather than reinterpreted.</div>

<h3>A bug found during this run</h3>
<p>Time of day is circular: 11:50 p.m. and 12:10 a.m. are twenty minutes apart, not twenty-three hours
and forty minutes. One of our calculations was treating them the wrong way, which inflated apparent
differences wherever a species' peak sat near midnight. Three separate places in the code carried the
error; five others were checked and found correct. All time-of-peak results were recomputed. One
consequence: after the fix, coyote time-of-peak has more noise than signal, so it carries no usable
information between arrays.</p>

<h3>Things we cannot verify, stated plainly</h3>
<ul>
<li>The true size of the local-versus-regional split, for the reasons above.</li>
<li>The distance over which similarity fades, for four species whose curves never flatten within the
continent.</li>
<li>Whether the 22 excluded cameras are the complete set of faulty ones. We can only confirm the ones
we found.</li>
<li>The structure of variation below 25 km, which our sampling design cannot resolve at all &mdash;
we can only see that it is there.</li>
</ul>
</section>""")

    # ---------- 8. WHERE NEXT
    S.append("""
<section id="next">
<h2><span class="num">8</span>Where this leaves the paper</h2>
<p class="lead">The result to lead with is the scale finding, supported by the species-by-species
verdict on geographic structure &mdash; not a map, and not a single percentage.</p>
<ol>
<li><b>Lead with the conserved-trait result.</b> Daily activity is largely a species-level property;
this holds across all 13 species including 8 added as a check.</li>
<li><b>Report which species are exceptions and how we know.</b> 7 of 13 show genuine large-scale
structure by a test that is precise enough to state; black bear and raccoon strongly enough to
predict.</li>
<li><b>Report the precision honestly.</b> The uncertainty for a new location exceeds the mapped
variation for every species. This is the sentence that keeps the paper defensible.</li>
<li><b>Do not report a variance percentage.</b> It is not identified, and a reviewer would find that.</li>
<li><b>Treat the local component as the open question.</b> Something real is happening below 25 km.
Answering it needs a different design &mdash; paired nearby sites, not a wider net.</li>
</ol>
</section>""")

    body = "".join(S)
    nav = ('<nav><a href="#goal">Overview</a><a href="#data">The data</a><a href="#season">Season</a>'
           '<a href="#measure">Measuring activity</a><a href="#scale">Scale of variation</a>'
           '<a href="#predict">Prediction</a><a href="#changes">What changed</a>'
           '<a href="#next">Where next</a></nav>')
    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Mammal daily activity across the United States - analysis report</title>
<style>{CSS}</style></head><body>
<header><h1>How mammal daily activity varies across the United States</h1>
<p>Analysis report &mdash; 13 species, 686 camera arrays, 1.22 million detections. Every figure is
clickable to enlarge.</p></header>
{nav}<main>{body}</main>
<footer>Self-contained file: all figures are embedded, nothing is loaded from the internet, and it can
be emailed or posted anywhere. Figures are also saved individually as artifacts.</footer>
<div id="lb"><img id="lbimg" alt="enlarged figure"></div>
<script>{JS}</script></body></html>"""
    with open(outfile, "w") as f:
        f.write(doc)
    return outfile, len(doc)
