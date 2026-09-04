
import html as _h, pandas as pd, numpy as np

def build_memo(FIGS, T, dfs, S, host, outfile="results_memo.html"):
    I = {k: _img64(v, host) for k, v in FIGS.items()}
    (vp, cv, rel, rsv, ess, f13, flex, EX, strong, MD, uac, uacon, H13, A, DN, CN, QC,
     aud3, dsi, excl) = dfs
    def fg(key, title, cap):
        if key not in I: return ""
        return (f'<figure><img src="{I[key]}" alt="{_h.escape(title)}">'
                f'<figcaption><b>{_h.escape(title)}.</b> {cap}</figcaption></figure>')
    def vb(v):
        c = {"yes":"v-yes","no":"v-no"}.get(str(v).lower(),"v-part")
        return f'<span class="verdict {c}">{_h.escape(str(v))}</span>'
    X = []

    X.append(f"""
<section id="scope">
<h2><span class="num">1</span>Scope and headline results</h2>
<p class="lead">This memo consolidates every analysis run on the continental camera-trap diel-activity
project, so a manuscript can be written from it without returning to the analysis session. Each result
carries the test that produced it and the checks it survived.</p>

<div class="grid">
 <div class="card"><div class="big">{S['n_species']}</div><div class="lbl">species</div></div>
 <div class="card"><div class="big">{S['n_arrays']}</div><div class="lbl">camera sites (arrays)</div></div>
 <div class="card"><div class="big">1.22M</div><div class="lbl">independent detections</div></div>
 <div class="card"><div class="big">{S['n_curves']:,}</div><div class="lbl">site-level activity curves</div></div>
</div>

<h3>Five results, in the order a paper would present them</h3>
<ol>
<li><b>Daily activity is largely a species-level trait.</b> Site-to-site spread within a species is
small relative to differences between species, and this holds across all 13 species.</li>
<li><b>Site-to-site differences are real and connect to measurable conditions.</b> {S['n_strong']} of
{S['n_tests']} species-by-condition-by-measurement relationships are significant both alone and with a
spatial trend added, same sign in both.</li>
<li><b>Seven of 13 species show large-scale geographic structure</b> in night activity by a permutation
test; six do not.</li>
<li><b>Predator and conspecific-density effects exist but are narrow.</b> {S['h4_surv']} of
{S['h4_tot']} predator relationships and {S['h5_surv']} of {S['h5_tot']} density relationships survive
a simulated-null test; all surviving predator effects are in white-tailed deer and raccoon.</li>
<li><b>Two species support a national map.</b> Black bear and coyote night activity predict
held-out regions better than a single nationwide schedule.</li>
</ol>

<div class="box bad"><h4>Constraint that applies to every mapped claim</h4>
The 95% interval on a prediction for an unsampled location is
{S['width_ratio'][0]:.2f} to {S['width_ratio'][1]:.2f} times the full range of variation across
sampled sites (median {S['width_ratio'][2]:.2f}). It never falls below 1. Any map figure needs its
uncertainty shown alongside it.</div>
</section>""")

    tf = f13.copy()
    tf["Species"] = tf.species; tf["Detections"] = tf.det; tf["Sites"] = tf.arrays25
    tf["Range band (lon)"] = [f"{a:.0f} to {b:.0f}" for a,b in zip(tf.range_lon_min, tf.range_lon_max)]
    X.append(f"""
<section id="data">
<h2><span class="num">2</span>Data</h2>
<p class="lead">Camera-trap records contributed by 60 projects, restricted to August through October and
the lower 48 states.</p>

{_tab(tf, ["Species","Detections","Sites","Range band (lon)"])}
<p class="note">Sites are camera arrays with at least 25 independent detections of that species. Range
bands exclude records outside each species' distribution; a Douglas's squirrel record appeared at
longitude -71.7 before screening.</p>

{fg("fig_sampling.png", "Sampling coverage",
 "Camera and site locations, sites per array, effort by month with the analysis window marked, "
 "detections and usable sites per species, and each species' position on the day-night axis.")}

{fg("camera_trap_diagnostics.png", "Data-quality checks",
 "Camera locations, deployment timelines by project, monthly effort, and gaps between consecutive "
 "deployments at one location.")}

<h3>Unit of analysis</h3>
<p>A site is a cluster of cameras within 25 km, formed by one global complete-linkage clustering of all
deployments at once. Maximum extent is 24.92 km and minimum separation between sites is 8.48 km.
299 of {S['n_arrays']} sites pool cameras from more than one project, which pools repeated annual
surveys of the same physical location so one location is not counted as several independent sites.</p>
<div class="box warn"><h4>Methods caveat to state explicitly</h4>
Pooling repeated annual surveys is correct when sites are needed as independent replicates, which is
the case here. It would be wrong for a question about change between years, where year-within-site is
the contrast of interest. State the choice and its scope in the methods.</div>

{fg("array_rebuild_map.png", "Site definition",
 "Site locations and the distribution of site extents under the definition used.")}

<h3>Timestamps</h3>
<p>Sunrise and sunset are computed in local clock time from each site's time zone in the IANA database,
including daylight saving; zone boundaries follow state lines and Arizona does not observe daylight
saving, so the zone is looked up by location. The calculation agrees with published sunrise and sunset
times to within 14 minutes at nine cities in four time zones.</p>
<p>Storage was verified per project using species with fixed schedules. Dataset-wide values are
chipmunk 4.9%, gray squirrel 5.7%, raccoon 91.2% and opossum 98.3% of activity at night, against
expectations near 5%, 5%, 90% and 95%. One project stored universal time and was shifted by five hours;
eleven projects were replaced with their original published archive; one was dropped as unrecoverable.
Twenty-two cameras with faulty clocks were excluded, two confirmed from photographs. No camera
timestamps were shifted to fit an expectation.</p>

{fg("november_assessment.png", "Why the window ends in October",
 "Monthly effort, the behaviour of strictly daytime species around the daylight-saving change, and the "
 "number of additional sites November would contribute.")}
<p>November was tested and excluded: strictly daytime species appear to become night-active immediately
after the clock change but not before it, day length accounts for about 3% of the jump, and November
would add one new site location.</p>
</section>""")
    return X, I, fg, vb
