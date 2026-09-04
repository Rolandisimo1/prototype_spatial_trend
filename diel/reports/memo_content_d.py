
def memo_close(X, dfs, S, T, FIGS, host):
    (vp, cv, rel, rsv, ess, f13, flex, EX, strong, MD, uac, uacon, H13, A, DN, CN, QC,
     aud3, dsi, excl) = dfs
    rows = []
    for nm, vid in sorted(T.items()):
        rows.append(dict(File=f'<a href="/artifacts/{vid}">{nm}</a>', Contents=DESC.get(nm, "")))
    tt = pd.DataFrame(rows)
    o = ["<table><thead><tr><th>File</th><th>Contents</th></tr></thead><tbody>"]
    for _, r in tt.iterrows():
        o.append(f"<tr><td>{r.File}</td><td style='text-align:left'>{_h.escape(r.Contents)}</td></tr>")
    o.append("</tbody></table>")
    files_tab = "".join(o)

    X.append(f"""
<section id="limits">
<h2><span class="num">8</span>Limitations</h2>
<ul>
<li>The analysis covers August through October. Schedules outside those months are untested, and the
window was chosen because survey effort outside it is too thin to separate geography from season.</li>
<li>The local-versus-regional variance split has a wide plausible range for 11 of 13 species and moves by
0.21 between two defensible site definitions.</li>
<li>Conditions are measured at 1 to 5 km resolution and explain a median of {S['ex_med']:.1f}% of the
site-to-site variation that is not counting noise.</li>
<li>Variation between sites closer than 25 km cannot be resolved, because sites are formed by clustering
cameras within 25 km.</li>
<li>Four species have variograms that do not flatten within the continent, so their correlation distance
and effective sample size are not estimated.</li>
<li>Activity level is excluded: fitted curves do not reproduce raw values for that measurement in 7
species.</li>
<li>Coyote time of peak has more measurement uncertainty than between-site variation.</li>
<li>Concentration, time of peak and activity breadth cannot support any claim involving detection rate,
following the simulated-null test.</li>
<li>All surviving predator effects are in the two species with the most sites, so the absence of effects
elsewhere may reflect statistical power rather than biology.</li>
<li>Density and predator effects add no out-of-sample accuracy (median change &minus;0.006), so they are
in-sample relationships rather than predictive ones.</li>
<li>Predator exposure from iNaturalist measures range overlap, not encounter rate.</li>
<li>Within-array temperature variation comes largely from different survey dates, so phenology cannot be
fully separated from temperature in the heat results.</li>
<li>Whether the 22 excluded cameras are the complete set of faulty ones cannot be established from the
data.</li>
</ul>

<div class="box warn"><h4>Claims retracted during the analysis, listed so they do not reappear</h4>
<ul>
<li><b>Lights and population push activity in opposite directions within a species.</b> Withdrawn: the
opposite-sign pattern is the signature of collinearity between two highly correlated predictors, and
lights flip sign between joint and single-predictor fits.</li>
<li><b>The local-versus-regional split reversed between analyses.</b> Withdrawn: that quantity is not
precisely estimable, so the apparent reversal was movement in an unidentified number.</li>
<li><b>Raccoon's activity map beats a nationwide schedule.</b> Withdrawn at the earlier stage when
bootstrapping over folds rather than events showed no combination with an interval excluding zero; the
current run's raccoon result is concentration predicted from map position.</li>
<li><b>Per-camera clock offsets fitted from activity patterns.</b> Withdrawn: searching 25 candidate
shifts improves fit for 44.5% of cameras known to be clean, against 60.7% of suspect ones, so the
procedure cannot identify individual cameras.</li>
</ul></div>
</section>

<section id="framing">
<h2><span class="num">9</span>Suggested framing for the manuscript</h2>
<ol>
<li><b>Lead with conservation of the daily schedule.</b> It holds across all 13 species, and the
activity-curve figure carries it directly.</li>
<li><b>Present the five predictions as a set, with the two failures reported as results.</b> H1 and H3
were built to predict opposite species rankings and neither ordering holds, which is a stronger
contribution than a list of confirmations.</li>
<li><b>Report the predator and density results with their null test.</b> The simulated-null test is what
makes the surviving effects credible and it disqualified three of five measurements, so it belongs in
the main text rather than the supplement.</li>
<li><b>Report the geographic-structure verdict per species, not a variance percentage.</b></li>
<li><b>Close with the two maps.</b> Black bear and coyote support a national surface; presenting them
after the scale result frames them as the exception the analysis identified.</li>
<li><b>Carry the sun-time effort correction as a methods contribution.</b> The bias it prevents would
have appeared as a latitudinal gradient.</li>
</ol>
</section>

<section id="files">
<h2><span class="num">10</span>Result files</h2>
<p class="lead">Every table behind this memo. Figures are embedded above and also saved individually.</p>
{files_tab}
</section>""")
    return X

DESC = {
 "all_hypothesis_results.csv":"H4 and H5 unified: effect, bootstrap interval, null-test verdict, artifact ratio, spatial shift, out-of-sample change",
 "hypothesis_tests_13sp.csv":"H1, H2, H3 across 13 species: 195 condition-by-measurement tests with and without a spatial trend",
 "density_null_test.csv":"H5 simulated-null test per species and measurement, with artifact ratios",
 "carnivore_null_test.csv":"H4 simulated-null test, 255 relationships, including effort-conditioned refits",
 "inat_taxonomic_qc.csv":"iNaturalist predator records: fetched, off-genus, captive and casual removed, clean cells",
 "midday_heat_13sp.csv":"H3 revised: heat effect on the share of activity around solar noon",
 "urban_vs_ag_concentration.csv":"H2 revised: urban against agriculture effect on concentration, paired per species",
 "urban_vs_ag_consistency.csv":"H2 second test: site-to-site spread within high-urban against high-farmland settings",
 "flexibility_13sp.csv":"Site-to-site spread per species, measurement-corrected, with distance from the day-night extremes and a logit-scale version",
 "explained_by_covariates.csv":"Share of apparent variation that is counting noise, and share of the rest explained by the six conditions",
 "strong_effects_13sp.csv":"All 56 condition effects significant alone and with a spatial trend",
 "variance_partition_rerun2.csv":"Variance split into counting noise, local and regional, per species and measurement, all four site-definition variants",
 "variance_partition_bootstrap_ci.csv":"Bootstrap intervals on the regional share",
 "variograms_rerun2.csv":"Semivariance by distance band, per species and measurement",
 "stage3_regional_structure_verdict.csv":"Permutation-test verdict on geographic structure per species",
 "stage3_circular_audit.csv":"Every code path touching time of peak, and whether it handled circularity correctly",
 "block_cv_summary_rerun2.csv":"All 104 block cross-validation tests with intervals and held-out distances",
 "effective_sample_size_rerun2.csv":"Effective sample size, interval width for a new location, and width against mapped range",
 "stage2_metric_reliability.csv":"Per-measurement reliability, fit-vs-raw correlation, and attainable ceiling",
 "prediction_raster.csv":"Predicted night activity per 25 km cell for black bear and coyote, with envelope and range flags",
 "species_final13.csv":"The 13 species with detections, usable sites and range bands",
 "array_audit_v3.csv":"All 686 sites with extent, deployment count and derivation",
 "data_source_inventory.csv":"60 contributing projects with citations and deployment counts",
 "excluded_deployments.csv":"The 22 excluded cameras with the evidence for each",
 "stage6_changes_vs_superseded.csv":"Every headline number against the previous analysis with attribution",
 "stage6_added8_vs_core5.csv":"Whether the 8 later-added species behave differently from the first 5",
}

def assemble_memo(X, outfile="results_memo.html"):
    nav = ('<nav><a href="#scope">Headline</a><a href="#data">Data</a><a href="#hypotheses">Predictions</a>'
           '<a href="#conditions">Conditions</a><a href="#scale">Spatial scale</a>'
           '<a href="#mapping">Mapping</a><a href="#methods">Methods</a><a href="#limits">Limitations</a>'
           '<a href="#framing">Framing</a><a href="#files">Files</a></nav>')
    doc = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Diel activity across the United States - consolidated results memo</title>
<style>{MEMO_CSS}</style></head><body>
<header><h1>Daily activity of 13 North American species: consolidated results</h1>
<p>Every analysis run on this project, with the test behind each result and the checks it survived.
Written to be handed to a manuscript-drafting session. Click any figure to enlarge.</p></header>
{nav}<main>{"".join(X)}</main>
<footer>Self-contained: all figures embedded, nothing loaded from the internet.</footer>
<div id="lb"><img id="lbimg" alt="enlarged figure"></div>
<script>{MEMO_JS}</script></body></html>"""
    with open(outfile,"w") as f: f.write(doc)
    return outfile, len(doc)
