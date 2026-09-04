
def memo_rest(X, fg, vb, dfs, S, T, host):
    (vp, cv, rel, rsv, ess, f13, flex, EX, strong, MD, uac, uacon, H13, A, DN, CN, QC,
     aud3, dsi, excl) = dfs
    ts = strong[(strong.p<0.05)&(strong.p_spatial<0.05)].copy()
    ts["abs_b"] = ts.beta_std.abs(); ts = ts.sort_values("abs_b", ascending=False).head(20)
    CL = {"pop_1km":"human population","ag_5km":"agriculture","tmax_hottest_month":"summer heat",
          "tcc_1km":"tree cover","rug_5km":"terrain"}
    ML = {"pct_noct":"% at night","pct_crep":"% at dawn/dusk","conc":"concentration"}
    ts["Condition"] = ts.covariate.map(CL); ts["Measurement"] = ts.metric.map(ML)
    ts["Species"]=ts.species; ts["Sites"]=ts.n; ts["Effect"]=ts.beta_std
    ts["With spatial trend"]=ts.beta_std_spatial

    te = EX.copy(); te["Species"]=te.species; te["Sites"]=te.n
    te["Counting noise (%)"]=te.pct_is_measurement.round(0).astype(int)
    te["Rest explained by habitat (%)"]=te.pct_real_explained.round(0).astype(int)

    tr = rsv.copy(); tr["Species"]=tr.species; tr["Sites"]=tr.n
    tr["Regional share"]=tr.frac_regional; tr["Plausible range"]=[f"{a:.2f} to {b:.2f}" for a,b in zip(tr.lo,tr.hi)]
    tr["p"]=tr.p_short_structure; tr["Geographic structure"]=tr.regional_structure

    X.append(f"""
<section id="conditions">
<h2><span class="num">4</span>Environmental conditions and their effects</h2>
<p class="lead">Six conditions were fitted per site: human population density, agriculture, tree cover,
summer maximum temperature, terrain ruggedness and elevation. Agriculture and urbanisation are kept as
separate axes throughout; no composite human-footprint index was constructed at any stage.</p>

<div class="box"><h4>Two covariate decisions worth carrying into the methods</h4>
<b>Agriculture and urbanisation are near-orthogonal in this sample</b> (maximum absolute correlation 0.17
across all pairs), so a composite index would average two nearly independent predictors into a term
measuring neither.
<br><br><b>Nighttime lights were tested and rejected as the urbanisation axis.</b> Lights flip sign
between joint and single-predictor fits in 4 of 5 species while population density keeps its sign in all
5. An earlier finding that lights and population push activity in opposite directions within a species
is most likely a collinearity artifact and should not appear in the paper.</div>

<h3>The largest effects that survive a spatial trend</h3>
{_tab(ts, ["Species","Condition","Measurement","Sites","Effect","With spatial trend"])}
<p class="note">{S['n_strong']} of {S['n_tests']} relationships are significant alone and with a spatial
trend added, same sign in both. Effect is in standard deviations of the measurement per standard
deviation of the condition. Human population appears 15 times, terrain 14, agriculture 12, summer heat 9,
tree cover 6.</p>

<h3>How much site-to-site variation the conditions account for</h3>
{_tab(te, ["Species","Sites","Counting noise (%)","Rest explained by habitat (%)"])}
<p class="note">Counting noise is computed from the standard error on each site's fitted curve. The
second column is the reduction in residual variance when all six conditions are fitted.</p>
<p>Counting noise accounts for 7% to 47% of apparent variation, largest where detections are sparse. Of
the variation that is not counting noise, the six conditions explain a median of {S['ex_med']:.1f}%,
from 34% (black bear) to none (wild turkey, red squirrel). All conditions are measured at 1 to 5 km
resolution; whether finer-resolution conditions would explain more is untested, so the residual is
unaccounted for rather than shown to be unexplainable.</p>
</section>

<section id="scale">
<h2><span class="num">5</span>Spatial scale of variation</h2>

{fg("fig_variance_simple.png", "Site-to-site differences, three questions separately",
 "How much sites differ within each species with the counting-noise component marked; strength of "
 "evidence that nearby sites resemble each other; share of real variation explained by the conditions.")}

{fg("fig_variograms_grouped.png", "Whether similarity fades with distance",
 "Difference between pairs of sites as a function of separation, grouped into species where the curve "
 "rises with distance and species where it is flat.")}

<h3>Geographic structure verdict per species</h3>
{_tab(tr, ["Species","Sites","Regional share","Plausible range","p","Geographic structure"],
      hl=lambda r: r["regional_structure"]=="yes")}
<p class="note">Verdict is a permutation test on short-distance semivariance. The regional share and its
plausible range are reported for completeness, not used for the verdict.</p>

<div class="box bad"><h4>Do not report the local-versus-regional split as a single percentage</h4>
The plausible range on the regional share spans most of [0,1] for 11 of 13 species, and the point
estimate moves by 0.21 on average between two defensible definitions of a site &mdash; black bear reads
0.74 with annual surveys pooled and 0.07 with one year per site. Use the permutation-test verdict, which
is precise enough to state, and report the share as a range if at all.</div>

<div class="box warn"><h4>A circular-arithmetic bug was found and fixed during the final run</h4>
Time of day is circular: 23:50 and 00:10 are 20 minutes apart. Three code paths differenced them as
23 hours 40 minutes &mdash; the variogram pair differences, the permutation statistic, and the partition
bootstrap. A unit test gives 1.03 corrected against 124.10 buggy. Five other paths were audited and
found correct. All time-of-peak results were recomputed. After the fix, coyote time of peak has more
measurement variance than between-site variance, so it carries no usable signal.</div>
</section>""")

    tcv = cv[cv.beats_null==True].sort_values("skill", ascending=False).copy()
    tcv["Species"]=tcv.species
    tcv["Measurement"]=tcv.metric.map({"pct_noct":"% at night","pct_crep":"% at dawn/dusk",
                                       "conc":"concentration","peak_h":"time of peak"}).fillna(tcv.metric)
    tcv["Predicted from"]=tcv.model.map({"covariate":"local conditions","position":"map position"}).fillna(tcv.model)
    tcv["Accuracy gain"]=tcv.skill
    tcv["Plausible range"]=[f"{a:.2f} to {b:.2f}" for a,b in zip(tcv.lo,tcv.hi)]
    tes = ess[ess.metric=="pct_noct"].copy()
    tes["Species"]=tes.species; tes["Sites"]=tes.n; tes["Effective sites"]=tes.n_eff.round(0)
    tes["Interval width (points)"]=tes.best_ci95_width_pp
    tes["Mapped range (points)"]=tes.mapped_range_p5_p95
    tes["Width / range"]=tes.width_vs_range

    X.append(f"""
<section id="mapping">
<h2><span class="num">6</span>Prediction and mapping</h2>
<p class="lead">Whole 400 km regions are held out, the model is fitted on the remainder, and the held-out
sites are predicted. The comparison is against using that species' average schedule everywhere. The
median held-out site sits {S['med_dist']:.0f} km from the nearest fitting site.</p>

{fg("fig_s4_block_cv.png", "Prediction accuracy for unsampled regions",
 "Accuracy gain per species and measurement with plausible ranges, error against distance from the "
 "nearest fitting site, and the regional blocks used.")}

<h3>The {S['cv_beats']} of {S['cv_tot']} combinations that beat a nationwide schedule</h3>
{_tab(tcv, ["Species","Measurement","Predicted from","Accuracy gain","Plausible range"])}

{fg("fig_prediction_raster.png", "National maps for the two species that support one",
 "Predicted night activity at 25 km resolution for black bear and coyote, masked to the fitted condition "
 "range and the species' known range. Open circles are the camera sites.")}

<div class="box good"><h4>Two maps are defensible</h4>
<b>Black bear</b>, predicted from summer heat and terrain ruggedness only. Fitting all six conditions to
72 sites raised block-tested accuracy from 0.21 to 0.26 but produced a surface where neighbouring cells
differed by 9.6 points on average against a total between-site spread of 18.1 points, and predicted
values from &minus;13% to 202%. The two-condition version differs by 2.3 points between neighbours and
stays inside the observed range.
<br><br><b>Coyote</b>, predicted from all six conditions, with a narrow predicted range of 57% to 86%
night activity.</div>

<h3>Precision of any prediction for a new location</h3>
{_tab(tes, ["Species","Sites","Effective sites","Interval width (points)","Mapped range (points)","Width / range"])}
<p class="note">Effective sites accounts for spatial autocorrelation. Width is the 95% interval on a
prediction for a new location from the best model for that species.</p>

{fg("fig_s5_precision.png", "Prediction uncertainty against the variation being mapped",
 "Interval width for a new location alongside the range of variation across sampled sites.")}
</section>

<section id="methods">
<h2><span class="num">7</span>Measurement and methods detail</h2>

<h3>Activity curves</h3>
<p>Detection times are converted to a scale with sunrise at 0 and sunset at 12, then counted in 48
half-hour intervals. Each site's counts are described by a negative-binomial model with two harmonic
cycles, one per day and one twice per day, giving five parameters: an overall level and two pairs
describing the position and amplitude of each cycle. Standard errors come from the delta method. All
{S['n_curves']:,} fits converged.</p>
<p>Independence filter: a 30-minute running anchor per camera and species.</p>

<div class="box warn"><h4>Sun-anchored time breaks the effort offset, and this needs stating</h4>
Rescaling so sunrise and sunset fall at fixed anchors compresses day and night by different amounts, so
a camera does not record for equal real time in each interval. The within-camera ratio reaches 1.64 on
full-day deployments. In ordinary clock time the same measure is flat to 1.000002. Every count model
carries an exact analytic per-interval effort offset. Ignoring it shifts a site's night activity by up to
10.6 points and the bias reverses sign across the equinox, so it would appear as a latitudinal gradient.
This consequence is not spelled out in the source literature for the transformation and may warrant a
short methods note.</div>

{fg("fig_s1_effort_map_timeline.png", "Effort and the sun-time correction",
 "Effort per site, deployment timelines, and the distribution of within-camera effort imbalance in "
 "sun-anchored against clock time.")}

<h3>Which measurements are usable</h3>
{_tab(rel.assign(Species=rel.species, Measurement=rel.metric, Reliability=rel.reliability,
                 **{"Fit vs raw":rel.r_obs, "Ceiling":rel.ceiling, "Gap":rel.gap})[
      ["Species","Measurement","Reliability","Fit vs raw","Ceiling","Gap"]].head(20))}
<p class="note">Reliability is the share of between-site variance that is real rather than counting
noise. Fit-vs-raw is the correlation between the value from the fitted curve and the same value computed
directly from raw counts. Ceiling is the highest correlation attainable given detection counts.</p>
<p><b>Activity level is excluded from all spatial analysis.</b> Its fit-vs-raw correlation falls below
its own reliability ceiling in 7 species, with reliability of exactly 0 for turkey, bear and opossum.</p>

{fg("fig_s2_fitted_vs_raw.png", "Fitted curves against raw counts",
 "Values from the fitted curves plotted against the same values computed from raw detection counts.")}

{fg("fig_curves_13species.png", "Every activity curve",
 "One panel per species. Heavy line is the species schedule pooled across sites; faint lines are "
 "individual sites; shaded band is night.")}

{fg("range_mask_comparison_panels.png", "Range masks",
 "Range polygon, iNaturalist occurrence grid and their union per species.")}
<p>Range masks combine the range polygon, an iNaturalist occurrence grid at two or more observations per
100 km cell, and any cell with a confirmed camera detection. The union never performs worse than the
polygon alone: for moose the polygon and iNaturalist grid agree at Jaccard 0.39, with iNaturalist adding
34% more cells in a coherent Southern Rockies block, and the gray squirrel polygon omits an established
Pacific Northwest population containing 103 detection-bearing deployments. The camera-detection term
rescues 22 cells across 5 species.</p>
</section>""")
    return X
