
def memo_hypotheses(X, fg, vb, dfs, S):
    (vp, cv, rel, rsv, ess, f13, flex, EX, strong, MD, uac, uacon, H13, A, DN, CN, QC,
     aud3, dsi, excl) = dfs
    h1 = H13[(H13.covariate=="pop_1km") & (H13.metric=="pct_noct")].sort_values("body_kg")
    th = h1.copy(); th["Species"]=th.species; th["Body mass (kg)"]=th.body_kg; th["Sites"]=th.n
    th["Effect"]=th.beta_std; th["p"]=th.p; th["With spatial trend"]=th.beta_std_spatial
    th["p (spatial)"]=th.p_spatial

    X.append(f"""
<section id="hypotheses">
<h2><span class="num">3</span>The five a priori predictions</h2>
<p class="lead">All five were committed to before testing. Two predict opposite species rankings, so
they are jointly falsifiable rather than a list of plausible expectations.</p>

<h3>H1 &mdash; Human pressure shifts activity into the night, most in large-bodied and hunted species</h3>
<p>Mechanism: large mammals are hunted and more conspicuous, so avoiding daylight near people should pay
off more for them than for an arboreal squirrel. Prediction: deer and bear respond most, squirrel least.</p>
<p><b>Verdict: {vb('not supported')}</b> for the species ordering. Seven of 13 species show a
relationship, but the direction is inconsistent and does not track body mass (rank correlation
{S['r_h3'][0]:+.2f}, p = {S['r_h3'][1]:.2f}). Chipmunk, opossum and coyote become more nocturnal where
people are dense; cottontail, red fox and mule deer become less so. Red fox has the largest effect in
the study at 0.63 standard deviations, opposite to the prediction.</p>
{_tab(th, ["Species","Body mass (kg)","Sites","Effect","p","With spatial trend","p (spatial)"],
      hl=lambda r: r["p"]<0.05)}
<p class="note">Effect is the change in night activity per standard deviation of human population
density, in standard deviations of that species' night activity. Highlighted rows are p &lt; 0.05.</p>

<h3>H2 &mdash; Agriculture and urbanisation act differently in kind</h3>
<p>The original prediction was that agriculture raises concentration more than it shifts nocturnality,
on the reasoning that crops are a predictable resource. That reasoning was revised during review: urban
food (bins, compost, feeders) is available year-round at consistent times, while crops are seasonal and
removed at harvest, so <b>urban</b> should be the more concentrating axis.</p>
<p><b>Verdict: {vb('supported in direction')}</b>. Urban density concentrates activity more than
farmland in {S['n_urb']} of 13 species. The paired difference across species is not significant
(p = {S['p_urb']:.2f}), so this is a tendency. Raccoon fits most cleanly, with urban concentrating
(+0.20, p &lt; 0.001) and farmland loosening (&minus;0.15, p = 0.001) &mdash; opposite signs on the two
axes within one species, which is the contrast the prediction was built to detect.</p>
<p>A second test of the same mechanism: if urban food is more consistent, urban sites should be more
consistent in behaviour. Among 8 species with enough sites in both settings, site-to-site spread is
smaller in high-urban settings for {S['n_con']} of 8 (paired p = {S['p_con']:.2f}).</p>

<h3>H3 &mdash; Thermal constraint drives midday avoidance</h3>
<p>The original prediction was that small, poorly-thermoregulating species avoid midday most. That was
revised on review: heat production scales with mass while heat loss scales with surface area, so a large
body sheds heat poorly and midday heat should cost it more. Prediction as revised: bear responds most,
squirrel least.</p>
<p><b>Verdict: {vb('mechanism supported, ordering not')}</b>. Measuring the share of activity in the
three hours around solar noon, {S['n_midday_neg']} of 13 species reduce midday activity where summers
are hotter. Two survive a spatial trend: black bear (&minus;2.1 points per standard deviation of heat,
p &lt; 0.001) and gray squirrel (&minus;2.4 points, p &lt; 0.001). Those sit at opposite ends of the
size range, so the ordering fails (rank correlation {S['r_md'][0]:+.2f}, p = {S['r_md'][1]:.2f}).
Black bear also has the largest heat effect on overall night activity in the study (+0.61 standard
deviations, p &lt; 0.001, unchanged by a spatial trend).</p>

<div class="box"><h4>H1 and H3 were designed to conflict, and both orderings failed</h4>
H1 predicted deer and bear would respond most to human pressure; H3 as originally stated predicted
squirrel would respond most to heat. Neither body-mass ordering holds. What does order site-to-site
variation is how far a species' average already sits from fully diurnal or fully nocturnal: rank
correlation {S['r_bound'][0]:+.2f} (p &lt; 0.001), against {S['r_body'][0]:+.2f} for body mass. Fitting
both jointly leaves body mass at p = 0.73 and distance from the extremes at p = 0.003.
<br><br>This relationship is partly statistical: a percentage bounded at 0 and 100 has mechanically less
room to vary near its bounds. On a logit scale, which removes that constraint, the correlation with body
mass is {S['r_logit'][0]:+.2f} (p = {S['r_logit'][1]:.2f}) &mdash; so the bound explanation is the more
defensible reading.</div>

{fg("fig_hypotheses.png", "H1, H2 and H3",
 "Urban against farmland effect on concentration paired within species; change in midday activity share "
 "under heat against body mass; site-to-site spread against distance from the day-night extremes.")}

<h3>H4 &mdash; Predator presence suppresses and shifts prey activity, with super-additive richness</h3>
<p>Following Moll et al. Predator exposure was measured four ways: large-carnivore presence on the same
cameras, large-carnivore detection rate, iNaturalist range overlap per carnivore species, and a predator
richness count with a squared term for the super-additive prediction.</p>
<p><b>Verdict: {vb('narrowly supported')}</b>. {S['h4_surv']} of {S['h4_tot']} relationships survive
both a bootstrap over sites and a simulated-null test. Every one is in white-tailed deer or northern
raccoon, the two species with the most sites.</p>
<ul>
<li><b>Deer respond to large carnivores as predicted.</b> Cougar on camera (&minus;0.33 standard
deviations), any large carnivore on camera (&minus;0.26), cougar range (&minus;0.21) and any
large-carnivore range (&minus;0.21) all reduce deer night activity.</li>
<li><b>Raccoon responds to mesopredator competitors, not to large carnivores.</b> Coyote detection rate
and black bear detection rate both reduce raccoon night activity (&minus;0.08 and &minus;0.20).</li>
<li><b>The super-additive richness prediction is {vb('not supported')}.</b> The squared richness term is
significant in 0 of 25 species-by-measurement combinations. The linear richness term survives in 2.</li>
<li><b>Richness is not driven by the most-observed carnivore.</b> Dropping black bear from the richness
count preserves the sign in 22 of 25 combinations, so the term is not a black-bear proxy &mdash; but no
richness-excluding-bear term is itself significant.</li>
</ul>

<h3>H5 &mdash; Density dependence broadens activity distributions</h3>
<p>Prediction: at high detection rate, competition pushes individuals off-peak, reducing concentration.</p>
<p><b>Verdict: {vb('partly supported')}</b>, with a large caveat. {S['h5_surv']} of {S['h5_tot']}
relationships survive the null test and bootstrap. Deer night activity falls with its own detection rate
(&minus;0.36 standard deviations) and deer dawn-dusk activity rises (+0.42); gray squirrel dawn-dusk
activity rises (+0.43). Bear shows the largest raw effect on night activity (&minus;12.4 points per log
unit) and survives the null test.</p>

<div class="box bad"><h4>The null test disqualified three of five measurements outright</h4>
Detection rate is the denominator of both the predictor and the response: a site with more photographs
yields a better-resolved activity curve as well as a higher rate. The test simulates photograph counts
from a <b>fixed</b> activity shape at each site's observed rate and effort, then re-runs the analysis;
whatever the simulation reproduces is arithmetic.
<br><br>Reproduced share of the observed effect: dawn-dusk 0.03, night activity 0.11, concentration 0.60,
time of peak 0.97, <b>activity breadth 1.92</b> &mdash; the simulation produces a larger effect than the
data. Concentration, time of peak and activity breadth cannot carry a density claim, which removes the
concentration measurement the prediction was originally stated in.</div>

{fg("fig_density_predators.png", "H4 and H5 after the simulated-null test",
 "Share of each measurement's apparent density effect reproduced by simulation from a fixed activity "
 "shape; surviving density effects; surviving predator effects.")}

{fg("fig_hypothesis1_density.png", "Density dependence, full detail",
 "Every species-by-measurement density relationship with its null distribution.")}

{fg("fig_hypothesis2_carnivores.png", "Predator effects, full detail",
 "Predator richness, carnivore presence and carnivore detection-rate terms per species and measurement.")}

{fg("fig_carnivore_null_test.png", "Predator null test",
 "Observed against simulated effects for all 255 predator relationships.")}

<h3>Predator layer construction</h3>
<p>Camera detections of large carnivores were too sparse to test each species individually, so an
iNaturalist-derived richness layer was built, counting large carnivore species per 100 km cell.
Taxonomic screening removed captive records and off-genus matches: iNaturalist name search does fuzzy
matching and returns unrelated taxa.</p>
{_tab(QC)}
<p class="note">Casual-grade records were dropped, which removes most captive zoo animals. Domestic dog
is a subspecies of grey wolf and would contaminate wolf presence; no dog records survived the
genus and grade filters. Red wolf has 4 presence cells, too few to test.</p>
</section>""")
    return X
