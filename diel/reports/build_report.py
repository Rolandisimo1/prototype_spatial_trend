"""Build the self-contained HTML results report with every figure embedded as base64."""
import base64, json, io
import numpy as np
import pandas as pd

FIGS = [
    ('fig1_sampling_coverage.png', 'Figure 1. Sampling coverage and effort',
     'Panel a places every site on the map, with dot size showing how many cameras it pools and dot colour '
     'showing total survey effort on a logarithmic scale. Panel b gives effort by survey year. Panel c gives, '
     'for each camera, the ratio of its busiest to its quietest half-hour bin in sun-anchored time.',
     'Sampling is heavily concentrated in the east and in recent years, so the dataset describes the eastern '
     'United States far better than the west. Panel c is the reason every model here carries an effort offset: '
     'if you convert clock time to sun time, a camera does not spend equal time in every bin. The median camera '
     'has 10 percent more effort in its busiest bin than its quietest, and some have three times more. Ignoring '
     'that would create a fake pattern that changes sign across the equinox.'),
    ('fig2_species_curves.png', 'Figure 2. Diel activity of all 13 species',
     'Each panel is one species. Thin grey lines are the fitted curves of individual sites; the thick coloured '
     'line is the species mean across sites. Time runs from sunrise through mid-day, sunset and mid-night back '
     'to sunrise, so the shaded band is night everywhere regardless of season or latitude. Curves are scaled to '
     'sum to 100 percent across the 48 half-hour bins. Colour marks the guild: diurnal is under 20 percent night, '
     'nocturnal is over 80 percent night.',
     'Species sort cleanly into daytime animals, night-time animals and a middle group. Within a species the '
     'individual site curves cluster tightly around the mean, which is the central finding of the whole analysis: '
     'a species keeps roughly the same daily rhythm across the country. The visible spread between sites is '
     'modest compared with the difference between species.'),
    ('fig3_variograms.png', 'Figure 3. Does similarity fade with distance',
     'Each line is one species. The horizontal axis is the distance between two sites on a logarithmic scale, '
     'and the vertical axis is how different their percent-nocturnal values are, scaled so that 1 equals the '
     'difference between two sites picked at random. A rising line means nearby sites resemble each other more '
     'than distant ones. Species are split by whether a permutation test found that pattern.',
     'Four of thirteen species show the classic pattern where nearby places are more alike. Nine do not: their '
     'lines are flat, meaning two sites 25 km apart differ about as much as two sites 2000 km apart. For those '
     'nine species there is no smooth geographic surface to map, because the variation is not organised by '
     'geography at all.'),
    ('fig4_variance_partition.png', 'Figure 4. Where the site-to-site variation comes from',
     'Three components of the between-site variance in percent nocturnal, one per panel: counting error, real '
     'differences at scales under 25 km, and real differences organised above 25 km. Filled circles use all '
     'sites with survey years pooled; open diamonds use one year per site. Species where a component cannot be '
     'estimated are marked.',
     'The three components are drawn separately because they are not a stacked whole that can be read off one '
     'bar. The gap between the circles and the diamonds in the third panel is the important feature: the '
     'regionally structured share moves substantially depending on whether repeated annual surveys of the same '
     'place are pooled. That share is not identified to better than a range, and the range spans 14 to 35 '
     'percent across the two site definitions and four choices of the local-versus-regional cutoff.'),
    ('fig5_block_cv_skill.png', 'Figure 5. Can the surfaces predict where you have not sampled',
     'For each species and each measure, the best of two competing models is shown: a smooth function of '
     'location, or a function of site covariates. The horizontal axis is predictive skill at held-out 400 km '
     'spatial blocks, measured against a null model of one nationwide value per species. Zero means no better '
     'than that single number. Bars are 95 percent intervals from resampling whole spatial blocks. The panel '
     'shows the better of the two models for each species and measure; the full set of 144 fits is in '
     'stage4_block_cv_summary.csv.',
     'Six of 144 species-by-measure-by-model fits predict held-out regions better than one national number. '
     'The median fit is slightly worse than that null. Coyote accounts for four of the six. This is the test '
     'that decides whether a national map is meaningful, and for most species and measures the answer is that '
     'a single number does the job as well as any surface.'),
    ('fig6_hypothesis_tests.png', 'Figure 6. The five a-priori predictions',
     'Panel a: effect of human population on percent nocturnal, by body mass, with 95 percent bootstrap '
     'intervals. Panel b: the urban effect on activity concentration minus the farmland effect. Panel c: effect '
     'of summer heat on the share of activity near solar noon, for the eight species with midday activity to '
     'lose. Panels d and e: for every effect that was significant, the share of it reproduced by a simulated '
     'null in which curve shape is held fixed at the national mean and only the site detection rate varies. '
     'Panel f: site-to-site spread against how far a species mean sits from the fully diurnal or fully '
     'nocturnal bound, on the raw percentage scale and on a bound-free logit scale.',
     'None of the five predictions holds in the general form it was posed. Human pressure does push some '
     'species into the night, notably black bear, coyote and raccoon, but red fox and gray squirrel move the '
     'other way, and the predicted ordering by body mass has an interval that includes zero. Urban and '
     'farmland effects differ in the predicted direction in seven of thirteen species, which is what a coin '
     'flip gives. Heat costs midday activity in two small species and in none of the large ones, which is the '
     'opposite of the heat-shedding prediction. For predators and for detection rate, the fixed-shape null '
     'reproduces most of the apparent effect: 74 and 82 percent of the median effect respectively. Panel f '
     'closes an older thread: the strong relationship between variability and body mass on the raw scale is an '
     'artefact of a bounded measure, and it disappears entirely on the logit scale.'),
    ('fig7_prediction_rasters.png', 'Figure 7. The only surfaces that earned a map',
     'Left column: predicted coyote night-time share and midday share across 25 km grid cells, masked to both '
     'the range of covariate values actually sampled and the species range. White dots in the top map are the '
     'fitted sites. Right column: the distribution of predicted values across the map, with the 95 percent '
     'interval for a single grid cell drawn to the same scale.',
     'These two surfaces are the only ones in the analysis that beat a single national number out of sample '
     'when the covariates are measured identically at the fitting sites and at the prediction cells. The east '
     'is predicted to be more nocturnal, the interior west more midday-active. The right column carries the '
     'main caveat: the honest uncertainty for any one place is about 50 percentage points wide for night-time '
     'share, while the whole map spans 27 points. The pattern is real in aggregate, and a prediction for one '
     'particular new location is not usable.'),
    ('fig8_example_curves.png', 'Figure 8. Fitted curves at contrasting locations',
     'Coyote activity at the two most and two least densely populated sites with at least 150 detections. '
     'Thin step lines are raw counts divided by the effort in each bin; thick lines are the fitted harmonic '
     'curves. Both low-population sites round to under 0.1 people per square kilometre.',
     'Where people are dense the coyote curve is almost flat through the night with a deep midday trough. '
     'Where people are sparse the same species keeps a substantial amount of daytime and dusk activity. The '
     'gap is about 18 percentage points of night-time share, which is the H1 effect made concrete at four '
     'real places.'),
    ('fig9_exclusion_effect.png', 'Figure 9. What removing the clock errors did',
     'Panel a plots, for every curve whose percent nocturnal moved more than half a point, the size of the '
     'move against its value before the exclusions, coloured by guild. Panel b gives the median move per guild.',
     'This is independent evidence that the 72 excluded deployments really did have broken clocks. A clock '
     'error puts detections at the wrong time of day. For a strictly daytime animal that shows up as night '
     'activity that cannot be real, and for a night animal it shows up as impossible daytime activity. After '
     'removing the affected deployments, 48 of 49 affected diurnal-species curves became less nocturnal and '
     'all 36 affected nocturnal-species curves became more nocturnal. A sign test on the diurnal group gives '
     'p = 2 in 10 trillion. The exclusions moved 105 of 1929 curves by more than half a point and left the '
     'other 1824 essentially untouched.'),
]

TABLES = [
    ('gate0_summary.csv', 'Reproduction gate',
     'Before applying any exclusion, the rebuilt pipeline was required to reproduce the previous run exactly '
     'from source data. It did: the same 1935 curves at the same 490 sites, the same species-by-site set, '
     'detections agreeing to 0.04 percent, and per-species curve correlations of 0.9998 or better.'),
    ('stage1_counts.csv', 'Stage 1 counts before and after the exclusions', ''),
    ('stage2_measure_verdicts.csv', 'Which activity measures are reliable enough to use',
     'A measure is carried forward only if its site-to-site differences are larger than its own measurement '
     'error in at least half the species. Activity level fails in 12 of 13 species and is excluded, the same '
     'verdict as the previous run.'),
    ('stage3_partition_anchor_sensitivity.csv', 'Variance partition under every defensible choice', ''),
    ('stage4_block_cv_summary.csv', 'Spatial block cross-validation, all 144 fits', ''),
    ('stage4_effective_sample_size.csv', 'Effective sample size and honest per-pixel uncertainty',
     'Sites are not independent, so the effective sample size is smaller than the site count. For percent '
     'nocturnal the median species has 99 sites but an effective sample of 14, an inflation factor of 5.2. '
     'The resulting 95 percent interval for a prediction at a new place is 43 percentage points wide, against '
     'an observed spread across sites of 35 points. The uncertainty is wider than the variation being mapped.'),
    ('stage5_hypothesis_results.csv', 'All hypothesis tests with nulls and spatial controls', ''),
    ('stage6_attribution.csv', 'Every headline number, before and after', ''),
    ('stage7_grid_covariate_validation.csv', 'Skill when covariates are measured the same way at fit and prediction',
     'Site covariates were aggregated from cameras, while the prediction grid carries its own layers. Human '
     'population agreed between the two at only r = 0.14, so fitting on one and predicting with the other '
     'would be predicting from a different quantity. Refitting on grid-extracted covariates leaves two of six '
     'combinations with genuine skill, and only those two were mapped.'),
]


def b64(path):
    with open(path, 'rb') as f:
        return base64.b64encode(f.read()).decode()


def table_html(df, maxrows=40):
    d = df.head(maxrows)
    h = '<table><thead><tr>' + ''.join(f'<th>{c}</th>' for c in d.columns) + '</tr></thead><tbody>'
    for _, r in d.iterrows():
        cells = []
        for v in r:
            if isinstance(v, float):
                cells.append('' if not np.isfinite(v) else (f'{v:.4g}'))
            else:
                cells.append(str(v) if v == v else '')
        h += '<tr>' + ''.join(f'<td>{c}</td>' for c in cells) + '</tr>'
    h += '</tbody></table>'
    if len(df) > maxrows:
        h += f'<p class="note">Showing {maxrows} of {len(df)} rows. The full table is the CSV artifact.</p>'
    return h


def build(summary_html, outfile='analysis_report_rerun3.html'):
    css = """
    body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',Helvetica,Arial,sans-serif;
         max-width:1080px;margin:0 auto;padding:36px 28px;color:#1a1a1a;line-height:1.62;font-size:15px}
    h1{font-size:26px;margin:0 0 6px;font-weight:650;letter-spacing:-0.2px}
    h2{font-size:19px;margin:44px 0 12px;padding-bottom:6px;border-bottom:2px solid #e3e3e3;font-weight:620}
    h3{font-size:15.5px;margin:26px 0 8px;font-weight:620}
    .sub{color:#666;font-size:13.5px;margin-bottom:26px}
    figure{margin:22px 0 30px}
    figure img{width:100%;border:1px solid #e0e0e0;border-radius:3px}
    figcaption{font-size:13px;color:#333;margin-top:10px}
    figcaption .tech{display:block;margin-bottom:8px}
    figcaption .plain{display:block;background:#f6f8fa;border-left:3px solid #8aa4c8;
                      padding:9px 12px;border-radius:2px;color:#25384f}
    figcaption b{font-weight:640}
    table{border-collapse:collapse;width:100%;font-size:11.5px;margin:12px 0}
    th,td{border:1px solid #e0e0e0;padding:4px 7px;text-align:left}
    th{background:#f2f4f6;font-weight:620}
    tbody tr:nth-child(even){background:#fafbfc}
    .note{font-size:12px;color:#777;margin:6px 0 0}
    .kv{background:#f6f8fa;border:1px solid #e3e6e8;border-radius:4px;padding:14px 18px;margin:18px 0}
    .kv b{font-weight:640}
    ul{margin:8px 0 8px 4px;padding-left:20px}
    li{margin:4px 0}
    .flag{background:#fff8e6;border-left:3px solid #d9a53b;padding:9px 12px;margin:14px 0;font-size:13.5px}
    """
    parts = [f'<!DOCTYPE html><html><head><meta charset="utf-8">'
             f'<title>Continental diel activity, rerun after clock-error exclusions</title>'
             f'<style>{css}</style></head><body>']
    parts.append('<h1>How mammal daily activity varies across the United States</h1>')
    parts.append('<p class="sub">Thirteen species, 1 August to 31 October, all years pooled. '
                 'Recomputed after excluding 72 deployments with clock errors confirmed against archived '
                 'photographs.</p>')
    parts.append(summary_html)
    parts.append('<h2>Figures</h2>')
    for fn, title, tech, plain in FIGS:
        parts.append(f'<figure><img src="data:image/png;base64,{b64(fn)}" alt="{title}">'
                     f'<figcaption><b>{title}.</b>'
                     f'<span class="tech">{tech}</span>'
                     f'<span class="plain"><b>What it means.</b> {plain}</span>'
                     f'</figcaption></figure>')
    parts.append('<h2>Tables</h2>')
    for fn, title, note in TABLES:
        df = pd.read_csv(fn)
        parts.append(f'<h3>{title}</h3>')
        if note:
            parts.append(f'<p>{note}</p>')
        parts.append(table_html(df))
    parts.append('</body></html>')
    html = ''.join(parts)
    with open(outfile, 'w') as f:
        f.write(html)
    return outfile, len(html)
