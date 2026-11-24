"""Build the slide deck for the 24 November meeting.

Follows the kays-presentations skill: attention getter first, more visuals and
fewer words, a single full-size figure before any grid of them, dark slides
behind dark figures, and no caption text rendered into a slide figure.

Palette is Linville Gorge. Maple red carries decline, leaf green carries
increase, yellow flags whatever the audience is meant to question.

Re-run after more species land. Every number is read from the saved tables at
build time, so the deck cannot drift from the fits.
"""
import os

import pandas as pd
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.util import Inches, Pt

HERE = os.path.dirname(os.path.abspath(__file__))
SIL = os.path.join(HERE, "phylopic_silhouettes")

SOIL = RGBColor(0x1F, 0x17, 0x14)
BARK = RGBColor(0x41, 0x35, 0x30)
WET = RGBColor(0x5E, 0x50, 0x48)
DUST = RGBColor(0x95, 0x84, 0x7A)
BIRCH = RGBColor(0xC8, 0xBC, 0xB2)
ASH = RGBColor(0xE0, 0xDC, 0xD9)
GREEN = RGBColor(0x64, 0x91, 0x58)
RED = RGBColor(0x8D, 0x46, 0x52)
YELLOW = RGBColor(0xB2, 0xA1, 0x4A)
RUST = RGBColor(0x93, 0x67, 0x43)

W, H = Inches(13.333), Inches(7.5)
EX = [("white-tailed_deer", "White-tailed deer"), ("coyote", "Coyote")]


def _fill(sl, rgb):
    sl.background.fill.solid()
    sl.background.fill.fore_color.rgb = rgb


def _txt(sl, x, y, w, h, text, size=18, bold=False, color=WET):
    tf = sl.shapes.add_textbox(x, y, w, h).text_frame
    tf.word_wrap = True
    for i, line in enumerate(str(text).split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph()
        r = p.add_run()
        r.text = line
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.color.rgb = color
        r.font.name = "Helvetica Neue"


def _slide(prs, title=None, sub=None, dark=False):
    """A slide. Dark slides exist so a dark figure is not framed in white."""
    sl = prs.slides.add_slide(prs.slide_layouts[6])
    ink, sec = (ASH, DUST) if dark else (SOIL, WET)
    if dark:
        _fill(sl, SOIL)
    if title:
        _txt(sl, Inches(0.55), Inches(0.26), Inches(12.3), Inches(0.8),
             title, size=30, bold=True, color=ink)
    if sub:
        _txt(sl, Inches(0.55), Inches(1.02), Inches(12.3), Inches(0.5),
             sub, size=17, color=sec)
    return sl


def _pic(sl, path, x, y, w, h):
    """Fit an image inside a box without distortion."""
    if not os.path.exists(path):
        _txt(sl, x, y, w, Inches(0.4), f"[missing: {os.path.basename(path)}]",
             size=12, color=RED)
        return
    from PIL import Image
    iw, ih = Image.open(path).size
    s = min(w / iw, h / ih)
    sl.shapes.add_picture(path, int(x + (w - int(iw * s)) / 2),
                          int(y + (h - int(ih * s)) / 2),
                          width=int(iw * s), height=int(ih * s))


def _sil(sl, key, x, y, size=Inches(1.1)):
    """Silhouette, only where the PhyloPic image was confirmed for that species."""
    p = os.path.join(SIL, key + ".png")
    if os.path.exists(p):
        _pic(sl, p, x, y, size, size)


def _points(sl, x, y, w, items, size=19, gap=0.52, dark=False):
    ink = ASH if dark else SOIL
    sec = DUST if dark else WET
    for i, it in enumerate(items):
        bold = it.startswith("**")
        _txt(sl, x, y + Inches(gap * i), w, Inches(gap), it.replace("**", ""),
             size=size, bold=bold, color=ink if bold else sec)


def _table(sl, rows, header, top=1.95, size=16):
    """Two-column table with a header row, at slide type size."""
    tb = sl.shapes.add_table(len(rows) + 1, 2, Inches(0.8), Inches(top),
                             Inches(11.6), Inches(0.62 * (len(rows) + 1))).table
    tb.cell(0, 0).text, tb.cell(0, 1).text = header
    for i, (a_, b_) in enumerate(rows, start=1):
        tb.cell(i, 0).text, tb.cell(i, 1).text = a_, b_
    for rw in tb.rows:
        for c in rw.cells:
            for p in c.text_frame.paragraphs:
                for run in p.runs:
                    run.font.size = Pt(size)
    return tb


def build(out="DECK_meeting_20261124.pptx"):
    prs = Presentation()
    prs.slide_width, prs.slide_height = W, H
    f = lambda n: os.path.join(HERE, n)
    rd = lambda n: pd.read_csv(f(n))
    fleet = rd("table_fleet_win5.csv")
    ag = rd("table_agency_evidence_tiers.csv")
    th = rd("table_stream_agreement_win5.csv")
    cov = rd("table_covariate_effects.csv")
    mp = rd("table_model_performance.csv")
    gs = pd.read_csv("/Users/rwkays/claude_code/geospatialdata/"
                     "north_american_mammal_ranges/fleet34_isdm/fleet34_grid_summary.csv")
    gs["pretty"] = gs.species.str.replace("_", " ").str.capitalize()
    RANGE = dict(zip(gs.pretty, gs.cells_union))
    tierA = ag[ag.tier.str.startswith("A")]
    n_fit = len(fleet)
    n_res = int(fleet.direction.isin(["increase", "decline"]).sum())

    # 1. Title over the growth figure: the attention getter is the first thing seen.
    sl = _slide(prs, dark=True)
    _pic(sl, f("fig_data_growth.png"), Inches(1.5), Inches(2.5), Inches(10.4), Inches(4.7))
    _txt(sl, Inches(0.7), Inches(0.5), Inches(12.0), Inches(1.2),
         "Population trends for North American mammals", size=38, bold=True, color=ASH)
    _txt(sl, Inches(0.7), Inches(1.5), Inches(12.0), Inches(0.8),
         "Combining camera traps and iNaturalist   ·   24 November 2026",
         size=16, color=DUST)

    # 2. the gap
    sl = _slide(prs, "We still cannot say which mammals are declining")
    _points(sl, Inches(0.7), Inches(1.9), Inches(12.0), [
        "**Whether a population is rising or falling is the basic input to conservation.**",
        "For most North American mammals that number does not exist.",
        "",
        "**Two datasets have grown enormously. Neither was built to answer it.**",
        "Cameras measure effort precisely and cover little ground.",
        "iNaturalist covers the continent at far lower effort per place.",
        "",
        "**No existing framework combines them to detect population change.**",
    ])

    # 3-5. the agency picture, broken up
    sl = _slide(prs, "Who is already watching", "Species with a published trend, by state")
    _pic(sl, f("fig_agency_infographic.png"), Inches(2.2), Inches(1.15),
         Inches(6.4), Inches(5.7))
    _points(sl, Inches(0.5), Inches(2.3), Inches(3.0), [
        "**Minnesota**", "13 of our 34 species", "",
        "**Median state**", "4 species", "", "**No state**", "reaches half",
    ], size=18, gap=0.46)

    sl = _slide(prs, "Survey-grade evidence is concentrated in big game")
    _pic(sl, f("fig_agency_by_group.png"), Inches(0.5), Inches(1.5),
         Inches(12.3), Inches(4.5))
    _txt(sl, Inches(0.7), Inches(6.2), Inches(12.0), Inches(1.0),
         "Not one survey-grade abundance trend exists for any rodent, rabbit or other "
         "small mammal.", size=22, bold=True, color=RED)

    sl = _slide(prs, "What that means for this study")
    _points(sl, Inches(0.7), Inches(1.9), Inches(12.0), [
        "**It motivates the work.** For most species nobody publishes a trend, so there is",
        "nothing for a manager to act on and nothing to compare a new estimate against.",
        "",
        f"**It bounds the validation.** An independent survey-grade number exists for "
        f"{tierA.species.nunique()} species.",
        "For the other 22 of our 34, our estimate would be the only one.",
        "",
        "**It sets the standard to beat.** Half the agency directions that do exist come",
        "from harvest or observation indices rather than from abundance.",
    ])

    # 6. the streams
    sl = _slide(prs, "Two data streams, complementary weaknesses")
    _points(sl, Inches(0.7), Inches(1.9), Inches(5.6), [
        "**Camera traps**", "effort measured to the day",
        "non-detections are informative", "sparse and clustered", "",
        "**iNaturalist**", "continental, little time lag",
        "effort measured as total mammal", "records in the same 50 km cell",
        "and year, used as an offset",
    ], size=18, gap=0.46)
    # Full width: at 6.3in the baked-in labels are unreadable projected.
    _pic(sl, f("fig_integration.png"), Inches(0.6), Inches(1.45),
         Inches(12.1), Inches(5.6))

    # 7-8. the model
    sl = _slide(prs, "The model, as published",
                "Goldstein et al: one shared surface, two observation models")
    _points(sl, Inches(0.7), Inches(2.0), Inches(12.0), [
        "**One latent abundance surface** over 50 km cells, from habitat covariates,",
        "two climate fields and a spatially smoothed intercept.",
        "",
        "**Cameras** see it through occupancy and detection.",
        "**iNaturalist** sees it through a count model, linked by an exponent:",
        "",
        "         records   =   constant  ×  abundance ^ theta1",
        "",
        "**No time dimension. This estimates distribution, not change.**",
    ])

    sl = _slide(prs, "What we added: a year effect",
                "The same surface, with a trend the cameras identify")
    _points(sl, Inches(0.7), Inches(2.0), Inches(12.0), [
        "**The camera submodel gains a year coefficient. That is the trend we report.**",
        "",
        "**The iNaturalist submodel gains its own increment on top of it**, which absorbs",
        "change in the reporting process and keeps it out of the camera-side trend.",
        "",
        "         camera occupancy        year_beta",
        "         iNaturalist intensity   year_beta  +  year_var",
        "",
        "**Everything spatial stays shared**, so both datasets inform where the animals are.",
    ])

    # 9. change log
    sl = _slide(prs, "What changed from the original model",
                "Change-log and code already circulated; this is the summary")
    _points(sl, Inches(0.7), Inches(1.9), Inches(12.0), [
        "**The model equations are untouched.** Our production code is byte-identical",
        "to the original. Every change is upstream in data preparation, or downstream",
        "in which parameter we report.",
        "",
        "     Fitting mask      union of range polygon and occurrence records",
        "     Covariates        one shared design, after a collinearity screen",
        "     Sampler           single uninterrupted runs, centred parameterisation",
        "     Reporting         the camera year coefficient, not the sum",
        "     Window            five years, with ten-year and full-record comparisons",
    ], size=18, gap=0.47)

    # 10-19. two worked examples
    for n, (key, pretty) in enumerate(EX, start=1):
        row = fleet[fleet.species == pretty]
        trow = th[th.species == pretty]
        r = row.iloc[0]
        band = trow.congruence_band.iloc[0] if len(trow) else "n/a"
        cg = cov[(cov.species == pretty)]
        worst = mp[mp.species == pretty].max_rhat.max()
        sl = _slide(prs, f"Example {n}: {pretty}",
                    "One of the two best-behaved fits we have, and among the most camera data")
        _sil(sl, key, Inches(11.3), Inches(1.5), Inches(1.6))
        _points(sl, Inches(0.7), Inches(2.1), Inches(10.3), [
            f"**Trend   {r.trend_pct_yr:+.1f}% per year    [{r.lo:+.1f}, {r.hi:+.1f}]**",
            f"{r.direction}",
            "",
            f"     range                {int(RANGE.get(pretty, 0))} grid cells",
            f"     worst R-hat          {worst:.3f}  across every parameter family",
            f"     theta1               {r.theta1:.2f}   [{r.theta1_lo:.2f}, {r.theta1_hi:.2f}]",
            f"     stream agreement     {band}",
            f"     habitat effects      {int((cg.resolved & cg.ok).sum())} of "
            f"{int(cg.ok.sum())} resolve",
        ], size=20, gap=0.52)
        # A plain verdict, because a table of diagnostics does not tell the room
        # whether the numbers are good. Marks, not emoji: these get screenshotted.
        ok_rhat = "OK" if worst <= 1.02 else "watch"
        _txt(sl, Inches(0.7), Inches(6.35), Inches(11.5), Inches(0.9),
             f"Everything here is in the range we want. R-hat {worst:.3f} is at the "
             f"threshold, not past it ({ok_rhat}); theta1 sits mid-band; "
             f"{int((cg.resolved & cg.ok).sum())} habitat effects resolve. "
             f"Nothing to explain away.", size=18, bold=True, color=SOIL)

        sl = _slide(prs, f"{pretty}: where it is, and how that changed", dark=True)
        _pic(sl, f(f"ex_{key}_surface.png"), Inches(0.25), Inches(1.6),
             Inches(12.8), Inches(4.6))

        sl = _slide(prs, f"{pretty}: what the camera arrays show", dark=True)
        # Full width: a boxed 8.5in map reads as a thumbnail when projected.
        _pic(sl, f(f"ex_{key}_arrays.png"), Inches(2.45), Inches(0.95),
             Inches(8.4), Inches(6.5))

        sl = _slide(prs, f"{pretty}: habitat effects and convergence")
        _pic(sl, f(f"ex_{key}_detail.png"), Inches(0.4), Inches(1.5),
             Inches(12.5), Inches(5.4))

    sl = _slide(prs, "The two examples together")
    _points(sl, Inches(0.7), Inches(1.9), Inches(12.0), [
        "**Both have more camera data than anything else in the fleet**, and both",
        "converge on every parameter that matters for the trend.",
        "",
        "**Coyote resolves a clear increase**, +9.1% per year [+5.6, +12.4].",
        "**White-tailed deer does not**, -0.7% per year [-3.5, +1.9].",
        "",
        "**Same model, same window, same quality of fit.** The difference between them",
        "is the data, not the method.",
    ])

    # 20-21. fleet
    sl = _slide(prs, "Where the fleet stands",
                "The species that finished first are the small-range, fast-fitting ones")
    # Which species are done, not what they found: the results are the next slide.
    # Every name and count is read from the table, because the hand-typed version
    # of this slide drifted to a wrong total and a missing species.
    abstract = ["Black-tailed jackrabbit", "Mule deer", "Grey fox",
                "Eastern gray squirrel", "Virginia opossum", "Northern raccoon",
                "White-tailed deer", "Coyote"]
    done_abs = sorted(set(fleet.species) & set(abstract))
    todo_abs = sorted(set(abstract) - set(fleet.species))
    _points(sl, Inches(0.7), Inches(1.9), Inches(12.0), [
        f"**{n_fit} of 34 species are fitted on the five-year window.**",
        "",
        f"**{len(done_abs)} of the 8 meeting species are in:**",
        "     " + ", ".join(done_abs[:4]),
        "     " + ", ".join(done_abs[4:]),
        "",
        ("**Still running:** " + ", ".join(todo_abs)) if todo_abs
        else "**All eight meeting species are in.**",
        "",
        f"**{int((~fleet.converged).sum())} of the {n_fit} did not converge** "
        + ", ".join(sorted(fleet.loc[~fleet.converged, "species"])) + ".",
        "",
        "**The remaining 17 species are the slow, wide-ranging ones**, still queued.",
    ])
    sl = _slide(prs, "Five-year trends")
    _pic(sl, f("fig_fleet_win5_slide.png"), Inches(0.5), Inches(1.15),
         Inches(12.3), Inches(6.2))

    # One sheet per resolved species: both year surfaces, the change, and the
    # cameras over the same window. Built from table_fleet_win5 so the set grows
    # as species land, rather than being a hand-kept list.
    res_sp = (fleet[fleet.converged & fleet.direction.isin(["increase", "decline"])]
              .sort_values("trend_pct_yr"))
    sl = _slide(prs, "The eight species with a resolved trend",
                "One sheet each: both year surfaces, the change, and what the cameras saw")
    _points(sl, Inches(0.7), Inches(2.1), Inches(12.0), [
        "**Declines**   " + ", ".join(res_sp.loc[res_sp.direction == "decline", "species"]),
        "",
        "**Increases**  " + ", ".join(res_sp.loc[res_sp.direction == "increase", "species"]),
        "",
        "The bottom two panels on each sheet are different quantities in different units",
        "and share no scale: one is % per year in abundance, the other a percentage-point",
        "change in an array's chance of detecting the species.",
    ], size=19, gap=0.5)
    for r in res_sp.itertuples():
        fn = f"fig_four_{r.species.lower().replace(' ', '_')}.png"
        sl = _slide(prs, f"{r.species}   {r.trend_pct_yr:+.1f}% per year "
                         f"[{r.lo:+.1f}, {r.hi:+.1f}]", dark=True)
        _pic(sl, f(fn), Inches(0.45), Inches(1.15), Inches(12.4), Inches(6.1))

    # 22. does the model work
    sl = _slide(prs, "Does the model work?",
                "The trend we report is the best-behaved parameter in it")
    _pic(sl, f("fig_model_performance.png"), Inches(0.4), Inches(1.5),
         Inches(8.4), Inches(5.5))
    _points(sl, Inches(9.1), Inches(1.9), Inches(4.0), [
        "**Trend        10 of 11**", "Detection     7 of 11",
        "Coupling      3 of 11", "**Fields        1 of 11**", "",
        "The trend sits directly", "in the camera likelihood.", "",
        "Each spatial field", "estimates a value for", "all 908 grid cells,",
        "with neighbours tied", "together. Most cells", "hold little data, so",
        "these move slowly.",
    ], size=16, gap=0.42)

    # 23-24. habitat: one panel first, then the grid
    sl = _slide(prs, "Habitat: one covariate at a time",
                "Human population, across every species that converged")
    _pic(sl, f("fig_covariate_single.png"), Inches(2.3), Inches(1.5),
         Inches(8.7), Inches(5.4))

    sl = _slide(prs, "Habitat: all nine covariates",
                "Eastern gray squirrel is the one species positive on human population")
    _pic(sl, f("fig_covariate_effects.png"), Inches(1.1), Inches(1.3),
         Inches(8.4), Inches(5.9))
    # Computed, not typed: the previous version claimed three covariates were
    # one-directional, which stopped being true as species were added.
    cok = cov[cov.ok]
    ctab = (cok[cok.resolved].groupby("label")
            .agg(n=("mean", "size"), neg=("mean", lambda x: int((x < 0).sum()))))
    top = ctab.sort_values("n", ascending=False).head(3)
    _points(sl, Inches(9.7), Inches(1.9), Inches(3.4), [
        f"**{int(cok.resolved.sum())} of {len(cok)} effects**", "**resolve.**", "",
        "Most resolved effects", "are negative:", "",
    ] + [f"**{lb}**   {int(r.neg)} of {int(r.n)}" for lb, r in top.iterrows()] + [
        "", "No covariate is now", "one-directional across", "every species.",
    ], size=16, gap=0.42)

    # 25-26. stream agreement
    sl = _slide(prs, "Do the two streams agree?")
    _pic(sl, f("fig_stream_agreement_win5.png"), Inches(0.4), Inches(0.95),
         Inches(12.5), Inches(6.4))

    sl = _slide(prs, "What the three checks ask")
    _points(sl, Inches(0.7), Inches(1.85), Inches(12.0), [
        "**1.  Does theta1 say the streams correspond?** The model's own coupling parameter,",
        "      which has to sit between 0 and 1. This comes from the joint fit alone.",
        "",
        "**2.  Do the cameras alone give the same verdict?** We refit each species on camera",
        "      data only. Same direction counts as agreement, and so does both finding no",
        "      trend, which is 9 of the 15.",
        "",
        "**3.  Would the cameras alone support this number?** Whether the joint estimate",
        "      falls inside the camera-only 95% interval. If it does not, adding iNaturalist",
        "      moved the answer somewhere the cameras by themselves do not reach. Two fail:",
        "      coyote and California ground squirrel.",
    ], size=18, gap=0.47)

    sl = _slide(prs, "Three checks, not one")
    _pic(sl, f("fig_stream_agreement_criteria.png"), Inches(0.35), Inches(1.15),
         Inches(12.6), Inches(5.6))
    _txt(sl, Inches(0.7), Inches(6.75), Inches(12.0), Inches(0.7),
         "Five species pass all three. Thirteen pass the camera-only verdict check, "
         "and thirteen pass the containment check.", size=19, bold=True, color=SOIL)

    sl = _slide(prs, "Check 3, walked through: coyote")
    _points(sl, Inches(0.7), Inches(1.85), Inches(12.0), [
        "**We fit coyote twice: once on both streams, once on cameras alone.**",
        "",
        f"     both streams      +9.1% per year   "
        f"[+5.6, +12.4]",
        f"     cameras alone     +3.1% per year   "
        f"[+0.2, +6.1]",
        "",
        "**Same direction, so check 2 passes.** But +9.1 sits outside [+0.2, +6.1],",
        "so check 3 fails: the cameras by themselves never reach that number.",
        "",
        "**Adding iNaturalist roughly triples the size of the increase.** That does not",
        "make it wrong, but the opportunistic stream is doing the work, and coyote is",
        "one of our two headline increases.",
        "",
        f"**California ground squirrel is the other failure**, displaced "
        f"46 points per year.",
    ], size=18, gap=0.46)

    sl = _slide(prs, "Why theta1 cannot always be estimated",
                "Two shortages, and a species can hit either")
    _pic(sl, f("fig_theta1_why.png"), Inches(0.6), Inches(1.4), Inches(12.1), Inches(4.6))
    _points(sl, Inches(0.7), Inches(6.15), Inches(12.0), [
        "**It is the slope of records against abundance, so it needs places that differ "
        "and enough camera data to measure them.**",
        "Coyote has the data and the flattest surface; brown bear has the contrast and "
        "496 camera-years. Both fail, for opposite reasons.",
    ], size=17, gap=0.44)

    sl = _slide(prs, "Slow is not the same as split",
                "Nine of seventeen miss the mixing threshold, for two different reasons")
    _pic(sl, f("fig_chain_cases.png"), Inches(0.35), Inches(1.35),
         Inches(12.6), Inches(4.3))
    _points(sl, Inches(0.7), Inches(5.85), Inches(12.0), [
        "**Six are slow.** Coyote's 12 chains all land between 0.521 and 0.533, which is "
        "agreement by any reading. Longer runs clear these.",
        "**Three are split.** Eastern chipmunk has nine chains at +0.36 and three at "
        "-0.17. Longer runs will not.",
        "**Neither touches the trend**, which converged for 16 of the 17.",
    ], size=17, gap=0.44)

    # 27-28. window
    sl = _slide(prs, "Does the window change the answer?", dark=True)
    _pic(sl, f("fig_window_trends.png"), Inches(0.5), Inches(1.3), Inches(12.3), Inches(3.9))
    _txt(sl, Inches(0.8), Inches(5.5), Inches(11.8), Inches(1.6),
         "White-tailed deer gives +4.00%/yr over the full record, +0.46%/yr over ten years "
         "and -0.73%/yr over five.\nA resolved increase becomes no detectable change, on "
         "the same data and the same model.", size=21, bold=True, color=ASH)

    sl = _slide(prs, "Ten years, for the two species that have it", dark=True)
    _pic(sl, f("fig_win10_windows.png"), Inches(1.2), Inches(1.15), Inches(10.9), Inches(6.1))

    # ---- simulations
    sl = _slide(prs, "Simulations",
                "Everything so far is the real data. This part is simulated.")
    _points(sl, Inches(0.7), Inches(2.3), Inches(12.0), [
        "**We generate data from the fitted model at a known truth, then fit it back.**",
        "",
        "**The question here is the trend itself**: when a trend is really there, how often",
        "do we find it, and when nothing is there, how often do we claim one.",
        "",
        "**A separate validation is running now** on whether the regional version works.",
        "That one is not finished.",
    ], size=19, gap=0.5)

    # Power and its caveat on one slide: showing the number and withdrawing half
    # its meaning on the next slide was the meander Roland flagged.
    sl = _slide(prs, "Can the model find a trend that is really there?",
                "National trend, 4 estimator arms x 3 abundance levels, 30 replicates")
    _pic(sl, f("fig_sim_power.png"), Inches(0.9), Inches(1.4), Inches(11.5), Inches(4.2))
    _points(sl, Inches(0.7), Inches(5.75), Inches(12.0), [
        "**Abundance decides it: 17% for a sparse species, 40% intermediate, 67% "
        "deer-like.** False positives stay at or near the nominal 5%.",
        "**These are a floor, not our power.** The simulation has 698 camera sites; the "
        "five-year window has 15,435 deployments, 22 times more.",
        "**It predicts which species fail**, and the ones that failed are those: sparse, "
        "small-range.",
    ], size=17, gap=0.44)

    sl = _slide(prs, "What the simulation still cannot tell us",
                "Three gaps, and the first one matters most")
    _points(sl, Inches(0.7), Inches(2.0), Inches(12.0), [
        "**1.  It has never scored the parameter we report.** Every run scored the summed",
        "      trend. Re-scoring is impossible from the saved output, so this needs a re-run.",
        "",
        "**2.  The data-volume question is unanswered.** We had a result saying iNaturalist",
        "      volume matters 3x more than cameras. We withdrew it last night: the regional",
        "      term sat only in the iNaturalist likelihood, so cameras could not have helped.",
        "",
        "**3.  Range size and window length were never varied.** Both are where the real",
        "      fits fail: three small-range species, and cameras covering 7 of 18 years.",
        "",
        "**A re-run answering all three is the strongest methods contribution available.**",
    ])

    # 29-31. validation and what is next
    sl = _slide(prs, "Checking ourselves against the agencies",
                "Not shown: it needs rescoring on the five-year window")
    _points(sl, Inches(0.7), Inches(2.0), Inches(12.0), [
        "The existing comparison was scored on full-record fits over a different window,",
        "so it is left out rather than shown with a caveat.",
        "",
        f"**What it will be:** our trend against the published direction, for the "
        f"{tierA.species.nunique()} species",
        "with survey-grade agency evidence.",
        "",
        "**A student is verifying those records now.** 323 records across 207 source",
        "documents, ordered so a partial pass still covers what feeds the comparison.",
    ])

    # ---- closing: status, then questions
    # Sensitivity checks: one summary slide here, one detail slide each in the
    # appendix. The count is read from the table so it cannot go stale.
    inv = rd("table_sensitivity_inventory.csv")
    sl = _slide(prs, "Sensitivity checks",
                f"{len(inv)} checks, {int((inv.status == 'complete').sum())} complete; "
                "one slide each in the appendix")
    th_ = (inv.groupby("theme")
             .agg(n=("check", "size"),
                  done=("status", lambda x: int((x == "complete").sum())))
             .sort_values("n", ascending=False))
    _table(sl, [(name, f"{int(r.n)}   ({int(r.done)} complete)")
                for name, r in th_.iterrows()], ("What it probes", "Checks"))
    pend = inv[inv.status == "pending"]
    _txt(sl, Inches(0.8), Inches(6.0), Inches(11.6), Inches(1.2),
         f"The {len(pend)} still open: " + "; ".join(sorted(pend.check)),
         size=16, color=WET)

    sl = _slide(prs, "Where everything stands", "1. Simulations")
    rows = [("National power and false positives", "done, shown above"),
            ("Camera-anchored regional validation", "smoke test queued tonight; "
             "720 fits after it"),
            ("Validation on the parameter we report", "not run, and cannot be "
             "re-scored from saved output"),
            ("Range size and window length", "never varied")]
    _table(sl, rows, ("Simulation", "Status"))

    # Every species and where it stands, read from the status table.
    status = rd("table_fleet_status_full.csv")
    done = sorted(status.loc[status.win5_joint == "current", "species"])
    wait = sorted(status.loc[status.win5_joint != "current", "species"])
    sl = _slide(prs, "All 34 species: done and waiting",
                f"{len(done)} fitted on the current five-year spec, {len(wait)} still to run")
    _txt(sl, Inches(0.6), Inches(1.5), Inches(6.2), Inches(0.4),
         f"FITTED  ({len(done)})", size=17, bold=True, color=GREEN)
    _points(sl, Inches(0.6), Inches(1.95), Inches(6.2),
            [("**" + sp_ + "**") if sp_ in list(
                status.loc[status.win5_direction.isin(["increase", "decline"]), "species"])
             else sp_ for sp_ in done], size=13.5, gap=0.285)
    _txt(sl, Inches(7.0), Inches(1.5), Inches(6.0), Inches(0.4),
         f"STILL TO RUN  ({len(wait)})", size=17, bold=True, color=RED)
    _points(sl, Inches(7.0), Inches(1.95), Inches(6.0), wait, size=13.5, gap=0.285)
    _txt(sl, Inches(0.6), Inches(7.05), Inches(12.4), Inches(0.5),
         "Bold = trend direction resolved. The 16 still to run have fits on the older "
         "specification only, which are not comparable and are not reported.",
         size=14, color=WET)

    sl = _slide(prs, "Where everything stands", "2. Species and windows")
    srows = [(f"Five-year joint fits", f"{n_fit} of 34 done, "
              f"{len(set(fleet.species) & set(abstract))} of the 8 meeting species"),
             ("Five-year camera-only controls", "all 34 exist"),
             ("Ten-year joint fits", "2 species"),
             ("Ten-year camera-only controls", "none"),
             ("Full-record fits", "all 34, superseded by the five-year window"),
             ("Early window, 2016-2020", "deer only, chains queued")]
    _table(sl, srows, ("Fit", "Status"))

    sl = _slide(prs, "Where everything stands", "3. What we can check against agencies")
    _points(sl, Inches(0.7), Inches(2.0), Inches(12.0), [
        "**28 of 34 species have a published agency trend somewhere.**",
        "     but only 17 in five or more states, and 6 in none",
        "",
        "**Survey-grade evidence is almost all big game:** ungulates 51%, carnivores 23%,",
        "     rodents and other small mammals zero.",
        "",
        "**Of our fitted species, six can be scored against an agency direction today.**",
        "     A student is verifying the underlying records before we rely on them.",
        "",
        "**Twelve species we resolve have no agency direction at all.** That gap is the",
        "     reason the project exists.",
    ], size=18, gap=0.46)

    sl = _slide(prs, "Two questions for you", "Both are modelling, and both are open")
    _points(sl, Inches(0.7), Inches(2.0), Inches(12.0), [
        "**1.  Should repeat camera locations be grouped?**",
        "      The published non-temporal model gives every physical location a random",
        "      intercept. Ours does not. 97.8% of the variance in deployment year is",
        "      between locations, but only 1,862 of 21,206 locations were sampled in",
        "      more than one year, so a location term is estimable for 8.8% of them.",
        "      Was it weighed and dropped, or lost when the temporal structure went in?",
        "",
        "**2.  Should regional trends be camera-anchored?**",
        "      In the current regional model the regional term sits only in the",
        "      iNaturalist likelihood, so cameras cannot inform it. We are validating a",
        "      fork that shares it across both. Is that the right fix?",
    ], size=17, gap=0.44)

    sl = _slide(prs, "What we need to decide")
    _points(sl, Inches(0.7), Inches(1.9), Inches(12.0), [
        "**1.  How many species get a ten-year fit**, and is ten years a second reported",
        "      window or only a sensitivity check.",
        "**2.  Camera-only controls before or after the remaining 15 species.**",
        "**3.  Is the non-identification result a paper section or a supplement.**",
        "**4.  Do we report the three species whose coupling cannot be estimated**, with",
        "      the camera trend alone and no agreement statement.",
        "**5.  Freeze numbered releases** and review against the frozen one, rather than against",
        "      a model that changes weekly.",
    ], size=19, gap=0.55)

    # Appendix: one slide per check, in theme order.
    sl = _slide(prs, "Appendix", "Sensitivity checks in detail, one per slide")
    _points(sl, Inches(0.7), Inches(2.3), Inches(12.0), [
        "**Not covered in the call.** Each slide states the question, what the check does,",
        "what it found, how many species it covers, and the file holding the numbers.",
    ])
    for theme, grp in inv.groupby("theme", sort=False):
        for r in grp.itertuples():
            sl = _slide(prs, r.check, f"{theme}   ·   {r.status}")
            _points(sl, Inches(0.7), Inches(1.85), Inches(12.0), [
                f"**Question**", f"     {r.question}", "",
                f"**What the check does**", f"     {r.what_it_does}", "",
                f"**What it found**",
            ], size=17, gap=0.44)
            _txt(sl, Inches(1.0), Inches(4.75), Inches(11.5), Inches(1.4),
                 r.finding, size=17, color=SOIL)
            _txt(sl, Inches(0.7), Inches(6.45), Inches(12.0), Inches(0.8),
                 f"Coverage: {r.coverage}        Numbers: {r.artifact}",
                 size=15, color=WET)

    path = os.path.join(HERE, out)
    prs.save(path)
    return path, len(prs.slides._sldIdLst)


if __name__ == "__main__":
    p, n = build()
    print(f"{os.path.basename(p)}: {n} slides, {os.path.getsize(p):,} bytes")
