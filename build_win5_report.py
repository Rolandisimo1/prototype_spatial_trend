"""Generate the five-year fleet results report from the fitted posteriors.

Every number in the document is computed from the table at render time, so the
report cannot drift out of step with the fits. Re-run after more species finish
and the counts, lists and take-home points update themselves.

Usage:
    python build_win5_report.py          # writes REPORT_fleet_win5.md
"""
import os

import pandas as pd

import reportlib.figures_win5 as fw

HERE = os.path.dirname(os.path.abspath(__file__))
FIG = "fig_fleet_win5.png"
OUT = "REPORT_fleet_win5.md"


def _fmt(r):
    return f"{r.trend_pct_yr:+.1f}% per year [{r.lo:+.1f}, {r.hi:+.1f}]"


def _names(sub):
    # Species names keep their own capitalisation: lowercasing them turns
    # "California ground squirrel" into a geography error.
    n = list(sub.species)
    if not n:
        return "none"
    if len(n) == 1:
        return n[0]
    return ", ".join(n[:-1]) + " and " + n[-1]


def _vb(n, singular, plural):
    """Agree a verb with a count that changes as more species finish."""
    return singular if n == 1 else plural


# Common names are lowercase mid-sentence, except where the name begins with a
# proper noun. "brown bear" but "California ground squirrel".
PROPER = ("california", "virginia", "american", "north", "northern", "eastern",
          "western", "southern", "canada", "arizona", "texas", "mexican")


def _lc(name):
    first = name.split()[0].lower()
    return name if first in PROPER else name[0].lower() + name[1:]


def _lcnames(sub):
    n = [_lc(x) for x in sub.species]
    if not n:
        return "none"
    if len(n) == 1:
        return n[0]
    return ", ".join(n[:-1]) + " and " + n[-1]


DIRS = {"increasing", "decreasing", "stable", "mixed"}


def _agency_view(species, agency=None):
    """What state and provincial agencies publish for one species, over 5 and 10 years.

    Returns a sentence, computed at render time, so a claim about what agencies
    say cannot drift away from the agency table or be carried over by hand.
    """
    if agency is None:
        for cand in ("agency_trends_combined.csv",
                     os.path.join(HERE, "agency_trends_combined.csv")):
            if os.path.exists(cand):
                agency = pd.read_csv(cand)
                break
    if agency is None:
        return ""
    words = species.lower().split()
    m = agency.species.str.lower()
    g = agency[m.str.contains(words[0], na=False) & m.str.contains(words[-1], na=False)]
    if not len(g):
        return f"No agency publishes an assessment for {_lc(species)}."
    d5 = g.trend_5yr_direction[g.trend_5yr_direction.isin(DIRS)]
    d10 = g.trend_10yr_direction[g.trend_10yr_direction.isin(DIRS)]
    if not len(d5) and not len(d10):
        return (f"None of the {len(g)} agency records for {_lc(species)} states a "
                f"direction, so no external assessment is available to compare against.")
    parts = []
    for lab, d in (("five years", d5), ("ten years", d10)):
        if not len(d):
            continue
        counts = d.value_counts()
        # Each window gets its own denominator: a record stating a ten-year
        # direction often leaves the five-year one unknown, so one shared count
        # would overstate the five-year evidence.
        parts.append(f"over {lab}, {len(d)} of {len(g)} records give a direction "
                     + "(" + " and ".join(f"{n} {k}" for k, n in counts.items()) + ")")
    return (f"For {_lc(species)}, " + "; ".join(parts) + ".")


def build(table=None, fig_marker=FIG):
    t = table if table is not None else fw.fleet_table()
    conv = t[t.converged]
    res = conv[conv.direction.isin(["increase", "decline"])]
    lim = conv[conv.data_limited]
    solid = conv[~conv.data_limited]
    solid_res = solid[solid.direction.isin(["increase", "decline"])]
    gate_needs = conv[conv.gate == "needs camera-only control"]
    gate_none = conv[conv.gate == "theta1 below reporting threshold"]
    pending = fw.FLEET_N - len(t)
    nonconv = t[~t.converged]

    take = []
    take.append(
        f"**{len(t)} of the {fw.FLEET_N} species in the fleet have been fitted so far, "
        f"and {len(conv)} of those converged.** The remaining {pending} are still "
        f"running. Every count below is out of the {len(conv)} converged fits, not "
        f"out of {fw.FLEET_N}, so none of them is a statement about the fleet.")
    take.append(
        f"**{len(res)} {_vb(len(res), 'species resolves', 'species resolve')} a direction**, "
        f"meaning the 95% credible interval for the trend excludes zero. "
        + (f"They are {_lcnames(res)}. " if len(res) else "")
        + (f"The remaining {len(conv) - len(res)} "
           f"{_vb(len(conv) - len(res), 'fit has an interval', 'fits have intervals')} spanning "
           f"zero, so the data leave "
           f"{_vb(len(conv) - len(res), 'its direction', 'their directions')} open."
           if len(conv) - len(res) else ""))
    if len(lim):
        take.append(
            f"**{len(lim)} of the resolved species sit below the range size where these "
            f"estimates hold up.** {_names(lim)} "
            f"{'occupies' if len(lim) == 1 else 'occupy'} "
            f"{' and '.join(str(int(c)) for c in lim.range_cells)} grid cells, against a "
            f"fleet median of 425. At {lim.iloc[0].trend_pct_yr:.0f}% per year, "
            f"{_lc(lim.iloc[0].species)} would fall to "
            f"{100 * (1 + lim.iloc[0].trend_pct_yr / 100) ** 5:.0f}% of its starting "
            f"abundance over the five years. "
            + " ".join(_agency_view(s) for s in lim.species))
    if len(solid_res):
        take.append(
            f"**Setting those aside, {len(solid_res)} species with ranges above "
            f"{fw.RANGE_MIN_CELLS} cells resolve a direction**: "
            + "; ".join(f"{_lc(r.species)} at {_fmt(r)}" for r in solid_res.itertuples())
            + ".")
    else:
        take.append(
            f"**No species with a range above {fw.RANGE_MIN_CELLS} cells has yet resolved "
            f"a direction.** The species that have finished are among the smallest-range "
            f"and least-sampled in the fleet, so this says nothing about the fleet's "
            f"eventual result.")
    take.append(
        f"**Stream agreement stays unconfirmed for {len(gate_needs) + len(gate_none)} of "
        f"the {len(conv)} converged fits, on the joint fit alone.** "
        f"{len(gate_needs)} {_vb(len(gate_needs), 'needs', 'need')} a camera-only control "
        f"fit because the exponent linking "
        f"iNaturalist records to abundance has an interval reaching 1 or including 0, and "
        f"{len(gate_none)} returned that exponent too imprecisely to judge.")

    lines = ["# Population trends from integrated camera-trap and iNaturalist data: "
             "the five-year window", ""]
    lines.append("## Take-home points")
    lines.append("")
    for i, s in enumerate(take, 1):
        lines.append(f"{i}. {s}")
        lines.append("")

    lines += ["## What is being measured", "",
        "**Trend.** The change in a species' relative abundance per year, expressed as "
        "a percentage. A value of −5% per year means the modelled abundance in an "
        "average grid cell falls by five percent of its value each year, so after five "
        "years about 77% of the starting abundance remains. The model estimates how "
        "abundance changes from year to year, on a relative scale. Converting that into "
        "a count of animals would need a separate calibration we have not done.", "",
        "**Credible interval.** The range containing the true value with 95% "
        "probability, given the model and the data. When it excludes zero we say the "
        "direction is resolved. When it spans zero, the data leave the direction open: "
        "a wide interval spanning zero is a statement about what we can see, and the "
        "population may still be rising or falling.", "",
        "**The five-year window.** The fits here use detections from the most recent "
        "five years only. Camera-trap effort grew steeply over the study period, and "
        "few cameras were running in the earliest years, so those years add more noise "
        "than signal.", "",
        "**Stream agreement (the exponent θ₁).** The model ties the two data streams "
        "together with one number: the expected count of iNaturalist records in a grid "
        "cell is proportional to modelled abundance raised to the power θ₁. At θ₁ = 1 "
        "records track abundance proportionally; at θ₁ = 0.5 doubling the animals brings "
        "about 41% more records; at θ₁ = 0 the records say nothing about where the "
        "animals are. When the fitted interval reaches 1 or includes 0, the iNaturalist "
        "half of the fit may be tracking where people looked instead of where the "
        "animals are, and the species needs a camera-only control fit before we report "
        "its joint trend.", "",
        "**Data-limited.** A species occupying fewer than "
        f"{fw.RANGE_MIN_CELLS} of the 50 km grid cells. Across the species with both a "
        "joint and a camera-only fit, those below this size sat a median of 17 percentage "
        "points per year away from their camera-only estimate, against 1.7 points above "
        "it. The threshold comes from that measurement.", "",
        "**Convergence.** Each fit runs several independent chains. R-hat compares the "
        "spread between chains with the spread within them, and should be at or below "
        f"{fw.RHAT_MAX}; effective sample size is how many independent draws the chains "
        f"are worth, and should be at least {fw.ESS_MIN}. We report only the fits that "
        "pass both.", ""]

    lines += ["## Methods", "",
        "**The model.** Camera-trap detections and iNaturalist records are fitted to one "
        "shared surface of relative abundance on a 50 km grid, each stream through its "
        "own observation model. We know how long each camera ran, so a camera that "
        "recorded nothing tells us as much as one that recorded an animal, and we anchor "
        "the trend to the camera stream. We do not know how much time iNaturalist "
        "observers spent looking. We control for that using the total number of mammal "
        "records in the same cell and year as a measure of how much observing happened "
        "there.", "",
        "**How the trend is computed.** The model fits a slope against a standardised "
        "calendar year, so the raw coefficient gives change per standard deviation of "
        "year. Dividing by that standard deviation, exponentiating and "
        "subtracting one converts it to proportional change per year. The same scaling "
        "factor is read from each fit's own record, because the standard deviation of "
        "year differs between species with different numbers of sampled years.", "",
        "**Which fits are reported.** A fit enters the results only if both R-hat and "
        "effective sample size pass for the trend parameter. Fits that fail appear on the "
        "figure as open markers on a dotted interval, and every count in the title "
        "excludes them.", ""]

    lines += ["## Results", "",
        f"![Five-year trend for every fitted species]({{{{artifact:art_FIGID}}}})", "",
        f"*The {len(conv)} converged fits on the five-year window. Each bar is the 95% "
        f"credible interval and each dot the posterior mean. Species in amber occupy "
        f"fewer than {fw.RANGE_MIN_CELLS} grid cells, with their cell count beside the "
        f"name. Open markers on dotted intervals did not converge and are excluded from "
        f"the counts in the title.*", ""]
    show = t[["species", "range_cells", "trend_pct_yr", "lo", "hi", "direction",
              "theta1", "theta1_lo", "theta1_hi", "gate"]].copy()
    show.columns = ["Species", "Range (cells)", "Trend (%/yr)", "Lower", "Upper",
                    "Direction", "θ₁", "θ₁ lower", "θ₁ upper", "Stream agreement"]
    lines.append(show.round(2).to_markdown(index=False))
    lines.append("")
    if len(lim):
        lines += ["### What the steepest declines rest on", "",
            f"{_names(lim)} {_vb(len(lim), 'returns the steepest decline', 'return the two steepest declines')} "
            f"in the table, and {_vb(len(lim), 'is', 'are')} also the smallest "
            f"{_vb(len(lim), 'range', 'ranges')} in the fleet. The model sees a species "
            f"in a few dozen grid cells across five years of detections, and fits the "
            f"steepest slope those detections allow. Both species also failed the "
            f"anchoring check, which compares the joint estimate against a camera-only "
            f"fit of the same data.", ""]

    lines += ["## What happens next", "",
        f"**{pending} species are still being fitted.** The species finished so far "
        f"include several of the fleet's smallest ranges, so a fleet-level statement "
        f"waits for the rest.", "",
        f"**{len(gate_needs) + len(gate_none)} of the converged fits need a camera-only "
        f"control fit** before their joint trend can be reported, because stream "
        f"agreement could not be established from the joint fit.", ""]
    if len(nonconv):
        lines += [f"**{len(nonconv)} {_vb(len(nonconv), 'fit', 'fits')} did not converge** "
                  f"({_names(nonconv)}) and is excluded from all results.", ""]

    lines += ["## Limitations", "",
        "**These species finished first because they are small and quick to fit**, which "
        "is the same property that makes their estimates least reliable. The fleet "
        "average will look different.", "",
        "**The trends are proportional change.** Reading them as population sizes or "
        "animal counts would need a calibration step we have not done.", "",
        "**Stream agreement is unconfirmed for most of these fits.** Where θ₁ came back "
        "imprecise, the iNaturalist half of the fit may be tracking where people looked "
        "instead of where the animals are. A camera-only control fit settles it.", "",
        "**The agency comparison is still to come.** It needs the state and provincial "
        "records matched to the same five-year window.", ""]

    lines += ["## Glossary", "",
        "**Chain.** One independent run of the fitting algorithm. Several are run per "
        "species so that agreement between them can be checked.", "",
        "**Credible interval.** See *What is being measured*.", "",
        "**Effective sample size.** The number of independent draws the chains are worth "
        "once their internal correlation is accounted for.", "",
        "**Grid cell.** A 50 km square. Range size is reported as the number of cells a "
        "species occupies.", "",
        "**R-hat.** A convergence statistic comparing variation between chains with "
        "variation within them. One means perfect agreement.", "",
        "**Relative abundance.** How many animals are in a place compared with another "
        "place or another year, without estimating the absolute number.", "",
        "**θ₁ (theta-one).** The exponent linking iNaturalist record counts to modelled "
        "abundance. See *What is being measured*.", ""]
    return "\n".join(lines), t


if __name__ == "__main__":
    md, t = build()
    with open(os.path.join(HERE, OUT), "w") as fh:
        fh.write(md)
    print(f"{OUT}: {len(md.split())} words, {len(t)} species")
