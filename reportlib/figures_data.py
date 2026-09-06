"""Data-context figures: survey effort over time, and the range-mask comparison."""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from . import inputs, style

SPECIES_KEY = {"bobcat": "Bobcat", "white-tailed_deer": "White-tailed deer",
               "moose": "Moose"}
SPECIES_COLOR = {"bobcat": "#b2182b", "white-tailed_deer": "#1b5e9c", "moose": "#4d9221"}
WINDOWS = {"10-year": 2016, "5-year": 2021}


def effort_table():
    """Camera deployments started per year, and iNaturalist records per year."""
    dep = pd.read_csv(os.path.join(os.path.dirname(inputs.REPO), "raw_cam_data",
                                   "combined_deployments_all.csv"), low_memory=False)
    start = pd.to_datetime(dep["start_date"], errors="coerce", format="mixed")
    cam = start.dt.year.value_counts().sort_index()
    cam = cam[(cam.index >= 2008) & (cam.index <= 2025)]

    yc = pd.read_csv(os.path.join(inputs.REPO, "inat_yearly_counts_fleet35_masked.csv"))
    inat = yc[yc.species.isin(SPECIES_KEY)].pivot(index="year", columns="species",
                                                  values="raw_count")
    return cam, inat


def fig_effort(out="fig_effort.png"):
    """Survey effort over the modelled period, for both data streams.

    This is the figure behind the reporting-window question: camera deployments
    are concentrated in the recent half of the record, so a trend fitted over
    the full 2008-2025 span leans on iNaturalist for its early years.
    """
    cam, inat = effort_table()
    fig, (axc, axi) = plt.subplots(2, 1, figsize=(7.0, 4.8), sharex=True)

    for ax, (wlab, wstart) in [(axc, ("10-year", 2016)), (axi, ("10-year", 2016))]:
        ax.axvspan(2015.5, 2025.5, color="#eef2f6", zorder=0)
    for ax in (axc, axi):
        ax.axvspan(2020.5, 2025.5, color="#dde6ee", zorder=0)

    axc.bar(cam.index, cam.values, color="#4d4d4d", width=0.72, zorder=3)
    axc.set_ylabel("Camera deployments\nbegun that year")
    axc.set_title("Camera effort is concentrated in the recent half of the record")
    share = cam[cam.index >= 2016].sum() / cam.sum()
    axc.annotate(f"{share:.0%} of all deployments\nbegan in 2016 or later",
                 xy=(2010.2, cam.max() * 0.86), fontsize=6.6, color="#22282e",
                 ha="left", va="top", linespacing=1.5)

    for sp, lab in SPECIES_KEY.items():
        if sp in inat.columns:
            axi.plot(inat.index, inat[sp], "-o", color=SPECIES_COLOR[sp], ms=3.0,
                     lw=1.3, label=lab, zorder=3)
    axi.set_yscale("log")
    axi.set_ylabel("iNaturalist records\nper year (log scale)")
    axi.set_xlabel("Year")
    axi.set_title("iNaturalist records grow steeply across the same period")
    axi.legend(frameon=False, fontsize=6.6, loc="upper left")
    axi.set_xlim(2007.4, 2025.6)
    axi.set_xticks(range(2008, 2026, 2))

    # Both window labels ride at the same height at the top of the bar panel.
    # An earlier draft dropped the 5-year label part-way down the axis, where it
    # landed on the 2021 bar.
    top = axc.get_ylim()[1]
    axc.text(2015.7, top * 0.985, " 10-year window", fontsize=6.2,
             color="#5f6a73", ha="left", va="top")
    axc.text(2020.7, top * 0.985, " 5-year", fontsize=6.2,
             color="#41586e", ha="left", va="top")
    fig.text(0.5, -0.02,
             "Shaded bands mark the 10-year (2016-2025) and 5-year (2021-2025) reporting windows. "
             "Camera deployment counts are fleet-wide, not\nper species; iNaturalist counts are "
             "per species after masking.",
             ha="center", va="top", fontsize=6.3, color="#5f6a73")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, (cam, inat)


# Short axis labels: the four full definitions are spelled out in the caption,
# and the long forms collided horizontally at this panel width.
