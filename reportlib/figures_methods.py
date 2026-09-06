"""Methods schematic for the real-data model: how the two data streams combine.

The counterpart to the simulation design figure. A reader who has never seen
this model should be able to tell, before meeting any numbers, what each data
stream contributes and what the model returns.
"""
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch

from . import style
from .conventions import verify_text_within_box

INK, SUB, RULE = "#22282e", "#5f6a73", "#c8cfd5"
INAT = "#e08214"      # opportunistic stream
CAM = "#2c6c9c"       # structured stream
SHARED = "#4d4d4d"

# Box geometry, reused by the containment checks so the two cannot drift.
X0, W0 = 2.0, 30.0     # data streams; width set by the longest panel title,
                       # which the containment check enforces
X1, W1 = 34.0, 30.0    # shared latent surface
X2, W2 = 69.0, 29.0    # what the model returns


def fig_integration(out="fig_integration.png"):
    """Camera surveys as the primary stream, iNaturalist as corroboration."""
    fig, ax = plt.subplots(figsize=(10.2, 4.3))
    ax.set_xlim(0, 100); ax.set_ylim(1.5, 55.5); ax.axis("off")
    panels = {}          # name -> (rect, [text artists])
    current = [None]

    def panel(x, w, y, h, n, title, colour=INK):
        current[0] = title
        panels[title] = ((x, x + w, y, y + h), [])
        ax.add_patch(Rectangle((x, y), w, h, facecolor="white", edgecolor=RULE,
                               linewidth=0.7, zorder=2))
        ax.add_patch(Rectangle((x + 1.0, y + h - 4.6), 3.0, 3.0, facecolor=colour,
                               edgecolor="none", zorder=3))
        ax.text(x + 2.5, y + h - 3.1, str(n), fontsize=6.4, color="white",
                ha="center", va="center", fontweight="bold", zorder=4)
        t = ax.text(x + 5.4, y + h - 3.1, title, fontsize=7.8, color=colour,
                    va="center", fontweight="bold", zorder=4)
        panels[title][1].append(t)
        return x + 2.2, y + h - 7.2

    def lines(x, y, items, colour=SUB, dy=3.3):
        for k, txt in enumerate(items):
            t = ax.text(x, y - k * dy, txt, fontsize=6.5, color=colour, va="center")
            if txt.strip():
                panels[current[0]][1].append(t)

    # 1 -- the primary stream, first because it is the standard the other is
    # checked against. 2 -- the corroborating stream.
    tx, ty = panel(X0, W0, 28.5, 26.0, 1, "Camera traps (primary)", CAM)
    lines(tx, ty, ["Structured surveys: known",
                   "  locations, known dates",
                   "Detection / non-detection",
                   "  at each site",
                   "Effort is measured, so an",
                   "  absence is real evidence"], CAM)
    tx, ty = panel(X0, W0, 2.0, 26.0, 2, "iNaturalist (corroborating)", INAT)
    lines(tx, ty, ["Public, opportunistic photos",
                   "Counts per 50 km cell per year",
                   "Effort controlled using total",
                   "  mammal records in the same",
                   "  cell and year",
                   "Fills in where cameras are few"], INAT)

    # 3 -- the shared latent surface, and the two trend structures fitted to it
    tx, ty = panel(X1, W1, 2.0, 52.5, 3, "One shared surface", SHARED)
    lines(tx, ty, [
        "Latent local abundance, on a",
        "  50 km grid, per year",
        "",
        "Explained by:",
        "  9 habitat covariates",
        "  a smoothed spatial field",
        "  spatially varying climate",
        "    responses",
        "",
        "Fitted twice, with two trend",
        "structures:",
        "  one national trend, or",
        "  one national trend plus a",
        "    deviation per ecoregion",
        "",
        "Both streams see this same",
        "surface, each through its",
        "own observation model.",
    ], SHARED, dy=2.48)

    # 4 -- outputs
    tx, ty = panel(X2, W2, 2.0, 52.5, 4, "What the model returns", INK)
    lines(tx, ty, [
        "Camera-anchored trend",
        "  what the cameras support",
        "  on their own",
        "",
        "iNaturalist-only increment",
        "  extra trend seen only by",
        "  iNaturalist",
        "",
        "Habitat associations",
        "  one per covariate",
        "",
        "Abundance surface, by year",
        "",
        "The reported trend is the two",
        "added together, so the split",
        "shows how far the cameras",
        "back it up.",
    ], INK, dy=2.60)

    for y in (42.0, 15.0):
        ax.add_patch(FancyArrowPatch((X0 + W0 + 0.4, y), (X1 - 0.6, y),
                                     arrowstyle="-|>", mutation_scale=9,
                                     color=SUB, linewidth=1.0, zorder=5))
    ax.add_patch(FancyArrowPatch((X1 + W1 + 0.4, 28.5), (X2 - 0.6, 28.5),
                                 arrowstyle="-|>", mutation_scale=9,
                                 color=SUB, linewidth=1.0, zorder=5))

    # Titles previously overflowed their boxes in the simulation schematic; the
    # same check is applied here rather than discovering it in the rendered PNG.
    for name, ((x0, x1, y0, y1), artists) in panels.items():
        verify_text_within_box(ax, artists, x0 + 0.5, x1 - 0.5, y0 + 0.5, y1 - 0.5,
                               what=f"panel {name!r}")

    # Caption sits just below the axes; a positive figure-fraction position left
    # a wide empty band because the panels already fill the axes.
    fig.text(0.5, -0.015,
             "The camera surveys are the standard: their effort is measured, so both presences "
             "and absences carry information. iNaturalist cannot stand alone -- effort is\n"
             "unknown and is controlled only indirectly -- but it covers the whole country and "
             "every year, so it extends the cameras' reach and provides an independent check.\n"
             "Because the model separates the camera-anchored trend from the iNaturalist-only "
             "increment, the report can state how much of each trend the cameras support.",
             ha="center", va="top", fontsize=6.4, color=SUB)
    fig.savefig(out, bbox_inches="tight", pad_inches=0.10, dpi=300)
    plt.close(fig)
    return out
