
import numpy as np, pandas as pd, matplotlib.pyplot as plt

def build_atlas2(D, SPECIES13, CUR5, outfile="atlas2_diel.png", nbins=48):
    """Atlas 2: per-species diel curves in sun time, plus array-level spread."""
    apply_figure_style(sizes=(8,7,6))
    order = D.groupby("common_name").night.mean().sort_values().index.tolist()
    nrow, ncol = 4, 4
    fig, axes = plt.subplots(nrow, ncol, figsize=(13.6, 10.4))
    edges = np.linspace(0, 24, nbins+1); ctr = (edges[:-1]+edges[1:])/2
    P = {}
    for k, sp in enumerate(order):
        ax = axes.ravel()[k]; P[sp] = ax
        g = D[D.common_name==sp]
        h, _ = np.histogram(g.suntime, bins=edges)
        dens = 100*h/h.sum()
        col = "#1f6f8b" if sp in CUR5 else "#e08a3c"
        # night shading: suntime 12-24 is night by construction
        ax.axvspan(12, 24, color="#2c3e50", alpha=.10, lw=0)
        # per-array curves, thin, to show spread
        per = g.groupby("final_array").size()
        for a_ in per[per>=50].index[:120]:
            hh, _ = np.histogram(g[g.final_array==a_].suntime, bins=edges)
            if hh.sum() >= 50:
                ax.plot(ctr, 100*hh/hh.sum(), lw=.35, color=col, alpha=.13, zorder=1)
        ax.plot(ctr, dens, lw=1.7, color=col, zorder=3)
        ax.set_xlim(0,24); ax.set_xticks([0,6,12,18,24])
        ax.set_xticklabels(["sunrise","midday","sunset","midnight","sunrise"], fontsize=5.4)
        ax.set_ylim(0, max(dens.max()*1.35, 4))
        ax.set_title(f"{sp}\n{100*g.night.mean():.1f}% nocturnal, {len(g):,} det, "
                     f"{int((per>=25).sum())} arrays", loc="left", fontsize=6.4)
        if k % ncol == 0: ax.set_ylabel("% of detections")
    for k in range(len(order), nrow*ncol):
        axes.ravel()[k].axis("off")
    ax_leg = axes.ravel()[len(order)]
    ax_leg.axis("off")
    ax_leg.annotate("Thick line: species mean curve\nThin lines: individual arrays "
                    "(\u2265 50 detections)\nShaded: night (sun time 12\u201324)\n\n"
                    "Blue = core 5, orange = added 8\n\nX axis is DOUBLE-ANCHORED sun time:\n"
                    "sunrise\u21920, sunset\u219212, so day and\nnight are each stretched to 12 h",
                    (.02,.94), xycoords="axes fraction", va="top", fontsize=6.2)
    fig.suptitle("Data atlas 2 \u2014 diel activity curves in sun-anchored time, per species and per array",
                 fontsize=9.8, y=.995)
    fig.tight_layout(rect=[0,0,1,.975])
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    return fig, P
