#!/usr/bin/env python3
"""Redraw fig_theta1_congruence_fleet.png from table_theta1_congruence.csv (2026-11-15).

Rows with estimable == False are drawn as a grey text label with no point or interval: their
chains do not mix, so no single theta1 value describes the posterior (Roland, 11/15: California
ground squirrel's published 2.24 [1.60, 3.16] and the centred refit's -0.01 [-3.88, 2.01] both come
from non-mixing chains). The original figure was made by a script that was never saved; this one
reproduces its layout.
"""
import csv, statistics
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
rows = list(csv.DictReader(open(HERE / "table_theta1_congruence.csv")))
BLUE, ORANGE, RED, GREY = "#1f6db5", "#c4621d", "#c0392b", "#8a8f98"

fig, ax = plt.subplots(figsize=(10.8, 8.6), dpi=220)
ax.axvspan(0, 1, color="#eef3f7", zorder=0)
ax.axvline(0, color=ORANGE, ls="--", lw=1.2)
ax.axvline(1, color=RED, ls="--", lw=1.2)
n = len(rows)
for i, r in enumerate(rows):
    y = n - 1 - i
    if r["estimable"] != "True":
        ax.text(0.02, y, "unestimable: chains do not mix (two modes)", color=GREY, va="center",
                ha="left", fontsize=8.5, style="italic")
        continue
    t, lo, hi = float(r["theta1"]), float(r["theta1_lo"]), float(r["theta1_hi"])
    c = ORANGE if r["ci_includes_0"] == "True" else (RED if r["theta1_ci_reaches_1"] == "True" else BLUE)
    ax.plot([lo, hi], [y, y], color=c, lw=2.2, solid_capstyle="butt")
    ax.plot(t, y, "o", color=c, ms=6.5)
    if r["camonly_supports"].startswith("NOT supported"):
        ax.plot(hi + 0.08, y, "x", color=RED, ms=7, mew=2)
ax.set_yticks(range(n))
ax.set_yticklabels([r["species"] for r in reversed(rows)])
for lab, r in zip(ax.get_yticklabels(), reversed(rows)):
    if r["estimable"] != "True":
        lab.set_color(GREY)
ax.set_ylim(-0.8, n - 0.2)
ax.set_xlim(-0.5, 2.15)
ax.set_xlabel("theta1   (exponent linking iNaturalist records to camera-side abundance)")
ax.set_title("Spatial congruence between the two data streams, by species", loc="left")
ax.text(-0.08, n * 0.5, "no spatial abundance signal", rotation=90, color=ORANGE, fontsize=8, va="center")
ax.text(1.03, n * 0.5, "proportional to abundance", rotation=90, color=RED, fontsize=8, va="center")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(handles=[
    Line2D([], [], color=BLUE, lw=2.5, label="partial correspondence (interval inside 0–1)"),
    Line2D([], [], color=ORANGE, lw=2.5, label="interval includes 0 — no correspondence not excluded"),
    Line2D([], [], color=RED, lw=2.5, label="interval reaches or exceeds 1 — over-response"),
    Line2D([], [], color=RED, marker="x", ls="", mew=2, label="camera-only fit does not support the joint trend"),
], loc="upper right", frameon=False, fontsize=8)
est = [float(r["theta1"]) for r in rows if r["estimable"] == "True"]
med = statistics.median(est)
fig.text(0.5, 0.012,
         "theta1 is the exponent in log(mu) = theta0 + theta1 log(sum lambda): a cell with twice the camera-side "
         "abundance yields 2^theta1 times the iNaturalist records.\n"
         f"Shaded band is the healthy range. Fleet median {med:.2f} ({len(est)} estimable species), so a tenfold "
         f"abundance difference typically shows up as only {10**med:.1f}x more records.\n"
         "Judge on the interval, not the mean: American badger's mean of 0.33 sits in the healthy band while its "
         "interval spans both thresholds. Grey rows are unestimable and carry no value.",
         ha="center", va="bottom", fontsize=7.5, color="#5f6670")
fig.tight_layout(rect=(0, 0.07, 1, 1))
fig.savefig(HERE / "fig_theta1_congruence_fleet.png")
print("median", round(med, 3), "n", len(est))
