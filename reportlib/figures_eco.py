"""Is the ecoregion parameterisation worth its cost?

The two parameterisations answer the same question with different assumptions.
The national model estimates one slope for the whole country. The ecoregion
model estimates that slope as the mean of eight regional deviations drawn from a
shared distribution, so the national number becomes a hyper-mean and inherits
the uncertainty in how much regions differ.

That trade has to be measured, not assumed. Cost is the loss of precision on the
national estimate. Benefit is regional inference -- but only if the regional
estimates are actually distinguishable from each other. A hierarchy that returns
eight copies of the national number has bought nothing but a wider interval.
"""
import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from . import inputs, style

SPECIES_ORDER = ["Bobcat", "White-tailed deer", "Moose"]
SPECIES_COLOR = {"Bobcat": "#b2182b", "White-tailed deer": "#1b5e9c", "Moose": "#4d9221"}

NATIONAL = {"Bobcat": ("pull", "bobcat_v2b_national_scalar_posterior.csv"),
            "White-tailed deer": ("pull", "white-tailed_deer_v2b_national_scalar_posterior.csv"),
            "Moose": ("repo", "moose_v2b_national_scalar_posterior_single.csv")}
ECO_GLOBAL = {"Bobcat": ("pull", "bobcat_v2b_ecoregion_global.csv"),
              "White-tailed deer": ("pull", "white-tailed_deer_v2b_ecoregion_global.csv"),
              "Moose": ("repo", "moose_v2b_ecoregion_global_single.csv")}
ECO_REGIONAL = {"Bobcat": ("pull", "bobcat_v2b_ecoregion_ecoregion_posterior.csv"),
                "White-tailed deer": ("pull",
                                      "white-tailed_deer_v2b_ecoregion_ecoregion_posterior.csv"),
                "Moose": ("repo", "moose_v2b_ecoregion_posterior_single.csv")}


def _dir(tag):
    return inputs.PULL if tag == "pull" else inputs.REPO


def ecoregion_value_table():
    """Cost and benefit of the regional layer, per species.

    Benefit is counted three ways, because "the model produced regional numbers"
    is not the same as "the regional numbers say anything":
      regions_clear           regions whose own interval excludes zero
      distinguishable_pairs   pairs of those regions whose intervals do not overlap
      opposing_national       regions that clearly move against the national direction
    """
    rows = []
    for sp in SPECIES_ORDER:
        nat = pd.read_csv(os.path.join(_dir(NATIONAL[sp][0]),
                                       NATIONAL[sp][1])).set_index("parameter")
        eco = pd.read_csv(os.path.join(_dir(ECO_GLOBAL[sp][0]),
                                       ECO_GLOBAL[sp][1])).set_index("parameter")
        reg = pd.read_csv(os.path.join(_dir(ECO_REGIONAL[sp][0]), ECO_REGIONAL[sp][1]))

        nat_mean = nat.loc["total_var_beta", "mean"]
        nw = nat.loc["total_var_beta", "q975"] - nat.loc["total_var_beta", "q025"]
        ew = eco.loc["total_var_beta", "q975"] - eco.loc["total_var_beta", "q025"]

        reg["clear"] = (reg.abs_trend_q025 > 0) | (reg.abs_trend_q975 < 0)
        clear = reg[reg.clear]
        pairs = sum(1 for i in range(len(clear)) for j in range(i + 1, len(clear))
                    if (clear.iloc[i].abs_trend_q975 < clear.iloc[j].abs_trend_q025
                        or clear.iloc[j].abs_trend_q975 < clear.iloc[i].abs_trend_q025))

        rows.append({
            "species": sp,
            "nat_mean": nat_mean, "nat_width": nw,
            "nat_clear": bool(nat.loc["total_var_beta", "q025"] > 0
                              or nat.loc["total_var_beta", "q975"] < 0),
            "eco_mean": eco.loc["total_var_beta", "mean"], "eco_width": ew,
            "eco_clear": bool(eco.loc["total_var_beta", "q025"] > 0
                              or eco.loc["total_var_beta", "q975"] < 0),
            "width_ratio": ew / nw,
            "sigma_region": eco.loc["sigma_region", "mean"],
            "sigma_q025": eco.loc["sigma_region", "q025"],
            "sigma_q975": eco.loc["sigma_region", "q975"],
            "n_regions": len(reg),
            "regions_clear": int(reg.clear.sum()),
            "regions_no_data": int(((reg.n_camera_sites == 0)
                                    & (reg.n_inat_cells == 0)).sum()),
            "distinguishable_pairs": pairs,
            "pairs_possible": len(clear) * (len(clear) - 1) // 2,
            "opposing_national": int((reg.clear & (np.sign(reg.abs_trend_mean)
                                                   != np.sign(nat_mean))).sum()),
        })
    return pd.DataFrame(rows)


def opposing_regions():
    """The specific regions that move against their species' national direction.

    These are the cases the national model cannot produce at all, so they are the
    concrete return on the regional layer.
    """
    out = []
    for sp in SPECIES_ORDER:
        nat_mean = pd.read_csv(os.path.join(_dir(NATIONAL[sp][0]), NATIONAL[sp][1])
                               ).set_index("parameter").loc["total_var_beta", "mean"]
        reg = pd.read_csv(os.path.join(_dir(ECO_REGIONAL[sp][0]), ECO_REGIONAL[sp][1]))
        reg["clear"] = (reg.abs_trend_q025 > 0) | (reg.abs_trend_q975 < 0)
        opp = reg[reg.clear & (np.sign(reg.abs_trend_mean) != np.sign(nat_mean))]
        for _, r in opp.iterrows():
            out.append({"species": sp, "region": r.region_name,
                        "regional_trend": r.abs_trend_mean,
                        "q025": r.abs_trend_q025, "q975": r.abs_trend_q975,
                        "national_trend": nat_mean,
                        "n_camera_sites": int(r.n_camera_sites),
                        "n_inat_cells": int(r.n_inat_cells)})
    return pd.DataFrame(out)


def fig_ecoregion_value(out="fig_ecoregion_value.png"):
    """Cost, whether regional structure exists, and what the regional layer returns."""
    t = ecoregion_value_table().set_index("species").loc[SPECIES_ORDER]
    y = np.arange(len(SPECIES_ORDER))[::-1]
    fig, axes = plt.subplots(1, 3, figsize=(10.8, 3.2))

    # Panel 1 -- cost. Each interval is drawn at its own estimate, so the panel
    # shows both the widening and the shift in the point estimate.
    ax = axes[0]
    for yi, sp in zip(y, SPECIES_ORDER):
        r = t.loc[sp]
        c = SPECIES_COLOR[sp]
        ax.plot([r.nat_mean - r.nat_width / 2, r.nat_mean + r.nat_width / 2],
                [yi + 0.17] * 2, color=c, lw=2.2, solid_capstyle="round")
        ax.plot(r.nat_mean, yi + 0.17, "o", color=c, ms=4.4)
        ax.plot([r.eco_mean - r.eco_width / 2, r.eco_mean + r.eco_width / 2],
                [yi - 0.17] * 2, color=c, lw=2.2, alpha=0.40, solid_capstyle="round")
        ax.plot(r.eco_mean, yi - 0.17, "o", color=c, ms=4.4, alpha=0.40)
        ax.annotate(f"{r.width_ratio:.1f}x wider", (r.eco_mean + r.eco_width / 2, yi - 0.17),
                    textcoords="offset points", xytext=(5, 0), fontsize=6.1,
                    color=c, va="center")
    ax.axvline(0, color="0.35", lw=0.8, ls="--")
    ax.set_yticks(y); ax.set_yticklabels(SPECIES_ORDER, fontsize=7)
    ax.set_xlabel("National trend (95% interval)")
    ax.set_title("Cost: national precision")
    ax.plot([], [], color="0.30", lw=2.2, label="national model")
    ax.plot([], [], color="0.30", lw=2.2, alpha=0.40, label="ecoregion model")
    ax.legend(frameon=False, fontsize=6.2, loc="lower left")
    ax.set_ylim(-0.85, len(SPECIES_ORDER) - 0.30)

    # Panel 2 -- is there regional structure at all?
    ax = axes[1]
    for yi, sp in zip(y, SPECIES_ORDER):
        r = t.loc[sp]
        ax.plot([r.sigma_q025, r.sigma_q975], [yi] * 2, color=SPECIES_COLOR[sp],
                lw=2.2, solid_capstyle="round")
        ax.plot(r.sigma_region, yi, "o", color=SPECIES_COLOR[sp], ms=5.0)
    ax.axvline(0, color="0.35", lw=0.8, ls="--")
    ax.set_yticks(y); ax.set_yticklabels([])
    ax.set_xlabel("Between-region spread\n(SD of regional deviations)")
    ax.set_title("Is there structure to find?")
    ax.set_xlim(-0.02, None)
    ax.set_ylim(-0.85, len(SPECIES_ORDER) - 0.30)

    # Panel 3 -- what the regional layer actually returns.
    ax = axes[2]
    w = 0.26
    metrics = [("regions_clear", "regions with a\nclear direction"),
               ("distinguishable_pairs", "region pairs that\ndiffer from each other"),
               ("opposing_national", "regions opposing\nthe national trend")]
    xm = np.arange(len(metrics))
    for k, sp in enumerate(SPECIES_ORDER):
        ax.bar(xm + (k - 1) * w, [t.loc[sp, m] for m, _ in metrics], w,
               color=SPECIES_COLOR[sp], zorder=3, label=sp)
    ax.set_xticks(xm); ax.set_xticklabels([lab for _, lab in metrics], fontsize=6.0)
    ax.set_ylabel("Count")
    ax.set_title("What it returns")
    ax.legend(frameon=False, fontsize=6.2, loc="upper right")

    fig.text(0.5, -0.11,
             "Left: 95% interval on the national trend under each parameterisation, each drawn at "
             "its own estimate. Middle: standard deviation of the regional\n"
             "deviations; an interval running down to zero means no regional structure was "
             "resolved. Right: regional inference actually delivered, counting only\n"
             "regions whose own interval excludes zero. All eight ecoregions are available to "
             "every species.",
             ha="center", va="top", fontsize=6.3, color="#5f6a73")
    fig.tight_layout()
    fig.savefig(out, bbox_inches="tight", pad_inches=0.06, dpi=300)
    plt.close(fig)
    return out, t
