
import numpy as np, pandas as pd, matplotlib.pyplot as plt

CUR5 = ["White-tailed Deer","Eastern Gray Squirrel","Northern Raccoon","Coyote","American Black Bear"]
BLUE, ORANGE = "#1f6f8b", "#e08a3c"

def build_atlas1(D, dv2, aud3, final13, outfile="atlas1_sampling.png"):
    """Atlas 1: sampling design, effort, species coverage. Panels held in a dict
    so a later edit targets a NAMED panel, never a positional index."""
    apply_figure_style(sizes=(8,7,6))
    l48 = dv2[dv2.longitude.between(-125,-66) & dv2.latitude.between(24,49.5)]
    fig = plt.figure(figsize=(13.5, 9.8))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.25,.95,.95], hspace=.52, wspace=.30)
    P = {}
    P["a"] = fig.add_subplot(gs[0,:2]); P["b"] = fig.add_subplot(gs[0,2])
    P["c"] = fig.add_subplot(gs[1,0]);  P["d"] = fig.add_subplot(gs[1,1:])
    P["e"] = fig.add_subplot(gs[2,0]);  P["f"] = fig.add_subplot(gs[2,1:])

    ax = P["a"]; ax.set_facecolor("#f4f4f2")
    draw_states(ax)
    ax.scatter(l48.longitude, l48.latitude, s=1.4, c="#b8c4cc", alpha=.55, lw=0,
               label=f"{len(l48):,} deployments")
    ac = aud3[aud3.mean_lon.between(-125,-66) & aud3.mean_lat.between(24,49.5)]
    s_ = ax.scatter(ac.mean_lon, ac.mean_lat, s=np.clip(ac.n_deployments*.55,5,190),
                    c=ac.extent_km, cmap="viridis", alpha=.85, lw=.35, edgecolor="white", zorder=3)
    cb = fig.colorbar(s_, ax=ax, fraction=.028, pad=.012)
    cb.set_label("Array extent (km)", fontsize=6.5); cb.ax.tick_params(labelsize=6)
    ax.set_xlabel("Longitude (\u00b0)"); ax.set_ylabel("Latitude (\u00b0)"); ax.set_aspect(1.18)
    ax.set_title(f"a  {len(ac)} arrays (bubble = deployments, colour = extent; 25 km cap)", loc="left")
    ax.legend(fontsize=6, frameon=False, loc="lower left")

    ax = P["b"]
    ax.hist(aud3.n_deployments, bins=np.logspace(0, np.log10(aud3.n_deployments.max()), 34),
            color="#2c6e8f", alpha=.85)
    ax.set_xscale("log"); ax.set_xlabel("Deployments per array"); ax.set_ylabel("Arrays")
    med = aud3.n_deployments.median()
    ax.axvline(med, color="#c0392b", ls="--", lw=1)
    ax.annotate(f"median {med:.0f}", (med*1.35, ax.get_ylim()[1]*.72), fontsize=6.2, color="#c0392b")
    ax.set_title("b  Array size is highly skewed", loc="left")

    ax = P["c"]
    vd = l48.dropna(subset=["sd","ed"]); mo=[]
    for s0, e0 in zip(vd.sd, vd.ed):
        if e0 < s0: continue
        mo.extend(pd.period_range(s0, e0, freq="M").month)
    mc = pd.Series(mo).value_counts().sort_index()
    ax.bar(mc.index, mc.values/1000, color=["#c0392b" if m in (8,9,10) else "#b8c4cc" for m in mc.index])
    ax.set_xticks(range(1,13)); ax.set_xticklabels(list("JFMAMJJASOND"), fontsize=6)
    ax.set_xlabel("Month"); ax.set_ylabel("Camera-months (thousands)")
    ax.set_title("c  Aug\u2013Oct window (red) holds most effort", loc="left")

    ax = P["d"]
    sp_n = D.groupby("common_name").size().sort_values()
    ax.barh(np.arange(len(sp_n)), sp_n.values,
            color=BLUE)
    ax.set_yticks(np.arange(len(sp_n))); ax.set_yticklabels(sp_n.index, fontsize=6)
    ax.set_xscale("log"); ax.set_xlim(right=sp_n.max()*3.4)
    ax.set_xlabel("Aug\u2013Oct detections (log scale)")
    ax.set_title("d  Detections per species", loc="left")
    for i, v in enumerate(sp_n.values): ax.annotate(f"{v:,}", (v*1.18, i), va="center", fontsize=5.6)

    ax = P["e"]
    f13 = final13.sort_values("arrays25")
    ax.barh(np.arange(len(f13)), f13.arrays25,
            color=BLUE)
    ax.set_yticks(np.arange(len(f13))); ax.set_yticklabels(f13.species, fontsize=5.6)
    ax.axvline(60, color="#555", ls=":", lw=.9)
    ax.set_ylim(-1.6, len(f13)-.4)
    ax.annotate("60-array floor", (60, -1.35), fontsize=5.8, color="#555", ha="center")
    ax.set_xlabel("Arrays with \u226525 detections")
    ax.set_title("e  Usable arrays per species", loc="left")

    ax = P["f"]
    nn = D.groupby("common_name").night.mean().sort_values()*100
    ax.axvspan(0,20, color="#f6d98a", alpha=.35); ax.axvspan(80,100, color="#2c3e50", alpha=.12)
    ax.scatter(nn.values, np.arange(len(nn)), s=42,
               c=BLUE, zorder=3)
    for i,(s,v) in enumerate(nn.items()):
        right = v > 70
        ax.annotate(s, (v + (-2.6 if right else 2.6), i), va="center",
                    ha="right" if right else "left", fontsize=5.8)
    ax.set_yticks([]); ax.set_xlim(-2, 104); ax.set_ylim(-1, len(nn))
    ax.set_xlabel("% of detections at night")
    ax.annotate("diurnal (<20% night)", (5, len(nn)-.45), fontsize=6, ha="center", color="#8a6d1f")
    ax.annotate("nocturnal", (95, len(nn)-.45), fontsize=6, ha="center", color="#2c3e50")
    ax.set_title("f  Diel niche coverage of the 13-species set", loc="left")

    fig.suptitle("Data atlas 1 \u2014 sampling design, effort, and species coverage "
                 "(Aug\u2013Oct, lower 48, corrected clock)", fontsize=9.8, y=.985)
    fig.savefig(outfile, dpi=200, bbox_inches="tight")
    return fig, P
