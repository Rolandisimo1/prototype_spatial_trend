"""Which camera-side estimator carries more information about the trend, per species.

Occupancy and Royle-Nichols are not competing estimators so much as instruments
with different working ranges. Occupancy reads presence/absence, so its
information about log(lambda) collapses once psi saturates. RN reads detection
FREQUENCY across replicate windows, so it keeps discriminating well past that
point -- but it is weaker at low abundance, where frequency carries little.

This module locates the crossover from data we already have, rather than
assuming an r or a cutoff:

  1. Per species, from the raw camera detection histories, take two statistics
     among DETECTING deployments -- mean proportion of windows with a detection,
     and the fraction of deployments detecting in EVERY window. Restricting to
     detecting deployments makes these range-independent: a site outside the
     species' range simply never detects, so it drops out rather than diluting
     the estimate.
  2. Solve for (lambda, r) by matching those two statistics under the RN
     observation model. Two moments, two unknowns. The heterogeneity in
     detection frequency is what separates lambda from r; the mean alone cannot.
  3. At that (lambda, r, J), compare Fisher information about log(lambda) under
     each observation model and recommend the larger.

CAVEAT, stated because it bounds the conclusion: step 2 fits (lambda, r) under
RN's own assumptions (Poisson N, independent per-individual detection). So this
answers "given RN's data-generating model, which observation model extracts more
information about the trend" -- not "is RN's model true". It is a comparison of
instruments under a common reference, and a species whose real detection process
violates RN badly is not diagnosed here.
"""
import numpy as np
from scipy.optimize import least_squares
from scipy.stats import poisson


def _pmf_y(lam, r, J):
    """P(y windows with detection | lambda, r, J) under Royle-Nichols."""
    Nmax = int(lam + 10 * np.sqrt(max(lam, 1.0)) + 30)
    N = np.arange(Nmax + 1)
    pN = poisson.pmf(N, lam)
    q = 1.0 - (1.0 - r) ** N                      # per-window detection given N
    y = np.arange(J + 1)
    from scipy.special import comb
    B = comb(J, y)[None, :] * q[:, None] ** y[None, :] * (1 - q[:, None]) ** (J - y)[None, :]
    return np.clip(pN @ B, 1e-300, None)


def _moments(lam, r, J):
    """Model-predicted (mean window proportion, fraction firing in all windows),
    both conditional on the deployment detecting at least once."""
    p = _pmf_y(lam, r, J)
    y = np.arange(J + 1)
    denom = p[1:].sum()
    return float((y[1:] / J * p[1:]).sum() / denom), float(p[J] / denom)


def fit_lambda_r(mean_prop, frac_all, J):
    """Method-of-moments (lambda, r) from the two observed detection statistics."""
    def resid(theta):
        lam, r = np.exp(theta[0]), 1 / (1 + np.exp(-theta[1]))
        mp, fa = _moments(lam, r, J)
        return [mp - mean_prop, fa - frac_all]
    # Bounded, not unconstrained: an unbounded search walks lambda to e^35 and
    # the Poisson support blows up before the residual is ever evaluated.
    lo = [np.log(1e-2), np.log(1e-4 / (1 - 1e-4))]
    hi = [np.log(500.0), np.log(0.99 / 0.01)]
    best = None
    for lam0 in (0.5, 2.0, 8.0, 30.0):            # multistart: the surface is not convex
        for r0 in (0.02, 0.1, 0.3, 0.6):
            x0 = [np.clip(np.log(lam0), *[lo[0], hi[0]]),
                  np.clip(np.log(r0 / (1 - r0)), *[lo[1], hi[1]])]
            s = least_squares(resid, x0, bounds=(lo, hi), method="trf", max_nfev=4000)
            if best is None or s.cost < best.cost:
                best = s
    return float(np.exp(best.x[0])), float(1 / (1 + np.exp(-best.x[1]))), float(best.cost)


def info_rn(lam, r, J, h=1e-4):
    """Fisher information about eta = log(lambda) from the RN detection history."""
    p0 = _pmf_y(lam, r, J)
    pu, pd_ = _pmf_y(lam * np.exp(h), r, J), _pmf_y(lam * np.exp(-h), r, J)
    return float((((pu - pd_) / (2 * h)) ** 2 / p0).sum())


def info_occ(lam, r, J, h=1e-4):
    """Fisher information about eta from the occupancy presence/absence signal.

    Occupancy reaches lambda only through psi = 1 - exp(-lambda); its per-visit
    detection p is a free nuisance parameter and carries no information about
    lambda, so p is held fixed here. p is set to the detection probability given
    presence implied by the same (lambda, r), which makes the two instruments
    comparable rather than differently calibrated.
    """
    def psi_star(l):
        psi = 1 - np.exp(-l)
        p = (1 - np.exp(-lam * r)) / (1 - np.exp(-lam))     # fixed at the anchor point
        return psi * (1 - (1 - p) ** J)
    s0 = psi_star(lam)
    dd = (psi_star(lam * np.exp(h)) - psi_star(lam * np.exp(-h))) / (2 * h)
    return float(dd ** 2 / (s0 * (1 - s0)))


def crossover_lambda(r, J, lo=0.05, hi=200.0):
    """Smallest lambda at which RN's information exceeds occupancy's; NaN if never."""
    grid = np.exp(np.linspace(np.log(lo), np.log(hi), 900))
    diff = np.array([info_rn(l, r, J) - info_occ(l, r, J) for l in grid])
    idx = np.where(diff > 0)[0]
    return float(grid[idx[0]]) if len(idx) else np.nan


def recommend(mean_prop, frac_all, J):
    """Per-species recommendation with the quantities behind it."""
    lam, r, cost = fit_lambda_r(mean_prop, frac_all, J)
    i_o, i_r = info_occ(lam, r, J), info_rn(lam, r, J)
    return dict(lambda_hat=lam, r_hat=r, fit_cost=cost, info_occ=i_o, info_rn=i_r,
                info_ratio_rn_over_occ=i_r / i_o if i_o > 0 else np.inf,
                recommended="RN" if i_r > i_o else "occupancy",
                crossover_lambda=crossover_lambda(r, J))


def moment_grid(J, n_lam=140, n_r=100):
    """Precompute (lambda, r) -> moments and both informations, for one J.

    The grid exists because two moments at small J do NOT identify (lambda, r):
    e.g. at J=3, (lambda=25, r=0.03) and (lambda=0.11, r=0.52) reproduce the same
    pair to machine precision. Point-fitting therefore lands on an arbitrary spot
    along a ridge. What we actually need is whether the RECOMMENDATION is stable
    across the whole ridge, which requires enumerating it rather than optimising.
    """
    lam = np.exp(np.linspace(np.log(0.05), np.log(200.0), n_lam))
    rr = 1 / (1 + np.exp(-np.linspace(np.log(1e-3 / 0.999), np.log(0.95 / 0.05), n_r)))
    L, R = np.meshgrid(lam, rr, indexing="ij")
    mp = np.empty_like(L); fa = np.empty_like(L)
    io = np.empty_like(L); ir = np.empty_like(L)
    for i in range(L.shape[0]):
        for j in range(L.shape[1]):
            mp[i, j], fa[i, j] = _moments(L[i, j], R[i, j], J)
            io[i, j] = info_occ(L[i, j], R[i, j], J)
            ir[i, j] = info_rn(L[i, j], R[i, j], J)
    return dict(J=J, lam=L, r=R, mean_prop=mp, frac_all=fa, info_occ=io, info_rn=ir)


def recommend_robust(mean_prop, frac_all, J, grid, tol_mp, tol_fa):
    """Verdict over every (lambda, r) consistent with the observed moments.

    Returns the range of the information ratio across the consistent set. The
    verdict is only 'RN' or 'occupancy' when the ratio stays on one side of 1
    everywhere on the ridge; otherwise it is 'undetermined', which is an honest
    answer rather than a coin flip on an arbitrary fitted point.
    """
    ok = (np.abs(grid["mean_prop"] - mean_prop) <= tol_mp) & \
         (np.abs(grid["frac_all"] - frac_all) <= tol_fa)
    n = int(ok.sum())
    if n == 0:
        return dict(n_consistent=0, verdict="no consistent (lambda, r)",
                    ratio_min=np.nan, ratio_max=np.nan,
                    lam_min=np.nan, lam_max=np.nan)
    ratio = grid["info_rn"][ok] / np.clip(grid["info_occ"][ok], 1e-300, None)
    lo, hi = float(ratio.min()), float(ratio.max())
    verdict = "RN" if lo > 1 else ("occupancy" if hi < 1 else "undetermined")
    return dict(n_consistent=n, verdict=verdict, ratio_min=lo, ratio_max=hi,
                lam_min=float(grid["lam"][ok].min()), lam_max=float(grid["lam"][ok].max()),
                r_min=float(grid["r"][ok].min()), r_max=float(grid["r"][ok].max()))


def fig_detection_structure(obs, grids, out="fig_estimator_choice.png"):
    """Observed detection-frequency structure against what RN can produce.

    x: mean proportion of windows with a detection, among detecting deployments.
    y: fraction of detecting deployments that detect in EVERY window.
    Both are conditional on detecting, which makes them range-independent.

    The shaded band is the region reachable by Royle-Nichols over ALL (lambda, r)
    at that window count. A species outside the band cannot be represented by RN
    at any parameter value -- which is a statement about the observation model,
    not about the species' abundance.
    """
    import matplotlib.pyplot as plt
    Js = sorted(obs.J.unique())
    fig, axes = plt.subplots(1, len(Js), figsize=(3.5 * len(Js), 3.5), sharey=True)
    axes = np.atleast_1d(axes)
    for ax, j in zip(axes, Js):
        g = grids[j]
        ax.scatter(g["mean_prop"].ravel(), g["frac_all"].ravel(), s=1.0,
                   color="#b9bcbe", alpha=0.35, edgecolors="none",
                   label="reachable by RN\n(all $\\lambda$, $r$)", rasterized=True)
        o = obs[obs.J == j]
        ax.scatter(o.mean_prop, o.frac_all_obs, s=np.clip(o.n_dep / 40, 6, 70),
                   color="#cb4335", edgecolors="white", linewidths=0.4, zorder=3,
                   label="observed species")
        for _, rr in o.nlargest(3, "n_dep").iterrows():
            ax.annotate(rr.species.replace("_", " "), (rr.mean_prop, rr.frac_all_obs),
                        textcoords="offset points", xytext=(-6, 6), fontsize=5.6,
                        ha="right", color="#22282e")
        ax.set_title(f"{j} ten-day windows", fontsize=8.4)
        ax.set_xlabel("Mean fraction of windows\nwith a detection")
        ax.set_xlim(0.15, 0.9); ax.set_ylim(-0.03, 0.72)
    axes[0].set_ylabel("Fraction of deployments detecting\nin every window")
    axes[0].legend(frameon=False, fontsize=6.0, loc="upper left")
    n_above = int((obs.excess > 0).sum()); n_below = int(((~obs.inside) & (obs.excess <= 0)).sum())
    fig.suptitle(f"Camera detections are more polarised than Royle-Nichols can produce "
                 f"({n_above} of {len(obs)} species-by-window cases)", fontsize=9.4, y=1.03)
    fig.text(0.5, -0.10,
             f"{n_above} of {len(obs)} cases sit ABOVE the reachable band: given how often they "
             "are detected on average, more deployments detect in every single window\nthan "
             "Royle-Nichols allows for any abundance and per-individual detectability. "
             f"{n_below} sit BELOW it -- all but one are species that never detected in every "
             "window\nagainst an RN ceiling of 0.003 or less, i.e. sampling zeros rather than "
             "evidence of the opposite pattern (white-tailed antelope squirrel at 5 windows is "
             "the\none real exception, 0.049 below). Point size is the number of deployments. "
             "Only deployments with exactly the panel's window count are used, since a shorter "
             "deployment\nfires in 'every window' more easily.",
             ha="center", va="top", fontsize=6.4, color="#5f6a73")
    fig.savefig(out, bbox_inches="tight", dpi=300)
    plt.close(fig)
    return out
