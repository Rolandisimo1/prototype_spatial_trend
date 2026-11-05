#!/usr/bin/env python3
"""Fleet convergence table (2026-11-01): one row per run from the per-run extraction CSVs.

usage: fleet_convergence.py <posteriors_dir> [out_csv]

For each of year_beta, year_var, total_var_beta, theta0, theta1: R-hat, bulk ESS, tail ESS (pooled over all chains),
and two counts of chains whose mean sits off the bulk:
  off_across  = chains with |chain mean - mean of chain means| > 2 x SD of chain means   (the brief's definition).
                With n chains a single chain can sit at most (n-1)/sqrt(n) SDs out: 3.18 for 12 chains, 1.15 for 3,
                so this is ALWAYS 0 for 3-chain runs, and several off-mode chains inflate the SD and mask each other.
  off_within  = chains with |chain mean - median of chain means| > 2 x median within-chain SD. Robust to how many
                chains are off and usable with 3 chains.
Prior-only parameters (camonly runs) are left blank: their R-hat describes prior sampling.
"""
import glob, os, sys
import pandas as pd

PARAMS = ["year_beta", "year_var", "total_var_beta", "theta0", "theta1"]


def main(d, out):
    rows = []
    for cf in sorted(glob.glob(os.path.join(d, "convergence_*.csv"))):
        run = os.path.basename(cf)[len("convergence_"):-4]
        c = pd.read_csv(cf).set_index("param")
        pf = os.path.join(d, "perchain_%s.csv" % run)
        pc = pd.read_csv(pf) if os.path.exists(pf) else None
        r0 = c.iloc[0]
        row = {k: r0[k] for k in ["run", "species", "mask", "window", "camonly", "affected_by_resume_defect",
                                  "n_chains", "n_draws"]}
        for p in PARAMS:
            if p not in c.index or bool(c.loc[p, "prior_only"]):
                for s in ["rhat", "ess_bulk", "ess_tail", "off_across", "off_within"]:
                    row["%s_%s" % (p, s)] = None
                continue
            row["%s_rhat" % p] = c.loc[p, "rhat"]
            row["%s_ess_bulk" % p] = c.loc[p, "ess_bulk"]
            row["%s_ess_tail" % p] = c.loc[p, "ess_tail"]
            if pc is None:
                row["%s_off_across" % p] = row["%s_off_within" % p] = None
                continue
            q = pc[pc.param == p]
            m, s = q["mean"], q["sd"]
            row["%s_off_across" % p] = int(((m - m.mean()).abs() > 2 * m.std()).sum()) if m.std() > 0 else 0
            row["%s_off_within" % p] = int(((m - m.median()).abs() > 2 * s.median()).sum())
        rows.append(row)
    f = pd.DataFrame(rows)
    f.to_csv(out, index=False)
    print("wrote %s: %d runs" % (out, len(f)))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(__doc__); sys.exit(1)
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else os.path.join(sys.argv[1], "fleet_convergence.csv"))
