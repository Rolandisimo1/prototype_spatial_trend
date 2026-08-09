# Pending Hazel extraction — moose full report (single remaining item)

All spatial CAR fields, mu/psi reconstruction, and R-hat-by-family items
from the original queued request have been resolved using data already on
disk (see moose_model_full_report.md v2). One item remains, blocked on
Hazel access:

**theta0 and theta1 full posterior summary, both models.**
Per Goldstein et al. (bioRxiv 2025.01.17.633640), theta1 is their primary
camera/iNat congruence metric. Our model's theta1 plays the identical
structural role (`log(mu) = theta0 + theta1*log(sum(lambda))` in
`calcIntensity_SVC`). R-hat is already extracted (moose_v1fix_*_rhat_by_family.csv,
iNat family: 1.021 national-scalar, 1.015 ecoregion) but posterior
mean/median/sd/CI is not.

Extraction script ready to run: extract_moose_v1fix_theta.R (adapted from
extract_moose_v1fix_posterior.R's load_chains() pattern). Both nodes are
directly monitored scalars in the same samples matrices already read for
the CAR field / R-hat extraction — no new Hazel connection type needed,
just execution.

Blocked this session: ssh:hazel target returned "Permission denied
(keyboard-interactive)" when probed. Needs Claude Code (has working Hazel
SSH access) or the user to run this script directly and report back the
output CSV (moose_v1fix_theta_posterior_combined.csv).
