# submit_guard

The only sanctioned way to submit single-shot iSDM chains on Hazel. It lives at `~/isdm/fleet_v2u/submit_guard/`.

## Why

Until 2026-09-19, every chain asked for 110 GB of memory. The highest peak any chain has ever used is 32.5 GB. Hazel bills memory at one unit per 4 GB, so a 110 GB chain counted as 28 units of usage, 27 of them memory. That drove our fair-share factor down to 0.025 on a 0–1 scale, and fair share is what made chains wait 30–40 h. The QOS CPU cap and cluster load were not the cause. When the three pending coyote chains were cut from 110 GB to 48 GB, all three started within a minute, about 15 h ahead of their scheduled start.

## Rules

| | Rule |
|---|---|
| Size class | small if cell50 ≤ 1.15 × the largest measured small-class cell50 (moose-sized, currently 526); otherwise large |
| Memory | small 16 GB, large 48 GB. The script refuses if the request is more than 2.0× the class's largest measured peak, less than 1.3× it (out-of-memory risk), or if cell50 is more than 1.15× the largest measured cell50 in the class |
| Wall | 2 × cell50 × (hours per cell50). The rate comes from completed chains of the same model kind: the same species if one exists, otherwise the slowest species |
| QOS | `normal` if the wall is ≤ 96 h, otherwise `long` (refuses above 240 h). Both have priority 0, so the QOS only changes the caps |

## Files

- `chain_calibration.csv`: one row per finished chain, with cell50, nsite, requested and peak memory, wall and state. It was seeded from sacct on 2026-09-19. Chunked v1fix runs contribute peak memory only, because their walls were split across chunks.
- `chain_submissions.csv` (on Hazel): one row per submitted or registered chain, recording requested memory and wall. `harvest` moves finished chains into the calibration table, and it runs automatically before every `plan` or `submit`.
- `bundle_meta.R`: reads ncell50 and nsite from a bundle. The result is cached per bundle in `<run_dir>/bundle_meta.json`, keyed on the bundle's size and mtime.

## Usage

```
python3 submit_chains.py plan   <run_dir> [--chains 1,2,3]   # dry run
python3 submit_chains.py submit <run_dir> [--chains 1,2,3]
python3 submit_chains.py harvest
python3 submit_chains.py register <run_dir> <chain> <jobid> <mem_gb> <wall_h> <qos>
```

## Pending

- Wrappers should refuse to start unless `SUBMIT_GUARD=1` (which the guard exports), so a raw `sbatch` can't bring back 110 GB. This has not yet been applied on Hazel.
- The large class spans roughly 700–3,300 cell50. Tier A species with mid-sized ranges will be over-requested at 48 GB. Add a middle class once measured points exist between 600 and 2,000 cell50.
