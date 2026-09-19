# submit_guard

The only sanctioned way to submit single-shot iSDM chains on Hazel. It lives at `~/isdm/fleet_v2u/submit_guard/`.

## Why

Until 2026-09-19, every chain asked for 110 GB of memory. The highest peak any chain has ever used is 32.5 GB. Hazel bills memory at one unit per 4 GB, so a 110 GB chain counted as 28 units of usage, 27 of them memory. That drove our fair-share factor down to 0.025 on a 0–1 scale, and fair share is what made chains wait 30–40 h. The QOS CPU cap and cluster load were not the cause.

### Evidence: the coyote chains, 2026-09-19

Coyote chains 864015–864017 were submitted at 110 GB on 9/17 at 15:00. After 41 h in the queue, the scheduler estimated they would start at 9/19 23:11, 9/20 03:45 and 9/20 06:28. At 08:14 on 9/19 their memory was lowered in place to 48 GB with `scontrol update`, and their wall was trimmed from 7 days to 142 h. The jobs were not resubmitted, so they kept their eligible time (9/17 15:00). **All three started at 08:14:23, within a minute of the change and 15–22 h ahead of the estimates**, on different nodes from the ones reserved for them. Nothing else changed. Over the same days, jobs requesting less than 100 GB started with a median wait of 0.02 h, while 110 GB chains waited a median of 35.8 h. Memory was the constraint, acting through our fair-share priority.

**When resizing a pending job, change it in place (`scontrol update`); don't cancel and resubmit.** Resubmitting resets the eligible time and throws away the queue age accumulated so far.

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

- **Applied 2026-09-19:** every `run_single*_sbatch.sh` exits with code 2 unless `SUBMIT_GUARD=1` is set, and the guard exports that variable, so a raw `sbatch` can't bring back 110 GB. The check covers the moose, bobcat, deer and coyote `_v2u_national_scalar` wrappers and the deer v2b count-arm wrapper. Each was backed up as `*.bak_20260919_preguard`. Any run dir built by copying one of these inherits the check.
- The large class spans roughly 700–3,300 cell50. Tier A species with mid-sized ranges will be over-requested at 48 GB. Add a middle class once measured points exist between 600 and 2,000 cell50.
