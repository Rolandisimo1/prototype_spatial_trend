Hazel update, Mon 12 Oct 2026, 20:30 EDT. Covers the work since 10 Oct. The support email was never sent.

**1. Needs a decision: the scheduler will not start any long fit before 21 Oct.**
- Slurm now gives a start estimate for our long jobs. All 24 show the same time, 2026-10-21 08:02, whatever their wall
  limit, from 1-16:00 (snowshoe hare, camera-only) to 5-21:00 (bobcat camdev).
- `scontrol show reservation` is still empty. The Hazel news and outage pages list nothing for October.
- 605 CPUs are idle. Short jobs do run in the end: the three diagnostics waited 34 h, then started on 11 Oct at
  07:29 and ran in 2–27 min.
- One shared start time for every wall limit suggests something we can't see holds the nodes until 21 Oct 08:02.
  That could be maintenance or a node set held back.
- Consequence: the bobcat camdev control (about 6 days a chain at the guard's 2x wall) can't finish before
  about 27 Oct. Neither can the remaining camera-only fits (badger, porcupine, elk, opossum, puma, skunk, snowshoe
  hare).
- Options:
  (a) Send the support email. The draft is at prototype_spatial_trend/draft_hazel_support_email.md and needs this
      21 Oct estimate added.
  (b) Accept the delay.
  Nothing has been cancelled or resubmitted.

**2. Western gray squirrel camdev test: the prediction held, with one caveat.**
| fit | camera %/yr [95% CI] | CI width (pp) | R-hat |
|---|---|---|---|
| joint (shared surfaces) | −13.22 [−20.51, −5.00] | 15.5 | 1.002 |
| camera-only | −2.68 [−12.42, +6.63] | 19.1 | 1.000 |
| camdev (camera-side deviation added) | −1.41 [−10.98, +7.53] | 18.5 | 1.001 |
- The trend moved from −13.2 to −1.4 %/yr, to the camera-only side, and the interval widened from 15.5 to 18.5 pp.
  Read as: the cameras can't tell whether western gray squirrels are declining. The diagnosis (iNaturalist-shaped
  shared surface, year_beta absorbing the mismatch) holds.
- Surface correlations across the 37 camera cells:
  - joint vs camera-only intercept: baseline 0.08.
  - camdev effective surface (shared + deviation) vs camera-only: 0.43.
  - camdev shared part vs camera-only: −0.05.
  - camdev shared vs shared + deviation: 0.58.
- **Caveat:** sigma_cam_dev = 1.86 [1.55, 2.00], under a dunif(0, 2) prior. The posterior presses against the upper
  bound, so the prior sets the deviation's size, not the data. The year_beta result may not depend on this, since
  the trend already sits next to the camera-only value. Still, the bound is binding. Widening it would be a model
  change, so it waits for Roland.
- Bobcat, the control, has not started (see 1). Per Roland's instruction nothing else camdev runs.

**3. Fisher: the trim hypothesis is wrong.**
- The build fails with the same dOcc_v dimension error whether y/yday are as-is (4673x10, max J 9) or trimmed to
  9 columns.
- Next step, queued (job 1220813): a read-only audit of every v2u/v2un bundle on the shape facts NIMBLE checks.
  These are J class/NA/zero/non-integer, y class and storage mode, yday location and dims, all-NA rows, and rows
  whose non-NA count differs from J.
- The same job answers the zero-night question: whether any fitted bundle has J = 0 or all-NA camera rows.
- Fisher stays unsubmitted until the audit names a cause.

**4. Numbers from the trend table (33 species).**
- **Elk refit** (corrected camera record): camera +4.25 %/yr [−1.13, +10.56], down from +6.59. iNat +2.44
  [−0.45, +5.71]. R-hat 1.007. The interval now includes zero. Elk's camera-only twin is queued (see 1). The
  clean-21 compression row will be recomputed when elk's camera-only twin finishes.
- **Tier C:**
  | species | camera %/yr | iNat %/yr | trend R-hat |
  |---|---|---|---|
  | groundhog | +2.70 [−0.56, +6.15] | +3.70 [+2.03, +5.62] | 1.003 |
  | grey wolf | −15.08 [−27.66, −1.29] | −6.79 [−14.44, +0.35] | 1.000 |
  | kit fox | −2.20 [−15.49, +9.80] | +4.67 [+0.10, +10.10] | 1.013 |
  | pronghorn | −16.13 [−23.82, −7.83] | +0.82 [+0.27, +1.42] | 1.018 |
  | rock squirrel | −14.49 [−24.62, −2.86] | −0.48 [−2.12, +1.16] | 1.002 |
  | white-tailed antelope squirrel | −11.61 [−21.33, −0.94] | +1.28 [+0.21, +2.48] | 1.005 |
  | golden-mantled ground squirrel | NOT converged | NOT converged | 1.69 |
  | brown bear | NOT converged | NOT converged | 2.89 |
  | fisher | no fit (see 3) | | |
- Pronghorn, rock squirrel and antelope squirrel show the western gray squirrel pattern: a steep camera decline
  beside a flat iNat trend. The same mechanism may be at work. None of the Tier C species has a camera-only twin
  yet; those were never requested.
- Fisher's camera-record flag (effectively 2019–2025; no 2009 data; zero detections in 2008, 2010 and 2014–2018)
  goes in when fisher has a fit.

**5. Decomposition (joint − camera-only shift).**
- 19 species summarised. Of those, 16 pass the convergence screen. Brown bear, peccary and California ground
  squirrel are excluded because their joint trend isn't converged.
- field_cor_intercept, the correlation between the joint and camera-only intercept surfaces: median 0.72 across
  the 16 (range 0.42–0.86 without western gray squirrel). Western gray squirrel is the outlier at 0.08, with 37
  camera cells.
- Residual sign: 10 negative, 6 positive. By the agreed rule this is scattered, not systematic.
- The trend table ran alongside the sweep and missed 13 of these 19 summaries, so its decomposition columns read
  "pending". A refresh is queued (1220808).

**6. Queued, not yet reported:**
- t1init_check (1220807): per-chain theta1 for porcupine, badger and brown bear, original fit vs inits
  0.3/0.6/0.9. If the split persists, the reported result is "integration not identified for this species".
- Table refresh (1220808).
- Bundle audit (1220813).
