# Open decisions on the integrated camera-trap and iNaturalist trend model

A register of every decision currently open, what each is waiting on, and what
turns on it. Current as of this entry. Twenty-five entries: one needs a decision
now, three are blocked on running jobs, three were specified and never came back,
four are settled and recorded here so they are not reopened, and the rest are open
without a blocker.

## Take-home points

1. **One decision is worth making now.** Whether to drop the camera-only control
   arm for the fifteen species whose full-record fit resolved no direction (W2),
   which is only useful before those jobs run. Reporting campaigns (M2) are now
   decided: exclude the affected cells.
2. **The headline window decision gates most of the rest.** Which of the
   five-year, ten-year or full-record trend becomes the reported result (W1)
   determines the agency comparison, the congruence gate and every trend figure.
   It is blocked on fits that need about nine more days.
3. **Three sensitivity checks were specified and never returned**, and all three
   are named as pending in the supplement.
4. **Six items are documentation or naming defects** that do not change a number
   but will mislead a reader, including one figure that disagrees with its own
   source table and should not circulate.
5. **One proposed test could settle the campaign question** without new fits, by
   asking whether the reporting step stops at the campaign's administrative
   boundary.

## Decisions needing a call now

### Ready to decide

W2 is not blocked; it expires once the jobs run. (M2, organised reporting campaigns, is now decided - exclude the cells; see the settled section.)

**W2. Drop the camera-only control arm for the 15 species whose full-record fit resolved no direction**

- Status: open - needs your call now. Decided by: Roland.
- Why it matters: Saves 5,409 chain-hours, about 2.6 days; decision is only useful before those jobs run
- Evidence: table_fleet_progress.csv


## Blocked on running work

### Waiting on the cluster

**W1. Which window becomes the headline trend: 5-year, 10-year or full record**

- Status: open - blocked. Decided by: Roland.
- Waiting on: Windowed fits still running, about 9 days at current slot count
- Why it matters: Everything downstream re-scores against it: agency comparison, congruence gate, every trend figure
- Evidence: table_trend_window_options.csv; table_camera_effort_by_year.csv

**V5. Re-score the agency concordance window-matched on both sides**

- Status: open - blocked. Decided by: Claude Code.
- Waiting on: Needs the windowed model fits
- Why it matters: Current comparison sets a full-record model trend against 5- and 10-year agency windows
- Evidence: table_model_vs_agency_concordance_by_window.csv

**F3. Eight species still lacking a camera-only control fit**

- Status: open - blocked. Decided by: Claude Code.
- Waiting on: Queued behind the windowed fleet
- Why it matters: Includes the three steep decliners whose trends are untested on both diagnostics
- Evidence: table_fleet_progress.csv


## Specified but never returned

### Queued or dispatched, no result delivered

**V2. Urbanisation interaction sensitivity check**

- Status: open - never returned. Decided by: Claude Code.
- Waiting on: Specified and queued, no result
- Why it matters: Named in the supplement as pending
- Evidence: table_sensitivity_inventory.csv

**V3. Camera deviation term sensitivity check**

- Status: open - never returned. Decided by: Claude Code.
- Waiting on: Specified and queued, no result
- Why it matters: Named in the supplement as pending
- Evidence: table_sensitivity_inventory.csv

**V4. Project fixed effect sensitivity check**

- Status: open - never returned. Decided by: Claude Code.
- Waiting on: Specified and queued, no result
- Why it matters: Named in the supplement as pending
- Evidence: table_sensitivity_inventory.csv


## Open, not blocked

### No dependency; these are waiting on attention rather than on data

**M1. Add a camera-location grouping term so repeat deployments at one site are not independent**

- Status: open - awaiting collaborators. Decided by: Arielle and Krishna.
- Waiting on: Their account of whether it was weighed and dropped, or lost when the temporal structure went in
- Why it matters: 97.8% of the year contrast is between locations, so the trend is identified by comparing places rather than revisits; intervals may be too narrow
- Evidence: REPORT_camera_location_effect.docx

**V1. Extend the fitted-occupancy saturation check beyond the three species it covers**

- Status: open - gap. Decided by: Roland.
- Why it matters: Reads existing posteriors, needs no new fits; currently the strongest saturation claim rests on 3 species
- Evidence: table_occupancy_saturation.csv

**V6. Spatial discontinuity test at campaign boundaries**

- Status: proposed - not run. Decided by: Roland to approve.
- Why it matters: The only test identified that could separate campaign reporting from local ecology given the camera gap; needs no new fits
- Evidence: this session

**S1. Confirm the metrics now score the reported trend parameter rather than the summed one**

- Status: open - needs confirmation. Decided by: Claude Code.
- Waiting on: Fix was specified; Tier 1 completed but the scoring change has not been confirmed in the delivered rows
- Why it matters: If unfixed, Tier 1 validates a parameter the paper does not report
- Evidence: table_simulation_coverage_gaps.csv

**S2. Whether the data-volume corner of Tier 1 can be reported**

- Status: open. Decided by: Roland.
- Waiting on: Up to 74% of replicates below an effective sample size of 100 in the sparse-camera corner
- Why it matters: Range-size and window power curves are reportable; the volume corner is not, and it is the axis the prospects claim needs
- Evidence: Tier 1 interim summary

**D1. Report what agency assessments actually measure in the manuscript**

- Status: open. Decided by: Roland.
- Why it matters: Method text exists for every record and appears nowhere in the manuscript; a minority rest on population estimates
- Evidence: agency_trends_combined.csv; MEMO_agency_assessment_caveats.docx

**D2. Record that one jackrabbit index pools two species, in the concordance notes**

- Status: open. Decided by: Claude Code.
- Why it matters: That record is the sole disagreement in the comparison and does not measure the species it is filed under
- Evidence: MEMO_agency_assessment_caveats.docx

**D3. Commit a canonical builder script for the concordance table**

- Status: open. Decided by: Claude Code.
- Why it matters: The published table has no builder anywhere; one dropped row prevents exact reproduction
- Evidence: table_model_vs_agency_concordance_by_window.csv

**D4. Rename the camera_anchored column to say what it computes**

- Status: open. Decided by: Claude Code.
- Why it matters: The name implies data sufficiency while the column tests interval containment; it has already caused one misreading
- Evidence: table_fleet_progress.csv

**D5. Fix fig_congruence.png, which disagrees with its own source table on bobcat**

- Status: open - do not circulate. Decided by: Claude Code.
- Why it matters: Flagged during the supplement work and still unfixed
- Evidence: fig_congruence.png

**D6. Correct any text asserting that all fits converged**

- Status: open. Decided by: Roland.
- Why it matters: Two of eight fits exceed an R-hat of 1.01 on the trend family; for bobcat the 1.1145 maximum belongs to a derived indicator we do not report, not to the trend parameter
- Evidence: table_convergence_and_trend.csv

**F1. Fisher: never fitted, and its R-hat check is still queued**

- Status: open. Decided by: Claude Code.
- Waiting on: Proposed to run off-cluster rather than queued behind our own chains
- Why it matters: 1 of 34 species has no trend at all
- Evidence: table_fleet_progress.csv

**F2. Five species whose fits did not converge**

- Status: open. Decided by: Roland.
- Why it matters: 28 of 34 converged; the remaining 5 carry no reportable trend
- Evidence: table_fleet_progress.csv

**P1. Freeze a numbered release and run expert review against it, rather than against the live model**

- Status: proposed. Decided by: Roland, with Ben.
- Why it matters: Ben cannot review a model that turns over weekly; a workshop on a live model spends its time on comprehension
- Evidence: this session


## Settled, recorded so they are not reopened

### These were closed on explicit reasoning and should not be relitigated

**M3. Let the record-to-abundance exponent vary freely by region**

- Status: decided - rejected. Decided by: Roland.
- Why it matters: Not identifiable against the abundance change we are measuring; closed on your objection
- Evidence: this session

**M4. Whether to correct the trend for drift in the stream-coupling exponent**

- Status: decided - no correction. Decided by: Roland.
- Why it matters: Reported as a diagnostic with three named exceptions, chiefly the opportunistic intensity surfaces
- Evidence: camera-inat-stream-agreement skill; table_theta_time_drift.csv

**M5. Adopt the centred parameterisation of the intensity intercept (decouple theta0 from theta1) fleet-wide**

- Status: decided - adopt, rolling out. Decided by: Roland.
- Why it matters: theta0 and theta1 traded off along an identifiability ridge, so the stream-coupling exponent theta1 failed to converge in 67 of 68 windowed joint fits and in the older 3-chain full-record fits, blocking every reportable joint trend and the congruence gate. Centring the abundance predictor log(sum lambda) by a fixed per-run constant C (the ridge slope, read from each run's existing per-chain summaries) decouples the two parameters. Verified a pure reparameterisation: likelihood and priors identical to 2e-13, pooled medians unchanged. Three-species win10 test: bobcat theta1 R-hat 1.24 to 1.004 (bulk ESS 37 to 1,889) and its apparent second mode dissolved, so that was the ridge and not real bimodality; white-tailed deer theta1 R-hat 1.045 to 1.020 (tail ESS 450 to 3,544), trend already converged; kit fox not fixed, with two genuine data-driven solutions, confirming the sparse hard-10 are data-limited rather than parameterisation-limited. Rolling out to the joint fleet at higher iterations to clear theta1, with a confirmatory simulation study of parameter recovery and coverage.
- Evidence: hazel_reparam_test_prompt.md; ctr_test_pull/ convergence summaries; camera-inat-stream-agreement skill

**M2. How to handle organised reporting campaigns: model term, exclude the cells, or limitation only**

- Status: decided - exclude the cells. Decided by: Roland, with Ben.
- Why it matters: Four strong-effect species had campaign-inflated iNaturalist reporting in specific cells (bobcat 47, coyote 5, California ground squirrel 3, northern raccoon 3; 58 cells total). A campaign model term is not identifiable, because no campaign footprint has a pre-campaign camera deployment to separate reporting from local ecology. The iNaturalist records in the campaign footprints are dropped for these species, camera data kept, and the exclusion is reported as a limitation. Applied in the bundle rebuild and folded into the centring rollout so each species is fit once.
- Evidence: table_campaign_strong_effects.csv; REPORT_campaign_effects.docx


## How to read the status values

- **open** means no decision has been made and nothing is preventing one.
- **open - needs your call** marks the subset where the evidence is assembled and
  the next step is a choice rather than more analysis.
- **open - blocked** means a running job must finish first.
- **open - never returned** means the work was specified and dispatched and no
  result arrived.
- **proposed** means the idea exists and has not been approved.
- **decided** means a conclusion was reached and the reasoning recorded.

## Glossary

**Camera-only control fit.** The same species refitted with the iNaturalist
records removed, used to test whether a joint trend depends on the opportunistic
stream.

**Congruence gate.** A screen on the exponent linking records to abundance, used
to decide whether a species' joint trend can be reported without a control fit.

**Full record.** The 2008 to 2025 fitting window, as against the ten-year and
five-year windows.

**Reporting campaign.** An organised project asking the public to report a named
species in a named region through iNaturalist.

**Selectivity.** How much faster a campaign's target species' records grew than
every other mammal's in the same cells, relative to cells outside the campaign.
