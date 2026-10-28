# Multiphase workplan — October 27, 2026

**Today's end goal, potentially a stretch: produce the primary/secondary
deformation plot comparing linear, Omiya, and joint nonlinear profiles at the
target outer amplitudes.** The deliverable is the six-curve overlay, signed
deformation panels for both clouds, and the extra joint-minus-Omiya change,
exported as PDF/PNG with the plotted data in CSV. Label this first figure a
preliminary fixed-amplitude comparison.

Organize the day's work toward that plot: preserve and match the weak seed,
validate its derivatives, then continue to the target secondary amplitude.
Prepare the linear/Omiya curves and plotting layout while the solves run, so
the matched target profiles can go straight into the figure. Convergence and
validation determine whether we reach the end goal today; the seed and
reference curves are intermediate milestones.

**Later October 27 priority:** the user requested completing derivative
verification autonomously before further continuation. See the
[verification record](multiphase-ad-verification-2026-10-27.md) for the current
jobs, criteria, and results. The plot remains the end goal.

## Recovery checkpoint after the n0003 outage

The shared working files, source snapshot, saved chat transcripts, and job
artifacts are accessible from n0001. Slurm confirms that the three verification
jobs ran on n0060.lr7 and completed successfully, independently of n0003:

| Job | Recovered result |
|---|---|
| 27333738, initial-state comparison | COMPLETED at 15:33:40 PDT; all 56 columns and both finite-difference directions pass |
| 27334006, weak seed | COMPLETED at 15:26:26 PDT; saved and reintegrated norm `5.306800040310494e-6` |
| 27334124, matched-state comparison | COMPLETED at 15:59:48 PDT; all 56 columns and both directions pass; worst column error `9.587805e-9 < 1e-8` |

**Assistant workstream:** the agreed local-Jacobian AD verification is complete
for the tested configuration and states; see the updated verification record.
The seed, both finite radial branches, full matrices, probe sweeps, and all 199
frozen-source files were checked. Review the next A2 continuation step from
`results/ad_seed_20261027_v1-seed.jls` before submitting it. No ramp is queued,
and the current Slurm queue is empty. The seed uses 0.1% of target A2.

**User workstream:** `scripts/reference_profiles.jl` retains the loading,
reference selection, explanatory input tree, and Omiya-primary reconstruction
blocks. On resumption, the user explicitly retained manual driver edits.
The original secondary test ran as **27337260**: matching passed, but direct
integration at A2 failed the `1e-10` sample-scaling check (`8.24126e-10`).
The revised block reconstructs the unit-amplitude solution and scales the saved
samples afterward. **27337856** passed its matching, sample-scaling, and
independent equation-scaling checks; use
`work/check_reference_secondary_20261027_v2.jl`. **27337858** passed both
linear solves and reconstructions in `work/check_reference_linears_20261027_v1.jl`.
All four reference profiles are saved across
`results/reference_secondary_20261027_v2.jls` and
`results/reference_linears_20261027_v1.jls`. The driver was not edited: present
the tested secondary block for the user to add, then the linear blocks.
Combined CSV/plot export and spatial refinement remain ahead; see the campaign
README for the checks and exact residuals. This is a fixed-amplitude reference
preparation, not a target joint solution or a spatial-convergence result.

The original plan and earlier checkpoints follow for context; this recovery
checkpoint supersedes their pending seed/derivative status.

Reviewed [yesterday's plan](multiphase-workplan-2026-10-26.md), the
[handoff](../HANDOFF.md), campaign logs, output files, and live Slurm accounting
on October 27 at approximately **13:08 PDT**. The queue for `clwelch` was empty.
Accounting from October 26 through this check lists the five jobs below; all
ended on October 26, with the last finishing at **16:50:10 PDT**. There was no
additional overnight run or queued continuation.

## 1. What yesterday accomplished and what remains

The source walkthrough and reviewed full-width Newton batching change are
complete. The new direct solve was assessed, and the primary reference
reconstruction passed its check. Neither the weak seed nor the target joint
solution was saved; the four reference curves and comparison figure remain
unfinished.

All links in this table point to the original job logs.

| Job | Final state / elapsed | Numerical outcome | Consequence for today |
|---|---|---|---|
| [27296945: v4 derivative gate](../runs/adiabatic/multiphase_two_cloud/logs/direct-profile-27296945.out) | FAILED / 3m10s | Directional finite-difference error `2.38e-4`, above `1e-5`; Newton never started. | Historical evidence; the removed preflight is not part of today's workflow. |
| [27304701: direct v5](../runs/adiabatic/multiphase_two_cloud/logs/direct-profile-27304701.out) | FAILED / 29m04s | Twelve corrections: scaled residual `1.90648e9 → 1.49475e8`; target `1e-5`. | No converged solution or profiles. Follow the weak-secondary approach supported by yesterday's controls. |
| [27305885: primary reconstruction](../runs/adiabatic/multiphase_two_cloud/logs/reference-primary-check-27305885.out) | COMPLETED / 19s | Scaled mismatch `5.58739e-11`; finite radial branches. | The next reference-driver block has been tested; reference export is still pending. |
| [27305905: weak-secondary seed](../runs/adiabatic/multiphase_two_cloud/logs/weak-secondary-27305905.out) | FAILED / 14m00s | Twelve corrections: weighted residual `1.32554e7 → 2.27420e-5`; target `1e-5`. | Promising but not converged. No accepted-state file exists, so the retry must start from the original guess. |
| [27306004: independent Jacobian check](../runs/adiabatic/multiphase_two_cloud/logs/weak-jacobian-27306004.out) | FAILED / 2m31s | Ordinary/optimized AD agree to `7.47785e-15` in direction 1. Finite-difference errors are `0.5103`, `0.3044`, `0.1636`; direction 2 was not reached. | Preserve the saved diagnostic; validation is incomplete, not evidence that the optimized AD method is wrong. |

The weak seed's last three accepted residuals were `124.01057`,
`0.001374344`, and `2.274197e-5`. It ended only **2.274 times** above tolerance,
but that alone does not guarantee that another correction will pass. The direct
and weak norms use different row weights and must not be compared as if they
were the same error measure. Neither is a percentage profile error.

The direct attempt has an input file and failure report. The independent
Jacobian diagnostic saved its original state, full optimized Jacobian, first
direction, and errors in
[`weak_jacobian_20261026_v2.jls`](../runs/adiabatic/multiphase_two_cloud/results/weak_jacobian_20261026_v2.jls).
The weak seed's final state cannot be recovered from this file: it contains the
**initial** guess. Its post-convergence derivative check was never reached.

## 2. First working block: a small, recoverable seed retry

Budget **20–30 minutes of active work**, followed by Slurm queue/runtime.
Work through the driver changes together in small tested blocks, as established
yesterday. The existing library callback already receives copies of the state
and residual before the iteration-limit check; no solver-library change is
needed for this proposal.

1. Review a minimal use of `accepted_state_callback` in
   [`weak_secondary_profile.jl`](../runs/adiabatic/multiphase_two_cloud/scripts/weak_secondary_profile.jl).
   Save the latest accepted packed state, iteration, weighted residual, problem
   settings, amplitudes, weights, and source provenance to an explicitly named
   accepted-state artifact. Write through a temporary file before replacing it.
   Keep it distinct from the converged solution. This is a small driver addition,
   not a restart framework or a restored preflight gate.
2. Propose **20 maximum corrections**, up from 12, with the **same `1e-5`
   tolerance** and existing stagnation checks. Keep the primary amplitude fixed
   at `0.036390480353373555` and the secondary at
   `1.0896058341192524e-7` (0.1% of its target). Preserve all spatial settings,
   the original coefficient scale, secondary row weights, and full-width AD.
3. Check the actual callback/save path on a cheap known-answer problem that
   reaches its iteration cap; load the artifact and verify its state and
   residual. Run the existing setup-only check. Freeze a fresh snapshot,
   proposed `weak_secondary_20261027_v1`, and verify its manifest before
   submission. Preserve yesterday's snapshots and outputs.
4. Use the existing one-CPU, 8-GB, one-hour Slurm allocation. Record the new job
   ID and prefix once submitted. Yesterday's 14 minutes supplies a runtime
   reference, not a convergence prediction. Report progress about once a minute
   while actively working, including what each change in residual means.

**Completion condition:** a freshly reintegrated weighted norm at most `1e-5`,
with the seed, both radial branches, and provenance saved before derivative
diagnostics. If the retry stalls, retain its accepted state and inspect the
component mismatches and accepted steps before selecting another calculation.
An accepted iterate above tolerance remains a diagnostic, not a joint prediction.

## 3. Derivative validation: choose probes that test a local change

Budget **20–30 minutes of active diagnosis**, with expensive evaluations on
Slurm. This can progress while the seed runs; the matched-state check depends
on obtaining that seed. Explain the proposed diagnostic before running it.

At yesterday's initial weak guess, `norm(J*d) = 5.8345761e17`, where `J` is the
weighted matching Jacobian and `d` is the saved first direction. Even the
smallest finite-difference step, `h=1e-8`, predicts a residual change
`h*norm(J*d)` about **440 times** the initial residual. Those probes did not
establish behavior in the local linear regime.

- Reuse the saved Jacobian/direction to select a smaller geometric sequence of
  probe sizes. For scale only, `h=1e-12` predicts about 4.4% of the original
  residual; it is a candidate probe, not a validated step size.
- Check admissibility, ensure perturbed coordinates actually change in floating
  point, and inspect error versus step size for a convergent range before
  cancellation or integration noise dominates. Near a matched root, do not
  choose the step solely as a fraction of the nearly zero residual.
- Evaluate and save **both planned directions**, including failures, before
  reporting the overall result. Compare optimized AD with ordinary-source AD
  and centered differences. Keep the existing relative criteria: `1e-8` for
  AD agreement and `1e-5` for a finite-difference step in each direction.
- Repeat at the saved matched seed. Agreement at the original guess alone does
  not complete the planned matched-state validation. If no useful finite-
  difference range emerges, record that limitation and review the diagnostic
  conditioning; do not relax its threshold to label it passed.

**Completion condition:** both directional comparisons at the matched seed
documented, or a specific remaining failure with the state and probe data
preserved. This tests derivatives of the discretized matching calculation;
spatial accuracy is a separate check.

## 4. Reference profiles while jobs run

Budget **45–60 minutes**, in the user's small-snippet workflow for
[`reference_profiles.jl`](../runs/adiabatic/multiphase_two_cloud/scripts/reference_profiles.jl).
The live driver currently loads inputs and selects the Omiya solutions. Start
with the tested `primary_profiles_at(primary_setup, primary_index)` block from
job 27305885, then reconstruct the averaged secondary and actual linear
eigenfunctions.

Produce the four linear/Omiya curves for both clouds at the target outer
coefficients `A1=0.036390480353373555`, `A2=1.0896058341192524e-4`, spin 0.99,
and alpha 0.3. These references prepare the final target-amplitude comparison;
a weak seed must not be overlaid as though it had the target secondary amplitude.
Use the existing `solve_linear_mode` machinery and retain the secondary's A2
scaling. Check radial-field and angular normalization, coordinate origin, and
both values and slopes where branches join.

The saved primary uses radial step 0.1; the averaged-secondary reconstruction
uses step 0.2. Record that difference and assess radial refinement before
interpreting a small curve separation. Preserve the original references.
Export the four curves and prepare the final overlay/difference-panel layout
with provenance while the joint solves run. A reference-only preview is an
intermediate artifact; the end-goal figure includes the matched joint curves.
Do not reuse the older plotter's small-amplitude/equal-energy rescaling as the
linear reference for this fixed-outer-coefficient comparison.

## 5. Reach the target amplitude and assemble the end-goal plot

Once the weak seed matches and its derivative checks pass, review the results
together before increasing A2. The next proposal is one small amplitude step
from the saved seed, followed by reintegration and saving its result before
choosing another step. Aim to continue to the target
`A2=1.0896058341192524e-4`, with A1 and spatial settings fixed. Choose successive
steps from the measured convergence behavior and preserve each matched result.

If continuation reaches the target, assemble the six-curve linear/Omiya/joint
overlay and signed difference panels defined in
[yesterday's comparison specification](multiphase-workplan-2026-10-26.md#3-get-to-the-mutual-nonlinear-deformation-plot).
Reserve **30–45 minutes of active work** for adding the joint curves, checking
branch joins and common-grid interpolation, and exporting CSV, PDF, and PNG.
Record amplitudes, discretization, matching norm, and source run IDs. Label it
a **preliminary fixed-amplitude comparison**.
Joint energies are outputs. Follow with radial step-halving and then separate
angular/phase, boundary/domain, and retained-mode checks before interpreting
small differences as resolved mutual deformation.

Reserve **20–30 minutes after the next convergence attempt** for the
[architecture review](../HANDOFF.md#important-revisit-multiphase-architecture):
trace ordinary source evaluation, the local Jacobian, and the ForwardDiff
adapter; identify duplicated mathematics, packing, or unnecessary wrappers.
Record a minimal simplification proposal and its decisive check. Keep the
user's learning-map boxes under their control and continue tested snippets for
library changes. The previously recorded evolution-test field mismatch is a
separate backlog item, not evidence about today's stationary solve.

**End-of-day target:** the exported linear/Omiya/joint deformation plot at the
target amplitudes, with its data and validation status recorded. This is the
day's end goal even if reaching it proves a stretch.

**If convergence or validation blocks the plot:** preserve the latest accepted
state, record exactly which check or amplitude step remains, and leave the
reference curves and plotting layout ready to receive the joint result.
Record the architecture review as well. Before leaving, update the
handoff/campaign record with exact job IDs, artifacts, and any jobs still
running. This planning review itself submitted no jobs and changed no solver
or driver code.
