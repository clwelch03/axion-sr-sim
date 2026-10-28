# Fixed-target profile investigation — October 28, 2026

The user requested outlining and beginning an investigation of direct solutions
at the plot amplitudes after the weak-secondary ramp continued to need roughly
eight Newton corrections for each 1.5625% amplitude increment. The goal remains
jointly solved primary/secondary profiles at A1=`0.036390480353373555` and
A2=`1.0896058341192524e-4`, followed by the linear/Omiya/joint comparison plot.

## Question and existing evidence

Does the current matching location amplify the reference guess so strongly that
Newton is unnecessarily difficult at the target amplitude? Moving the join
could help, or could merely move the difficulty to the other shooting branch.
We have not established that weak-amplitude continuation is required.

The October 26 direct attempt used the reference primary and secondary boundary
data and zero mixed-component coefficients. At r*=5, the horizon and outer
field norms were approximately 0.05055 and 129715.5, dominated by the secondary
l=6 inward branch. The initial scaled mismatch was `1.906482699496033e9`.
Twelve Newton corrections did not converge. Halving radial spacing barely
changed the initial mismatch. These observations motivate a shooting diagnosis;
they do not prove a cause, absence of a target solution, or numerical correctness
of every part of the finite-amplitude calculation.

## First experiment: change only the matching location

Use the original saved target-amplitude input
`runs/adiabatic/multiphase_two_cloud/results/direct_fig7_20261026_v5-input.jls`
and verified source snapshot `work/amplitude_continuation_20261028_v1` in that
campaign. Keep amplitudes, frequencies, boundary guesses, mode content,
quadrature, radial endpoints, and radial spacing fixed. Test r*=5, 10, 20, 40,
and 60. The r*=5 result must reproduce the recorded baseline within `1e-8`
relative tolerance.

At each join:

1. Save both integrated branches, field/slope norms, relative value/slope
   jumps, the scaled matching norm, and the modes dominating the mismatch.
2. Build the full 56-column matching Jacobian with the verified local-source
   AD specialization and one full-width ForwardDiff batch.
3. Record its condition number in native packed coordinates and after dividing
   rows by their largest absolute entry, then columns by theirs. Save those
   scales. These are coordinate-dependent diagnostics, not convergence tests.
4. Compute the proposed Newton correction without applying it; record its size,
   admissibility, and linear-system residual. Save results after each join so
   a later failure does not erase earlier evidence.

Standalone scripts: `work/target_match_scan_20261028_v1.jl` and adjacent `.sh`.
The setup assertions and all five problem reconstructions were checked before
submission. The experiment records input, script, and source-manifest hashes.
It does not change the production library or continuation driver.

Job **27380795** COMPLETED on n0072.lr7 in **5m19s**, October 28,
12:01:41–12:07:00 PDT, using one CPU and 8 GB. Log:
`logs/target-match-scan-27380795.out` in the campaign. Output prefix:
`results/target_match_scan_20261028_v1`. Existing amplitude continuation job
**27378594** remains separate.

## Scan results

All five profile integrations and full-width Jacobians completed. The r*=5
baseline reproduces the recorded residual exactly. R denotes the vector of
radial field coefficients; Q=dR/dr* denotes the slopes. Relative jumps use the
larger norm of the two endpoint vectors, separately for R and Q.

| Join r* | Horizon norm R | Outer norm R | Relative R jump | Relative Q jump | Scaled residual | Native condition | Equilibrated condition | Full step admissible |
|---|---:|---:|---:|---:|---:|---:|---:|---|
| 5 | 0.05055 | 129715.5 | 1.000 | 1.000 | 1.90648e9 | 7.53e15 | 1.68e4 | no |
| 10 | 0.08505 | 807.953 | 1.000 | 1.000 | 9.47348e6 | 3.00e13 | 1.58e3 | no |
| 20 | 0.12218 | 5.33903 | 1.00044 | 0.99999 | 51008.4 | 1.60e14 | 2.00e5 | no |
| 40 | 0.10425 | 0.10648 | 0.41652 | 0.96710 | 413.240 | 3.01e14 | 4.61e3 | yes |
| 60 | 0.06573 | 0.06340 | 0.09335 | 0.29657 | 56.5614 | 1.95e16 | 2.99e3 | yes |

The largest mismatch is the secondary l=6 mode through r*=40; at r*=60 it
is the driven `(l,q)=(4,(2,-1))` mode. Jacobian passes took 47.6–63.9 s.
Linear-system relative residuals ranged from `1.90e-15` to `4.54e-13`.
None of these initial guesses satisfies the `1e-5` matching tolerance.

The changed join greatly reduces the outer branch's growth before matching.
At r*=40 and 60, the initial branch field norms are comparable and full Newton
steps satisfy the frequency/branch admissibility check. This supports testing
a direct solve there. Native conditioning remains severe, especially at 60;
row/column equilibration is only a diagnostic here and is not used by Newton.
The r*=20 result illustrates that a smaller residual need not mean better
conditioning. The scan does not isolate a physical mechanism or prove convergence.

## Bounded direct solve at r*=60

Selected r*=60 for the smallest initial mismatch, comparable branch sizes,
admissible full Newton step, and a moderate condition estimate after the stated
equilibration. Job **27381429** was submitted using
`work/target_match_newton_20261028_v1.jl` and adjacent `.sh`, with one CPU,
8 GB, a 30-minute limit, **six Newton corrections**, full-width local-source
AD, and unchanged `1e-5` tolerance. The original target amplitudes and reference
guess are retained. Setup assertions and shell syntax passed before submission.

Input: `results/target_match_scan_20261028_v1-match60-profile.jls`.
Output prefix: `results/target_match60_20261028_v1`.
Log: `logs/target-match-newton-27381429.out`.
The script saves setup provenance and every accepted Newton state. It saves a
complete matched profile only after fresh integration passes tolerance and
both radial branches are finite. A correction-limit failure preserves the last
accepted state and an error report. No production library was changed.

**Outcome:** job **27381429 FAILED** at the six-correction limit after **13m32s**,
ending October 28 at **12:27:05 PDT**. The accepted scaled residuals were
`56.5614 → 54.1791 → 54.1119 → 54.0477 → 53.9864 → 53.9279 → 53.8721`:
only a **4.75%** decrease, still far above `1e-5`. The first correction evaluated
six trial residuals; each later correction evaluated nine before accepting a
decrease. Jacobian passes took 56.4–76.6 s, and trial integrations about 8 s.
Large proposed corrections produced residuals around 6000, so admissibility
alone did not make them useful Newton steps. The last accepted state and failure
report are preserved; no matched target profile was produced.

This test shows that changing the join to 60 improves the starting mismatch but
does not by itself give useful convergence for this guess. Increasing the
iteration cap is not the recommended next experiment. Before changing the
boundary-value formulation, compare the saved Newton direction with an explicitly
equilibrated linear solve and inspect the components responsible for poor descent.
That separates numerical error in the linear correction from a poor nonlinear
prediction. Residual weighting can also affect line-search acceptance and should
be distinguished from scaling the linear solve. Shorter shooting segments remain
a candidate if sensitivity across the long intervals is the limiting issue.

At the same review, continuation **27378594** remained RUNNING after about
1h07m. It had saved seven new matched points, the latest at
A2=`5.069450807458185e-6` (**4.65255%** of target), with freshly integrated
weighted norm `5.828310199789505e-6`. Attempt 8 was underway.

## Why the first Newton step fails: saved-matrix check

In response to the user's request for the cause, the standalone
`work/target_match60_linear_check_20261028_v1.jl` re-solved the saved initial
56-by-56 system with explicit row/column equilibration and with 256-bit
arithmetic. It ran no radial integrations. Log:
`logs/target-match60-linear-check-20261028-v1.out`.

Both alternatives agree with the original correction to `9.46e-15` relative
norm; the maximum difference per coordinate divided by `max(abs(correction),1)`
is below `7e-12`. The original relative linear-system residual is `4.54e-13`.
This rules out inaccurate solution of the saved linear system as the main
reason that the first step fails. Higher precision used the same saved Float64
matrix: it does not validate the Jacobian entries or radial discretization.

The linearized equations predict cancellation of the residual, but integrating
the full proposed step gives `6118.8183` versus initial `56.5614`. All initial
steps from full size through 1/256 are admissible. The log's sixth trial is
therefore the **1/32 step**, the first to decrease the residual to `54.1791`.
Subsequent corrections required nine trial integrations each, implying at least
eight halvings before acceptance. This establishes poor finite-step prediction
as the immediate numerical obstacle; agreement of linear solvers does not prove
that the underlying cause is entirely physical nonlinearity.

The largest physical boundary-coefficient changes in that first correction are
in mixed components initialized to zero:

- Horizon `(l,q)=(0,(2,-1))`: change `-0.0204086 + 0.0199907im`.
- Outer `(l,q)=(3,(-1,2))`: change `0.0363277 + 0.0918369im`.
- Outer `(l,q)=(0,(2,-1))`: change `0.0142871 + 0.0090578im`.

Different components have different asymptotic radial factors, so these
coefficients alone do not compare local field strengths or energies. Their
dominance in the proposed correction makes the missing driven-field profiles a
specific initial-guess hypothesis. A discriminating follow-up is to reconstruct
the existing primary-secondary driven solutions at the target amplitudes and
test their use as mixed-component seeds, or isolate those correction blocks in
trial integrations. The local derivative along the actual target Newton direction
also remains to be checked; the October 27 full derivative verification covered
weak-secondary states. No new heavy solve was submitted for this matrix check.

## Adding the driven fields to the initial guess

The user explicitly requested adding the driven fields. Job **27384564** was
submitted using `work/driven_field_seed_20261028_v1.jl` and adjacent `.sh`,
with one CPU, 8 GB, and a 30-minute limit. It uses the verified 201-file
`amplitude_continuation_20261028_v1` snapshot and the original frozen point-72
primary/secondary references. Live library files have concurrent changes, so
the experiment continues to load the frozen code. Setup, mode-order, and
full-script syntax checks passed before submission.

The m=0 modes are l=0,2,4 and the m=3 modes l=3,5, matching the joint truncation.
Solve the existing forced equations with unit secondary amplitude, check their
relative matching residuals and power identities, then scale boundary
coefficients by A2 for m=0 and A2 squared for m=3. The primary/secondary
coefficients and complex frequencies remain those of the original guess.

There is a boundary-convention difference to expose explicitly: the established
m=3 solver constructs outgoing data at r*=5000 (order 16) and propagates it to
300. The current joint solver applies an order-7 series directly at 300. The
reference campaign already documents poor convergence of this near-boundary
outgoing series. The experiment retains the standard far-leg forced solution
and a second forced solution using `far_position=300, far_order=7`, which is
compatible with the current joint boundary prescription. The latter is only
an initial-guess control, not a validation of that outgoing boundary.

At real carrier frequencies, boundary values and slopes for m=0 and the near
m=3 variant must match the joint boundary construction to `1e-12` relative.
Record the far/near outer discrepancy and the changes caused by the joint
complex frequency combinations. Sharing a leading asymptotic coefficient does
not imply identical finite-radius boundary data.

At fixed target amplitudes and join r*=60, save and compare five initial guesses:
zero driven fields, m=0 only, near-boundary m=3 only, both near-boundary fields,
and m=0 plus the far-leg m=3 reference. Record the total matching norm and its
primary, secondary, m=0, and m=3 contributions. The zero case must reproduce
the saved 56.5614 baseline. Each saved guess can be supplied to the existing
bounded-Newton diagnostic without editing the production solver.

Log: `logs/driven-field-seed-27384564.out`. Output prefix:
`results/driven_field_seed_20261028_v1`, with setup provenance, forced solutions,
a summary CSV, and separate initial-profile files. No production library or
reference driver was changed.

**Construction outcome:** job **27384564 COMPLETED** on n0003.lr7 in **1m33s**,
12:54:38–12:56:11 PDT on October 28. Relative forced matching residuals are
`1.36e-15` (m=0), `6.38e-15` (far m=3), and `1.84e-15` (near m=3).
Their column-scaled condition numbers are 1538.6, 29.87, and 24.20. The m=0
and standard far-leg m=3 power ratios are `0.9999999253711056` and
`1.0000000001966904`, passing the unchanged `1e-6` check.

At real carrier frequencies, all expected matching boundary constructions agree
exactly. For identical far-leading coefficients, the joint near-series m=3
outer state differs from the far-propagated state by **0.283333** in relative
norm. This comparison fixes the coefficients; it does not measure the difference
between the boundary subspaces after freely refitting coefficients. Changing
to the joint complex frequency combinations changes the m=0 boundary states
by up to `1.31e-4` and m=3 states by up to `4.27e-4` relative.

All five initial guesses are saved with finite branches. The total and block
norms below use the same coefficient scale; block norms combine field and slope
mismatches. They are not percentage profile errors.

| Added fields | Total norm | Primary | Secondary | m=0 | m=3 |
|---|---:|---:|---:|---:|---:|
| None | 56.5614 | 24.0759 | 33.8779 | 38.0646 | 4.78717 |
| m=0 only | 3529.9916 | 3493.8218 | 434.9336 | 254.3749 | 13.2701 |
| Near-series m=3 only | 57.0283 | 24.6753 | 33.7887 | 38.7334 | 1.18132 |
| Both, near-series m=3 coefficients | 3529.2624 | 3492.8481 | 435.5649 | 256.3772 | 16.0870 |
| Both, far-leg m=3 coefficients | 3529.0943 | 3492.6665 | 435.6328 | 256.4717 | 15.2849 |

Adding m=3 alone reduces its own mismatch but does not improve the total.
Adding m=0 greatly increases the primary mismatch in the joint integration.
The far/near radiative coefficient choice barely changes that overall result.
Thus the missing fields alone do not provide a smaller-residual starting point;
this weakens the simplest version of the initial-guess hypothesis. It does not
by itself establish how Newton will behave from the new guess.

**Controlled Newton comparison:** submitted job **27384946** using the existing
`target_match_newton_20261028_v1.sh` and unchanged diagnostic code. Input:
`results/driven_field_seed_20261028_v1-both_near.jls`. Output prefix:
`results/target_match60_driven_20261028_v1`. Log:
`logs/target-match-newton-27384946.out`. It uses six corrections, full-width
local AD, join 60, full plot amplitudes, and `1e-5` tolerance, matching the
zero-driven test. Accepted states are checkpointed. The near-series forced
variant keeps the joint boundary prescription fixed for this comparison.

**Driven-seed Newton outcome:** **27384946 FAILED** at the six-correction cap
after **8m36s**, 12:58:26–13:07:02 PDT on n0013.lr7. Accepted norms:
`3529.26235 → 2267.89812 → 2147.61360 → 2072.02371 → 2022.38703 →
1896.03094 → 1692.49867`. This is about a 52% decrease, still far above
`1e-5`. Trial residual evaluations per correction were **4,3,4,3,3,1**.
The final correction reduced the norm by 10.7% and accepted its first evaluated
trial; the log does not separately count inadmissible trials that skipped
integration. Jacobian passes took 51.5–66.6 s. The final accepted state and
failure report are preserved; no matched target artifact exists.

The driven-seed case is making more fractional progress and needs fewer trial
integrations than the zero-driven test. It has not entered rapid convergence,
and its absolute residual remains much larger. A further bounded solve could
resume the saved accepted state; no extension was submitted during this review.
The earlier recommendation against merely increasing the cap referred to the
nearly stalled zero-driven attempt and should not be read as a result for this
different trajectory.

At this review continuation **27378594** was still RUNNING after about 2h01m,
with **13** newly saved matched points. Latest A2=`5.394443404258194e-6`,
**4.95082%** of target, weighted norm `6.00678166743477e-6`. Attempt 9 took
14 corrections and halved the next increment to **0.78125%**; attempts 10–13
each took six corrections. Attempt 14 was underway.

## Requested twelve-correction extension

The user explicitly requested twelve more corrections for the driven-field test.
Submitted job **27388374** using standalone
`work/resume_target_match_newton_20261028_v1.jl` and adjacent `.sh`.
It loads the saved correction-6 state from
`results/target_match60_driven_20261028_v1-accepted.jls`, SHA-256
`f680a44bb983ede293e36575f27300ffe4472d4a063ba2061611f0278fa4476f`, with stored
norm `1692.4986686897735`. The iteration cap is **12 additional**, up to
**18 total**; accepted-state logs/checkpoints and any matched solution preserve
the cumulative count. Convergence or numerical failure can stop it earlier.

The frozen library, target amplitudes, join 60, boundary prescription, residual
weights, full-width local AD, and `1e-5` tolerance are unchanged. Setup-only
execution checked exact state identity, the saved mismatch, all settings, and
the total correction limit. Full-script and shell syntax passed, and all 201
snapshot files verified. The actual first residual evaluation must reproduce
the stored mismatch to `1e-10` relative before proceeding.

Output prefix: `results/target_match60_driven_20261028_v2`.
Log: `logs/target-match-resume-27388374.out`.
Resources: one CPU, 8 GB, one-hour limit. Prior input and results are preserved.
The script saves every accepted state and requires fresh matching and finite
radial branches before saving a complete matched profile.
At launch, **27388374** is RUNNING on n0060.lr7 from **13:27:03 PDT**. The
first fresh residual exactly reproduced the stored `1692.4986686897735`;
correction 6 is checkpointed under the new prefix and correction 7 is beginning.

## Decisions after the scan

A smaller initial residual alone is insufficient. Inspect both branches and
Jacobian diagnostics, then choose a promising location for a bounded Newton
attempt at the same target amplitudes, with accepted-state checkpoints and the
existing `1e-5` matching tolerance. Preserve the baseline for comparison.

If a target solve converges, freshly integrate it and test sensitivity to radial
spacing and the join before interpreting profile differences in the plot.
Matching at another location changes where numerical error is assessed; it does
not by itself establish spatial accuracy.

If no join gives a useful solve, use the branch growth and sensitivity data to
design shorter shooting segments or a solver for the full radial profile.
Do not infer a physical obstruction from failure of this particular shooting
formulation. Efficient time evolution remains a separate benchmark including
profile corrections and charge-response calculations.
