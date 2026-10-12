# QCD two-mode grid, October 12, 2026

**Disconnect checkpoint: stopped; no production array is running.** The user
asked to leave a stable state before disconnecting. The production array
26875559 was automatically canceled by its failed pilot dependency. The
summary job finished. The remaining diagnostic 26875867 was explicitly
canceled before disconnect. No automatic retry or future submission is set.
Inputs, source hashes, successful initial points, and failure logs are saved.

This campaign follows the requested 10 by 10 logarithmic grid in initial BH
mass (2–100 solar masses) and axion rest energy (1e-14–1e-11 eV). The final
requested exclusion is **initial alpha > 0.45**, superseding the earlier
alpha <= 1 proposal. There are **83 included and 17 excluded cells**;
the largest retained alpha is 0.4376165685906593. The requested physical age
limit is **10 million years**.

## Inputs and interpretation

- Initial spin is 0.99. The QCD mapping is `fa_GeV = 5.7e6 / axion_eV`.
- Initial physical cloud energies divided by initial BH mass are fixed at
  `2.1796991094505996e-6` (211) and `2.1796991094505997e-9` (322), matching
  the original QCD run. Code energies therefore vary with fa.
- `grid.csv` records all 100 cells. `manifest.json` records case ordering,
  stopping rules, revision, and SHA-256 hashes. The campaign imports the shared
  repository `src/` directly and checks its recorded hashes before execution.
  Source edits require an explicit campaign update before any new submission.
- Forty cells satisfy both alpha < 0.02 and fewer than two estimated primary
  energy e-foldings in 10 million years. These receive explicitly labeled
  **small-alpha linear age screens**, using both hydrogenic growth rates and
  evolving BH mass and angular momentum. They omit mode-coupling fluxes and
  are not presented as nonlinear solutions or rigorous exclusion bounds.
- The other 43 cells use the existing nonlinear 211 primary, averaged 322
  secondary, and driven m=0,3 flux solvers with coupled BH/cloud balance.
  The outer boundary scales as at least `45/alpha0^2`, with step 0.25 M.
  Each case builds 23 stationary points on its own primary-spin-down pilot.
  A stiff time integrator restores the horizon frequency factor using the
  actual evolving spin. The pilot only chooses stationary solve locations.

The integration stops at the age limit, stationary energy boundary, a spin
offset of 0.002 or relative alpha offset of 1e-4 from its local stationary
track, a cloud/BH mass ratio of 0.1, or a secondary cubic/binding-scale ratio
of 0.01. A track-limit stop can precede the secondary-dominated epoch; it does
**not** establish the eventual maximum over 10 million years. Failed solves
and stopped evolutions remain explicit outcomes.

At each stationary point, the bound-cloud diagnostics maximize over a polar
grid including the equator and sample the radial profiles. Absolute values
of the different m blocks are summed to bound their relative phases on that
grid. Adding the separate 211 and 322 maxima gives a conservative envelope
estimate, rather than their measured coincident maximum. It excludes the
driven fields. This campaign inherits the original cosine-potential,
instantaneous-Kerr approximation, without GW losses or self-gravity.

## Execution and checks

The new time driver reproduces the saved 20-solar-mass reference trajectory
with a maximum spin difference of 1.71e-8. The analytic 211 field normalization,
fixed physical seed fractions, and energy-balance regression pass.
Each numerical trajectory also repeats time integration with tighter tolerances.
These checks do not establish radial, angular, or stationary-table convergence
throughout the new grid.

Initial jobs:

| Job | Purpose | Status when submitted |
|---|---|---|
| 26875539 | Initial stationary checks, indices 6, 9, 63 | Two concurrent checks |
| 26875559 | 83-cell array, maximum four concurrent cases | Gated on success of 26875539; cancel on invalid dependency |
| 26875675 | Summary CSV and heatmap | After any termination of 26875559 |
| 26875867 | Low-alpha absolute-residual and radial-refinement diagnostic | Investigates the failed initial check |

Initial checks at alpha=0.149663 and 0.437617 passed. At alpha=0.0149663,
the linear shooting solve failed to reduce its absolute matching residual
below 0.000279892. The dependency gate therefore prevents the first production
array from running. This is a numerical failure, not evidence about the field
amplitude. The source and failure logs are retained.

The canceled diagnostic completed its step=0.25 calculation before stopping:
with a diagnostic-only tolerance of 1e-3, it found
`omega = 0.014965856988351504 + 6.518757974394914e-19 im` and a matching residual
of `8.75923e-6` relative to the matching state. This suggests an absolute
normalization/conditioning issue but does not resolve it. The half-spacing
comparison was not completed. **Production tolerances were not loosened.**

## Resume here

1. Confirm with the user that the survey should resume. The current stopping
   point was requested before disconnect, not successful completion.
2. Read the low-alpha pilot and diagnostic logs. Complete the half-spacing
   comparison and establish a justified relative residual criterion, or improve
   the shooting conditioning. Do not just lower acceptance thresholds until a
   solver returns. Check physical fluxes and profile normalization too.
3. Preserve these inputs and results if solver changes are needed; prepare a
   new campaign manifest and new output prefixes using the shared `src/`.
   Source hashes intentionally reject changes to submitted code.
4. Re-run small-, reference-, and upper-alpha checks. Before a broad release,
   also complete an end-to-end numerical cell: only initial stationary points
   have passed so far. The reference time-driver regression passes separately.
5. Submit a fresh production array gated on successful checks, then a collector.
   The original dependency chain is canceled and must not be described as queued.

The batch is exploratory: 40 low-alpha cells use labeled linear age screens;
43 request the coupled stationary treatment. The one-dimensional primary
track may terminate before later secondary evolution. If comprehensive maxima
through 10 million years are required for those cells, extending the track/model
is still work to do. Existing heatmaps contain only two local age-screen results
and must not be used as a completed grid study.

Two inexpensive age-screen cases (indices 0 and 78) were executed locally
before array release to check output generation. Their maximum envelope
estimates are 5.25e-14 and 3.11e-9, respectively, in units of fa; they remain
explicitly labeled linear estimates.

## Commands and outputs

The duplicated `source/` tree was removed on October 12 at the user's request.
All 62 source files were byte-for-byte identical to their counterparts in
`src/` before removal. New campaign scripts, diagnostics, inputs, and results
remain here. Preparation records hashes and does not create another copy.

The original `source_sha256` map and `source_fingerprint` remain in the manifest
as historical provenance and checkpoint identity. `runtime_sha256` records
repository-relative paths for the shared solvers and updated launchers and is
what execution now verifies. Preserving the checkpoint identity is justified
by the identical solver files and the import-path-only cleanup; it does not
authorize resuming with modified numerical code.

Use the BHsuper Python from the repository root. Preparation is one-time and
refuses to overwrite a frozen manifest:

```bash
python runs/adiabatic/qcd_grid_20261012/scripts/campaign.py prepare
python runs/adiabatic/qcd_grid_20261012/scripts/check_campaign.py
python runs/adiabatic/qcd_grid_20261012/scripts/campaign.py smoke 9
python runs/adiabatic/qcd_grid_20261012/scripts/campaign.py run 9
python runs/adiabatic/qcd_grid_20261012/scripts/collect.py
```

Actual numeric solves must run through Slurm. `scripts/submit.sh` requests
one CPU, 24 GB, and at most three days per cell, with four cells concurrent.
Each cell has its own output directory and writer lock. Completed stationary
points resume only when source and numerical inputs match. A failed or
wall-time-limited run must not be interpreted as a completed survey cell.

Per-case files are in `results/mXX_aYY/`: target tables, stationary CSV/JLS
checkpoints, `trajectory.csv`, and `summary.json`. The collector preserves all
100 cell statuses in `results/summary.csv` and makes
`plots/qcd_grid_20261012.{png,pdf}` in the repository plot directory.
