# Joint two-cloud response and evolution

## Missing linear-secondary retries (October 29)

Array **27425079** retries five secondary references: alpha=0.3/spin=0.9,
then alpha=0.2 at spins 0.99/0.95/0.9/0.8. The separately recovered spin-0.95,
alpha-0.3 pair remains in `results/omiya_spin095_linear_20261029_v2/` below.
There are no new primary solves. Outputs use a fresh directory and preserve
the original grid, failed attempts, and comparison references.

The opt-in precision rule follows the spin-0.95 recovery. A coarse solve at
`1e-5` locates the root; with its boundary coefficients held fixed, the driver
measures changes in the residual vector at frequency plus/minus one real-frequency
ULP. A restarted final solve and fresh profile reconstruction must both pass
`max(1e-8, 2*max(probe_changes))`. Coarse tolerance alone never accepts a
reference. Each report retains the requested `1e-8`, actual acceptance tolerance,
coarse root/residual, probe sizes and changes, and the existing profile/energy/flux
diagnostics. These absolute residuals use unit first outer coefficient. This
precision allowance does not establish spatial convergence or tiny growth rates.
The original strict grid mode remains the default.

Frozen source/config: `work/linear_grid_retry_20261029_v1`; results:
`results/linear_grid_retry_20261029_v1`; logs: `logs/linear-grid-27425079_TASK.out`.
All source/config hashes and Python/Julia/shell syntax passed before submission.
At most four cases run concurrently, each with one CPU, 8 GiB, and 20 minutes.
**Outcome:** all five tasks COMPLETED in 23–26 seconds. Fresh profile matching,
positive free-field energies, and all solution/profile/CSV outputs were verified.

| alpha | spin | Matching residual | Precision-adjusted tolerance |
|---|---|---|---|
| 0.3 | 0.9 | 1.47906e-8 | 1.06769e-7 |
| 0.2 | 0.99 | 1.49161e-6 | 1.54520e-5 |
| 0.2 | 0.95 | 3.61701e-7 | 1.50881e-5 |
| 0.2 | 0.9 | 6.04092e-7 | 1.44533e-5 |
| 0.2 | 0.8 | 1.35327e-6 | 1.40253e-5 |

All twelve reference pairs are now available across the preserved exports:
all primaries and six secondaries in `results/linear_grid_20261029_v1`,
the alpha=0.3/spin=0.95 secondary in
`results/omiya_spin095_linear_20261029_v2/linear_reference.jls`, and these five
secondaries in `results/linear_grid_retry_20261029_v1`. The spin-0.95 combined
export is energy-normalized (its metadata gives amplitude factors); the grid
exports retain unit outer coefficient and record unit-energy factors. Preserve
that distinction when loading seeds. Reproduce the retry with a fresh name:

```bash
python runs/adiabatic/multiphase_two_cloud/scripts/prepare_linear_reference_grid.py UNIQUE_NAME --retry-missing-secondaries
sbatch --parsable --array=1-5%4 \
    runs/adiabatic/multiphase_two_cloud/scripts/submit_linear_reference_grid.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/UNIQUE_NAME" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/UNIQUE_NAME"
```

## Linear references across twelve backgrounds (October 29)

Slurm array **27418189** attempts the two existing linear eigenmode solves for
every alpha `[0.2,0.3,0.4]` and spin `[0.8,0.9,0.95,0.99]` combination. Each
task saves m=1 with ell `[1,3,5]` and m=2 with ell `[2,4,6]`. At submission
verification, four tasks were RUNNING and eight PENDING. Each task has one CPU,
8 GiB, and a 20-minute limit; up to four run concurrently.

**Final first-pass outcome:** all 12 tasks finished in 23–37 seconds each.
All 12 primary profiles and six secondary profiles passed. Tasks 1, 4, and
9–12 COMPLETED; tasks 2, 3, and 5–8 FAILED solely on the secondary Newton
descent check. At alpha=0.3, spin=0.95/0.9, the stalled residuals are
`1.16608e-8`/`1.47904e-8`. At alpha=0.2, spins 0.99/0.95/0.9/0.8, they are
`1.00375e-6`, `1.33584e-7`, `2.01530e-7`, and `1.07736e-6`. The unchanged
threshold is `1e-8`. All alpha=0.4 pairs passed, and the alpha=0.3, spin=0.99
frequencies exactly reproduce the previous reference. These failures need a
focused residual-scaling/precision diagnostic before accepting or rerunning
the missing modes; no automatic tolerance relaxation was made.

The driver calls `solve_linear_mode` without changing the library. Initial
frequency estimates use principal number n=m+1. Absolute matching tolerance is
`1e-8`, with at most 20 Newton iterations. Radial steps are 0.1 for m=1 and 0.2
for m=2; the matching point is 5 and inner endpoint -100. The outer endpoint is
675 for alpha=0.2 and 300 otherwise. Angular quadrature uses 16 polar and 31
azimuthal nodes; boundary series retain the library defaults (orders 3/7).
These are initial reference settings, not a spatial or endpoint convergence
study. The saved flux-balance diagnostic is not an acceptance criterion and
does not certify accuracy of tiny imaginary frequencies.

The frozen source and 12-case TOML are under `work/linear_grid_20261029_v1`;
the source, project, and config are covered by `source-manifest.sha256`.
Outputs are under `results/linear_grid_20261029_v1`, with one background prefix
and separate primary/secondary solution JLS, profile JLS/CSV, and metadata TOML
files. Logs are `logs/linear-grid-27418189_TASK.out`. A mode failure is saved
without preventing the other mode's attempt; the pair summary passes only if
both modes pass. Existing output prefixes are never overwritten.

Stored profiles have first outer coefficient one. Each report records the
quadratic free-field energy E and unit-energy amplitude factor `1/sqrt(E)`;
the arrays themselves are not rescaled. Reports also retain complex frequency,
matching residuals, signed finite-boundary fluxes, superradiant margin, and
provenance. Both sides of the matching point remain separate in the CSV.
No cloud energy or BH mass is chosen by this dimensionless linear batch.

Shell/Python/Julia syntax checks passed before submission. Reproduce using a
fresh name and the same Julia environment as the other campaign jobs:

```bash
python3.11 runs/adiabatic/multiphase_two_cloud/scripts/prepare_linear_reference_grid.py UNIQUE_NAME
sbatch --parsable --array=1-12%4 \
    runs/adiabatic/multiphase_two_cloud/scripts/submit_linear_reference_grid.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/UNIQUE_NAME" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/UNIQUE_NAME"
```

## Independent energy and background sweep (October 29)

Slurm array **27409247** contains all 60 combinations of alpha `[0.2,0.3,0.4]`,
spin `[0.8,0.9,0.95,0.99]`, and energy factors
`[(1,1),(0.5,1),(1,0.5),(1.5,1),(1,1.5)]` relative to
`E0=650.5419114886128 F_a^2 M`. The two energies are independent inputs.
They use the comparison convention described below; full interacting-field E/J
and signed horizon/outer fluxes are saved separately. These are stationary
profile attempts, not time evolution.

The frozen snapshot is `work/parameter_energy_20261029_v1`, with the complete
case list in `inputs/sweep.toml` and last night's accepted equal-energy solution
in `inputs/initial.jls`. A SHA-256 manifest covers the source, project, and inputs.
The array runs at most four cases at once, each with one CPU, 24 GiB, and a
four-hour limit. Each case takes bounded steps in alpha/spin/log energy, adjusts
both boundary amplitudes, and solves the full radial equations after changes.
Intermediate differential tolerance is `1e-5`; final tolerance is `1e-6`.
Final energy errors must be below `1e-4`, with quadrature agreement below `5e-6`.
Every computed radial result and accepted continuation stage is checkpointed.
Failure stops that case without accepting or deleting its saved intermediate data.

The provisional domain is `[-100,675]` for alpha=0.2 and `[-100,300]` otherwise;
the extension supplies an initial guess from the old outer series and then
re-solves the extended domain. These choices do not establish endpoint or
angular convergence. No sign restriction is imposed on frequency imaginary parts.
Some requested configurations may fail to converge or lie beyond the branch
reachable from the supplied seed. Every case starts from the preserved seed.

Outputs live under `results/parameter_energy_20261029_v1/`, with one unique
case prefix, request/history TOMLs, solve/stage JLS files, and final JLS/TOML/CSV
files only for accepted targets. The final CSV contains all 14 retained modes.
Use a saved accepted stage as the seed with a fresh output directory to restart.
The baseline configuration reuses the original result after its energy check.

Validation job **27409219** checks the 60 configurations and exercises unequal
energy targeting on a small cosine problem with changed alpha/spin/domain.
It was RUNNING when the production array was submitted. The user explicitly
requested proceeding without waiting; the array has no validation dependency.
Validation outputs: `work/parameter_energy_20261029_v1_validation/`.
Logs: `logs/parameters-27409247_TASK.out` and
`logs/parameters-27409219_4294967294.out`. Submission status: array PENDING.

First status check: the reused baseline (task 1) completed, four energy-variation
tasks (2–5) are running with their initial radial solves accepted, and 55 tasks
are queued. Validation **27409219** failed after 2m54s on `KeyError: 2`, after
153 successful setup checks and one successful small radial-solve check.
The energy helper built an azimuthal table only through the primary's maximum m;
the small test needs secondary m=2 even though it retains only primary m=1.
The live helper now includes both clouds when sizing the table. The production
14-mode table already extends to m=5 and is unaffected. Its frozen source remains
unchanged. Validation **27409399** was submitted from
`work/parameter_energy_20261029_validation_v2`, a copy of the original snapshot
with only the helper fix and two regression assertions changed, recorded in its
manifest/provenance. Results: `work/parameter_energy_20261029_validation_v2_results`;
log: `logs/parameters-27409399_4294967294.out`.

At the next check, tasks 2–5 had passed 3–4 radial solves each; the accepted
baseline and 55 queued cases were unchanged. Validation **27409399** failed
after 4m26s: the final radial refinement moved the primary's energy by
`4.43e-4` relative, beyond the `1e-4` acceptance threshold, despite successful
coarse energy targeting. The final acceptance check correctly rejected it.
The live driver now repeats energy targeting after refinement, retaining
differential tolerance `1e-6` and using quadrature step 0.5. Validation
**27410452** uses `work/parameter_energy_20261029_validation_v3` and writes
to `work/parameter_energy_20261029_validation_v3_results`; the log is
`logs/parameters-27410452_4294967294.out`. The original production snapshot
remains unchanged. A production case with the same final mismatch will fail
explicitly and can be corrected from its preserved radial checkpoint.

Validation **27410452 COMPLETED** in 4m47s, with `validation.toml` reporting
passed. Both refined energy targets pass after one additional correction.
At ~50 minutes, production tasks 2–5 had 7–8 accepted radial solves each,
but only the baseline was final; no production failures were recorded.
Most radial solves take 5–7 minutes. The half-energy targets use three
intermediate stages and the 1.5x targets use two, with a fresh Broyden estimate
per stage. The live driver now skips the redundant same-background/same-amplitude
solve at each stage boundary. This small change is syntax-checked and postdates
the successful validation snapshot. Existing running/queued jobs retain the
original frozen driver. Last night's single refinement was 12m00s, and its
complete equal-energy target adjustment was 41m53s.

```bash
python3.11 runs/adiabatic/multiphase_two_cloud/scripts/prepare_parameter_energy_sweep.py UNIQUE_NAME \
    --energy-pair 1 1 --energy-pair 0.5 1 --energy-pair 1 0.5 \
    --energy-pair 1.5 1 --energy-pair 1 1.5
sbatch --parsable --array=1-60%4 \
    runs/adiabatic/multiphase_two_cloud/scripts/submit_parameter_energy_sweep.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/UNIQUE_NAME" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/UNIQUE_NAME"
```

## Omiya comparison for the first four energy variants (October 29)

Production tasks 2–5 completed and passed at alpha=0.3, spin=0.99. The comparison
uses our existing multimode primary/secondary formalism for the Omiya curves:
six nonlinear cosine primary modes and three secondary modes on the averaged
primary background. Job **27417583 COMPLETED** in **3m13s** and produced all four
references in `results/omiya_energy_20261029_v1/`.

The primary is solved at each distinct E1 target; the half- and 1.5x-energy
profiles are new nonlinear solves seeded from the old continuation branch.
Its E0 solution is reconstructed and rechecked. The secondary is solved on
each resulting primary, then its linear equation permits normalization to E2.
No nonlinear primary or joint profile is rescaled. Measured Omiya energy errors
are at most `8.23e-7` relative. Primary matching uses tolerance `1e-9`; the
secondary retains the existing frequency-ULP tolerance rule, with the actual
residual and tolerance saved per case. Both clouds use the same energy convention,
alpha, spin, radial endpoints, and mode support as the previous reference.
The two E1=E0 exports reproduce the earlier complex profiles and slopes within
`3.6e-15` relative after the requested secondary normalization. These checks
do not establish spatial/angular/endpoint convergence.

Frozen source and inputs: `work/omiya_energy_20261029_v1`; the SHA-256 manifest
includes the seed branches and all four accepted joint-case metadata files.
Each primary solve and unit-secondary solution is checkpointed. The final JLS
preserves all retained modes; CSV exports the two dominant modes on the same
0.25-spaced grid as the joint results, using cubic Hermite interpolation of R/Q.
Log: `logs/omiya-energies-27417583.out`.

The repository-level `plots/parameter_energy_four_cases_omiya_20261029_v1.png`
and `.pdf` show joint (solid), Omiya (dashed), and linear (dotted) curves with
matched energy targets. Lower panels show
`100*(abs(R_joint)-abs(R_omiya))/max_displayed(abs(R_linear))`.
Individual-case figures, plotted-data CSV, and provenance JSON share this prefix.
The earlier linear-only figures are preserved. Reproduce with
`scripts/prepare_omiya_energy_references.py`,
`scripts/submit_omiya_energy_references.sh`, then
`scripts/plot_parameter_energy_cases.py --omiya-directory OUTPUT_DIRECTORY`;
each script's usage lists its positional arguments and requires fresh outputs.

## Spin 0.95 comparisons (October 29)

Production tasks **6 and 7 COMPLETED**, accepting alpha=0.3, spin=0.95 with
energy pairs `(1,1)` and `(0.5,1)` in units of E0. Their comparisons are in
`plots/parameter_energy_spin095_completed_20261029_v1.{png,pdf}` at repository
level, with individual figures, a plotted-data CSV, and input hashes alongside.
The line styles and lower-panel definition match the first four comparisons.

`omiya_background_references.jl` transports the existing multimode Omiya
solution from spin 0.99 in increments no larger than 0.02, solving the primary
and correcting its energy at every stage. It then solves the secondary on the
final primary and normalizes its quadratic energy to E2. Each intermediate
primary is checkpointed. Alpha and radial endpoints must match the seed;
this driver does not yet transport alpha or extend the domain.

The equal-energy Omiya export passed in job **27418199**, saved in
`results/omiya_spin095_20261029_v1/`, with maximum energy error `4.04e-7`.
That job subsequently failed during the linear-secondary export: its unit-amplitude
absolute residual stopped at `1.17e-8`, just above the requested `1e-8`.
The accepted Omiya files were preserved. Linear-only job **27419139 COMPLETED**
in 33s using `work/omiya_spin095_linear_20261029_v2`, writing
`results/omiya_spin095_linear_20261029_v2/linear_reference.{jls,toml,csv}`.
The linear-secondary residual changes by `5.48e-8` for one representable-frequency
step; its final residual is `1.64e-8`. The linear driver now measures this floor
near the root and uses `max(1e-8, 2*floor)`, recording both values and the actual
residual. Both linear energies are reintegrated after normalization and checked
within `1e-12` relative.

The half-primary-energy Omiya job **27419349 COMPLETED** in 6m17s, writing
`results/omiya_spin095_half_20261029_v1/` from the matching frozen work directory.
Its maximum energy error is `1.80e-7`. This job uses `--omiya-only`, reusing the
new spin-0.95 linear references for plotting. The submission wrapper also accepts
`--background` (both references) and `--linear-only`. Each snapshot includes the
accepted joint metadata, seeds, source, and SHA-256 manifest. The plotting driver
accepts repeated `--omiya-directory` arguments to combine separate reference
batches and rejects ambiguous case matches. Spatial/angular/endpoint convergence
remains a separate check.

## Equal-energy joint comparison (October 28)

At the user's request, job **27401326** adjusts the two outer amplitudes until
both cloud energies match the reference `650.5419114886128 F_a^2 M`, re-solving
the nonlinear radial equations after every adjustment. The initial joint
solution is `results/src_target_20261028_v2.jls`; its measured cloud energies
are **562.6763156010196** and **545.072913465489**, respectively **13.51%** and
**16.21%** below the references. Thus the fixed-outer-coefficient overlay also
compares different cloud populations.

The comparison reuses the Omiya reference energy convention: full isolated
cosine energy of all pure-primary harmonics for E1, and quadratic energy of
the secondary on the phase-averaged primary for E2. These are normalization
coordinates, not a unique split of the full interacting-field energy. Mixed
driven modes are excluded from these two labels; the full joint E and J are
reported separately. The energy evaluator reproduces both saved Omiya energies
within `3e-11` absolute at radial quadrature spacing 0.5.

The outer iteration uses a two-variable Broyden update in log amplitude/log
energy, initialized with the amplitude-squared scaling derivative. Every step
is capped at 0.2 in log amplitude and followed by a full radial solve and
independent residual/defect assessment. Targeting uses differential tolerance
`1e-5`; the final solution is refined to `1e-6`. Final energy errors must be
below `1e-4` relative, and quadrature spacings 1 and 0.5 must agree within
`5e-6`. Eight solves and a two-hour batch limit bound the attempt.

Frozen source and inputs: `work/equal_energy_20261028_v1`, verified by its SHA-256
manifest. Output prefix: `results/equal_energy_20261028_v1`; `-history.toml`
records each energy measurement, and `-stepN.jls` preserves corrected solutions.
Log: `logs/equal-energy-27401326.out`. The job **COMPLETED** in **41m53s**.
Accepted energies are `[650.5436220051258,650.5431992413581]`, within
`2.63e-6` relative of the target. Boundary residual is `1.82e-13`, independent
defect `3.38e-7`, and final tolerance `1e-6`. The accepted amplitudes are
`[0.04216226110599378,0.00012807425525547365]`. The resulting plot is
`plots/multiphase_deformation_20261028_v2.png` at repository level.

```bash
python3 runs/adiabatic/multiphase_two_cloud/scripts/prepare_equal_energy_target.py UNIQUE_NAME
sbatch runs/adiabatic/multiphase_two_cloud/scripts/submit_equal_energy_target.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/UNIQUE_NAME" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/UNIQUE_NAME"
```

The final CSV contains energy-normalized references already. Plot it without
`--energy-references`; that flag is only for the historical fixed-amplitude
export. The plot driver rejects double normalization and checks the measured
joint energies before accepting an equal-energy label. Its historical plot
output remains pixel-identical to the original.

## Current-source collocation comparison (October 28)

Job **27398811 COMPLETED** in 12m00s, reproducing the refined 14-mode cosine
target through the new `src/adiabatic/multiphase/load_collocation.jl` backend.
It holds spin 0.99, alpha 0.3, the 16×31×31 polar/phase grid, the interval
`[-100,300]`, boundary orders 3/7, and outer coefficients
`A1=0.036390480353373555`, `A2=0.00010896058341192524` fixed.

The 881-point solution passes independent assessment: boundary residual
`1.68835e-13`, parameter drift `9.69e-21`, and off-mesh defect `3.79897e-7`
against requested tolerance `1e-6`. Its largest frequency change from the
earlier refined checkpoint is `3.21794e-13`; the largest per-mode relative
profile change, including slopes, is `5.49e-7`. Solve work took 664.7 seconds,
including 127.7 seconds of compilation. Peak batch memory was about 8.7 GiB.

Source snapshot: `work/src_target_20261028_v2`; output prefix:
`results/src_target_20261028_v2`; log: `logs/src-target-27398811.out`.
The JLS preserves the problem, mesh, scaled states, physical samples, and
diagnostics; TOML records validation/settings; CSV contains the six physical
profiles on one grid with their complex R and Q values. The snapshot manifest
hash-checks source and all inputs. Earlier job 27398712 stopped before solving
on an old-checkpoint field-name mismatch; v2 validates the saved diagnostics
directly. Its failed output prefix is preserved.

The preliminary comparison uses the energy-matched `reference_curves_20261028_v2`
linear/Omiya curves, each with reference energy `650.5419114886128` in code
units. Linear amplitude factors relative to the old fixed-amplitude references
are `0.7255014992542931` and `0.7541792626953592`. Both R and Q are scaled.
The joint curves retain the Omiya outer coefficients; their individual energies
have not been independently retuned. No nonlinear profile is rescaled by hand.
The reference conversion agrees with all 12,008 original reference rows to
`2.78e-17` absolute. Angular, domain/boundary, and retained-mode systematics
remain separate checks.

Artifacts:

- [PNG comparison](../../../plots/multiphase_deformation_20261028_v1.png) and
  [PDF comparison](../../../plots/multiphase_deformation_20261028_v1.pdf).
- [Six complex curves and slopes](../../../plots/multiphase_deformation_20261028_v1-curves.csv).
- [Plotted data and signed changes](../../../plots/multiphase_deformation_20261028_v1-plot-data.csv).
- [Normalization, diagnostics, and source/input hashes](../../../plots/multiphase_deformation_20261028_v1-plot.json).

The upper panel plots six magnitudes without individual peak normalization.
For each cloud, lower panels use `D_X=(abs(R_X)-abs(R_linear))/max(abs(R_linear))`,
with the linear peak measured over the displayed interval `[-100,200]`.
The extra curve is `D_joint-D_omiya`. Plot data retains the full `[-100,300]`
interval and stores fractional differences; the displayed lower axes use percent.

To reproduce a run, choose a fresh snapshot/output name:

```bash
python3 runs/adiabatic/multiphase_two_cloud/scripts/prepare_collocation_target.py src_target_NEW
sbatch runs/adiabatic/multiphase_two_cloud/scripts/submit_collocation_target.sh \
  "$PWD/runs/adiabatic/multiphase_two_cloud/work/src_target_NEW" \
  "$PWD/runs/adiabatic/multiphase_two_cloud/results/src_target_NEW"
```

To redraw the completed result, choose a fresh figure prefix; existing figures
are preserved:

```bash
MPLCONFIGDIR=/tmp/axion_collocation_mpl PYTHONPATH=src \
  /global/home/users/clwelch/.conda/envs/BHsuper/bin/python \
  runs/adiabatic/multiphase_two_cloud/scripts/plot_collocation_comparison.py \
  runs/adiabatic/multiphase_two_cloud/results/src_target_20261028_v2 \
  plots/multiphase_deformation_NEW \
  --energy-references runs/adiabatic/multiphase_two_cloud/results/reference_curves_20261028_v2.toml
```

## Start here

This campaign tests implicit amplitude responses and a short E/J-balanced
trajectory of the joint four-component cosine model on a fixed Kerr background.
The October 15 pilot completed two trajectories and the derivative checks.
On October 19, a single primary-amplitude doubling run completed as job
**27100648**; see the run record below. It is separate from the pilot's two
time resolutions.
It is a finite-domain numerical pilot. It does not establish spatial or
Fourier-support convergence, boundary-tail accuracy, or the physical validity
of the frozen-profile slow-growth approximation.

The library is in `src/adiabatic/multiphase/response.jl` and `evolution.jl`.
`core/sensitivity.jl` supplies the general multiple-amplitude derivative.
The user authorized direct implementation on October 15, with explanatory
comments and docstring skeletons left for them to fill in.

## Amplitude continuation

Current source layout: the run script calls `continue_two_phase_amplitude` from
the multiphase library. Shared stepping and I/O are in
`src/adiabatic/core/continuation.jl` and `continuation_io.jl`; the model predictor,
correction, and radial validation are in `src/adiabatic/multiphase/profiles.jl`.
The script accepts an optional third argument, `OUTPUT_PREFIX`, for a fresh run
or restart. The two-argument batch command and per-attempt checkpoint schemas
are preserved. Frozen snapshots described below retain their original source.

## Model and numerical settings

- Kerr spin `0.99`, scalar mass parameter `alpha=0.42`, fixed background.
- Independent clouds have azimuthal integers `(1,2)` and initial outer
  normalization coefficients `(1e-4,1e-5)`; these are not cloud energies.
- Retained `(l,q)` pairs: `(1,(1,0))`, `(2,(0,1))`, `(0,(2,-1))`, `(3,(-1,2))`.
- Cosine potential, 8 polar nodes and 9 nodes for each independent phase.
- Tortoise endpoints `(-30,80)`, matching radius `5`, maximum radial step `0.4`.
- Horizon/outer series orders `3/7`. Newton scaled residual tolerance `1e-9`.
- Reference coefficient and residual scale `1e-5`, held fixed when amplitudes vary.

The trial vector contains both complex frequencies, all horizon coefficients,
and the free outer coefficients. Real and imaginary parts are packed into a
real vector `u`; the two positive amplitudes are `A=(A1,A2)`. The response solves
`G_u u_A = -G_A` with one factorization for both columns. The observable
response includes explicit and induced dependence: `dO/dA = O_A + O_u u_A`,
where `O=(E,J)`.

The stress tensor uses the real parts of the complex matching frequencies as
carrier frequencies. It omits slow profile time derivatives. Agreement with
finite differences verifies differentiation of this specified numerical
family, not the physical slow-growth approximation.

## Evolution and validation

The balance equations are `(dO/dA) dA/dt = -(F_E,F_J)`, with horizon plus outer
losses in each row. Horizon losses negate the signed inner radial current.
Columns are scaled by amplitudes and rows by their largest coefficient before
reporting the balance condition number; the code rejects ill-conditioned
systems instead of using a pseudoinverse.

Evolution uses explicit midpoint steps in `log(A)`. Every midpoint and accepted
endpoint gets an independent Newton correction, initialized using the previous
implicit profile response. Boundary losses are integrated separately with the
midpoint rule and compared with changes in total E and J.

The driver checks both response columns against independently corrected
neighbors at relative amplitude offsets 1% and 0.5%. It also changes the fixed
coefficient scale by a factor of seven and checks invariance of physical E/J
derivatives. It evolves with 4 and 8 time steps to a final time chosen for an
approximately 3% logarithmic change in the faster cloud.

Analytic unit tests check a nonlinear two-parameter family with explicit
observable dependence, singular response systems, E/J unit invariance, and
second-order convergence of amplitudes and integrated charge budgets.

## Commands and preservation

From the repository root, create a new source snapshot and output prefix:

```bash
python3 runs/adiabatic/multiphase_two_cloud/scripts/freeze_source.py SNAPSHOT_NAME
sbatch runs/adiabatic/multiphase_two_cloud/scripts/submit.sh \
    /absolute/path/to/work/SNAPSHOT_NAME /absolute/path/to/results/PREFIX
```

The submission wrapper uses `lr7`, runs the complete Julia unit suite, then
runs the campaign driver. `AXION_RUN_TESTS=0` skips the suite for a subsequent
spatial check using already-tested code. `AXION_RADIAL_STEP=0.2` selects a radial
refinement; use a distinct prefix. Source snapshots include a Git revision,
working-tree patch/status, and SHA-256 hashes. The submission wrapper verifies
those hashes before running. Snapshots live under ignored `work/`; preserve the
snapshot used for a campaign until its results have been superseded explicitly.

Each accepted endpoint updates `PREFIX-steps4.jls` or `PREFIX-steps8.jls` and
its CSV table. A resumed run continues from its last accepted endpoint. A
prefix is restricted to the same source manifest and spatial settings. Never
run two writers for the same prefix. An abrupt job kill can leave a `.lock`
directory; check that its job has ended before removing that stale lock.

The initial solution, neighbor-check table, scale check, and final comparison
summary are also saved separately. Times are `t/M`; E, J, and their losses have
the overall `Fa^2` divided out, matching the observable library.

Plot a completed run using the project's Python environment:

```bash
PYTHONPATH=src python runs/adiabatic/multiphase_two_cloud/scripts/plot_evolution.py \
    runs/adiabatic/multiphase_two_cloud/results/PREFIX \
    plots/multiphase_two_cloud_evolution.png
```

## October 15 run record

- Initial pilot: job **26977745**, submitted on `lr6`, moved while pending to
  `lr7`, running on `n0060.lr7`.
- Snapshot: `work/response_20261015_v1`; base revision `70f445a` plus the saved
  uncommitted implementation. All 176 snapshot files were hash-verified.
- Output prefix: `results/cosine_pilot_20261015`.
- The initial snapshot predates the automatic hash-verification/resume-configuration
  guards in the current driver. Its hashes were verified separately; resume that
  prefix only with its original snapshot. The source copy succeeded, but its
  provenance step initially encountered the host's Python 3.6 `text=True`
  incompatibility. Provenance was finalized separately without changing frozen
  source files. The current snapshot script uses `universal_newlines=True`.
- Job 26977745 passed all **985** Julia unit checks (13 loader, 350 core,
  243 single-cloud, 16 primary/secondary, 363 multiphase). It then stopped in
  the linear seed solve at residual `1.25747e-10`, above the requested `1e-11`.
  No trajectory was produced by that attempt.
- Corrected pilot: job **26978242** on `lr7`, snapshot
  `work/response_20261015_v2`, prefix `results/cosine_pilot_20261015_v2`.
  It uses the established `1e-9` matching tolerance and the completed unit
  suite (`AXION_RUN_TESTS=0`). It includes automatic snapshot verification and
  resume-configuration guards. Its 178 files were verified at job startup.
- Job 26978242 **completed successfully** in 4 minutes 47 seconds, with peak
  batch RSS 822012 KiB. The current Julia source/test files match its snapshot.
  Re-running its frozen driver after completion successfully loaded the saved
  results and reproduced the comparison without repeating trajectory steps.

### Completed pilot results

| Quantity | Result |
|---|---:|
| Final time, t/M | 199406.62115563318 |
| A1, initial → final | 1e-4 → 1.0304542260926823e-4 |
| A2, initial → final | 1e-5 → 1.000012690345167e-5 |
| First amplitude change | +3.045422609% |
| Second amplitude change | +0.0012690345% |
| Largest relative E/J neighbor-derivative discrepancy | 9.495e-9 |
| Largest relative state-response discrepancy | 2.064e-9 |
| Largest relative endpoint amplitude difference, 4 vs 8 steps | 7.056e-11 |
| Energy budget error / initial E, 8 steps | 4.836e-10 |
| Angular-momentum budget error / initial J, 8 steps | 2.459e-10 |
| Energy/angular budget improvement on halving dt | 3.99977 / 3.99959 |
| Scaled E/J condition range, saved endpoints | 1140.47–1210.96 |
| Largest saved endpoint matching residual | 1.965e-10 |

The eight-step budget errors are 2.084e-6 of the accumulated energy change and
1.885e-6 of the accumulated angular-momentum change. This second normalization
matters because the charges change by only a small fraction of their initial
values. Both charge budgets show the expected second-order time convergence.
The extremely close amplitude trajectories measure time-step sensitivity on
this fixed numerical family; they are not absolute physical accuracy estimates.

The frozen-scale check changes the coefficient scale by seven and finds a
relative matrix-norm change of 6.40e-20 in the E/J derivative matrix. The
matching Jacobian condition is about 6.31e7 in its chosen state coordinates;
this is distinct from the scaled E/J balance condition reported above.

Files to inspect:

- [Eight-step trajectory](results/cosine_pilot_20261015_v2-steps8.csv)
- [Four-step trajectory](results/cosine_pilot_20261015_v2-steps4.csv)
- [Neighbor derivatives](results/cosine_pilot_20261015_v2-neighbors.csv)
- [Evolution and balance plot](../../../plots/multiphase_two_cloud_evolution.png)

The associated JLS files preserve the initial problem, response matrices,
accepted endpoint histories, and summary. These numerical checks validate the
implementation on the stated finite grid. The next physical checks are radial
and domain refinement, phase/angular and Fourier-support convergence, tail
accounting, and a justified slow-growth approximation.

## October 19 primary-amplitude doubling

The "Review adiabatic workplan" chat requested one run to double the primary
amplitude and compare with linear growth. That chat ended during preparation;
no doubling job had been submitted. Job **27100648** now runs on `lr7`, using
snapshot `work/doubling_20261019_v2` (183 hash-verified files) and output prefix
`results/cosine_doubling_20261019_v2`. It completed successfully in 17m52s
(exit 0), reaching A1=2e-4 at t/M=4607986.0447. The secondary grew 0.029340%.
The primary ended 0.016464% above the linear exponential and 0.010753% below
the initial-E/J-rate exponential. Relative budget errors are 1.958e-7 in E
and 9.954e-8 in J, normalized by their initial charges.

The first attempt, job **27100618**, used snapshot/prefix ending in `v1`.
All 47 balance, midpoint, and amplitude-target checks passed, but the linear
reference calculation received unsupported boundary-order keywords and stopped
before any evolution. The corrected `v2` driver passes the supported radial
settings and explicitly checks the shared series orders 3/7. Both attempts
are preserved.

The corrected job also passed all 47 startup checks, corrected the initial
joint profile to residual `1.826e-11`, and completed both linear-reference
solves. It has entered the evolution and saved its initial endpoint. The
reference amplitude growth rates are `1.50387294386e-7` and
`6.37113618587e-11`; the corresponding initial E/J rates are
`1.50446358432e-7` and `6.36395540283e-11`. These already differ slightly at
the starting profile, which motivates keeping both exponential comparisons.

The run starts from the pilot's saved **initial** joint cosine profile, with
`A1=1e-4`, `A2=1e-5`, and unchanged spatial settings. The absolute stopping
target is `A1=2e-4`, with relative target tolerance `1e-9`. The initial balance
rate estimates a doubling time `t/M=4607271.2413`. The nominal time step is
that estimate divided by 64; steps shorten near the measured target. Every
midpoint and endpoint is independently profile-corrected. A safety time limit
of twice the initial doubling estimate and explicit stop reasons prevent an
unbounded run if the target cannot be reached. Only one trajectory is run;
there is no doubling time-step comparison in this campaign.

The submission first runs the balance/midpoint and amplitude-target tests,
then the evolution and plot. It requests one CPU, 4 GiB, and two hours.
Reproduce with a fresh snapshot and prefix, or resume the existing prefix
with its original snapshot after confirming that no writer remains:

```bash
python3 runs/adiabatic/multiphase_two_cloud/scripts/freeze_source.py \
    doubling_20261019_v2 \
    runs/adiabatic/multiphase_two_cloud/results/cosine_pilot_20261015_v2-initial.jls
sbatch runs/adiabatic/multiphase_two_cloud/scripts/submit_doubling.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/doubling_20261019_v2" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/cosine_doubling_20261019_v2"
```

The snapshot command refuses to overwrite an existing snapshot. For a resume,
run only the submission command. The existing configuration/locking rules
also apply to this driver. Abrupt interruption requires checking the job's
state before removing its stale prefix lock.

Saved outputs use the prefix above:

- `-steps64.csv` and `.jls`: every accepted endpoint, frequencies, rates,
  fluxes, conditioning, residuals, and cumulative charge budgets.
- `-linear-reference.csv`: the two independently solved linear frequencies
  on the same spatial grid, plus the initial E/J logarithmic amplitude rates.
- `-summary.jls`: measured doubling time and endpoint comparison/budget errors.
- `-comparison.png`, `.pdf`, and `.csv`: joint amplitudes compared with
  `A_a(0) exp(Im(omega_a_linear) t)` and with the initial E/J rate held fixed.
  The amplitude exponent is `Im(omega)`, not the energy exponent `2Im(omega)`.
- `logs/doubling-27100648.out`: progress and startup test results.

Regenerate the comparison with the initial dimensionless BH spin and coupling
in the title (both remain fixed throughout this run):

```bash
PYTHONPATH=src python runs/adiabatic/multiphase_two_cloud/scripts/plot_doubling.py \
    runs/adiabatic/multiphase_two_cloud/results/cosine_doubling_20261019_v2 \
    --spin 0.99 --alpha 0.42
```

No physical BH mass is assigned; the evolution uses mass-scaled units.

## Omiya-approximation comparison (October 19)

The user requested a comparison of the joint evolution with linear growth and
the existing Omiya nonlinear-primary/linear-secondary approximation. The old
Figure 6 trajectory has different initial cloud populations; its coefficient
table starts at A1=0.0011390625, above this entire doubling run. Consequently
`scripts/omiya_comparison.jl` recomputes a small branch with the existing
`primary_secondary/` solvers and evolves its energy/angular-momentum equations.
It does not extrapolate the saved Figure 6 tables.

The comparison starts at the same outer normalization coefficients
`A1=1e-4, A2=1e-5`, spin 0.99, and alpha 0.42. Matching outer coefficients
does not force the different approximations to have identical total charges.
It uses the joint pilot's domain `(-30,80)`, match position 5, one angular
function in each of the 211, 322, driven m=0 and driven m=3 components,
8 polar/9 azimuthal nodes, and boundary orders 3/7. Primary radial step 0.1
and secondary/driven step 0.2 place all secondary RK4 stages on saved primary
samples. The outgoing m=3 series is imposed at 80 to match the joint pilot's
prescription; this differs from the old Figure 6 reproduction's far-boundary
propagation. Thus this is an Omiya-approximation comparison on the pilot's
finite domain and support, not a new full-domain Figure 6 reproduction or a
converged measurement of mutual deformation. The two radial discretizations
also differ from the joint run's maximum step 0.4.

The new branch covers `A1=0.9e-4` through `2.2e-4` and explicitly includes
the initial amplitude. The evolution uses the established `two_cloud_rates`
with all four flux channels. It evolves E1 and E2, then recovers amplitudes
using the solved maps `A1^2/E1` and `E2/A2^2`, including the latter's change
with the primary. Times are t/M and coincide with the saved joint endpoints.
Checks compare 4/8 RK4 substeps per output interval and full/thinned branch
tables. A separate analytic test of the energy/amplitude maps, varying
secondary normalization, physical time, and domain guards passes 18 checks:

```bash
julia --project=. tests/adiabatic/primary_secondary/test_omiya_comparison_driver.jl
```

Job **27103035** uses frozen snapshot `work/omiya_comparison_20261019_v2` and
output prefix `results/omiya_comparison_20261019_v2`. The earlier job 27103000
stopped at driver loading because a file-header docstring attached to a
`using` statement; it produced no branch points. Both snapshots are preserved.
The corrected driver's analytic test passed before submission.

Job 27103035 **completed successfully in 1m01s** (exit 0, peak batch RSS
599304 KiB). All 18 branch points were solved. The largest relative amplitude
change under time-step halving was 4.66e-15, and thinning the branch table
changed amplitudes by at most 4.11e-14. These check this interpolation/time
integration, not the physical accuracy of the finite-domain approximation.

At the joint run's final time t/M=4607986.0447, the Omiya reference has
A1=1.9996657138e-4 and A2=1.0002939355e-5. The joint endpoint is 0.0167171%
higher in A1 and 0.0000539919% lower in A2. The reference's primary energy
grows from 0.000211641 to 0.000846279 in Fa^2 M units, while its secondary
stays near 0.0631. The joint total energy grows from 0.0632982 to 0.0639701.
This dilute pilot is far below the O(100–1000) cloud energies of the paper's
nonlinear turnover, despite using the same alpha and spin. Small differences
between these curves have not been isolated from discretization, finite-domain
balance, and slow-growth prescription differences.

The output prefix owns `-branch.csv`/`.jls`, `-trajectory.csv`, and
`-summary.jls`. Checkpoints record the frozen source manifest and joint CSV
hash, and only resume the same configuration. Reproduce with a fresh snapshot
and prefix; resume only after confirming that no writer remains:

```bash
python3 runs/adiabatic/multiphase_two_cloud/scripts/freeze_source.py \
    omiya_comparison_20261019_v2 \
    runs/adiabatic/multiphase_two_cloud/results/cosine_pilot_20261015_v2-initial.jls
sbatch runs/adiabatic/multiphase_two_cloud/scripts/submit_omiya_comparison.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/omiya_comparison_20261019_v2" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/cosine_doubling_20261019_v2" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/omiya_comparison_20261019_v2"
PYTHONPATH=src python runs/adiabatic/multiphase_two_cloud/scripts/plot_doubling.py \
    runs/adiabatic/multiphase_two_cloud/results/cosine_doubling_20261019_v2 \
    --spin 0.99 --alpha 0.42 \
    --omiya runs/adiabatic/multiphase_two_cloud/results/omiya_comparison_20261019_v2-trajectory.csv
```

With `--omiya`, the top panels compare all three amplitudes, and the bottom
panels show joint/linear and Omiya/linear amplitude ratios with axes zoomed
around one. All y axes show actual values without additive offset notation;
no subtraction of one is plotted.
The exported comparison CSV also contains joint/Omiya fractional differences.
The previous plot and CSV are preserved under the doubling prefix's
`-initial-rate-comparison` suffix.

### What follows this run

1. Measure departures from both exponentials. A difference from the fixed
   initial E/J rate isolates rate variation along this numerical family;
   neither comparison alone attributes that difference to physical nonlinearity.
2. Use selected saved profiles to check radial/domain and angular/Fourier
   support sensitivity at the size of any measured effect. The earlier
   short-run time-step result does not establish long-run or spatial accuracy.
3. Establish the physical slow-growth formulation before interpreting extended
   evolution. The current complex matching frequencies and frozen real-carrier
   observables still specify a numerical family with this open approximation.
4. The next two-cloud physics milestone is a converged change in the primary
   caused by a populated secondary. Doubling a dilute primary may show little
   of this effect; use the result to choose further continuation amplitudes.
   Phase-action closure for a third independent cloud, evolving BH parameters,
   self-gravity, and GW losses remain later work.

## Direct Fig. 7 profile attempt (October 23)

Prepared `scripts/direct_profile.jl` for one joint cosine solve at alpha=0.3,
spin=0.99, seeded from point 72 of the saved Omiya primary/secondary branches.
That reference has E1=E2=650.541911489 Fa^2 M, A1=0.03639048035337356, and
A2=0.0001089605834119252. The joint solve fixes those outer coefficients;
it does **not** constrain the joint cloud energies to equal the reference.
There is no amplitude scan or time evolution.

The retained components are q=(1,0) with l=[1,3,5], q=(3,0) with l=[3,5],
q=(5,0) with l=[5], q=(0,1) with l=[2,4,6], q=(2,-1) with l=[0,2,4], and
q=(-1,2) with l=[3,5]. The primary and secondary boundary guesses are copied
from their saved solutions, with the secondary scaled by A2; the two driven
blocks start at zero and are solved for jointly. Radial settings are
(-100,300), match=5, step=0.1, boundary orders 3/7; quadrature is 16x31x31.
This is a preliminary retained-support calculation, not a convergence result.

The scaled matching tolerance is explicitly 1e-5, with at most 12 Newton
iterations and no automatic relaxation or continuation. The saved reference
secondary has a residual floor of a few 1e-6 per unit outer coefficient.
The Slurm script allocates four hours; this is a cap, not a runtime estimate.

Prepare a fresh named snapshot, then submit from the repository root:

```bash
python3 runs/adiabatic/multiphase_two_cloud/scripts/prepare_direct_profile.py direct_fig7_20261023_v1
sbatch runs/adiabatic/multiphase_two_cloud/scripts/submit_direct_profile.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/direct_fig7_20261023_v1" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/direct_fig7_20261023_v1"
```

The preparation copies and hashes the three reference inputs along with the
source. It does not submit a job. Check input selection, coefficient mapping,
packing, and boundary assembly without integrating:

```bash
julia --startup-file=no --project=. runs/adiabatic/multiphase_two_cloud/scripts/direct_profile.jl \
    --check runs/adiabatic/multimode_secondary/results
```

These input checks passed. No target joint solve has been run by this preparation.
The output prefix owns `-input.jls`, a successful `-solution.jls` with full radial
histories, and `-profiles.csv` with complex/absolute R11 and R22 on both branches.
A failed attempt records `-failure.txt`; existing prefixes are never overwritten.
The log is `logs/direct-profile-JOBID.out`. A successful solve supplies the two
new profile curves; the final comparison still needs explicit normalization
labels and numerical-convergence checks.

The live driver was subsequently rebuilt with the user as sequential commented
sections. Running it without arguments now performs setup checks only; the
`--check INPUT_DIRECTORY` and `INPUT_DIRECTORY OUTPUT_PREFIX` interfaces remain
available. The latter runs one Newton solve and saves the input, complete
solution, and profile CSV. The `direct_fig7_20261023_v1` snapshot preserves the
earlier draft, so create a fresh snapshot name to run the rebuilt driver.

Submission **27241322** failed after four seconds, before Julia started: the
shell received whitespace as its snapshot-directory argument and `cd` failed.
The submitted command referred to v1. The wrapper now checks its argument count
and directory explicitly. Snapshot `direct_fig7_20261023_v2` preserves the rebuilt
driver and the same reference inputs for the next attempt; preparing it did not
submit another job.

Submission **27241421** uses v2. Snapshot verification and setup passed; its
initial scaled matching residual was **1.9064826994960344e9**, first observed
after about ten minutes of elapsed job time. This is the initial mismatch,
not an accepted correction or a converged result.

The live driver now times every residual evaluation and flushes each message.
Messages distinguish the initial evaluation, ForwardDiff evaluations used to
build each Jacobian, trial residuals, and accepted Newton iterates. Each finished
evaluation reports its duration and scaled residual norm. The driver also
announces input saving and the final integration for radial output. It uses the
same residual, admissibility check, and Newton settings as
`correct_two_phase_profile`; only the driver adds the logging. A single
integration can still take minutes between its start and finish messages.
Already-running jobs and the v1/v2 snapshots keep their original logging.
Setup checks and 17 focused logging checks passed: scalar and 56-variable
analytic problems followed exactly the same residual-evaluation sequence and
returned the same Newton result with and without logging, including damped
trials and multiple ForwardDiff passes; evaluation errors still propagated.
These checks do not run or validate the target physical solve. Prepare a new
snapshot to use the extra logging; job 27241421 was not restarted.

Submission **27241725** runs the added logging from the verified 193-file
snapshot `work/direct_fig7_20261023_v3`, using the submit wrapper inside that
snapshot. Its output prefix is `results/direct_fig7_20261023_v3` and its log is
`logs/direct-profile-27241725.out`. Slurm reported RUNNING after submission.
The physical inputs, numerical settings, and four-hour limit are unchanged;
job 27241421 was left running with its separate v2 outputs.

Job 27241725 reported the initial residual after **537.04 seconds** of evaluation,
again `1.9064826994960344e9`, and began its first Jacobian evaluation at 537.73 s
on the matching timer. This is consistent with the substantially larger source
calculation: the earlier pilot used 4 radial functions, an 8x9x9 grid, and 276
RK4 steps across both branches; this attempt uses 14, 16x31x31, and 4000.
The projection-loop work count per residual increases by about 1204 times.

A small warmed-source timing check is preserved in
`work/measure_direct_profile_source_20261023.jl`. It uses the v3 source and
synthetic coefficients at one radius, not a physical solution or a full solve.
Median times over five samples on the interactive host were 0.493 ms at the
old dimensions and 33.244 ms at the current dimensions; multiplying by RK4
stage counts estimates 0.544 s and 531.903 s of source work per residual.
With 12 derivative directions, the current source took 42.711 ms. The installed
ForwardDiff selects width 12 for 56 unknowns, giving five residual passes per
Jacobian; a source-only extrapolation is 683.379 s per pass and 3416.896 s per
Jacobian. Actual derivative integrations and Newton convergence remain to be
measured; these timings are estimates, not a promised completion time. Live
Slurm accounting showed CPU time close to elapsed time for both jobs, with
less than 1 GB reported peak resident memory at that check.

The first actual ForwardDiff residual pass then finished in **832.78 s**;
pass 2 began at 1370.69 s on the matching timer. Its unchanged residual is
expected because derivative passes evaluate at the same Newton state. Five
passes at this observed cost would take about 69 minutes for the Jacobian,
before any trial residual evaluations. Later passes may be faster after
compilation; no Newton correction had yet been accepted at this checkpoint.

Pass 2 completed in **839.90 s**, and pass 3 began at 2210.59 s. The first two
passes therefore have comparable costs; the long wait is not just compilation.

### Source performance investigation (October 23)

At the user's request, standalone alternatives were tested in
`work/cached_source_experiment_20261023.jl`, loading the frozen v3 library.
No production library or running snapshot was modified. All **91 checks**
passed: source agreement on two grids (including unequal, shifted phase grids),
dilute and nonlinear amplitudes, all three potentials, analytic quartic source
coefficients, mixed/empty reconstruction inputs, and derivatives with respect
to all coefficient components, radius, spin, and alpha.

The reconstruction's runtime field-type selection prevents efficient
specialization of its inner loop. A separate helper receiving the allocated
field array preserves the existing type-selection procedure and gave exactly
equal source values and Jacobians in the checked cases. Precomputing real and
imaginary angular/phase basis matrices was also tested; the cached version's
source Jacobians agreed to about 1.4e-15 relative error, with 3.69 MB of cache
storage. These are discrete source checks, not full radial or profile solves.

Latest warmed medians over five samples at the target dimensions, with one
BLAS thread on the interactive host:

| Source evaluation | Original | Reconstruction helper | Cached basis |
|---|---:|---:|---:|
| Ordinary coefficients | 27.545 ms | 10.299 ms | 0.892 ms |
| 12 derivative directions | 64.142 ms | 19.174 ms | 20.813 ms |
| Ordinary allocated bytes/call | 14,040,840 | 376,360 | 616,888 |
| Derivative allocated bytes/call | 80,486,344 | 4,819,680 | 8,005,472 |

The small helper split was the initial suggested implementation step: it improves
the expensive derivative source by about 3.3 times in this sample without
changing arithmetic order. Ordinary cached-source speedups varied from about
31 to 42 times across samples, but derivative speedups were only about 3 times;
those ordinary timings must not be presented as whole-solve speedups. Library
integration and a full matching comparison remain to be done together.

The user then requested an order-of-magnitude improvement. A second standalone
experiment, `work/batched_source_experiment_20261023.jl`, caches the combined
real basis and compares two derivative strategies. Packing automatic
derivative directions into BLAS matrix columns helps, but the larger gain comes
from forming the local source Jacobian with respect to the 28 real field
coefficients once, then applying it to every incoming direction. Scalar
automatic differentiation of `potential_remainder` preserves its small-field
polynomial branch. Background derivatives with respect to radius, spin, and
alpha are included in the experiment as well.

All 84 grouped checks passed for source values and Jacobians on both grids,
two amplitudes, all three potentials, and chunk widths 12 and 31. Twelve
additional assertions compared arbitrary incoming directions at widths 12,
28, and 56. Warmed five-sample medians on the same interactive host, using
one BLAS thread, were:

| Source evaluation | Original | Batched basis products | Local source Jacobian |
|---|---:|---:|---:|
| Ordinary | 32.008 ms | 0.872 ms | Uses ordinary cached path |
| 12 derivative directions | 44.970 ms | 4.614 ms | 3.301 ms |
| 28 derivative directions | 85.584 ms | 12.270 ms | 4.026 ms |
| 56 derivative directions | 176.087 ms | 29.012 ms | 3.378 ms |

The current five width-12 passes require about 224.851 ms of source work per
radial stage in this benchmark. One width-56 local-Jacobian pass requires
3.378 ms: about **66.6 times faster for that source work**. This is not yet a
measurement of the full matching Jacobian or Newton solve. It keeps the same
angular quadrature, retained modes, and potential; no library or running
snapshot has been changed.

The follow-up `work/short_radial_jacobian_experiment_20261023.jl` loads the
saved v3 input and compares derivatives through the actual boundary helpers,
linear radial operators, and RK4 integration, with the candidate source and
width 56 versus the original source and width 12. Each branch covers only
0.2 in tortoise radius from its actual endpoint, using the target step limit.
The endpoint difference is a numerical test observable, not the physical
matching residual. All **58 checks** passed, including separate comparisons
of all 56 Jacobian columns at relative tolerance 1e-9. Combined with the 96
source checks, this gives 154 passing assertions for the second experiment.

Warmed five-sample medians were **5.849 s** for the original short radial
Jacobian and **0.105 s** for the candidate: **55.8 times faster**. Allocations
fell from about 8.06 GB to 111 MB per Jacobian (cumulative allocated bytes,
not peak resident memory). This establishes a substantial speedup through
the short integration pipeline, not the full radial domain or Newton solve.
The source tests include nonlinear fields; the short endpoint test uses the
saved starting profile, whose boundary fields are dilute. Full-domain
matching accuracy, runtime, and Newton convergence still need validation
before adopting this prototype in the library.

### Guided implementation and job checkpoint (October 23, 17:25 PDT)

The [October 23 work summary](../../../docs/multiphase-work-2026-10-23.md)
records the normalization decisions, driver construction, optimization
experiments, and exact implementation resume point. The live `projection.jl`
now contains the optional `source_cache` adapter, `cached_two_phase_source`,
and the ordinary matrix `local_chain_source` method. The cache builder,
optimized derivative method, radial/matching forwarding, and Newton batching
option remain pending. Seventeen existing nonlinear-source checks and 21
adapter checks with the experimental helper passed; the ordinary matrix
snippet passed 36 value/derivative checks. These edits do not affect frozen
v2/v3 or their running jobs.

At 17:25 PDT, job **27241421** was RUNNING at 2h13m07s and **27241725** was
RUNNING at 1h56m38s. Both logs record accepted Newton iteration 1 at scaled
residual **6.832061250638738e8**, down 64.2% from
**1.9064826994960344e9**, but still far above the requested 1e-5 tolerance.
Only the v2/v3 `-input.jls` checkpoints exist; neither job has saved a solution
or profile CSV.

The v3 first Jacobian ran from matching-timer 537.73 to 4776.00 seconds,
about **70.6 minutes**. Its five passes took 832.78, 839.90, 835.22, 871.10,
and 859.09 seconds. The full-step trial took 554.21 seconds and increased the
residual to 2.253960839827775e19, so it was rejected. The next, half-step trial
took 547.12 seconds and was accepted at timer 5877.72 seconds. Correction 2's
first Jacobian pass then took 865.25 seconds; pass 2 had started at this check.
No target convergence or final profile overlay has been established. Both
jobs were left running under their four-hour limits.

### Final status of the October 23 attempts (checked October 26)

Both target jobs reached their four-hour time limits on October 23:

| Job | State | Elapsed | End time (PDT) | Last accepted scaled residual |
|---|---|---|---|---:|
| 27241421 (v2) | TIMEOUT | 4h00m19s | October 23, 19:13:04 | 4.6531599593434143e8 |
| 27241725 (v3) | TIMEOUT | 4h00m26s | October 23, 19:29:40 | 4.6531599593434143e8 |

Both logs record accepted Newton residuals
`1.9064826994960344e9 -> 6.832061250638738e8 -> 4.6531599593434143e8`.
This is a 75.6% reduction from the initial mismatch, but is still far above
the requested `1e-5`. The `TIMEOUT` state and explicit Slurm time-limit
messages establish termination; the accounting `ExitCode=0:0` does not mean
that these solves converged.

The detailed v3 second Jacobian ran from timer 5877.72 to 10156.81 seconds
(71.3 minutes). Four evaluated trials then gave residuals
`1.2565082495518079e11`, `4.948017400002172e9`,
`1.9504016356406722e9`, and `4.6531599593434143e8`; only the last was accepted.
Those trial evaluations took about 35.4 minutes in total. Accepted iteration
2 was logged at timer 12279.90 seconds (204.7 minutes). Correction 3's first
two Jacobian passes completed in 868.82 and 857.38 seconds, and pass 3 had
started when Slurm stopped the job. The v2 log has the same two accepted
corrections but no per-pass timing.

An October 26 file check found only the v2/v3 `-input.jls` files. Neither
prefix has a saved solution, profile CSV, or accepted-iterate checkpoint.
The driver only saves its initial state and successful final solution, so
the two accepted numerical states cannot be resumed from its saved artifacts.
The logs retain their norms, not their frequency/boundary-coordinate vectors.

The working source still has the cache adapter, dictionary wrapper, and
ordinary matrix source installed on October 23. The optimized derivative
method and cache builder remain experimental; neither timed-out job used
them. Today's priorities are to finish the source optimization incrementally,
connect it to the full matching Jacobian, add accepted-state checkpointing,
and measure full-domain agreement/runtime before a fresh frozen solve.
The final Fig. 7 overlay remains dependent on a converged joint profile;
speedup alone does not establish matching or physical discretization accuracy.

### Source cache type and local derivative checkpoint (October 26)

Guided implementation now includes `build_two_phase_source_cache` and
`local_source_jacobian`. The user then authorized introducing an immutable
`TwoPhaseSourceCache` container for `pairs`, `basis`, `weights`, and
`cosine_theta_squared`. Its constructor preserves the supplied storage types;
the builder returns this type and the three cached evaluators accept it
explicitly. Arrays are treated as read-only after construction.

All 87 focused checks passed against the actual library methods: 30 cache,
field-reconstruction, and source comparisons; 40 local-Jacobian comparisons
against generic automatic differentiation; and the 17 existing nonlinear-source
checks. Repeated opt-in loading also succeeded. The checks reused the fixtures
and assertions in the October 26 scratch scripts without loading their
prototype function definitions. The fast derivative-specific dispatch,
radial/matching cache forwarding, and full-width Newton batching remain pending.
No full-domain benchmark or new solve was performed for this type change.

### Source method cleanup (October 26; superseded interface)

This intermediate interface and scratch runner were superseded by the
mandatory-cache change below.

The source calculation now has two methods named `nonlinear_two_phase_source`:
the dictionary interface, with an optional cache, and the packed-vector method
using a `TwoPhaseSourceCache`. The dictionary method handles packing and
unpacking directly; the separate `cached_two_phase_source` wrapper was removed,
and `local_chain_source` became the vector method. `build_two_phase_source_cache`
retains its name, and the cache's ordered mode labels are stored in `mode_labels`.

The scratch runner `work/source_method_cleanup_check_20261026.jl` reuses earlier
fixture/assertion sections without loading their prototype implementations.
Run it from the repository root with:

```bash
julia --project=. runs/adiabatic/multiphase_two_cloud/work/source_method_cleanup_check_20261026.jl
```

All 111 checks passed: 30 cache/reconstruction/source comparisons, 40 local
Jacobian comparisons, 21 dictionary-adapter checks, 17 existing nonlinear-source
checks, and three checks of fallback for an input mode outside the cache.
Repeated loading also succeeded. These establish local source and derivative
agreement, including mode ordering, missing inputs, and cache fallback; they
do not establish full-domain runtime or matching convergence. The optimized
ForwardDiff method, cache forwarding, Newton batching, and accepted-state
checkpointing remain the next implementation steps.

### Required precomputed source basis (October 26)

The production source now requires `TwoPhaseSourceBasis`, which stores fixed
angular/phase basis and quadrature data, not past field or source evaluations.
The sole source interface is
`nonlinear_two_phase_source(parameters, source_basis; potential=:cosine)`.
The dictionary interface and direct-quadrature fallback were removed.
`TwoPhaseProblem` builds the basis data once by default, and the radial RHS
passes its packed coefficient vector to the source evaluator at each stage.
Matching and response evaluations reuse that same object. A basis with the
wrong radial mode order is rejected; changes to mode/phase/quadrature settings
require rebuilding it.

The multiphase test command was:

```bash
julia --project=. -e 'using LinearAlgebra; BLAS.set_num_threads(1); include("tests/adiabatic/multiphase/runtests.jl")'
```

It passed 388 checks through source, radial integration, Newton matching,
observables, and profile coordinates, including 38 new source/basis checks
against explicit direct quadrature and its derivatives. It then stopped on
two existing errors in the first evolution testset: the test expects
`amplitude_rates`, while `two_cloud_balance_rates` returns `log_rates`.
The committed function was checked independently and returns the same fields;
the evolution implementation and tests were not changed. Later evolution tests
were not reached by this run.

The user then requested the descriptive names `TwoPhaseSourceBasis`,
`build_two_phase_source_basis`, and `source_basis` throughout the library and
current documentation. After that rename, 147 focused source, RHS, matching,
and profile-coordinate checks passed, and repeated loading succeeded.

The direct-quadrature reference now lives in `test_source_basis.jl`. The tests
establish source/derivative agreement and short integration/matching behavior;
no full-domain benchmark or production solve was run. The optimized derivative
method, full-width Newton batching, and accepted-state checkpoints remain pending.

### v4 check and driver rollback (October 26)

User-submitted job **27296945**, frozen in `work/direct_fig7_20261026_v4`,
matched the old initial residual in 14.63 s and built the full-width Jacobian
in 92.51 s. The first directional check reached 2.38e-4 relative agreement,
above the assistant-added 1e-5 threshold; the job stopped before Newton.
See `logs/direct-profile-27296945.out`. The discrepancy remains unresolved.

At the user's request, the live direct and submit drivers were restored to
their prior versions, removing full-width batching, accepted-state saving,
and the preflight invocation; the live preflight script was deleted. The
library optimization and Newton options remain. Preserve v4 and its log as
historical evidence. Resume toward a convergence test one reviewed block at
a time; further driver additions must be walked through with the user.

### Reviewed v5 preparation (October 26)

After the source walkthrough, the user added
`jacobian_chunk_size=length(initial)` to the existing `newton_solve` call in
`logged_profile_newton`. This is the only change to the restored driver.
Five checks against the actual function passed on the 56-variable analytic
problem `residual(x) = 2x - target`: the known root, residual tolerance, one
Newton correction, exactly one derivative evaluation carrying 56 directions,
and preservation of the input vector. The setup-only command above also
passed with point 72, the same prescribed amplitudes, 14 radial modes,
16x31x31 quadrature, step 0.1, endpoints -100/300, match 5, boundary orders
3/7, matching tolerance 1e-5, and at most 12 corrections. The focused
`test_source_basis.jl` checks passed 38/38; these do not resolve v4's
full-domain derivative discrepancy.

Prepared `work/direct_fig7_20261026_v5` using the existing prepare script.
All 194 manifest files verify, the frozen driver matches the checked driver,
and all three reference-input hashes agree with v4. The user submitted this
snapshot as job **27304701**. The command from the repository root was:

```bash
sbatch runs/adiabatic/multiphase_two_cloud/work/direct_fig7_20261026_v5/runs/adiabatic/multiphase_two_cloud/scripts/submit_direct_profile.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/direct_fig7_20261026_v5" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/direct_fig7_20261026_v5"
```

This is a direct Newton convergence diagnostic using the existing four-hour
allocation. It has no preflight gate or accepted-iterate saving. Success
still requires the reintegrated scaled residual to be at most 1e-5 and the
saved solution/profile files; neither submission nor a residual decrease
establishes convergence.

Slurm confirms job 27304701 uses the frozen v5 wrapper, snapshot, and output
prefix shown above. Its first observed state was **PENDING**, with no log or
numerical outputs yet. The expected log is `logs/direct-profile-27304701.out`.
It subsequently started and verified all 194 snapshot files. The initial
scaled residual was **1.906482699496033e9** in **12.44 s**, identical to v4.
Correction 1's first Jacobian evaluation began at matching-timer 13.45 s.
At a subsequent **18-minute** Slurm check it remained RUNNING. The inspected
log recorded accepted scaled residuals `1.90648e9 -> 6.83206e8 -> 4.28028e8
-> 3.12652e8 -> 2.82101e8 -> 2.82100e8 -> 2.82073e8 -> 2.80571e8`;
correction 8 had started. The first seven Jacobians took 51–70 s each.
The factor-6.80 decrease still leaves the residual `2.81e13` times the target.

**Final v5 outcome:** job 27304701 FAILED with exit 1 after **29m04s**,
reporting `Newton iteration limit reached` after 12 accepted corrections.
The final scaled residual is **1.4947468124727777e8**, versus target 1e-5.
Every Jacobian used one pass, taking 51.41–69.71 s. Selected accepted residuals
and elapsed times from the matching timer are:

| Correction | Elapsed seconds | Scaled residual |
| --- | ---: | ---: |
| 0 | 12.63 | 1.906482699496033e9 |
| 1 | 105.75 | 6.832061138791364e8 |
| 3 | 317.89 | 3.126519191605745e8 |
| 7 | 932.05 | 2.805706872654804e8 |
| 9 | 1220.28 | 1.5729072959620476e8 |
| 12 | 1721.05 | 1.4947468124727777e8 |

The complete trace is in the job log. Only `-input.jls` and `-failure.txt`
exist for the v5 output prefix: no converged solution, profile CSV, or saved
accepted iterate is available. This establishes faster Jacobian construction
and failure of this direct attempt, not derivative correctness or convergence.
Preserve the attempt and review the diagnosis before selecting another solve.

Residual interpretation: the logged norm is
`norm([R_horizon-R_outer; Q_horizon-Q_outer]) / coefficient_scale`, with
`Q=dR/dr_*`, M=1, and fixed `coefficient_scale=A2=1.0896058341192524e-4`.
Thus correction 7's unscaled mismatch is about 30571.1, versus a target of
`1.08961e-9`. This normalization is not a percentage profile error.
A diagnostic reintegration of `direct_fig7_20261026_v5-input.jls` with the
frozen v5 source and `save_radial_result=true` reproduced the initial norm:
unscaled 207731.4672. At the matching point, the initial field norms are
0.0505503 from the horizon and 129715.5172 from the outer boundary; the field
mismatch divided by the larger field norm is approximately one. The largest
field mismatch is the secondary `(l,q)=(6,(0,1))`, about 129709.6888. Slope
norms are 0.00708174 and 162253.6504 respectively. These measurements establish
a severe initial joining failure, not its cause or the current iterate's
component breakdown: accepted intermediate profiles are not saved.

Further diagnosis (October 26): single residual evaluations using frozen v5
source, the saved input, and original reference inputs give the following
**unscaled** matching norms. No Newton corrections were run in these controls.

| Equations / boundary guess | Matching norm |
| --- | ---: |
| Original primary | 2.0332764801e-12 |
| Original averaged secondary, unit outer amplitude | 3.1578141557e-6 |
| Full joint system, secondary boundary coefficients zeroed | 2.1184413690e-12 |
| Full joint system, secondary coefficients multiplied by 0.1 | 900.4051881 |
| Full joint system, intended initial guess, step 0.1 | 207731.4672 |
| Same intended guess, step 0.05 | 207731.7767 |

The original secondary norm becomes `3.4408e-10` at the prescribed A2.
The halved-step joint norm differs by `1.49e-6` relative, so coarse radial
spacing does not explain the enormous initial residual. Along the full joint
guess's inward integration, secondary `abs(R_62)` grows from `4.96571e-5`
at tortoise radius 100 to `0.0378649` at 40, `5.31502` at 20, `807.767`
at 10, and `129709.689` at 5. The original reference matches, and the full
joint calculation matches in the zero-secondary limit, while the finite
secondary changes the equations and produces a rapidly growing inward branch
with the old boundary data. This supports an unsuitable direct guess amplified
by inward shooting as the immediate diagnosis. It does not isolate which
additional coupling dominates, rule out all implementation errors, or establish
that amplitude continuation will converge. The full-domain Jacobian discrepancy
remains a separate unresolved issue for Newton's progress.

### Weak-secondary seed and derivative validation (October 26)

Following the direct failure and controls above, the user authorized a matched
weak-secondary seed, a Jacobian check there, and small secondary-amplitude
continuation steps after those checks pass. The standalone
`scripts/weak_secondary_profile.jl` reuses the saved v5 problem; all Julia
library sources were checked byte-for-byte against v5 before reuse.
The primary amplitude remains `0.036390480353373555`; the secondary starts
at `1.0896058341192524e-7`, or 0.1% of its target. All spatial settings remain
those of v5. Newton still allows 12 corrections with full-width ForwardDiff.

Coordinates retain v5's coefficient scale. Residual rows belonging to components
with nonzero second-phase index are additionally divided by 0.001, so their
physical mismatch is normalized by the weak secondary amplitude. Primary-only
rows keep v5's scale. The stopping norm is `1e-5`; at the target amplitude this
row weighting becomes exactly the original residual. A converged profile is
reintegrated and saved before derivative validation, including both radial
branches and source provenance.

The diagnostic compares the optimized full Jacobian with two ordinary-source
AD directional derivatives. A diagnostic-only ForwardDiff tag invokes the
existing generic source method, bypassing the optimized local-Jacobian adapter
without changing the library. It also evaluates centered differences at steps
`1e-6`, `1e-7`, and `1e-8`. Required agreement is `1e-8` for the two AD routes
and `1e-5` for at least one finite-difference step in each direction. The check
records failures as well as successes; no tolerance is automatically relaxed.

Before submission, the setup check passed; the derivative-check routine passed
a nonlinear two-variable known-answer check; the diagnostic tag agreed with
the optimized local source derivative; and all 38 direct-quadrature source
checks passed. The submission shell also passed syntax checking. The existing
freeze utility saved the v5 input as `inputs/initial.jls`; all 195 snapshot files
verified. Submitted as job **27305905**, confirmed RUNNING:

```bash
sbatch runs/adiabatic/multiphase_two_cloud/work/weak_secondary_20261026_v1/runs/adiabatic/multiphase_two_cloud/scripts/submit_weak_secondary.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/weak_secondary_20261026_v1" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/weak_secondary_20261026_v1"
```

The allocation is one CPU, 8 GB, one hour. Log:
`logs/weak-secondary-27305905.out`. A successful seed is written to
`results/weak_secondary_20261026_v1.jls`, with derivative diagnostics in
`results/weak_secondary_20261026_v1-jacobian.jls`. Continuation has not yet run.

After four seed corrections showed slow progress, the user selected an immediate
independent Jacobian comparison at the original weak-secondary guess while the
seed attempt finishes. Job **27306004** runs `submit_weak_jacobian.sh` from
snapshot `work/weak_jacobian_20261026_v2`, with all 196 files verified. It reuses
`weak_secondary_setup` and `check_matching_jacobian`, runs no Newton solve,
and writes `results/weak_jacobian_20261026_v2.jls`. The one-CPU allocation is
capped at 30 minutes; log `logs/weak-jacobian-27306004.out`. The unused v1
snapshot predates a correction to the diagnostic's printed state description;
it was not submitted. The user requested explanations and updates roughly
every minute, with a decision point before any amplitude ramp.

**Stopping point at the user's departure request:** Slurm confirms both jobs
finished. No continuation, retry, or further diagnostic was submitted.

- **Seed job 27305905:** FAILED (exit 1) after 14m00s at 12 corrections.
  The weighted residual decreased from `1.3255350665e7` to
  `2.274196817352566e-5`, 2.274 times the `1e-5` target. Its last three accepted
  values were `124.0105717`, `0.001374343783`, and `2.274196817e-5`.
  This is close but does not satisfy the requested tolerance. No final seed,
  radial profile, or accepted-state checkpoint was saved; the post-convergence
  derivative check was not reached. The input snapshot and full log remain.
- **Jacobian job 27306004:** FAILED (exit 1) after 2m31s at the first direction's
  finite-difference check. Ordinary-source and optimized AD agreed to
  `7.4778485345e-15` relative. Centered-difference errors were
  `0.51032445`, `0.30443178`, and `0.16355185` at steps `1e-6`, `1e-7`,
  and `1e-8`. The second direction was not evaluated. The diagnostic saved
  the original weak state, full optimized Jacobian, first direction, and these
  errors in `results/weak_jacobian_20261026_v2.jls`.

For the saved first direction, `norm(J*direction)=5.8345761248e17`. Multiplying
by even the smallest probe step gives a predicted mismatch change 440.17 times
the initial residual `1.3255350665e7`. These parameter steps are therefore not
demonstrably small for this sensitive shooting map; finite-difference failure
does not identify an optimized-derivative bug. The independent AD agreement
is strong evidence for that direction only, not full spatial convergence or
all Jacobian directions.

On return, review a minimal accepted-state save using the existing Newton
callback and a modestly larger iteration cap before rerunning the weak seed.
The nearly matched state cannot be resumed from this attempt. Keep the
residual tolerance unchanged, and use the saved Jacobian to choose smaller,
properly scaled difference probes. No such change or rerun has yet been made.

Reference preparation: the older `plot_fig7_analogue.py` uses the smallest
saved-amplitude profiles as linear approximations and rescales them to equal
energy. Today's fixed-outer-coefficient comparison instead needs actual
`solve_linear_mode` profiles normalized to A1/A2. Reuse the point-72 primary
and averaged-secondary reconstruction routines, keeping the secondary's A2
scaling; do not reuse the old equal-energy linear rescaling.

A read-only audit of the frozen inputs confirms spin 0.99, alpha 0.3, the
averaged secondary prescription, and primary radial step 0.1. The secondary
routine takes twice that spacing, so its saved profiles use step 0.2. At
point 72 the saved primary residual norm divided by A1 is 5.58739e-11; the
unit-outer-amplitude secondary residual norm is 3.15781e-6. The secondary
checkpoint also contains an actual linear solution with saved residual norm
6.15553e-9. These are saved-solve diagnostics, not fresh reintegration or
refinement results; reference-grid accuracy remains to be checked.

The reference driver is now being written by the user in small reviewed,
tested blocks, including the driver code itself. The next proposed block,
`primary_profiles_at(primary_setup, primary_index)`, passed Slurm check
**27305885** using `work/check_reference_primary_20261026_v1.jl`.
The log `logs/reference-primary-check-27305885.out` reports scaled mismatch
5.587385657966086e-11, residual change from the saved solve
8.671068890112971e-29, and finite branch arrays of sizes 12x1051 and 12x2951.
The six primary modes supply six R and six Q rows; both branches end at
the match. This reproduces the saved primary on its original step-0.1 grid;
it does not establish radial refinement. No reference curves were exported.

### Full derivative verification and matched weak seed (October 27)

The user authorized completing verification autonomously after the guided
construction of `src/adiabatic/multiphase/autodiff_vs_optimized_test.jl`.
See [the verification record](../../../docs/multiphase-ad-verification-2026-10-27.md)
for acceptance criteria, tests, artifacts, and current results.

Frozen snapshot `work/ad_verification_20261027_v1` verifies 199 files and
contains the original v5 input at `inputs/initial.jls`. No numerical library
algorithm was changed; the user's `two_phase_real_matching_residual` rename
is included. The 55 focused local checks passed, including independent direct
quadrature and a short full-matrix matching regression.

- **27333738:** COMPLETED after 33m44s; full 56-column ordinary/optimized AD
  comparison at the initial weak guess and both finite-difference directions
  pass (worst column relative error `8.87909e-9 < 1e-8`); prefix
  `results/ad_initial_20261027_v1`.
- **27334006:** weak-seed retry, COMPLETED after 22m45s. Correction 13 reached
  weighted residual `5.306800040310494e-6`, reproduced by fresh reintegration.
  Saved `results/ad_seed_20261027_v1-seed.jls` contains the matched profile and
  both radial branches; `-accepted.jls` and `-setup.jls` preserve the last state
  and settings/provenance. A1 is fixed and A2 remains 0.1% of its target.
- **27334124:** COMPLETED after 33m19s at 15:59:48 PDT; matched-state
  verification reproduces the seed norm, and all 56 columns and both
  finite-difference directions pass (worst column `9.58781e-9 < 1e-8`); prefix
  `results/ad_matched_20261027_v1`.

All three completed with exit `0:0` on n0060.lr7. Their shared artifacts survived
the n0003 SSH outage. Recovery verified the manifest, loaded both completed
AD/finite-difference records, checked the seed and its provenance hash against
the matched-stage setup, and exported the matched named-column CSV. The agreed
derivative verification is complete at these tested states; target-amplitude
convergence and spatial accuracy remain separate checks.

All use `scripts/submit_ad_verification.sh` from the frozen snapshot, with
arguments `SNAPSHOT STAGE OUTPUT_PREFIX [SEED_JLS]`. Stages are `initial`,
`seed`, and `matched`; the last requires the saved seed path. Comparisons have
one CPU, 8 GB, and a two-hour limit; the seed used a one-hour limit and 20
maximum corrections at unchanged tolerance `1e-5`. The matched job used
`--dependency=afterok:27334006 --kill-on-invalid-dep=yes`. Logs are
`logs/ad-verify-JOBID.out`. No amplitude continuation was submitted by that
verification work; the later launch follows.

### Amplitude continuation launch (October 27)

The user built `scripts/continue_amplitude.jl` through commented, tested paste
blocks, then explicitly authorized creating and running the Slurm wrapper.
The driver preserves A1, the spatial settings, and the fixed coefficient scale.
It starts from the matched weak seed and targets A2=`1.0896058341192524e-4`.
The initial fractional increase is 10%. Successful corrections taking at most
four iterations grow the next increment by 50%; twelve or more halve it.
Fractional increments are bounded by `1e-4` and `0.25`; the last step is clipped
to the exact target. Recognized numerical Newton failures halve the actual
attempted increment and retry from the last matched point. Falling below the
minimum increment stops the job. Fresh-integration assertion failures also stop
the job. The norm tolerance remains `1e-5`, with 20 corrections per attempt and
full-width AD.

Frozen snapshot `work/amplitude_continuation_20261027_v1` verifies **201 files**.
All core/multiphase Julia library files match `ad_verification_20261027_v1`
byte-for-byte. Its `inputs/initial.jls` is the matched weak seed, SHA-256
`42c0d13e388bd7f1efe0f771741603aad2a478f7966dac82f3fae8aa7511fbb8`.
The shell syntax check passed. Earlier controlled tests covered loop behavior,
retries, target stopping, weight propagation, and actual file saves; these used
simulated solves and do not establish physical continuation convergence.

Submitted as job **27337026**, RUNNING on **n0060.lr7** from
**17:54:21 PDT**, one CPU, 8 GB, four-hour limit (21:54:21 PDT):

```bash
sbatch runs/adiabatic/multiphase_two_cloud/work/amplitude_continuation_20261027_v1/runs/adiabatic/multiphase_two_cloud/scripts/submit_amplitude_continuation.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/amplitude_continuation_20261027_v1" \
    "$PWD/runs/adiabatic/multiphase_two_cloud"
```

Log: `logs/amplitude-27337026.out`. Each attempt writes
`results/amplitude_step_20261027_v1-attemptN-accepted.jls` for its latest Newton
iterate. A successfully reintegrated point is saved separately as
`results/amplitude_step_20261027_v1-attemptN.jls`, including the profile,
amplitudes, residual weights, and both radial branches. These matched files
share the seed schema. Restarting requires choosing the matched input and a
fresh output prefix in a reviewed driver edit; the driver does not automatically
resume existing output files or recover its previous adaptive step size.

At the launch check, attempt 1 uses A2=`1.1985664175311778e-7`. Its weighted
residual falls from `2354.7637801303267` to `0.009104341308015403` after correction
1. An accepted-state checkpoint exists; no new matched point is saved yet.
The target-amplitude comparison plot remains pending.

### Continuation timeout and attempt 25 recovery (October 28)

Slurm accounting confirms job **27337026** ended with **TIMEOUT** after
4h00m22s on October 27 at **21:54:43 PDT**. It saved 24 matched points, with
no numerical-failure retries. Attempt 24 reached A2=`4.4781006633078346e-6`
and freshly integrated weighted norm `6.742071675579288e-6`. Attempt 25's
checkpoint reached correction 8 and norm `2.153421086390111e-6` just before
the time limit, leaving only its accepted-state file.

At the user's request, recovery job **27377949** ran on **n0003.lr7** on
October 28 from **11:06:56 to 11:07:29 PDT** and **COMPLETED** (exit 0).
It used the original `amplitude_continuation_20261027_v1` source snapshot,
verified all **201 files**, and freshly integrated the saved state without
changing its amplitudes or taking any Newton corrections.

- A1=`0.036390480353373555`, A2=`4.54807098617202e-6`, **4.17405%** of
  the target A2=`1.0896058341192524e-4`.
- Fresh weighted norm: **`2.153421086390111e-6 < 1e-5`**. The fresh residual
  vector is exactly equal to the checkpoint's stored vector.
- Both branches are finite: horizon states `(28, 1051)` on r*=-100 to 5;
  outer states `(28, 2951)` on r*=300 to 5. Matching includes values and slopes.
- Saved complete seed: `results/amplitude_step_20261027_v1-attempt25.jls`.
  It uses the existing matched-point schema and retains correction count 8.
  Serialization readback checked the state, amplitudes, weights, residual norm,
  and both radial state arrays before moving the temporary file into place.
- The original `-accepted.jls` is unchanged, SHA-256
  `2f277dffc9b68e79d81fb6b47d992bc9da8c9471bcec4ce2b0301d101a8396ac`.
  Recovered output SHA-256:
  `b94963278abc201fe374650e87277e10f60678de0e1b6b9b728ab9286694e18a`.

Recovery code and submission wrapper are
`work/recover_amplitude_attempt25_20261028_v1.jl` and the adjacent `.sh`.
Log: `logs/recover-attempt25-27377949.out`; it records the checkpoint, script,
source-manifest, and output hashes. The command was:

```bash
sbatch runs/adiabatic/multiphase_two_cloud/work/recover_amplitude_attempt25_20261028_v1.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/amplitude_continuation_20261027_v1" \
    "$PWD/runs/adiabatic/multiphase_two_cloud"
```

This recovery completes attempt 25 at the existing spatial settings; it does
not establish spatial convergence or a target-amplitude solution. No further
continuation was submitted. Restart from the recovered seed with a fresh output
prefix and the existing **1.5625%** increment (8 corrections leave that increment
unchanged). The current driver would reset to 10% and reuse the old output
names, so make the small restart edits through the user's reviewed snippets.
The goal remains both main profiles at the exact plot amplitudes.

### Continuation restart from recovered attempt 25 (October 28)

After recovery the user explicitly requested further continuation. The only
behavioral changes to `scripts/continue_amplitude.jl` are the initial fractional
step `0.015625` and fresh prefix `results/amplitude_step_20261028_v1-attemptN`.
The target A2, fixed A1, residual weights, spatial settings, `1e-5` tolerance,
20-correction limit, full-width AD, and adaptive rules remain unchanged.
Attempt numbering starts at 1 for this new prefix.

Frozen snapshot `work/amplitude_continuation_20261028_v1` verifies **201 files**.
Its core/multiphase library matches `amplitude_continuation_20261027_v1`
byte-for-byte. `inputs/initial.jls` is the recovered attempt-25 matched point,
SHA-256 `b94963278abc201fe374650e87277e10f60678de0e1b6b9b728ab9286694e18a`.
Shell syntax and a setup-only execution of the frozen Julia driver passed,
including the seed's amplitudes, iteration count, matching norm, and finite
radial branches. The driver comparison confirmed only the two restart edits
and addition of a terminal newline; the new output prefix was unused.

Submitted at **11:20:16 PDT** as job **27378594**, subsequently **RUNNING** on
**n0064.lr7**, with one CPU, 8 GB, and a four-hour limit:

```bash
sbatch runs/adiabatic/multiphase_two_cloud/work/amplitude_continuation_20261028_v1/runs/adiabatic/multiphase_two_cloud/scripts/submit_amplitude_continuation.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/amplitude_continuation_20261028_v1" \
    "$PWD/runs/adiabatic/multiphase_two_cloud"
```

Log: `logs/amplitude-27378594.out`. The first new point targets
A2=`4.6191345953309576e-6` from recovered A2=`4.54807098617202e-6`.
Its initial weighted residual is `606644.1247245938`; this launch check does
not establish a newly matched point. Each accepted Newton state is checkpointed,
and each successfully reintegrated matched profile is saved under the new prefix.
The end goal remains A2=`1.0896058341192524e-4` for both main profiles in the plot.

**Performance follow-up discussed, not yet run:** the user is concerned about
the cost implied for time evolution. Each midpoint evolution step currently
requires two corrected samples and their charge-response calculations. Evolution
already predicts profiles with an implicit amplitude derivative; continuation
currently scales boundary coefficients and carries over frequencies. Compare
predictors at the same amplitude increment and tolerance, then benchmark a
complete evolution sample/step. Repeatedly shrinking continuation increments
alone does not establish that a long finite-amplitude trajectory is practical.

### Fixed-target matching-location investigation (October 28)

The user authorized outlining and beginning an investigation of direct solves
at the plot amplitudes. The [plan](../../../docs/multiphase-target-solve-investigation-2026-10-28.md)
specifies the measurements and decisions. Job **27380795** compares joins
r*=5, 10, 20, 40, 60 at fixed A1=`0.036390480353373555`,
A2=`1.0896058341192524e-4`, with the original reference-based guess from
`results/direct_fig7_20261026_v5-input.jls`. All other numerical settings stay
fixed. It uses verified source snapshot `amplitude_continuation_20261028_v1`.

The standalone `work/target_match_scan_20261028_v1.jl` saves both radial
branches, field/slope norms and relative jumps, dominant mismatch modes,
full-width local-AD Jacobians, native and explicitly equilibrated condition
estimates, and proposed Newton corrections. It applies no Newton corrections.
Input, script, and snapshot-manifest hashes are recorded in its setup artifact.
The baseline must reproduce the original r*=5 norm within `1e-8` relative.
Syntax, target-input assertions, and reconstruction of all five problem settings
passed before submission. No production driver or library was changed.

```bash
sbatch runs/adiabatic/multiphase_two_cloud/work/target_match_scan_20261028_v1.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/amplitude_continuation_20261028_v1" \
    "$PWD/runs/adiabatic/multiphase_two_cloud"
```

Resources: one CPU, 8 GB, 30-minute limit. Log:
`logs/target-match-scan-27380795.out`. Output prefix:
`results/target_match_scan_20261028_v1`, with `-setup.jls`, `-summary.csv`,
and separate `-matchN-profile.jls` / `-matchN-jacobian.jls` files. It was
PENDING at submission; continuation **27378594** remains separate.

**Scan outcome:** COMPLETED on n0072.lr7 in **5m19s**, from **12:01:41 to
12:07:00 PDT**. All five integrations and Jacobians completed, with the r*=5
baseline reproduced exactly. Scaled residuals for r*=5, 10, 20, 40, 60 were
`1.90648e9`, `9.47348e6`, `51008.4`, `413.240`, and `56.5614`.
At 60 the horizon/outer field norms are 0.06573/0.06340 and relative
field/slope jumps are 0.09335/0.29657. Full Newton steps satisfy admissibility
at 40 and 60. Native/equilibrated Jacobian condition estimates at 60 are
`1.95e16`/`2.99e3`; severe native scaling remains a concern. At 20 the scaled
condition estimate worsens despite its smaller residual. See the investigation
plan and saved CSV for all values. No initial guess passes the `1e-5` tolerance.

**Bounded direct test:** selected r*=60 and submitted job **27381429**, RUNNING
on n0064.lr7 from **12:13:33 PDT**, with one CPU, 8 GB, and a 30-minute limit.
It uses the same full target amplitudes, original guess, numerical library,
residual normalization, full-width local AD, and tolerance, with **six corrections**
maximum. Each accepted state is saved. The standalone Newton diagnostic's
setup checks and shell syntax passed before submission.

```bash
sbatch runs/adiabatic/multiphase_two_cloud/work/target_match_newton_20261028_v1.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/amplitude_continuation_20261028_v1" \
    "$PWD/runs/adiabatic/multiphase_two_cloud" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/target_match_scan_20261028_v1-match60-profile.jls" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/target_match60_20261028_v1"
```

Log: `logs/target-match-newton-27381429.out`. Output prefix:
`results/target_match60_20261028_v1`, including `-setup.jls`, `-accepted.jls`,
and either a freshly reintegrated matched `.jls` or a `-failure.txt` report.
The changed join is promising for this guess but neither the scan nor launch
establishes a target solution or spatial accuracy.

**Direct-test outcome:** job **27381429 FAILED** at its six-correction limit
after **13m32s**, ending **12:27:05 PDT** on October 28. Accepted residuals:
`56.5614343, 54.1790533, 54.1118685, 54.0476949, 53.9864134, 53.9279074,
53.8720631`, versus the `1e-5` target. The decrease is only **4.75%**.
The first correction needed six trial residual evaluations; all five later
corrections needed nine. Jacobian passes took 56.4–76.6 s, with trial integrations
around 8 s each. The final `-accepted.jls` and `-failure.txt` exist; no matched
`.jls` was saved. Moving the join helped the starting residual but did not give
useful convergence in this attempt. See the investigation plan for the proposed
linear-solve scaling and component diagnostics before a formulation change.

At this review, continuation **27378594** remained RUNNING after roughly 1h07m,
with seven newly saved matched points. The latest,
`results/amplitude_step_20261028_v1-attempt7.jls`, has
A2=`5.069450807458185e-6` (4.65255% of target) and weighted norm
`5.828310199789505e-6`; attempt 8 targets `5.148660976324719e-6`.

**Saved-matrix cause check:** `work/target_match60_linear_check_20261028_v1.jl`
compares the saved first correction with an explicitly equilibrated solve and
a 256-bit solve of the same saved Float64 matrix. Both agree with the original
direction to `9.46e-15` relative norm; maximum coordinate-scaled differences are
below `7e-12`. Log: `logs/target-match60-linear-check-20261028-v1.out`.
No radial integrations were run. The first full Newton trial raises the norm
from 56.56 to 6118.82 despite accurately solving the linearized equations; its
first accepted damping is 1/32. Largest coefficient changes occur in initially
zero mixed driven components. This motivates testing driven-field initial guesses
and the actual target-state directional derivative; it does not establish the
complete physical/numerical cause. See the investigation record for details.

### Driven-field initial guesses (October 28)

The user requested adding the existing forced m=0 and m=3 solutions to the
joint initial guess. Standalone `work/driven_field_seed_20261028_v1.jl` uses
the frozen `amplitude_continuation_20261028_v1` library and point-72 reference
inputs from `work/direct_fig7_20261026_v5/inputs`. It scales m=0 coefficients
by A2 and m=3 by A2 squared, retains the original fundamental coefficients
and frequencies, and evaluates joint matching at the full plot amplitudes and
r*=60. Mode sets are l=0,2,4 for m=0 and l=3,5 for m=3.

```bash
sbatch runs/adiabatic/multiphase_two_cloud/work/driven_field_seed_20261028_v1.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/amplitude_continuation_20261028_v1" \
    "$PWD/runs/adiabatic/multiphase_two_cloud"
```

Job **27384564 COMPLETED**, 12:54:38–12:56:11 PDT, **1m33s** on n0003.lr7.
Setup, mode-order, syntax, and all 201 frozen-source checks passed. The forced
matching residual/source ratios are below `7e-15`; m=0 and far-leg m=3 power
ratios agree with one within `8e-8`. All saved joint trial branches are finite.
Log: `logs/driven-field-seed-27384564.out`. Output prefix:
`results/driven_field_seed_20261028_v1`, including `-setup.jls`, `-forced.jls`,
`-summary.csv`, and `-zero.jls`, `-m0_only.jls`, `-m3_only_near.jls`,
`-both_near.jls`, `-both_far.jls`.

The standard m=3 reference propagates outgoing data from r*=5000 at order 16;
the joint solver uses an order-7 series directly at 300. A second forced solve
with `far_position=300, far_order=7` matches the current joint boundary
prescription. It is a controlled seed variant, not a validated outgoing boundary.
The near/far outer states differ by 28.3% for identical leading coefficients;
this does not compare their subspaces after freely refitting coefficients.

Scaled joint residuals: zero fields **56.5614**; m=0 only **3529.9916**;
m=3 only **57.0283**; both near-series fields **3529.2624**; both with far-leg
m=3 coefficients **3529.0943**. The m=0 addition drives the increase, predominantly
in the primary mismatch. See the investigation record for block norms and
boundary checks. Adding the fields did not improve the initial total residual.

A controlled six-correction Newton comparison was submitted as **27384946**:

```bash
sbatch runs/adiabatic/multiphase_two_cloud/work/target_match_newton_20261028_v1.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/amplitude_continuation_20261028_v1" \
    "$PWD/runs/adiabatic/multiphase_two_cloud" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/driven_field_seed_20261028_v1-both_near.jls" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/target_match60_driven_20261028_v1"
```

It uses the same library, full target amplitudes, join 60, boundary prescription,
and `1e-5` tolerance as the earlier zero-driven test. Log:
`logs/target-match-newton-27384946.out`. Every accepted state is checkpointed;
the script saves a matched profile only after successful fresh reintegration.
No production library or reference driver was changed.

**Driven-seed Newton outcome:** **27384946 FAILED** at its six-correction cap
after **8m36s** on n0013.lr7, ending **13:07:02 PDT**. Accepted residuals:
`3529.26235, 2267.89812, 2147.61360, 2072.02371, 2022.38703, 1896.03094,
1692.49867`, about a 52% decrease versus the `1e-5` target. Trial integration
counts were 4,3,4,3,3,1. The last correction decreased the residual by 10.7%
and accepted its first evaluated trial; inadmissible trials are not logged.
`results/target_match60_driven_20261028_v1-accepted.jls` and `-failure.txt` are
saved, with no matched `.jls`. No target extension was submitted at this review.
The checkpoint allows a further bounded comparison without restarting from
the original driven guess.

At the same review, continuation **27378594** remained RUNNING after about
2h01m and had saved 13 new matched points. Latest
`results/amplitude_step_20261028_v1-attempt13.jls` has
A2=`5.394443404258194e-6` (4.95082% of target), weighted norm
`6.00678166743477e-6`. The current amplitude increment is 0.78125%; attempt 14
targets `5.4365874933539615e-6`.

### Driven-field Newton extension: twelve additional corrections (October 28)

At the user's explicit request, job **27388374** resumes the saved sixth
correction from `results/target_match60_driven_20261028_v1-accepted.jls` and
allows **twelve additional corrections**, up to **18 total**. The stored norm
is `1692.4986686897735`; checkpoint SHA-256:
`f680a44bb983ede293e36575f27300ffe4472d4a063ba2061611f0278fa4476f`.
The target amplitudes, join 60, boundary prescription, residual weights,
full-width local AD, and `1e-5` tolerance are unchanged. It stops earlier if it
converges or a numerical failure occurs. The saved total iteration count and
accepted-state logs include the original six corrections.

```bash
sbatch runs/adiabatic/multiphase_two_cloud/work/resume_target_match_newton_20261028_v1.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/amplitude_continuation_20261028_v1" \
    "$PWD/runs/adiabatic/multiphase_two_cloud" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/target_match60_driven_20261028_v1-accepted.jls" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/target_match60_driven_20261028_v2"
```

The standalone restart script and shell passed syntax/setup checks, including
exact state identity, stored norm, unchanged settings, and total correction
limit. All 201 frozen-source files verified. The first fresh residual must
agree with the saved mismatch to `1e-10` relative. Resources: one CPU, 8 GB,
one-hour limit. Log: `logs/target-match-resume-27388374.out`. New output prefix:
`results/target_match60_driven_20261028_v2`, with setup, accepted-state, failure,
and possible matched-profile artifacts. The original checkpoint is preserved.
The job is RUNNING on n0060.lr7 from **13:27:03 PDT**. Its first fresh
integration exactly reproduced `1692.4986686897735`, and cumulative correction
6 was saved under the new prefix before beginning correction 7.

### Reference reconstruction snippets tested (October 27)

The user requested testing the secondary block and then the two linear
references, and explicitly asked to make the driver edits themselves. The live
`scripts/reference_profiles.jl` was left unchanged at primary reconstruction.

The first test, `work/check_reference_secondary_20261027_v1.jl`, ran as
**27337260** and FAILED after 23 s. Its normalized matching norm was
`2.8370280671313985e-6 < 1e-5`, and its unit-amplitude reconstruction exactly
reproduced the saved solution residual. However, integrating with A2-scaled
boundary coefficients gave an inward-branch scaling error `8.241260253527706e-10`,
above the proposed `1e-10` check; the horizon-branch error was `1.264e-15`.
The original script and failure log are preserved.

The revised candidate, `work/check_reference_secondary_20261027_v2.jl`,
reconstructs the unit-amplitude secondary and multiplies its saved states and
residual by A2. This uses linearity on the fixed Omiya primary background; the
nonlinear primary is not rescaled. Job **27337856** COMPLETED, exit `0:0`, in
18 s. The scaled mismatch is `3.157814155718417e-6`, and direct subtraction of
the scaled endpoint states gives `3.157814155718426e-6`, both below the unchanged
`1e-5` tolerance. Sample-scaling errors are at most `1.835e-16`, below the
unchanged `1e-10` criterion. Separate evaluations of the averaged-secondary RHS
at seven samples on each radial branch commute with amplitude scaling to at
most `1.567e-16`, against `1e-12`. Both branches are finite.

Job **27337858** COMPLETED, exit `0:0`, in 21 s, running
`work/check_reference_linears_20261027_v1.jl`. It solves the actual linear
eigenmodes with the existing solver, retaining `(l,m)=(1,1),(3,1),(5,1)` for the
primary and `(2,2),(4,2),(6,2)` for the secondary. The fundamental m is passed
explicitly to radial reconstruction, so the secondary uses omega2, not 2omega2.
Both unit-outer-coefficient solves meet `1e-8`; unit reconstruction reproduces
each saved solve residual exactly. Integration with A1/A2-scaled boundary data
also passes the `1e-5` matching and `1e-10` profile-scaling checks:

| Reference | Frequency | Unit solve norm | Scaled reconstruction norm |
|---|---|---|---|
| Linear primary | `0.29630433912074067 + 2.680675437876823e-8im` | `4.712206398423109e-13` | `3.7386254061443036e-13` |
| Linear secondary | `0.2984522683211546 + 1.7965646025702976e-12im` | `6.523329776273922e-9` | `6.78907457596775e-9` |

The maximum linear-profile scaling error is `8.81233e-13`. All reference curves
use spin 0.99, alpha 0.3, fixed A1=`0.036390480353373555`, fixed
A2=`1.0896058341192524e-4`, endpoints -100/300, match 5, and boundary orders 3/7.
The primary profiles use step 0.1; the secondary profiles use step 0.2. This
matches the Omiya reconstruction grids but does not establish spatial accuracy.

Saved results under this campaign:

- `results/reference_secondary_20261027_v2.jls`: both Omiya profiles, the unit
  secondary, settings, matching norms, input/check hashes, and job ID.
- `results/reference_linears_20261027_v1.jls`: both linear solutions and radial
  profiles, mode lists, settings, input/check hashes, and job ID.
- `logs/reference-secondary-check-27337260.out`,
  `logs/reference-secondary-v2-27337856.out`, and
  `logs/reference-linears-check-27337858.out`: complete diagnostic output.
- `work/reference_snippet_checks_20261027_v1-provenance.json`: 68 source,
  dependency, input, and check-file hashes, verified unchanged after the jobs.

All three used one CPU and 8 GB with a 15-minute allocation. The two successful
jobs started at 18:15:44 PDT. Next, the user adds the tested secondary block and
then the linear blocks to the reference driver. The four sampled reference
profiles exist; combined CSV/plot export and spatial refinement remain pending.

### Four reference curves exported through the shared branch interface (October 28)

The user authorized completing the reference workflow directly. The driver now
reads both primary and secondary checkpoints with `read_branch`, selects their
`ContinuationPoint`s with `branch_point`, and reconstructs from the shared
`TwoPhaseProfile` fields. The secondary adapter preserves its primary-point
reference, actual per-point tolerances, and linear seed. The linear matcher now
accepts `save_profile=true`, returning the same radial-result layout as the
nonlinear reconstructions. No new reconstruction framework is required.

The first export (`v1`, job **27397837**, 32 seconds) reproduced the October 27
radial arrays exactly, but held the linear outer coefficients fixed. The user
correctly identified that this was the wrong normalization for Fig. 7, whose
caption specifies equal cloud energies. Those outputs are preserved as history.

The corrected export is **`reference_curves_20261028_v2`**. Job **27398471
COMPLETED**, exit `0:0`, in **42 seconds**, using the 84-file frozen source copy
`work/reference_curves_20261028_v2`. Its manifest verifies successfully.
Before normalization, both radial state arrays for every curve agree exactly
with the October 27 archives. Normalized matching residuals remain:

| Curve | Matching norm |
|---|---:|
| Omiya primary | `5.587385657966086e-11` |
| Omiya secondary | `3.157814155718417e-6` |
| Linear primary | `3.7386254061443036e-13` |
| Linear secondary | `6.78907457596775e-9` |

Outputs:

- [PNG](figures/reference_curves_20261028_v2.png) and
  [PDF](figures/reference_curves_20261028_v2.pdf): four fundamental radial curves,
  shown over `-100 <= r*/M <= 200`.
- [CSV](results/reference_curves_20261028_v2.csv): 12,008 samples over the full
  domain `[-100, 300]`, with complex coefficients, magnitudes, mode labels, and
  integration-branch labels. Both samples at the matching radius are retained.
- `results/reference_curves_20261028_v2.jls`: all radial modes and slopes,
  linear solutions, and metadata.
- [TOML provenance](results/reference_curves_20261028_v2.toml): energies before
  and after scaling, amplitude factors, actual outer amplitudes for each curve,
  settings, point IDs, input/source hashes, residuals, and archive comparisons.
- `logs/reference-profiles-27398471.out`: reconstruction and normalization results.

All four curves now have energy **650.5419114886128 F_a^2 M**, using saved
primary/secondary point 72 at spin 0.99 and alpha 0.3. This is the saved nonlinear
point nearest the paper's approximate 660; the nonlinear samples are unchanged.
The Omiya primary uses its cosine-potential energy, and the secondary uses its
quadratic energy on the averaged primary background. Each linear energy uses
the free massive-field potential U(phi)=phi^2/2. Since E(cR)=c^2 E(R), the driver
multiplies every linear radial coefficient, slope, and matching residual by
sqrt(reference_energy / measured_linear_energy), then measures the energy again.

| Linear cloud | Energy before scaling | Amplitude factor | Corrected peak | Peak r*/M |
|---|---:|---:|---:|---:|
| l=m=1 | 1235.9435998974182 | 0.7255014992542931 | 0.11519751877506104 | 24.5 |
| l=m=2 | 1143.7368349843898 | 0.7541792626953592 | 0.030068575045116026 | 71.8 |

All four reintegrated energies agree with the target to roundoff. The corrected
linear outer amplitudes are `0.026401348054956413` and `8.217581246046196e-5`;
the Omiya outer amplitudes remain `0.036390480353373555` and
`1.0896058341192524e-4`. The resulting plot agrees visually with the supplied
Fig. 7. No digitized pointwise comparison or spatial convergence has been
established. Coupled multiphase target curves remain separate work.

The targeted branch, linear-sample, loader, and existing continuation checks
passed (98 assertions). Checks cover reordered/missing/duplicate IDs, preservation
of secondary tolerances and seeds, serialization, and identical linear residuals
with and without saved radial samples. The live adapter also handles an empty
secondary history; that small addition after the frozen run was tested separately.
CSV ordering, complex magnitudes, and both matching-point samples were checked,
and the resulting figure was inspected.

The normalization correction passed **63 focused assertions** covering energy,
potentials, branch access, and linear samples. These include analytic Kerr
free-field energies for m=1 and m=2, amplitude-squared energy scaling, and cubic
quadrature with even or odd interval counts in either direction. The existing
shared Simpson quadrature now also handles the secondary's odd interval counts.
Every v2 CSV sample was checked against v1: Omiya values are identical, linear
values have exactly the recorded scaling within roundoff, and complex magnitudes,
ordering, finite values, and both matching-point samples are correct.

Reproduce the numerical run from the repo root, choosing a fresh output prefix:

```bash
sbatch runs/adiabatic/multiphase_two_cloud/scripts/submit_reference_profiles.sh \
    "$PWD/runs/adiabatic/multiphase_two_cloud/work/reference_curves_20261028_v2" \
    "$PWD/runs/adiabatic/multiphase_two_cloud" \
    "$PWD/runs/adiabatic/multiphase_two_cloud/results/reference_curves_recheck"
```

Render the exported curves with the BHsuper Python environment:

```bash
MPLCONFIGDIR=/tmp/axion_reference_mpl PYTHONPATH=src \
    /global/home/users/clwelch/.conda/envs/BHsuper/bin/python \
    runs/adiabatic/multiphase_two_cloud/scripts/plot_reference_profiles.py \
    runs/adiabatic/multiphase_two_cloud/results/reference_curves_20261028_v2 \
    runs/adiabatic/multiphase_two_cloud/figures/reference_curves_20261028_v2
```
