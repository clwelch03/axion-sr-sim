# Joint two-cloud response and evolution

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
