# Joint two-cloud response and evolution

## Start here

This campaign tests implicit amplitude responses and a short E/J-balanced
trajectory of the joint four-component cosine model on a fixed Kerr background.
The October 15 pilot completed two trajectories and the derivative checks.
It is a finite-domain numerical pilot. It does not establish spatial or
Fourier-support convergence, boundary-tail accuracy, or the physical validity
of the frozen-profile slow-growth approximation.

The library is in `src/adiabatic/multiphase/response.jl` and `evolution.jl`.
`core/sensitivity.jl` supplies the general multiple-amplitude derivative.
The user authorized direct implementation on October 15, with explanatory
comments and docstring skeletons left for them to fill in.

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
