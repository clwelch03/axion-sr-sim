# Reproducing adiabatic axion-cloud evolution

Reference: Omiya, Takahashi & Tanaka,
[2201.04382v3 (PDF)](https://arxiv.org/pdf/2201.04382v3).
See [adiabatic-design.md](adiabatic-design.md) for equations, conventions,
module responsibilities, and the planned nonlinear extension.

## Status

2026-09-21: the Julia solver includes linear matching, cosine nonlinear source
projection, nonlinear matching with a real Jacobian, and a small-amplitude
continuation driver. The local saved continuation contains seven solutions
from A0=1e-4 to 0.0011390625. Radial profiles can now be retained for inspection.
Energy/flux integrals, the other two potentials in Fig. 2, continuation into
the strongly nonlinear regime, and slow amplitude evolution remain unfinished.

## Fig. 2: current ingredient — retaining radial profiles

The target is both panels at spin=0.99, alpha=0.42, with quartic, sextic,
and cosine potentials. We are working through the ingredients incrementally.
The saved small-amplitude branch is not yet a reproduction of these curves.

`integrate_radial(...; save_profile=true)` now returns `positions` and `states`.
Each column of `states` contains `[R; Q]`, where Q=dR/dr_*. The default still
returns only the final state, avoiding profile storage during Newton solves.
Calling `nonlinear_matching_residual(...; save_profile=true)` with an already
solved frequency and boundary amplitudes returns `from_horizon`, `from_outer`,
and `residual`. Both branches run from their boundary to the matching point;
the outer branch is therefore in decreasing coordinate order. Both matching
states are retained without averaging. Request this after solving, not as an
option to `solve_nonlinear_mode`, which requires a residual vector.

At the existing settings (step=0.05, inner=-100, outer=100, match=5), replaying
the first and last locally saved solutions reproduces their residuals exactly:

| A0 | norm of R mismatch / A0 | norm of Q mismatch / A0 |
|---|---:|---:|
| 0.0001 | 2.31e-12 | 1.80e-12 |
| 0.0011390625 | 1.17e-14 | 1.56e-14 |

Analytic integration tests check every retained sample in both directions.
This validates profile retention and matching continuity at these settings;
it does not establish accuracy of energy or flux integrals.

Next, derive and implement the energy and signed boundary fluxes with the
paper's angular normalization. The decisive small-amplitude check is
F_tot approximately -2 Im(omega) E, with E proportional to A0 squared.
Then introduce the quartic and sextic potentials consistently in both the
source and energy, extend continuation, and check radial/angular truncation
and endpoints before comparing full curves. Plot F_tot versus E in amplitude
order so the quartic branch's turning point is preserved.

## Completed checkpoint: the linear radial problem

The implementation starts with PDF (16), projects onto spherical harmonics,
and solves the linear limit of (21) using boundaries based on (22)–(23).
Units are M=1, spin a/M, alpha=M mu, with time convention exp(-i omega t).
It includes Kerr angular mixing through the projected cos-squared coupling.

For spin=0.99, alpha=0.42, and the fundamental radial mode with angular
l=[1,3,5], m=1, the checked result is

```
M*omega = 0.408805564968254 + 1.503864692665e-7im
```

The initial solve converged in four Newton iterations with residual about
1.3e-14. Halving the integration spacing from 0.05 to 0.025, extending endpoints,
and moving the matching point shifted the real frequency by less than 1e-13
and the imaginary frequency by less than 2e-16.
The single-l=1 version was independently checked against `common.leaver` with
its matching leading angular correction; that reference is not an exact
target for the full l=[1,3,5] block.

The last numerical test run passed 208 Julia assertions. Run from the repo root:

```bash
julia --project=. tests/adiabatic/runtests.jl
```

## Fig. 13 comparison

The adiabatic solver supplies both radial overtones at 15 values of alpha
from 0.30 through 0.44. All 30 points converged. The imaginary starting guess
is zero; Leaver real frequencies select the radial branches. At the two ends
of the sweep, halving the step and increasing the outer boundary from 15 to
20 decay lengths changed the growth rates by less than 4e-9 fractionally.
This is a linear validation plot, not a reproduced nonlinear evolution.

```bash
PYTHONPATH=src python -m dolan.fig13_growth_rates --include-adiabatic --grid-spacing-n 6
```

Pass `--julia /path/to/julia` if needed. The plot uses open diamonds for
adiabatic results and the M/6 time-domain runs; M/6 is now the plotting default.
[The PNG](../plots/fig13_growth_rates.png) is included as a repository artifact.
Generated frequencies and diagnostics are written to
`plots/fig13_adiabatic_points.csv`. Add `--skip-eigenvalue` to omit the ARPACK
comparison when redrawing.

## Following checkpoints

1. Compute cloud energy and horizon/outer fluxes, then evolve the amplitude
   using E'(A0) dA0/dt = -F_tot (PDF equation 24). Check the small-amplitude limit.
2. Add the quartic and sextic potentials and extend amplitude continuation.
3. Reproduce the energy/flux curves and evolution at the initial parameter
   point, then study parameter and harmonic-truncation dependence.

## Code boundaries

Adiabatic profile solving, continuation, and energy/flux evolution belong here.
Use `common.leaver` as an independent benchmark. Keep `dolan.qbs_eigenmode`
with its time-domain operator. Share nonlinear angular routines only after
checking coordinate, radial-field, and harmonic normalization conversions;
the current Dolan source is not assumed to be directly interchangeable.

The existing `runs/` data stays in place. Future adiabatic run recipes and
outputs should use `runs/adiabatic/`.
