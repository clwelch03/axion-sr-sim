# Reproducing adiabatic axion-cloud evolution

Reference: Omiya, Takahashi & Tanaka,
[2201.04382v3 (PDF)](https://arxiv.org/pdf/2201.04382v3).
See [adiabatic-design.md](adiabatic-design.md) for equations, conventions,
module responsibilities, and the planned nonlinear extension.

## Status

2026-09-18: the Julia linear radial matching solver is implemented and checked.
It includes spherical harmonics, angular couplings, stable tortoise geometry,
boundary series, radial RK4 integration, and a normalized complex-frequency
solve. Nonlinear profiles, continuation, energy/flux integrals, and slow
amplitude evolution are not implemented yet.

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

1. Retain and inspect the radial profiles; check continuity at the matching point.
2. Reconstruct the real field and project the nonlinear remainder. Resolve the
   source convention noted in the design doc; test cubic scaling and the linear limit.
3. Add amplitude continuation and a real-and-imaginary Newton system, including
   the radiative harmonics.
4. Compute cloud energy and horizon/outer fluxes, then evolve the amplitude
   using E'(A0) dA0/dt = -F_tot (PDF equation 24). Check the small-amplitude limit.
5. Reproduce the energy/flux curves and evolution at the initial parameter
   point, then study parameter and harmonic-truncation dependence.

## Code boundaries

Adiabatic profile solving, continuation, and energy/flux evolution belong here.
Use `common.leaver` as an independent benchmark. Keep `dolan.qbs_eigenmode`
with its time-domain operator. Share nonlinear angular routines only after
checking coordinate, radial-field, and harmonic normalization conversions;
the current Dolan source is not assumed to be directly interchangeable.

The existing `runs/` data stays in place. Future adiabatic run recipes and
outputs should use `runs/adiabatic/`.
