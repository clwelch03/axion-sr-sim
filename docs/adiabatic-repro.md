# Reproducing adiabatic axion-cloud evolution

Reference: Omiya, Takahashi & Tanaka,
[2201.04382v3](https://arxiv.org/html/2201.04382v3).

## Status

2026-09-17: package layout established. `src/adiabatic/` is a placeholder;
no adiabatic profile solver or reproduced nonlinear result exists yet.
The existing time-domain solver and its nonlinear extensions live in
`src/dolan/`. Independent benchmarks live in `src/common/`.

## First checkpoint: the linear radial problem

Start from the paper's separated linear equations (10)–(13), in
Boyer–Lindquist coordinates. Use units M=1, spin a/M, alpha=M mu,
and the time convention exp(-i omega t), so positive Im(omega) means growth.
Solve for a bound mode with an ingoing horizon boundary and a decaying
outer boundary, retaining the spheroidal angular eigenvalue.

The initial target is l=m=1, no radial overtone, spin=0.99, alpha=0.42.
Table 1 gives M omega_R=0.4088 and M omega_I=1.504e-7.
Compare both components against `common.leaver` at runtime; the rounded
table values are a sanity check. Vary radial endpoints and numerical accuracy
to distinguish boundary error from integration error, especially in omega_I.

Before coding, derive the radial residual and boundary expansions together.
The first implementation should expose those equations and this one numerical
checkpoint. Keep amplitude and phase conventions explicit so the later
nonlinear branch can grow from this solution.

## Following checkpoints

1. Translate the linear solution into the spherical-harmonic representation
   used by the nonlinear radial system; verify harmonic and field normalization.
2. Implement the nonlinear profile problem and amplitude continuation using
   Section 3 and Appendix A, retaining the published complex-frequency
   prescription and boundary conditions.
3. Compute cloud energy and horizon/outer fluxes, then evolve the amplitude
   using E'(A0) dA0/dt = -F_tot (equation 32). Check the small-amplitude limit.
4. Reproduce the energy/flux curves and evolution at the initial parameter
   point, then study parameter and harmonic-truncation dependence.

## Code boundaries

Adiabatic profile solving, continuation, and energy/flux evolution belong here.
Use `common.leaver` as an independent benchmark. Keep `dolan.qbs_eigenmode`
with its time-domain operator. Share nonlinear angular routines only after
checking coordinate, radial-field, and harmonic normalization conversions;
the current Dolan source is not assumed to be directly interchangeable.

The existing `runs/` data stays in place. Future adiabatic run recipes and
outputs should use `runs/adiabatic/`.
