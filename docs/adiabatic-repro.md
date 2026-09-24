# Reproducing adiabatic axion-cloud evolution

Reference: Omiya, Takahashi & Tanaka,
[2201.04382v3 (PDF)](https://arxiv.org/pdf/2201.04382v3).
See [adiabatic-design.md](adiabatic-design.md) for equations, conventions,
module responsibilities, nonlinear matching, and observables.

## Status — 2026-09-24

The Julia implementation includes linear and nonlinear matching, cosine,
quartic, and sextic source/energy potentials, continuation predictors, retained
radial profiles, and energy integration. Campaign drivers calculate signed
finite-radius horizon and outer fluxes. Newton corrections support centered
finite differences and ForwardDiff Jacobians on the real-packed residual.

All three preliminary Fig. 2 branches have saved data at spin=0.99, alpha=0.42.
The [campaign README](../runs/adiabatic/fig2_20260922/README.md) records table
sizes, current defaults, restart recipes, and plotting commands. Its
[history](../runs/adiabatic/fig2_20260922/HISTORY.md) records the numerical
experiments. The earlier seven-point
[small-amplitude branch](../runs/adiabatic/small_amplitude/README.md) remains a
separate reference.

These are stationary families, not yet a slow-time evolution. Full
high-amplitude radial, boundary, quadrature, and harmonic-truncation convergence
has not been established. Predictor and Jacobian comparisons establish local
solve behavior, not discretization accuracy.

## Figure 2 observables and data

`integrate_radial(...; save_profile=true)` returns `positions` and `states`,
with one `[R; Q]` column per position and `Q=dR/dr_*`. After a solve,
`nonlinear_matching_residual(...; save_profile=true)` reconstructs both branches
and their mismatch. The horizon branch runs outward and the outer branch runs
inward; both matching states are retained without averaging. Ordinary Newton
residual evaluations do not retain profiles.

`integrate_energy_branch` reverses descending profiles and uses composite
Simpson quadrature, requiring uniform spacing and an even interval count.
Sum the integrals of both branches for the cloud energy. Pass the same
`potential` selection to source evaluation, profile reconstruction, and energy.

The campaign flux convention is `total_flux = horizon_flux + outer_flux`.
The horizon contribution is minus the outward current at the inner endpoint;
a negative total flux corresponds to cloud energy gain. The stored
`balance_ratio = -total_flux / (2imag(omega) * energy)` is a small-amplitude
check, not an identity to impose on nonlinear solutions.

CSV columns include amplitude, energy, `normalized_energy=energy/A0^2`, real
and imaginary frequency, horizon/outer/total flux, balance ratio, sampled peak
field, scaled matching residual, and Newton iterations. Plot `energy`, not
`normalized_energy`, for Fig. 2. Connect flux-versus-energy points in amplitude
order to preserve the quartic energy turning point. The plotting recipe merges
production, refined, and gap tables, retaining production values at shared
amplitudes.

```bash
python runs/adiabatic/fig2_20260922/scripts/plot_branch.py
```

Outputs: `plots/preliminary_fig2_20260922.png` and `.pdf`.

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

The suite covers angular/radial operators, boundaries, integration, matching,
nonlinear sources, predictors, energy, and polynomial potentials. Run from the repo root:

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
The generated PNG is `plots/fig13_growth_rates.png`.
Generated frequencies and diagnostics are written to
`plots/fig13_adiabatic_points.csv`. Add `--skip-eigenvalue` to omit the ARPACK
comparison when redrawing.

## Remaining checkpoints

1. Check high-amplitude radial spacing, finite boundaries, angular quadrature,
   and retained-mode dependence. Extend or refine branches where needed for
   the intended Figure 2 comparison.
2. Implement slow amplitude evolution using `E′(A0) dA0/dt = -F_tot`, checking
   the small-amplitude limit and behavior near energy turning points.
3. Compare converged curves and evolution with the target figure before
   claiming a completed reproduction; then vary physical parameters.

## Code boundaries

Adiabatic profile solving, continuation, and energy/flux evolution belong here.
Use `common.leaver` as an independent benchmark. Keep `dolan.qbs_eigenmode`
with its time-domain operator. Share nonlinear angular routines only after
checking coordinate, radial-field, and harmonic normalization conversions;
the current Dolan source is not assumed to be directly interchangeable.

Run recipes and outputs are grouped under `runs/adiabatic/`; see the
[campaign index](../runs/README.md) for the layout and preservation policy.
