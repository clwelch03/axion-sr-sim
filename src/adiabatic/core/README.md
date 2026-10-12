# Shared adiabatic building blocks

Load `src/adiabatic/core/load.jl` for mathematics independent of a cloud model.

| Files | Responsibility |
|---|---|
| `geometry.jl`, `angular.jl` | Kerr/tortoise geometry and angular bases/quadrature |
| `potentials.jl` | Cosine, quartic, and sextic potentials and stable nonlinear remainders |
| `radial_operator.jl` | Radial coefficients for an explicitly supplied frequency and angular indices |
| `boundary_series.jl`, `boundaries.jl` | Horizon and bound/outgoing series and per-component endpoint values |
| `integration.jl`, `newton.jl` | Radial RK4, damped Newton, finite/automatic Jacobians, complex packing |
| `quadrature.jl` | Uniform Simpson quadrature and cumulative trapezoidal integration |
| `sensitivity.jl` | Implicit derivatives of a residual-defined branch and its observables |
| `units.jl` | Physical normalization and reference BH units |

Frequency assignment from `m / fundamental_m`, nonlinear field reconstruction,
and cloud-energy balance belong to the models, not this folder. The present
boundary helpers retain their original supported domains, including the
positive-frequency outer branch and exclusion of the exact mass threshold.
The refactor does not generalize those domains for future multiphase components.

Check this layer independently:

```bash
julia --project=. tests/adiabatic/core/runtests.jl
```
