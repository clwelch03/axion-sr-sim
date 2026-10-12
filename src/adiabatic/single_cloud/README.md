# One nonlinear cloud

This model has one fundamental phase. Its additional `(l,m)` components are
angular terms and harmonics of that phase, not independently populated clouds.
Load `src/adiabatic/single_cloud/load.jl`; it loads core automatically.

| Files | Responsibility |
|---|---|
| `profiles.jl` | Harmonic frequency/boundary assembly, matching residuals, linear and nonlinear mode solves |
| `projection.jl` | Real-field reconstruction, angular projection, and nonlinear radial source |
| `observables.jl` | Energy density, integrated energy, and boundary fluxes |
| `continuation.jl`, `continuation_predictor.jl` | Branch setup, checkpoint I/O, continuation, and guesses |
| `interpolation.jl`, `branch_fitting.jl` | Tabulated branch interpolation and local energy fits |
| `evolution.jl` | Slow amplitude evolution from energy balance |
| `branch_observables.jl` | Correct saved profiles and compute energy/flux amplitude derivatives |
| `linear_benchmarks.jl`, `fig13_points.py` | Batch linear eigenvalue checks and the Python Fig. 13 bridge |
| `drivers/` | Command-line entry points and single-cloud profile tools |

Library users include the loader and call functions directly; they do not need
to include a driver. Use `drivers/continue_branch.jl` for cosine, quartic, or
sextic continuation. `drivers/branch_energy_derivatives.jl` processes saved
branches. See the [design](../../../docs/adiabatic-design.md) and
[amplitude campaign](../../../runs/adiabatic/amplitude_evolution/README.md)
for equations and reproducible inputs.

```bash
julia --project=. tests/adiabatic/single_cloud/runtests.jl
```
