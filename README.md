# AxionCloud.jl

Multiphase adiabatic evolution of a superradiant axion cloud around a Kerr
black hole, with self-interactions treated non-perturbatively. This is the
Julia implementation of the method described in `docs/multiphase.tex`.

## Layout

| File | Contents |
| --- | --- |
| `src/kerr.jl` | Kerr background (`M = 1`), horizon data, tortoise coordinate and its inverse |
| `src/potential.jl` | Free, cosine, quartic and chiral (pion-integrated) axion potentials via `h(Φ) = Φ − U'(Φ)` |
| `src/modes.jl` | Bound modes, Fourier components `q`, retained polar indices, channel bookkeeping |
| `src/angular.jl` | Normalised polar harmonics, Gauss–Legendre grid, `x²` mixing matrix |
| `src/field.jl` | Field reconstruction on the `(x, ξ)` grid, phase projection, nonlinear source `N_{l,q}` |
| `src/series.jl` | Frobenius horizon expansion and large-`r` expansion, solved order by order from the operator |
| `src/radial.jl` | Mesh, problem definition, first-order radial system in `r*` |
| `src/bvp.jl` | Unknown layout, MIRK4 residual, structured sparse Jacobian, Newton with line search |
| `src/initial_guess.jl` | Hydrogenic seed for the Newton iteration |
| `src/charges.jl` | Phase-averaged `E`, `J`, boundary fluxes, amplitude derivatives via the implicit-function theorem |
| `src/calibration.jl` | Amplitudes for prescribed per-cloud energies |
| `src/evolution.jl` | Tsit5 drivers on a fixed background and with black-hole spindown (background derivatives, metric work) |

## Conventions

* Units `M = 1` for the radial problem; energies in `f_a² M`, angular momenta in `f_a² M²`.
  The evolution with spindown uses the *initial* mass `M₀` as reference.
* Radial functions take the horizon offset `s = r − r₊` rather than `r`, which is exact at
  `r* = −80` where `s ≈ 1e-11`.
* The amplitude `A_i` of a bound mode is the leading outer coefficient of its fundamental
  channel: `R ≈ A e^{ikr} r^ν (1 + …)`.
* Frequencies and leading boundary coefficients are a separate parameter block of the Newton
  system (bordered sparse matrix) rather than constant ODE states; the equations are identical.
* The Jacobian is analytic (`rhs_derivatives`; the nonlinear source's Jacobian is a phase-grid
  contraction), and warm-started solves run a chord iteration on the previous factorisation.
* Meshes are uniform in `r*` by default and can be stretched geometrically beyond `stretch_from`
  so that the outer boundary sits far out; the charge integrals use non-uniform Simpson weights.
  `RadialMesh(background, α, mode_set)` sizes the domain from the parameter point: the inner
  endpoint from the horizon gap, the outer one from the slowest bound component's containment and
  the propagating components' asymptotic requirement (`domain_estimate` reports the ingredients).
* Black-hole spindown couples through `(f_a/M_Pl)²` (non-reduced Planck mass), passed to
  `evolve_with_spindown` as `decay_constant_ratio = f_a/M_Pl`.
* Each independent phase carries a wave action `N_i = Σ_q q_i n_q` (`phase_actions`), with
  `J = Σ m_i N_i`. The evolution drivers take `balance=:energy` (default, energy and angular momentum,
  two modes) or `balance=:action` (one action balance per phase, any number of modes, no metric work).
  `metric_work` is the action-based `(∂E)|_N` work, independent of the inner endpoint.
* The evolution state is `ln A` by default (`log_amplitudes`), with the first step set from the
  initial e-folding time; read amplitudes with `amplitude_trajectory(result)`.

## Quick start

```julia
using AxionCloud
background = KerrBackground(0.9)
modes = ModeSet([BoundMode(2, 1, 1), BoundMode(3, 2, 2)])
problem = RadialProblem(background, CosinePotential(), 0.3, modes)
solution = solve_radial(problem, [1e-6, 1e-7])
mode_frequencies(solution), cloud_charges(solution), boundary_losses(solution)
calibrated = amplitudes_for_energy(problem, [650.0, 650.0])      # amplitudes with E₂₁₁ = E₃₂₂ = 650 f_a²M
result = evolve_fixed_background(problem, [1e-6, 1e-7], (0.0, 1e8); verbose=true)
trajectory = amplitude_trajectory(result)                         # times, amplitudes[time, mode]
phase_actions(solution), amplitude_rates(solution; balance=:action)
spun = evolve_with_spindown(problem, [1e-6, 1e-7], (0.0, 1e8), 1e-2)   # f_a/M_Pl = 0.01
far = RadialMesh(background; rstar_min=-100.0, rstar_max=2000.0, step=0.4, stretch_from=150.0)
sized = RadialMesh(background, 0.3, modes; max_radius=2000.0)        # endpoints from the physics at this point
```

Diagnostic and convergence scripts live in `validation/` (logs in `validation/logs/`).

Progress, numerical findings and open questions are tracked in `NOTES.md`.

Run the tests with `julia --project -e 'using Pkg; Pkg.test()'` (Perlmutter: `module load julia/1.11.7` first).
The plot script lives in `scripts/` with its own environment: `julia --project=scripts scripts/profile_comparison.jl`.
