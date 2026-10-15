# Mutually deforming independent clouds — foundations

This is the home for the proposed general multiphase model. The matching
residual and Newton wrapper are implemented and tested against the library.
Joint-field local and integrated E/J observables are also implemented.
Two-amplitude profile responses, action closure, and evolution are not yet
implemented. The next requested step is direct implicit differentiation of
the matching residual, followed by independently re-solved neighbor checks.

`load.jl` loads `core/`, cloud/component labels in `components.jl`, linear
operators, boundary assembly, and the full nonlinear radial RHS in `profiles.jl`,
and two-phase field reconstruction, projection, and nonlinear-source evaluation
in `projection.jl`, without importing the other models.
The root adiabatic loader leaves this library opt-in.
`observables.jl` currently requires an explicit include after this loader.

Run the unit checks with
`julia --project=. tests/adiabatic/multiphase/runtests.jl`.
The radial and source tests load the single-cloud library for comparison only. Current
checks validate label bookkeeping and linear assembly, not profile convergence.
An explicit two-phase quadrature benchmark checks library reconstruction by
separating equal-angular-pattern components. The function
`reconstruct_two_phase_field_at_radii` evaluates one radius per call.
The reusable `project_two_phase_field` passes coefficient recovery, absent-component,
conjugate, and analytic sine-field checks. `nonlinear_two_phase_source` passes
five quartic checks: two exact driven-source formulas, two forbidden outputs,
and a zero field. Sources also agree with the single-cloud implementation for
quartic, sextic, and cosine potentials when only q=(1,0) is populated, including
the generated third harmonic; an output with nonzero q_2 vanishes. These are
fixed-grid checks. `test_source_validation.jl` adds the complete quartic table,
phase-origin covariance and averaged-potential invariance, and independent
polar/phase convergence for a fixed two-cloud cosine example. All 159 registered
checks pass. This does not establish profile or retained-support convergence.

`nonlinear_two_phase_rhs!` is implemented in `profiles.jl`; its 21 checks in
`test_nonlinear_rhs.jl` are registered and pass against the library.
The implemented `two_phase_matching_residual` passes 24 registered checks in
`test_matching_residual.jl`, including integration histories and radial
step refinement. The implemented `solve_two_phase_mode` passes 20 registered checks in
`test_nonlinear_matching.jl`, including a four-component quartic matching solve
in three Newton corrections at modest resolution. The October 15 library run
reproduces the candidate result. This
does not yet establish a converged physical cloud family or evolution.

Keep the profile solver here. The first two-cloud evolution will use total
joint-field E/J balance at fixed BH parameters, including full profile-response
derivatives. Phase actions are deferred to the extension to three independent
clouds; a consistent slow-growth profile approximation remains required.

`observables.jl` evaluates one scalar stress tensor from the joint field and
its coordinate derivatives. `two_phase_observables_at_radius` returns E/J
per unit tortoise radius and signed `radial_energy_flux` and
`radial_angular_momentum_flux`, positive toward increasing radius.
`two_phase_integrated_observables` integrates those densities over the supplied
finite domain and returns horizon/outer losses, negating the inner currents.
Its top-level helpers are `ordered_radial_branch` and
`integrate_two_phase_branch_observables`. A radial branch means spatial
profile samples, not time evolution or a family of amplitude values.

These evaluators require real carrier frequencies and omit slow profile time
derivatives; they do not interpret the matching solver's complex trial
frequencies as a physical growth law. No tail correction or converged physical
profile family is established. The 164 observable and field-derivative checks
are now registered in the runner, with an explicit include of `observables.jl`;
the local-observables test uses the current `radial_*` field names. All 159
existing checks and the 164 added checks passed in separate October 15 runs
(323 combined). Library-loader inclusion remains pending. See the
[handoff](../../../HANDOFF.md#current-multiphase-checkpoint-and-resume-point-october-15)
for the validation counts, integration gaps, and exact resume point.

See the [radial RK4 design](../../../docs/multiphase-rk4-design.md) for diagrams
of one radial step and the joint RHS, with state ordering and a
function map.

See the [implementation plan](../../../docs/multiphase-implementation-plan.md)
for the remaining boundary, projection, and closure work. The first physical
target is two mutually deforming clouds; three-cloud evolution will test the
additional closure equations.
