# Mutually deforming independent clouds — foundations

This is the home for the proposed general multiphase model. The matching
residual and Newton wrapper are implemented and tested against the library.
No action closure or evolution is implemented yet.

`load.jl` loads `core/`, cloud/component labels in `components.jl`, linear
operators, boundary assembly, and the full nonlinear radial RHS in `profiles.jl`,
and two-phase field reconstruction, projection, and nonlinear-source evaluation
in `projection.jl`, without importing the other models.
The root adiabatic loader leaves this library opt-in.

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

See the [radial RK4 design](../../../docs/multiphase-rk4-design.md) for diagrams
of one radial step and the joint RHS, with state ordering and a
function map.

See the [implementation plan](../../../docs/multiphase-implementation-plan.md)
for the remaining boundary, projection, and closure work. The first physical
target is two mutually deforming clouds; three-cloud evolution will test the
additional closure equations.
