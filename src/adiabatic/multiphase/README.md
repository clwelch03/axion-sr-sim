# Mutually deforming independent clouds — foundations

This is the home for the proposed general multiphase model. No multiphase
profile solver, action closure, or evolution is implemented yet.

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
polar/phase convergence for a fixed two-cloud cosine example. All 115 registered
checks pass. This does not establish profile or retained-support convergence.

`nonlinear_two_phase_rhs!` is implemented in `profiles.jl`; its 21 checks in
`test_nonlinear_rhs.jl` are registered and pass against the library. Short
radial integration and step-refinement checks are next.

See the [radial RK4 design](../../../docs/multiphase-rk4-design.md) for diagrams
of one radial step and the joint RHS, with state ordering and a
function map.

See the [implementation plan](../../../docs/multiphase-implementation-plan.md)
for the remaining boundary, projection, and closure work. The first physical
target is two mutually deforming clouds; three-cloud evolution will test the
additional closure equations.
