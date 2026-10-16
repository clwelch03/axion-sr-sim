# Joint two-cloud profiles, response, and evolution

This is the home for the proposed general multiphase model. The matching
residual and Newton wrapper are implemented and tested against the library.
Joint-field local and integrated E/J observables are also implemented.
Two-amplitude profile responses and a short E/J-balanced evolution are now
implemented in `response.jl` and `evolution.jl`. Analytic checks, matched-cloud
neighbor checks, and the 4/8-step cosine evolution pilot pass. General
action closure remains deferred to three independent clouds.

`load.jl` loads `core/`, cloud/component labels in `components.jl`, linear
operators, boundary assembly, and the full nonlinear radial RHS in `profiles.jl`,
and two-phase field reconstruction, projection, and nonlinear-source evaluation
in `projection.jl`, without importing the other models.
The root adiabatic loader leaves this library opt-in.
The loader now also includes `observables.jl`, `response.jl`, and `evolution.jl`.

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
polar/phase convergence for a fixed two-cloud cosine example. All 159 original
foundation checks pass within the expanded 363-check multiphase suite. This does not establish profile or retained-support convergence.

`nonlinear_two_phase_rhs!` is implemented in `profiles.jl`; its 21 checks in
`test_nonlinear_rhs.jl` are registered and pass against the library.
The implemented `two_phase_matching_residual` passes 24 registered checks in
`test_matching_residual.jl`, including integration histories and radial
step refinement. The implemented `solve_two_phase_mode` passes 20 registered checks in
`test_nonlinear_matching.jl`, including a four-component quartic matching solve
in three Newton corrections at modest resolution. The October 15 library run
reproduces the candidate result. This
does not yet establish a converged physical cloud family or evolution.

Keep the profile solver here. The two-cloud evolution uses total joint-field
E/J balance at fixed BH parameters, including full profile-response derivatives.
Phase actions are deferred to the extension to three independent clouds; a consistent slow-growth profile approximation remains required.

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
are registered in the runner and load observables through the library loader;
the local-observables test uses the current `radial_*` field names. All 159
existing checks and the 164 added checks passed in separate October 15 runs
(323 combined before the new response/evolution checks). See the
[handoff](../../../HANDOFF.md#current-multiphase-checkpoint-and-resume-point-october-15)
for the validation counts, integration gaps, and exact resume point.

See the [radial RK4 design](../../../docs/multiphase-rk4-design.md) for diagrams
of one radial step and the joint RHS, with state ordering and a
function map.

See the [implementation plan](../../../docs/multiphase-implementation-plan.md)
for the remaining boundary, projection, and closure work. The first physical
target is two mutually deforming clouds; three-cloud evolution will test the
additional closure equations.

## Implicit response and short evolution

`two_phase_response_problem` fixes the background, component support, grids,
boundary branches, and coefficient scale. `pack_two_phase_solution` converts a
matched solution to real Newton coordinates. `two_phase_real_residual` exposes
both amplitudes explicitly; `correct_two_phase_profile` corrects a trial at
fixed amplitudes. `two_phase_profile_response` differentiates the matching
system and then total E/J, including changes in profiles and frequencies.

`two_cloud_balance_rates` solves the two global balance equations and reports
conditioning after physically meaningful column/row scaling.
`two_phase_evolution_sample` combines a corrected profile, response, and rates.
`evolve_two_phase` uses that evaluator for explicit midpoint steps in log
amplitudes, correcting profiles at every stage and recording integrated losses.
New docstrings intentionally contain skeletons for the user to complete.

See the [campaign](../../../runs/adiabatic/multiphase_two_cloud/README.md) for
frozen-source commands, neighbor checks, trajectory files, and validation limits.

The completed pilot reaches t/M=199406.62 with amplitude changes of +3.04542%
and +0.00126903%. Neighbor E/J derivatives agree within 9.50e-9 relative;
charge-budget errors improve fourfold when the time step is halved. This tests
the discretized finite-domain family, not spatial convergence or the physical
slow-growth closure. The full Julia unit suite passes all 985 checks.
