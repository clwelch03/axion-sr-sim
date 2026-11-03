# Joint two-cloud profiles, response, and evolution

This is the home for the proposed general multiphase model. The matching
residual and Newton wrapper are implemented and tested against the library.
Joint-field local and integrated E/J observables are also implemented.
Two-amplitude profile responses and a short E/J-balanced evolution are now
implemented in `response.jl` and `evolution.jl`. Analytic checks, matched-cloud
neighbor checks, and the 4/8-step cosine evolution pilot pass. General
action closure remains deferred to three independent clouds.

`load.jl` loads `core/`, cloud/component labels in `components.jl`, linear
operators, boundary states, and the full nonlinear radial RHS in `radial.jl`,
profile types, coordinate conversions, and matching in `profiles.jl`,
and two-phase field reconstruction, projection, and nonlinear-source evaluation
in `projection.jl`, without importing the other models.
The root adiabatic loader leaves this library opt-in.
The loader also includes `observables.jl`, `averaged_observables.jl`,
`response.jl`, and `evolution.jl`.

Amplitude continuation is available as
`continue_two_phase_amplitude(seed, output_prefix; target_amplitude, kwargs...)`.
Its predictor and profile correction live in `profiles.jl`, alongside matching;
adaptive stepping and checkpoint handling live in `core/continuation.jl` and
`core/continuation_io.jl`. The campaign's `scripts/continue_amplitude.jl` only
selects inputs, output prefix, and run settings. Matched seeds and per-attempt
Newton recovery files retain their existing formats.

Run the unit checks with
`julia --project=. tests/adiabatic/multiphase/runtests.jl`.
The radial and source tests load the single-cloud library for comparison only. Current
checks validate label bookkeeping and linear assembly, not profile convergence.
An explicit two-phase quadrature benchmark checks library reconstruction by
separating equal-angular-pattern components. The function
`reconstruct_two_phase_field_at_radius` evaluates one radius per call.
The reusable `two_phase_field_component_coefficients` passes coefficient recovery, absent-component,
conjugate, and analytic sine-field checks. `nonlinear_two_phase_source` passes
five quartic checks: two exact driven-source formulas, two forbidden outputs,
and a zero field. Sources also agree with the single-cloud implementation for
quartic, sextic, and cosine potentials when only q=(1,0) is populated, including
the generated third harmonic; an output with nonzero q_2 vanishes. These are
fixed-grid checks. `test_source_validation.jl` adds the complete quartic table,
phase-origin covariance and averaged-potential invariance, and independent
polar/phase convergence for a fixed two-cloud cosine example. All 159 original
foundation checks pass within the expanded 363-check multiphase suite. This does not establish profile or retained-support convergence.

`nonlinear_two_phase_rhs!` is implemented in `radial.jl`; its 21 checks in
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

Response and evolution samples now default to `observable_method=:averaged`.
`averaged_observables.jl` contracts the quadratic stress-tensor terms using
phase orthogonality, retaining interference among angular modes with the
same Fourier label. The full nonlinear potential still uses the original
polar/phase quadrature. Its first derivatives use the equivalent local chain
rule, `U'(Phi) = Phi - potential_remainder(Phi)`. Geometry and angular data
are cached only within a fixed-problem response calculation; profile values
and frequencies are always reevaluated and differentiated. AD chunk sizes,
radial matching, boundary conditions, integration and balance equations are
unchanged.

Use `two_phase_evolution_sample(problem, amplitudes, guess;
observable_method=:direct)` to select the original tensor/phase sums.
The same keyword is available on `two_phase_profile_response` and
`two_phase_trial_observables`; the original low-level observable functions
remain available. The averaged path rejects phase grids that alias quadratic
products, zero labels and separately stored conjugate/duplicate labels.
Keep all angular harmonics of one Fourier label in one component. Cached
backgrounds, grids and bases must be treated as read-only.

The averaging kernel accepts arbitrary-length phase labels and is checked
against an explicit three-phase tensor sum. The nonlinear potential quadrature
still grows with the number of phases, and the matching/evolution system
remains two-cloud; this change does not supply additional balance equations.
See the [benchmark record](../../../runs/adiabatic/multiphase_two_cloud/README.md#november-2-analytic-phase-averaging)
for measured complete-step timings and rate comparisons.

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

The response/evolution driver directly constructs a `TwoPhaseProblem` containing the
fixed background, component support, grids, boundary branches, and coefficient
scale. `pack_two_phase_profile` converts a physical `TwoPhaseProfile` to
real Newton coordinates; `unpack_two_phase_profile` reconstructs it from those
coordinates, prescribed amplitudes, and the problem. Physical matching accepts
`(profile, problem; save_radial_result=false)`. `two_phase_real_residual` exposes
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
slow-growth closure. The historical pilot source passed all 985 Julia unit
checks. On November 3 the complete current suite passed 1,399 assertions,
including all 555 multiphase checks, after removal of a stale regression for
the retired standalone AD diagnostic. See the campaign's analytic
phase-averaging record for the preserved historical checks and timings.
