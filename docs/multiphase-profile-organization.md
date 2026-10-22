# Proposal: organize two-phase profiles and solver coordinates

Proposed on October 20, 2026. This document recommends a refactor; it does not
record an implemented or tested interface. Function signatures below describe
the intended API, rather than code ready to paste.

October 21 status: the physical-profile/coordinate migration, caller and test
updates, type annotations, and live-driver cleanup are now implemented.
The actual matching keyword is `save_radial_result`, replacing this proposal's
`save_profile`. The old unpacker and matching wrapper were removed. Boundary
assembly, admissibility unification, and adopting `TwoPhaseMatchedSolution`
remain deferred. See [HANDOFF.md](../HANDOFF.md) for current verification; the
recommendations below preserve the original proposal.

Use `TwoPhaseProblem` for the fixed setup and `TwoPhaseProfile` for the changing
frequencies and physical boundary coefficients. Make the physical matching
interface `two_phase_matching_residual(profile, problem; save_profile=false)`.
Centralize the conversion between a physical profile and the numerical
coordinates used by Newton and implicit differentiation.

This consolidates the current parallel ways of passing profile data: three
separate vectors, a tuple returned by unpacking, and a `TwoPhaseProfile`. Keep
the fixed setup separate because one problem is reused for many trial profiles
and amplitudes.

The recommendation follows the current responsibilities in
[profiles.jl](../src/adiabatic/multiphase/profiles.jl),
[response.jl](../src/adiabatic/multiphase/response.jl), and
[evolution.jl](../src/adiabatic/multiphase/evolution.jl). Names in the working
tree are currently being edited; `TwoPhaseResponseProblem` and
`TwoPhaseProfileProblem` should converge on the single proposed name
`TwoPhaseProblem`.

| Representation | Responsibility | Contents |
|---|---|---|
| `TwoPhaseProblem` | Define the fixed profile calculation and its numerical coordinates | Background, retained components, angular data, potential, integration settings, coefficient scale, coefficient ordering, and reference branches |
| `TwoPhaseProfile` | Describe one physical trial or matched profile | Two complex frequencies, all physical horizon coefficients, and all physical outer coefficients |
| Packed state plus prescribed amplitudes | Supply independent numerical inputs to Newton and differentiation | A real vector of unknowns, with the two prescribed outer amplitudes supplied separately |
| `TwoPhaseMatchedSolution` | Record a converged numerical solve | Packed state, prescribed amplitudes, scaled residual norm, and iteration count |

There are two physical concepts here: the problem and a particular profile.
The packed vector is a numerical representation of that profile, and the
matched solution adds the outcome of a solve. Neither constructing a profile
nor constructing a result object establishes convergence by itself.

`TwoPhaseProfile` continues to mean the boundary data that specify a radial
trial. It does not contain sampled radial functions. The histories returned
when `save_profile=true` remain integration outputs, with their existing
`from_horizon`, `from_outer`, and `residual` fields.

The problem should retain the simple field types already chosen:
`Float64` for `spin`, `alpha`, and `coefficient_scale`; `Symbol` for `potential`;
`Vector{Int}` for `cloud_m` and the free coefficient indices; and `BitVector`
for `initial_bound`. Keep type parameters for `components`, `couplings`,
`polar_values`, `grid`, the two phase-node collections, and `integration`.
These choices fix the background calculations to double precision while
allowing the profile's numeric type to vary during differentiation.

Keep the solver metadata in this same problem for now. The radial evaluator
uses only its physical and integration settings; the packing and solver
routines also use the scale, coefficient indices, and branch restrictions.
That division does not yet justify additional configuration structs.

Put `TwoPhaseProblem`, `TwoPhaseProfile`, and their coordinate conversion
helpers near the beginning of `profiles.jl`, before signatures that reference
them. Put `TwoPhaseMatchedSolution` there as well if both profile solvers use
it. The existing loader already loads `profiles.jl` before `response.jl`.
There is no need for a separate `types.jl` file or a new Julia module.

| Proposed interface | Contract |
|---|---|
| `pack_two_phase_profile(profile, problem)` | Select the free unknowns, scale boundary coefficients, and return their real packing |
| `unpack_two_phase_profile(packed_state, amplitudes, problem)` | Restore physical coefficients, insert prescribed amplitudes, and return a `TwoPhaseProfile` |
| `two_phase_matching_residual(profile, problem; save_profile=false)` | Evaluate physical radial matching, optionally retaining histories |
| `two_phase_real_residual(packed_state, amplitudes, problem)` | Unpack, evaluate matching, pack the complex mismatch, and divide by the fixed coefficient scale |
| `correct_two_phase_profile(problem, amplitudes, packed_state_guess; ...)` | Solve the real matching equations and return a `TwoPhaseMatchedSolution` |
| `two_phase_profile_response(problem, solution; ...)` | Differentiate the converged matching equations and the observables |

The unpacking interface above accepts the real-packed vector. It replaces
the current repeated combination of `unpack_complex` and
`unpack_two_phase_trial`. The existing unpacking arithmetic can be retained
internally while its output is changed to a `TwoPhaseProfile`. The generic
`pack_complex` and `unpack_complex` utilities stay in `core/newton.jl`.

Preserve the current coordinate convention. If there are N retained radial
modes, the complex unknown vector contains two frequencies, N horizon
coefficients, and N minus two free outer coefficients. There are therefore
2N complex unknowns, or 4N real unknowns. Boundary coefficients are divided
by `coefficient_scale`; frequencies are unscaled. All real parts precede all
imaginary parts. The matching mismatch has the same real dimension.

The prescribed amplitudes still need their own argument. A physical profile
contains them among its outer coefficients, but they are excluded from the
packed unknown vector. Newton holds them fixed, and response calculations
differentiate with respect to them independently. Combining them into the
Newton unknown vector would change the mathematical problem.

Use one shared rule for identifying the prescribed coefficient positions:
the first retained angular coefficient of q=(1,0), followed by the first of
q=(0,1), in cloud order. Derive the free positions from that rule and use it
consistently in setup, packing, and unpacking. Taking the complement of the
free indices alone does not establish cloud order when components are
reordered. For a profile whose fixed outer coefficients equal the supplied
amplitudes, packing followed by unpacking should recover the physical profile.
Supplying different amplitudes deliberately replaces those two coefficients.
Require both fundamentals when setting up a two-amplitude solve. The physical
matching evaluator should still support the reduced component sets used by
the existing single-cloud-limit tests, which do not need that solve setup.

Preserve numeric promotion in unpacking. When differentiation changes either
the packed state or the amplitudes, all three profile arrays must support the
resulting numeric type. In particular, amplitude differentiation must retain
ForwardDiff's derivative-carrying values even when the packed state contains
ordinary floating-point numbers. Do not convert those arrays to
`ComplexF64` inside the conversion or matching functions.

The observable adapter should unpack once, pass the profile to matching with
saved histories, and obtain the real carrier frequencies from
`profile.frequencies`. Keep the existing finite-domain integration and real
carrier convention. This refactor changes data organization, not the physical
time-derivative prescription.

Keep a single admissibility implementation with a documented input
representation. If `two_phase_profile_admissible` continues to receive a
packed vector, both Newton entry points should call it in that representation;
the older complex-state solver needs to pack its candidate first. Separately
decide the intended frequency restrictions before consolidating their checks:
the current standalone solver and the edited response predicate differ on
component positivity and exact-threshold rejection. A structural refactor
should not silently settle that behavior. Reintroducing the deleted second
finiteness check is not part of this proposal.

I recommend keeping `TwoPhaseMatchedSolution`'s current numerical fields for
this refactor and making the corrector actually return that type. The
evolution predictor currently computes a new guess from `solution.state`
and `response.state_derivative`; both use the same packed coordinates. A
physical profile can be reconstructed from the result when needed, without
storing another copy of it in the result.

This refines the earlier suggestion to make the result contain only a profile
and diagnostics. That alternative is possible, but it would also require
changing how the generic evolution routine obtains its packed state and
amplitudes. It is unnecessary to achieve the immediate consolidation of
physical interfaces. Keep the coefficient scale and ordering fixed across
each response calculation and its predictors; a derivative in one coordinate
system cannot be added to a state packed in another.

Implement the proposal in small steps:

1. Establish the common names and move the type definitions into `profiles.jl`.
   Update annotations, constructors, and docstrings together. Check that the
   opt-in loader works, including repeated loading.
2. Add the profile-returning unpacking interface and centralize coefficient
   ordering. Reuse the coordinate tests for reordered components, scale
   changes, and derivatives with respect to both numerical inputs.
3. Change matching to accept `profile, problem`. Update the standalone solver,
   response adapter, matching tests, and saved-history tests. Keep integration
   defaults in one setup location and populate `problem.integration` completely
   before evaluation; `save_profile` remains a per-call keyword.
4. Route the physical observable adapter through the same profile. Retain a
   thin packed matching adapter only if it usefully serves both residual and
   observable evaluation; its sole job should be unpacking and forwarding.
5. Return the existing matched-result type from the response corrector and
   update its consumers. The standalone solver may retain its established
   return shape temporarily, but should use the same conversion and physical
   matching functions internally.
6. Update the live campaign driver and documentation. Preserve frozen source
   snapshots and saved numerical results. If saved checkpoints contain the old
   problem type or a named tuple, handle their conversion explicitly before
   using them with the new typed interfaces.

Validation should establish unchanged numerical behavior at each step:
packing round trips with matching prescribed amplitudes; preservation of
amplitude derivatives; agreement with the existing single-cloud matching
reference; identical residual-only and saved-history results; and a corrected
joint profile with unchanged observables and responses within existing
tolerances. Once the numerical result type is used, include the existing
evolution tests because they consume `solution.state` and
`solution.amplitudes`. These are proposed checks, not results of this document.

The first useful stopping point is small: one fixed problem object, one
physical profile object, and an explicit pair of conversion functions used by
both solvers. Broader result-format changes can be evaluated after that
organization is working.
