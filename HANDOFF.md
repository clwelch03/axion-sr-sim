# Project handoff

Updated 2026-10-21. This is the authoritative project status: what works,
open caveats, and next steps. Campaign READMEs own exact commands, saved-result
paths, and detailed validation records; design documents own equations and
implementation decisions.

Active work: reproducing Omiya, Takahashi, Tanaka & Yoshino,
*Impact of multiple modes…* ([2211.01949v1](docs/papers/2211.01949v1.pdf)),
Sec. III — a nonlinear l = m = 1 primary with a linear l = m = 2 secondary on it,
evolved by energy and angular-momentum balance (Eqs. 62–63). The single-cloud
adiabatic work (OTT 2201.04382) is the primary's input; see
[adiabatic-repro.md](docs/adiabatic-repro.md).

Read first: [multimode-adiabatic-design.md](docs/multimode-adiabatic-design.md)
(equations, decisions, component map) and the
[secondary campaign README](runs/adiabatic/multimode_secondary/README.md)
(commands, results, checks). Cite equation numbers from the PDF in
`docs/papers/`, not the arXiv HTML.

## Where things stand

### Resume here: multiphase organization refactor (October 21)

The physical-profile and solver-coordinate migration is complete for this pass.
The user authorized the solver migration and the remaining annotation, test,
driver, naming, and documentation cleanup. The broader docstring overhaul is
unfinished. Ordinarily the user prefers to make library edits themselves with
small guided steps; this task-specific authorization does not replace that
preference. Preserve the user's learning-map checkboxes.

Current interfaces and responsibilities:

- `TwoPhaseProblem` stores fixed settings; `TwoPhaseProfile` stores complex
  frequencies and physical boundary coefficients. Neither stores sampled
  radial functions. Both are defined in `profiles.jl`.
- `pack_two_phase_profile` and `unpack_two_phase_profile` convert between a
  physical profile and 4N real coordinates (2N complex unknowns). Frequencies
  are unscaled; boundary coefficients use the fixed problem scale. The two
  prescribed amplitudes are separate inputs to unpacking. Numeric promotion
  preserves derivatives from either the packed state or amplitudes.
- `two_phase_matching_residual(profile, problem; save_radial_result=false)`
  evaluates physical matching, including reduced single-cloud component sets.
  Saving returns `from_horizon`, `from_outer`, and `residual`. The integration
  calls translate the flag to the core keyword `save_profile`.
- `solve_two_phase_mode` constructs a problem and packs its complex Newton
  candidates before unpacking profiles. Integration defaults, the Newton
  backend, scaling, admissibility behavior, and return fields are preserved.
- Response residuals and observables use the same unpacker and matching
  interface. `two_phase_trial_matching` and `unpack_two_phase_trial` are removed.
  Library and live-driver definitions accepting a problem now annotate it as
  `TwoPhaseProblem`; physical profile arguments use `TwoPhaseProfile`.
- `two_phase_integrated_observables` calls its integration bundle
  `radial_matching_result`; this bundle is not a `TwoPhaseProfile`.
- The unpacker, standalone solver, and live response/evolution driver use
  `two_phase_coefficient_indices`. Preserve the descriptive field names
  `fixed_outer_coefficient_indices` and `free_outer_coefficient_indices`.
  The driver's admissibility check now receives its already-packed guess
  directly.
- `radial.jl` contains `component_tortoise_operator`, `component_boundary_state`,
  and `nonlinear_two_phase_rhs!`. Keep these helpers independent of the profile
  types. The opt-in loader includes this file before `profiles.jl`.

Validation against the edited working tree:

- All **97 focused checks** passed in a fresh Julia process: 24 physical
  matching checks, 49 integrated-observables checks, and 24 coordinate checks.
  They cover single-cloud limits for three potentials, ordering with unequal
  angular blocks, radial RK4 refinement, saved versus plain residuals,
  manufactured charges, boundary signs, packing/scaling, exact Jacobians with
  respect to state and amplitudes, and invalid inputs. Command:

  ```bash
  julia --startup-file=no --project=. -e 'using Test; include("src/adiabatic/multiphase/load.jl"); include("tests/adiabatic/multiphase/test_matching_residual.jl"); include("tests/adiabatic/multiphase/test_integrated_observables.jl"); include("tests/adiabatic/multiphase/test_response_coordinates.jl")'
  ```

- The preceding solver migration passed all **20 nonlinear matching checks**,
  including a reordered-component solve. Newton took three iterations to a
  scaled residual norm of approximately 1.22e-10 on the test's modest grid.
- The live response/evolution driver loaded without running its campaign;
  two setup checks confirmed the prescribed and free coefficient indices.
- These are algorithm/interface checks, not physical spatial convergence.
  The full suite and campaign were not rerun; saved campaign results below
  remain historical evidence for their frozen source.

Deferred work and known caveats:

- Continue the user's docstring walkthrough using the
  [learning map](docs/multiphase-learning-map.md). The
  [organization proposal](docs/multiphase-profile-organization.md) includes
  later work beyond this completed interface migration.
- Boundary-assembly cleanup and making the corrector return
  `TwoPhaseMatchedSolution` remain deferred. The corrector still returns a
  named tuple with packed state, amplitudes, and diagnostics.
- Resolve admissibility semantics before sharing predicates: the standalone
  solver rejects nonpositive component frequencies and exact thresholds;
  the response predicate checks fundamental bounds and branch classification.
- `two_cloud_balance_rates` returns `log_rates`, `scaled_condition`, and
  `balance_residual`, but `test_evolution.jl` still expects `amplitude_rates`.
  That pre-existing mismatch is outside this pass and remains unresolved.
- Preserve existing working-tree changes, frozen source snapshots, and saved
  results. Saved objects require explicit compatibility review before reuse
  with new interfaces; do not resume a historical campaign under edited source.

### Current multiphase checkpoint and resume point (October 15)

The joint profile solver, finite-domain E/J observables, two-amplitude implicit
responses, and a short-evolution integrator are implemented. On October 15 the
user authorized direct implementation of all three response/evolution steps,
with readable comments and docstring skeletons left for them to complete.
The matched-cloud derivative validation and short cosine evolution are complete;
results, files, and commands are in the
[joint two-cloud campaign](runs/adiabatic/multiphase_two_cloud/README.md).

- `multiphase/load.jl` now includes observables, `response.jl`, and `evolution.jl`.
  It remains opt-in; the other models do not import the joint model.
- The two positive amplitudes fix the first retained outer coefficient in
  q=(1,0) and q=(0,1). E and J remain outputs, not amplitude coordinates.
  `unpack_two_phase_trial` is shared by the original solver and response code.
- `response.jl` supplies a real-packed matching residual, a Newton corrector,
  and total E/J derivatives. The coefficient/residual scale is fixed at the
  reference amplitude rather than differentiated. `core/sensitivity.jl`
  factors the state Jacobian once for both amplitude directions and applies
  the full observable chain rule. Existing scalar-amplitude APIs remain intact.
- The observable adapter explicitly takes real parts of the complex matching
  frequencies as carriers, omitting slow profile time derivatives. This
  specifies the numerical family; it does not validate the physical
  slow-growth approximation. Charges are integrated over the supplied finite
  radial domain without tail corrections. Horizon losses negate inner currents.
- `evolution.jl` solves the 2x2 E/J balance system, reports a dimensionless
  condition number after row/column scaling, and rejects ill-conditioned
  matrices. Explicit midpoint steps evolve log amplitudes. Every stage is
  independently Newton-corrected from an implicit-response predictor.
  Independently accumulated boundary losses supply numerical charge budgets.
- All 985 Julia unit checks passed in Slurm: 13 loader, 350 core, 243
  single-cloud, 16 primary/secondary, and 363 multiphase. These include 14
  vector-response, 13 unpacking, 8 balance/scaling, and 19 midpoint/budget
  checks. The campaign then validates both E/J response columns against
  independently corrected neighbors at 1% and 0.5% offsets, checks fixed-scale
  invariance, and compares 4-step and 8-step cosine trajectories.
- Completed pilot job 26978242 uses snapshot `response_20261015_v2`, matching
  tolerance 1e-9, and prefix `cosine_pilot_20261015_v2`. Both derivative columns
  agree with independently corrected neighbors to at worst 9.50e-9 relative
  in E/J and 2.07e-9 in the state response. Changing the fixed coefficient scale
  by seven preserves the E/J derivatives.
- Four-step and eight-step cosine trajectories reach t/M=199406.6212.
  A1 grows 3.04542%, A2 grows 0.00126903%. Endpoint amplitudes differ by at most
  7.06e-11 relative when halving dt. Eight-step budget errors are 4.84e-10 of
  initial E and 2.46e-10 of initial J (2.08e-6 and 1.88e-6 of their respective
  accumulated changes), improving by a factor of four from the coarse run.
  The scaled balance condition stays between 1140 and 1211; saved endpoint
  matching residuals are below 1.97e-10. These are finite-grid numerical checks.
- The preceding job 26977745 passed all 985 unit checks, then stopped in a
  linear seed solve at residual 1.26e-10 above its requested 1e-11 tolerance.
  The corrected pilot uses the established 1e-9 matching tolerance. Both source
  snapshots and logs are preserved; see the campaign README for exact recipes.

The next scientific checkpoints remain spatial and retained-support
convergence, finite-boundary/tail sensitivity, and a consistent slow-growth
formulation. An accurate numerical E/J budget alone does not establish these.
Phase actions remain deferred to the three-independent-cloud extension.

The detailed implementation evidence and other project tracks follow.

- **Multiphase foundations (October 15):** independent-cloud labels, Fourier
  component construction, single-component linear radial assembly, and
  two-phase field reconstruction, projection, nonlinear-source evaluation, and
  full radial RHS assembly are present in the opt-in multiphase library.
  All 159 multiphase checks pass: 12 frequency/azimuthal-label
  checks, 7 radial-operator checks, 8 reconstruction/reference-quadrature checks,
  5 library-projector checks, 5 quartic-source checks, 12 single-cloud-limit checks,
  14 quartic-table checks, 18 phase-origin checks, 13 cosine-convergence checks,
  21 full-RHS checks, 24 matching/integration checks, and 20 Newton matching checks.
  The existing 10 loader checks also passed at
  the linear-assembly checkpoint. The radial checks compare to the single-cloud
  assembly and an analytic Schwarzschild m=0 result, and check signed-frequency/
  azimuthal symmetry and angular ordering. The phase benchmark separates two
  coefficients sharing l=m=3 to 3.11e-16 absolute error; an azimuth-only
  projection returns their sum. The field now comes from the library's
  `reconstruct_two_phase_field_at_radii`, which evaluates one radius per call;
  projection is checked both with an explicit reference sum and with the library's
  `project_two_phase_field`. Library checks also cover an absent component,
  conjugate recovery, and an analytic sine field on unequal phase grids.
  The source checks compare the two driven components to exact integrated
  quartic formulas and verify angular-parity, absent-phase, and zero-field zeros.
  With only q=(1,0) populated, the fundamental and generated third-harmonic
  sources agree with the existing single-cloud implementation for quartic,
  sextic, and cosine potentials; a q_2=2 output vanishes. The full quartic table,
  independent and combined phase shifts, and phase-averaged potential invariance
  now pass. A separate two-cloud cosine example converges under independent
  polar/phase refinement: 32 nodes in each refined direction agree with a 64^3
  reference to at most 1.24e-13 absolute source error. A warmed 32^3 evaluation
  including table setup took about 0.020 seconds and allocated 4.98 MB on this
  host (allocation total, not peak RSS); one real field occupies 256 KiB.
  This establishes quadrature convergence for the tested fixed support, not
  convergence of solved profiles or a universal grid choice.
  The RHS is implemented in `profiles.jl`; all 21 checks in
  `test_nonlinear_rhs.jl` are registered and pass against the actual library.
  These cover the single-cloud limit for all potentials, ordering, state
  preservation, dimensions, driven-source signs/scaling, and zero-field slopes.
  The state stores all R values then all Q values,
  ordered by component and then l; the source is computed once from all R.
  The implemented `two_phase_matching_residual` now connects this RHS to the shared
  integrator, assembles boundary blocks into global [R; Q] order, and returns
  both profile and slope mismatches. All 24 registered checks in
  `test_matching_residual.jl` pass against the library: full single-cloud integration histories for
  all potentials, saved/plain output agreement, joint ordering, input checks,
  and radial refinement in both directions (error-reduction ratios 17.60 and
  14.28 for step sizes 0.4/0.2/0.1, consistent with fourth-order convergence).
  The inserted `solve_two_phase_mode` passes 20 registered checks in
  `test_nonlinear_matching.jl` against the actual library. It fixes the first retained
  l coefficient in each fundamental q block, scales coefficients/residuals by
  the smaller positive outer normalization, and solves for two complex trial
  frequencies plus free boundary coefficients. Frequencies remain on their
  initial positive-frequency bound/outgoing branches.
  The four-component quartic trial at spin=0.99, alpha=0.42, outer normalizations
  [1e-4,1e-5], l=[1,2,0,3], inner=-30, outer=80, match=5, radial step=0.4,
  and quadrature 8x9x9 converges in 3 corrections (scaled residual 13.3874 ->
  1.2248e-10). Both driven outer coefficients are nonzero. This is a matched
  numerical prototype, not an endpoint/angular/radial convergence result or
  a validated slow-growth prescription. The October 15 library run reproduces
  the candidate result; all 159 foundation checks pass (now part of the larger suite).
  A positive-frequency
  boundary wrapper is present; its prototype checks still need adding to the
  repository suite. A converged physical joint-profile family and the action closure
  remain pending. The October 15 sequencing decision keeps total joint-field
  E/J balance for the first two-cloud evolution at fixed BH parameters, with
  full amplitude derivatives and matrix-conditioning checks. Phase actions
  are deferred until the extension to three independent clouds; the profile
  slow-growth approximation still needs to be established. Keep the profile
  solver in `multiphase/`; next validate boundaries and matched-profile
  convergence before physical continuation. Run
  `julia --project=. tests/adiabatic/multiphase/runtests.jl`.
  The first observables increment is
  `two_phase_potential_energy_density` for the joint potential contribution
  per unit tortoise radius. Its 10 checks in `test_potential_energy.jl` pass
  against the inserted helper, loaded explicitly: single-cloud agreement for all potentials,
  phase-origin invariance, an analytic quartic interaction term, and zero/input
  checks. The helper is in `multiphase/observables.jl`; loader inclusion and test
  registration are now complete. These checks are additional to the original 159.
  The complete local tensor-based observables candidate is described below.
  This potential-only check requires no matched profile or
  choice of time-derivative convention for complex trial frequencies.
  The user selected a shared stress-tensor architecture instead of implementing
  each energy contribution separately. Its first candidate helpers,
  `build_polar_theta_derivatives` and `reconstruct_two_phase_field_derivatives`,
  pass 32 checks in `test_field_derivatives.jl`. They reconstruct the real field
  and its Boyer-Lindquist coordinate derivatives (t,r,theta,varphi), with
  Q=dR/dr_* converted to dR/dr. Frozen-field time derivatives require explicitly
  supplied finite real carrier frequencies; complex trial frequencies are
  rejected, and slow amplitude/profile time derivatives are not included.
  Angular derivatives agree with automatic differentiation; all four field
  derivatives agree with an independently differentiated field to <4.5e-16.
  Signed/zero m, repeated angular patterns across distinct q, multiple l in
  one q, input preservation, and zero fields are checked. The helpers have
  been inserted into `multiphase/projection.jl`, including the user-requested
  interface that builds the polar derivative table inside reconstruction from
  `grid`. All 32 checks pass and are registered in the suite.
  The derivative checks do not validate a slow-growth model.
  Two further candidates, `kerr_metric` (for `core/geometry.jl`) and
  `scalar_stress_tensor` (for `multiphase/observables.jl`), pass 115 checks:
  86 metric checks in `core/test_kerr_metric.jl` and 29 tensor/contraction checks
  in `multiphase/test_stress_tensor.jl`. The metric accepts r_* as its position
  input but returns Boyer-Lindquist components in (t,r,theta,varphi) order,
  with signature (-+++). The pointwise tensor returns T_mu_nu/Fa^2.
  Raising its first index gives energy density -T^t_t, angular-momentum density
  T^t_phi, and outward currents -T^r_t and T^r_phi. Test quadrature includes the
  rho^2 spatial measure, 2pi azimuthal factor, phase average, and dr/dr_* only
  for charge densities. Energy agrees with the single-cloud routine for real
  carriers and all three potentials to 3.1e-15 relative error. Joint-field
  fluxes and Schwarzschild angular momentum agree with analytic expressions.
  `scalar_stress_tensor` and `kerr_metric` are now inserted; the 115 checks
  pass against the actual library and are registered. Observables load through
  the multiphase loader. Local checks
  do not establish radial-integrated charges, endpoint convergence, or evolution.
  The inserted `two_phase_observables_at_radius` now passes all 44
  checks in `test_observables_at_radius.jl` against the actual library. It returns energy and
  angular momentum per unit tortoise radius plus both outward currents from
  one quadrature loop. The single-cloud energy agrees to 2.0e-15 relative error
  for three spins and all potentials. Checks also cover joint analytic fluxes,
  Schwarzschild angular momentum, flux reversal, phase-origin invariance,
  exact quartic cross energy, zero fields, input preservation/validation,
  and automatic derivatives of a prescribed local field scaling against finite
  differences. That last check is not a matched-profile amplitude response.
  Observables loader inclusion and test registration are now complete.
  All observables use explicit real carrier frequencies.
  The `two_phase_integrated_observables` candidate passes 49 checks in
  `test_integrated_observables.jl`. It accepts saved matching histories,
  orients each branch toward increasing r_*, requires uniform spacing and a
  common match position, reuses `integrate_uniform` for total E and J, and
  returns horizon/outer energy and angular-momentum losses. Horizon losses
  reverse the inward endpoint's outward-current sign. Checks recover exact
  cubic E/J integrals, including odd interval counts and reversed storage,
  the existing single-cloud integrated energy for all potentials, endpoint
  signs, component reordering, and saved-history compatibility. Nonfinite
  profile samples or local observables are rejected. The first short-domain
  trial fixture overflowed the observables (finite profiles grew to 4.9e247);
  reducing its boundary coefficients by 1e-8 gave a finite dilute fixture.
  Manufactured snapshots and trial histories test bookkeeping, not physical
  matching or convergence. The radial wrapper is now inserted into
  `multiphase/observables.jl`; the local flux fields are named
  `radial_energy_flux` and `radial_angular_momentum_flux` (positive toward
  increasing radius). At the user's request, a tested refactor moves the
  nested helpers to top-level `ordered_radial_branch` and
  `integrate_two_phase_branch_observables`, passing their inputs explicitly.
  All 49 checks passed with that refactor, which is now inserted; comments
  and the radial-branch docstring are also present. The tests are registered.
  Outputs cover only the supplied
  finite radial domain; tails and radial/boundary/truncation convergence remain
  separate. Matched-profile amplitude responses are implemented as described
  in the current checkpoint above.
- **Adiabatic layout refactor (October 12):** branch
  `refactor/adiabatic-layout` groups source and tests by `core/`,
  `single_cloud/`, and `primary_secondary/`. At that checkpoint `multiphase/`
  had only a scope README and core-only loader; the October 15 work above
  supersedes that historical status. Shared operators,
  potentials, Newton/integration routines, quadrature, and sensitivities were
  separated from single-cloud assembly. Reusable workflow functions now load
  from model libraries; model drivers contain command execution. The root
  Julia loader remains available. See the [source guide](src/adiabatic/README.md)
  for current paths and imports, and the
  [implementation plan](docs/multiphase-implementation-plan.md) for the next
  physics milestones. Validation: 509 grouped Julia unit checks plus 10 loader
  checks, 5 Python checks, and the 44-check coupled-balance benchmark pass.
  A nonlinear profile, energy, frequencies, and both complete saved-table
  Fig. 6 trajectories agree exactly before/after the refactor. The QCD grid's
  independent reference-trajectory check retains its 1.71e-8 maximum spin
  difference. The Python wheel builds and contains all 51 Julia source/driver
  files. Numerical campaign artifacts and frozen manifests retain their
  locations and contents; source hashes intentionally reject moved/changed
  implementations. This validates the refactor, not new numerical convergence
  or multiphase physics.
- **QCD grid stopped before disconnect (October 12):** the requested 10x10 log
  grid spans 2–100 solar masses and 1e-14–1e-11 eV, with the final exclusion
  `alpha0 > 0.45`, leaving 83 cells (maximum 0.437617). The user selected a
  10-million-year age cap. Inputs, shared-source hashes, scripts, and restart
  instructions are in the [grid campaign](runs/adiabatic/qcd_grid_20261012/README.md).
  Initial stationary checks at alpha=0.149663 and 0.437617 passed; the
  alpha=0.0149663 check failed at an absolute shooting-residual floor. Array
  26875559 was automatically canceled before starting any production cases.
  The remaining diagnostic 26875867 was canceled at the user's disconnect
  request. No survey jobs or automatic retries remain active. Do not resubmit
  unchanged: the low-alpha conditioning and an end-to-end numerical cell need
  checking first. Forty planned cells are explicitly labeled linear age
  screens, not nonlinear simulations; only two screens have been run locally.
  The new time driver reproduces the saved reference spin curve to 1.71e-8.
- **First coupled QCD spin-down trajectory:** at the user's explicit request,
  a preliminary run now evolves a 20-solar-mass BH from chi=0.99 with
  `m_a=1e-12 eV`, using a nonlinear 211 primary, a 322 secondary on its averaged
  background, and driven m=0,3 scalar channels. The QCD relation fixes
  `f_a=5.7e18 GeV`; small classical seeds are specified in the
  [campaign README](runs/adiabatic/qcd_spindown/README.md).
  A fresh 21-point stationary track at alpha approximately 0.150–0.137
  supplies the leading coupled balance. It reaches chi=0.75 in 33.245 days,
  versus 29.789 days for two linear nonrelativistic clouds with the same seeds;
  the one-year spins are 0.508953 and 0.509951. PNG/PDF plots, trajectory CSVs,
  radial checks, and provenance are saved in that campaign. All 80 stationary
  comparisons pass; time refinement changes spin by at most 1.15e-8.
  Initial secondary horizon-flux radial sensitivity is 0.47%, while the
  primary's is below 5.5e-7 relative. This is a
  trajectory-specific interpolation with refreshed horizon factors, not a
  general background table. GW losses, self-gravity, and the subleading metric
  work/background angular response remain absent. The cloud reaches about
  9.2% of the remaining BH mass; timing is preliminary under that approximation.
  An October 12 field-amplitude audit restores `Phi=fa*phi` along this saved
  trajectory: sampled primary `|Phi|/fa` grows from `5.38e-6` to `9.21e-4`
  (171 times), or `3.06e13` to `5.25e15 GeV`. The cubic/mass ratio stays below
  `1.42e-7`; comparison with the binding scale gives below `2.96e-5`.
  Eight post-processing checks pass; the campaign README records the plot,
  CSV, normalization, and limits of this primary-only amplitude diagnostic.
- **Multimode adiabatic evolution:** all six Sec. III components are implemented
  at $M\mu = 0.42$ and $0.3$, $a/M = 0.99$: averaged primary background,
  secondary eigenproblem and energy, driven $m = 0$ and $m = 3$ fluxes, and
  evolution from energy and angular-momentum balance. The nonperturbative
  Fig. 6 curves agree with the paper to reading precision at $M\mu = 0.42$;
  at $0.3$ the primary peak is about 10% lower and plateaus 5–10% higher.
  Fig. 7 profile peaks at $M\mu = 0.3$ agree to reading precision. See the
  [secondary campaign](runs/adiabatic/multimode_secondary/README.md) for the
  saved tables, comparisons, and checks. Pieces 5–6 are committed in `5b1f322`.
- **Single-cloud adiabatic work:** cosine, quartic, and sextic stationary
  branches are preserved in the [Figure 2 campaign](runs/adiabatic/fig2_20260922/README.md).
  The [amplitude campaign](runs/adiabatic/amplitude_evolution/README.md) has a
  merged 426-point derivative table and a saved cosine trajectory through
  $t/M = 4\times10^7$. Its approach to $A_0 \simeq 12.96268$ is a result on
  the saved spatial grid, subject to the high-amplitude caveat below.
- **Stationary background continuation:** ten primary-only solves at fixed
  $A_0=0.0011390625$ around $(\alpha,\chi)=(0.42,0.99)$ pass twenty observable
  comparisons. Offsets of $\pm0.001$ in alpha/spin converge in two or three
  Newton corrections; two routes to the joint downward offset agree in growth
  rate to $4.01\times10^{-10}$ relative. Halving the radial step changes growth
  by at most $1.91\times10^{-8}$ relative. See the
  [background-continuation campaign](runs/adiabatic/background_continuation/README.md).
  The same ten solves and twenty comparisons also pass at
  $A_0=2.051758324382913$, where the primary peak field is about 1.446 and
  $M\omega_I$ is 6.137 times its dilute value. Half-step changes in imaginary
  frequency and horizon flux are below $1.1\times10^{-8}$ relative; the small
  outward flux changes by $6.8\times10^{-5}$ relative. These are local primary
  stationary checkpoints, not a coupled BH trajectory; angular/endpoint checks
  remain separate work.
- **Secondary and driven background response:** the averaged secondary and
  $m=0,3$ fields now solve on all ten nonlinear-primary configurations above.
  The initial 133/135 passing checks include route agreement, secondary energy
  balance, and driven power identities. Two coarse-grid $m=3$ horizon-flux
  comparisons miss the $10^{-4}$ relative tolerance. Re-solving baseline and
  the joint offset at primary step 0.025 passes all 28 follow-up checks, with
  changes in that flux reduced from about $1.6\times10^{-4}$ to at most
  $5.10\times10^{-6}$ relative. Original failures are preserved. See the
  [response results and refinement](runs/adiabatic/background_continuation/README.md#secondary-and-driven-fields-on-the-finite-amplitude-backgrounds).
  This establishes local scalar response at fixed $A_0$, not BH–cloud evolution.
- **Coupled-balance prototype:** the leading small-cloud-mass scalar-only
  balance is derived in [the balance note](docs/coupled-bh-cloud-balance.md),
  with BH mass/angular momentum and both cloud energies in fixed initial-BH
  units. The campaign's `check_coupled_balance.jl` passes 44 checks against
  the existing nonlinear RHS, dilute channel counting, unit conversion, and
  a changing-background occupation-conservation example. The latter exposes
  the metric-work term required beyond the leading approximation. Actual
  nonlinear metric work and cloud angular-momentum background derivatives
  remain unevaluated; that checkpoint checked instantaneous rates only.
  The separate QCD campaign above now integrates the leading law. See the
  [checkpoint recipe and illustrative RHS](runs/adiabatic/background_continuation/README.md#coupled-balance-checkpoint--october-8-2026).
- **Dolan time-domain work:** the coupled-tower solver, eigenvalue/Leaver
  comparison, and frequency-filter analysis have a
  [linear validation record](docs/dolan-repro.md) and
  [Fig. 13 campaign](runs/dolan/fig13/README.md). The nonlinear source and
  synchronized RK4 have [separate checks](docs/nonlinear-progress.md);
  those do not establish accurate collapse evolution.
- **Self-gravity benchmark:** a standalone Newtonian Poisson check on an
  analytic hydrogenic $211$ density recovers the monopole/quadrupole shift
  $\delta\omega/\mu=-(237/1280)\alpha^2 M_c/M_{\rm BH}$. All 16 checks pass,
  including radial refinement, angular/outer-radius checks, and an independent
  exact radial integral. See the [derivation and command](docs/adiabatic-self-gravity-review.md#211-benchmark-checkpoint--october-7-2026).
  Saved nonlinear profiles and scalar evolution still have no gravity correction.
- **Dilute self-gravity screening:** published hydrogenic self/cross shifts
  have been evaluated on both saved two-cloud energy trajectories. At
  $f_a=10^{16}$ GeV, peak $|\delta\omega_1|/\mu$ is $7.70\times10^{-6}$
  (alpha 0.42) and $5.80\times10^{-6}$ (0.30); the outward scalar threshold
  margin changes by at most 0.0219% and 0.0483%. At $10^{17}$ GeV these
  changes become 2.19% and 4.83%. See the
  [campaign recipe and limits](runs/adiabatic/multimode_secondary/README.md#dilute-self-gravity-estimates--october-7-2026).
  This is an analytic-coefficient estimate, not the nonlinear frequency
  response or a bound on scalar flux errors. BH evolution remains separate.
- **Extensions:** beyond the preliminary coupled trajectory above, GW losses, self-gravity feedback,
  mutual nonlinear cloud response, additional independent clouds, and a PDE
  handoff remain proposals. The [October 6 research note](docs/further_work_1006.md)
  discusses their ordering and scientific motivation.

## Implementation decisions

The [multimode design](docs/multimode-adiabatic-design.md#decisions-and-their-reasons)
records the averaged reading of Eq. (61), restored $\mu^2$, secondary domain
and normalization, and the energy variables used for the balance equations.
The campaign records the $M\mu = 0.3$ tolerance setting, phase test, resonance
scan, and convergence evidence. Use its averaged secondary results;
`cosine_secondary_fixed_time.*` is retained only as a superseded experiment.

## Known caveats

- The Mμ = 0.42 primary branch is **not step-converged above A₁ ≈ 9**: its
  step-0.1 solutions re-evaluated at step 0.05 have scaled residual 4e-5 at
  A₁ = 9.1 rising to 4e-2 at 13. This bears on the cosine flux zero near
  A₁ ≈ 12.96. It does not affect Fig. 6 at Mμ = 0.42, which stays below
  E₁ ≈ 1000 (A₁ ≲ 0.8, step-converged to ~1e-7).
- Unchecked: secondary l-truncation (l ≤ 6) and outer boundary at large A₁.
  The paper reports neither.
- Full high-amplitude angular-quadrature and mode-truncation convergence of
  the primary remains unestablished. Interpolation and time-step checks on
  saved tables do not supply those spatial checks.
- Pieces 4–5: the outer-radius convergence of the $m = 0$ field is unchecked, and
  so is each coefficient against an independent evaluation of the perturbative
  $F_0$ and $F_3$ (Eqs. (38)–(39), (A11), (A14), e.g. by a Green's function).
  The power laws match; the prefactors have not been compared.
- The physical effect of the Sec. III approximations, including neglect of
  secondary self-interaction and its deformation of the primary, is unchecked.
- The saved two-cloud trajectories include scalar fluxes only. GW power,
  torque, and radiation reaction have not been implemented; physical energy
  normalization alone does not add them.

## Next steps

- **Active multiphase work:** the response, neighbor checks, and short cosine
  evolution are complete. Next assess radial/domain/support convergence and
  the physical slow-growth closure before interpreting a long trajectory.
  The campaign preserves both time resolutions and all accepted endpoint data.
- **Coupled-BH preparation:** local primary continuation in alpha and spin
  now works at both dilute amplitude and $A_0\simeq2.052$; the secondary and
  driven channels are checked on the nonlinear backgrounds, with the targeted
  refinement described above. Use primary step 0.05 or finer for the tested
  small $m=3$ horizon-flux accuracy requirement.
  The leading balance is now derived and its algebra checked. Its functions
  now live in `src/adiabatic/primary_secondary/coupled_balance.jl`, with units and
  equations documented in the source. The user stopped the time-integrator
  draft to require small, explained implementation steps. Stationary-table
  job 26757461 was cancelled; `runs/adiabatic/coupled_local/` contains paused
  drafts, not a validated trajectory. Review the balance piece together before
  resuming that incremental teaching track. Later on October 8 the user
  explicitly authorized completing the separate QCD spin-down campaign directly;
  its commented driver and Python integrator now provide a first limited
  trajectory with refreshed stationary data, fixed reference units, and explicit
  seeds/$f_a$. That exception does not remove the usual incremental workflow.
  The fixed outer amplitude used in the earlier stationary check is not
  an evolution prescription. Beyond the leading approximation, evaluate
  metric work, background angular response, and consistent gravitational
  energy/boundary accounting as described in the balance note.
- **Perturbative coefficients.** The small-$E_1$ limits of the spline
  coefficients `m0_flux_coefficient` and `m3_flux_coefficient` are the paper's
  $F_0$ and $F_3$ (Eqs. (38)–(39)). An independent evaluation of those (for
  example by the Green's-function formula) would check the prefactors of
  pieces 4–5, which so far are checked only through power laws and the power
  identity. The same coefficients would also give the perturbative curves of
  Fig. 6.
- **$M\mu = 0.3$ details.** Our $E_1$ peak looks ~10% below the figure and the
  plateaus 5–10% above; this is at the limit of reading a log axis. Not
  investigated.
- **Single-cloud saturation:** refine the high-amplitude spatial solve before
  interpreting the saved cosine flux zero or saturation trajectory physically.
- **Proposed next research direction:** physical normalization of the saved
  two-cloud trajectories, with GW power/torque estimates and a frozen-profile
  self-gravity diagnostic, followed by a limited coupled-BH calculation
  including GW losses and a matched perturbative
  comparison. See the [October 6 proposal](docs/further_work_1006.md) for the
  rationale and conditional extensions; these are not completed capabilities.

## Working with the user

Follow [CLAUDE.md](CLAUDE.md) for incremental development, notation,
validation, and Slurm conventions. The
[piece 4 guide](docs/multimode-piece4-guide.md) is the model for explanations:
define each object, describe the check procedure, and state what it establishes.
Use LaTeX for math in Markdown documentation.

The September 24 handoff is in Git history. The
[September 28 branch-work record](docs/adiabatic-branch-work-2026-09-28.md)
preserves the following stage's experiments; its next-step list is historical.
