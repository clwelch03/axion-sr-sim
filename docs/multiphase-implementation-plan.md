# Multiphase adiabatic evolution: implementation and organization plan

Updated 2026-10-12. The model-based source reorganization is implemented on
`refactor/adiabatic-layout`; the multiphase physics remains a development plan. The action-flux closure below still
requires a derivation for our open Kerr problem. [HANDOFF.md](../HANDOFF.md)
remains the authority on current capabilities.

## Objective and initial scope

Extend the existing nonlinear-primary/test-secondary framework to several
bound clouds whose profiles can respond to their joint amplitudes. Preserve
the computational benefit of averaging rapid oscillations: the initial slow
evolution system should evolve amplitudes or actions, without tracking fast
phases in physical time.

The first nonlinear target is **two mutually deforming clouds**. The first
test requiring additional evolution equations is **three independent clouds**.
Begin on a fixed Kerr background. Coupling this extension to BH mass and spin
evolution comes after the fixed-background closure passes its checks; changing
the background introduces additional response/work terms, discussed in the
[coupled balance note](coupled-bh-cloud-balance.md).

The existing [multimode design](multimode-adiabatic-design.md) records the
Sec. III reproduction and the original extension sketch. The
[literature review](multimode-literature-review-2026-09-28.md) and
[October 6 proposal](further_work_1006.md) supply research context.

## Recommended repository organization

**Keep this repository, but organize `src/adiabatic/` around physical models
before adding the new implementation.** Use a feature branch, preferably in
a separate Git worktree. The intended source groups are `core/`,
`single_cloud/`, `primary_secondary/`, and `multiphase/`.

Before the October 12 refactor, the structure mixed two organizing principles: `stationary/`,
`observables/`, and `amplitude_evolution/` name numerical tasks, whereas
`multimode/` names a physical approximation. Understanding the
single-cloud model required traversing several sibling folders, and some
apparently general files assumed one fundamental frequency. The refactor
separates those model assumptions from shared mathematics.

The organizing rule should be: **choose the physical model first; find its
profiles, observables, and evolution together.** Keep only demonstrated shared
building blocks in `core/`. These repository choices serve different purposes
and can be used together:

| Choice | Recommendation | Reason |
|---|---|---|
| New repository | Do not create one at this stage | The new method shares geometry, angular bases, potentials, radial operators, boundary series, Newton solves, dependencies, and reference results with this project. A separate repository would duplicate those components or force premature packaging. |
| Feature branch | First use `refactor/adiabatic-layout`, then `feat/multiphase-foundations` from the validated layout | Separates changes to file ownership and imports from new physics. Merge validated increments instead of maintaining one indefinitely diverging research branch. |
| Git worktree | Prefer a separate worktree while the current solver is also being developed or used | A worktree is another checkout of the same repository, on its own branch. It avoids switching the shared checkout underneath other work or run scripts. It is not a separate project or fork. |
| Source folders | Use `core/`, `single_cloud/`, `primary_secondary/`, and `multiphase/` | Gives each approximation a readable home. `primary_secondary/` states the asymmetric approximation more precisely than `multimode/`; two mutually deforming clouds belong in `multiphase/`. |
| New campaign folder | Create `runs/adiabatic/multiphase_two_cloud/` when the first solves are ready | Keeps new seeds, results, diagnostics, and provenance distinct from the established secondary campaign. |

The layout refactor uses branch `refactor/adiabatic-layout` in the existing
checkout; a second worktree was not created. Future implementation should
start from a recorded baseline revision
containing the required shared helpers. Existing uncommitted changes in another
checkout do not automatically appear in a new worktree; incorporate only the
changes actually needed, through explicit reviewable commits. Do not stash,
discard, or move unrelated work merely to create the new checkout.

A worktree isolates checked-out source files, but it does not automatically
isolate external output paths. Give each campaign a distinct checkpoint prefix.
Record the source revision and any uncommitted source changes for runs; use a
stable source checkout or snapshot for long jobs. Follow the existing
[campaign preservation rules](../runs/README.md).

### Current layout and future additions

The [source guide](../src/adiabatic/README.md) owns the current file map and
commands. Source and tests now follow the same model groups:

```text
src/adiabatic/
  README.md
  load.jl                # implemented models, preserving the public entry point
  core/                  # reusable geometry, operators, boundaries, numerics
  single_cloud/          # profiles, observables, continuation, evolution, drivers
  primary_secondary/     # existing asymmetric approximation and its drivers
  multiphase/            # scope README and core-only loader; physics still planned
tests/adiabatic/
  runtests.jl
  core/
  single_cloud/
  primary_secondary/
```

For the new physics, add `components.jl`, `projection.jl`, `profiles.jl`,
`actions.jl`, `fluxes.jl`, and `evolution.jl` under `multiphase/` as the milestones
are implemented, with matching tests and model-specific drivers. Create
`runs/adiabatic/multiphase_two_cloud/` when physical runs are ready, using the
existing campaign structure. No new physics stubs or campaign outputs are
needed merely to reserve the model's home.

Reuse the existing Julia project and dependencies. `single_cloud/` depends on
`core/`; `primary_secondary/` depends on `core/` and the single-cloud machinery
that constructs its primary. The multiphase library should depend on `core/`;
its seed-generation and comparison drivers may explicitly use the other
models. No existing model should depend on experimental multiphase code.

Give each group a loader with an explicit dependency order and a strategy that
avoids repeated includes. Preserve the root loader's existing public behavior
during migration; keep the new multiphase loader opt-in until validated. A
directory is an ownership boundary, not automatically a new Julia package or
module. Namespace/package changes are a separate decision and need not accompany
this reorganization.

### Migration map from the previous layout

The refactor extracted shared functions from mixed files. Moving all of
`stationary/` into `core/` would have incorrectly classified the one-phase
ansatz as general machinery. This table records the old ownership and its
replacement; exact current filenames are in the source guide.

| Previous location | Current ownership |
|---|---|
| `angular.jl`, `tortoise.jl`, `observables/physical_normalization.jl` | Shared `core/` basis, geometry, and units |
| Potential functions at the start of `stationary/nonlinear.jl` | `core/potentials.jl` |
| Field reconstruction, projection, and harmonic RHS in `stationary/nonlinear.jl` | `single_cloud/projection.jl` and profile assembly |
| Explicit-frequency coefficients in `stationary/radial.jl`, boundary series and per-block boundary data | `core/` |
| Assembly that assigns frequency from `m / fundamental_m` in radial and boundary routines | `single_cloud/` |
| Generic RK integration and Newton routines in `stationary/integration.jl` and `stationary/matching.jl` | `core/integration.jl` and `core/newton.jl` |
| Matching residuals, mode solves, and continuation predictors tied to the single-cloud ansatz | `single_cloud/profiles.jl` and `continuation.jl` |
| `observables/energy.jl` | `single_cloud/observables.jl`; extract shared low-level currents only when another model needs them |
| Generic residual/observable differentiation in `amplitude_evolution/branch_sensitivity.jl` | `core/sensitivity.jl` |
| Remaining `amplitude_evolution/` interpolation, fitting, and amplitude evolution | `single_cloud/` |
| Existing `multimode/` library and current `qcd_spindown.py` | `primary_secondary/` |
| Existing global `drivers/` | Each model's `drivers/`; reusable solve/setup functions currently imported from drivers belong in that model's library |
| `plot_profiles.py` and campaign-specific plotting | Single-cloud tooling or the owning campaign, according to actual reuse |

The core must not select which clouds are populated, prescribe the
primary/secondary approximation, or contain special cases for a particular
evolution campaign. Keep generalization driven by real callers rather than
creating a large framework in advance. Small files with a clear purpose remain
useful; the aim is coherent ownership, not merging everything into giant files.

### Layout migration validation

1. Record a baseline and run the existing cheap checks. Inventory Julia include
   paths, Python imports, campaign scripts, and documented commands that depend
   on the old locations.
2. Extract the genuinely shared routines, then move each model's code in
   reviewable changes. Preserve function signatures, equations, normalization,
   serialized data layouts, and numerical settings. Separate algorithm changes
   from this refactor.
3. Update loaders, tests, live callers, and documentation together. Use thin
   compatibility entry points where old commands or imports need to remain
   usable; compatibility files must forward to the single implementation. Check
   Python package discovery when placing code in new subpackages. Preserve
   historical provenance and frozen campaign records instead of rewriting them
   as if the new paths had always existed.
4. Check the existing unit suite and a cheap representative profile/flux/rate
   calculation against the recorded baseline. Validate driver entry points and
   Python imports as well as Julia loading; expensive convergence campaigns are
   unnecessary for a behavior-preserving move unless discrepancies appear.
5. Add a short root README and one per model, naming its physical approximation,
   inputs, outputs, dependencies, and entry point. Merge the layout change before
   adding new multiphase physics.

Campaign directories and expensive saved results can stay at their current
paths. A source reorganization does not require migrating numerical artifacts.

## Notation and what is averaged

The following is notation for the proposed implementation. Let $K$ be the
number of independent bound clouds and $a=1,\ldots,K$ label them. Define
$\mathbf A=(A_1,\ldots,A_K)$ as their slowly evolving amplitude coordinates;
their precise normalization is part of the profile-solver design. They are
not assumed to equal energies or occupations.

Let $r,\theta,\varphi$ be spatial coordinates, with $\theta$ the polar angle
and $\varphi$ the azimuth. The scalar field is $\phi$, a different symbol from
the azimuth. Define $\vartheta_a$ as the temporal oscillation phase of cloud
$a$, $m_a$ as its azimuthal integer, and

$$
\psi_a=\vartheta_a-m_a\varphi,
\qquad \dot\vartheta_a=\omega_a(\mathbf A),
\qquad
\phi=\Phi(r,\theta;\boldsymbol\psi;\mathbf A).
$$

Here $\omega_a$ is a real carrier frequency and $\Phi$ is periodic in each
phase. The phases are coordinates used to construct and average a nonlinear
configuration. Phase quadrature inside a profile solve does not require
integrating these rapid phases in physical time.

Label each Fourier component by an integer vector
$\mathbf q=(q_1,\ldots,q_K)$, so that

$$
\Omega_{\mathbf q}=\mathbf q\cdot\boldsymbol\omega,
\qquad m_{\mathbf q}=\mathbf q\cdot\mathbf m.
$$

A component has phase factor $e^{-i\mathbf q\cdot\boldsymbol\psi}$ and an
angular expansion with indices $l\geq |m_{\mathbf q}|$. Independent clouds,
Fourier components, and angular basis functions are distinct objects. Driven
components need not acquire independent populations when their response can
be determined by the occupied clouds.

## Implementation milestones

### 1. Derive the closure and establish its approximation order

For $K>2$, the total energy and angular-momentum balances supply only two
constraints on $K$ amplitude rates. The missing information must come from
mode-resolved dynamics.

The candidate construction defines phase actions
$I_a=\partial\overline L/\partial\omega_a$, where $\overline L$ is the
spatially integrated, phase-averaged conservative Lagrangian. The partial
derivative is taken with the profile variables held fixed before evaluation
on the consistent family. Their subsequent amplitude derivatives include
the profile response. The proposed evolution structure is

$$
\sum_b D_{ab}\dot A_b=S_a,
\qquad D_{ab}=\frac{\partial I_a}{\partial A_b},
\qquad
S_a=-\sum_{\mathbf q}q_a\frac{F_{E,\mathbf q}}{\Omega_{\mathbf q}}.
$$

Here $F_{E,\mathbf q}$ is the signed energy loss through the horizon and
outer boundary for a resolved channel. Positive means loss from the cloud;
a superradiant horizon contribution is negative. Count each real-field
conjugate pair once with a documented normalization. The quotient assumes
nonzero frequency and separately resolved channels; zero-frequency components
and coherent degeneracies need their own treatment.

Derive these equations using a multiple-timescale expansion or appropriate
solvability projections, retaining boundary terms. Specify how slow growth
enters the profile problem. Arbitrary-amplitude profiles cannot generally be
exactly periodic real-frequency solutions with unbalanced open-boundary fluxes.
Simply assigning complex frequencies to every combination does not establish
the required nonlinear modulation equations. Do not count growth twice through
both an imaginary frequency and the same horizon flux.

On a consistent conservative family at fixed BH parameters, test the expected
identities

$$
dE=\sum_a\omega_a\,dI_a,\qquad J=\sum_a m_a I_a,
$$

where $E$ and $J$ are the total cloud energy and angular momentum. These are
conditional identities to establish for the chosen approximation, not exact
additional conserved charges of the dissipative real scalar field. Do not
substitute $I_a=E_a/\omega_a$ at finite amplitude without a derivation.

**Checkpoint:** recover the single-cloud amplitude equation and dilute
two-cloud transfer equations. Compare to the current finite-amplitude model
only at matching approximation order. The shared representation and projector
can be developed while this derivation is underway; the general production
evolution must wait for a defensible closure.

### 2. Introduce explicit component labels

Represent each independent cloud by its linear-limit quantum numbers and
normalization convention. Represent each Fourier block by $\mathbf q$ and its
retained $l$ values; derive its frequency and azimuthal integer from the cloud
vectors.

For the existing two-cloud example:

| Component | $\mathbf q$ | Frequency | Azimuthal integer for $(m_1,m_2)=(1,2)$ |
|---|---|---|---|
| Primary fundamental | $(1,0)$ | $\omega_1$ | 1 |
| Secondary fundamental | $(0,1)$ | $\omega_2$ | 2 |
| Driven dissipative component | $(2,-1)$ | $2\omega_1-\omega_2$ | 0 |
| Driven radiative component | $(-1,2)$ | $2\omega_2-\omega_1$ | 3 |

The assembled stationary code currently derives frequency from
`frequency * m / fundamental_m`. Add explicit-frequency block assembly using
the lower-level helpers in [radial.jl](../src/adiabatic/core/radial_operator.jl)
and [boundaries.jl](../src/adiabatic/core/boundaries.jl). Define and check
the conventions for zero/signed $m$, negative-frequency conjugates, and
bound/outgoing boundary selection. Do not silently merge distinct components
merely because they have the same $m$.

**Checkpoint:** reproduce existing radial operators and boundary states with
the new representation, including the driven $m=0$ block.

### 3. Generalize nonlinear projection

Reconstruct the real field from the retained components on a polar/phase
quadrature grid, evaluate the existing cosine/quartic/sextic potential
functions, and project the nonlinear remainder onto each retained component.
The current fixed-time azimuthal projection cannot distinguish every
independent-frequency component. The new representation must preserve the
existing factor-of-two and angular normalization conventions.

Start with two phase dimensions and a small explicitly chosen component set.
Compare quartic projections against direct algebraic products. Then test the
cosine potential and refine phase quadrature independently of the Fourier
truncation and angular resolution. Measure sources in selected omitted
components to guide expansion of the retained set.

**Checkpoint:** recover the single-phase projection, exact quartic selection
rules and coefficients, and the known dilute driven sources. Verify that
changing phase origins transforms coefficients consistently without changing
averaged observables.

Tensor phase grids grow rapidly with $K$. Measure the prototype's cost before
choosing a more scalable projection method; do not promise arbitrary cloud
counts from a two-phase implementation.

### 4. Construct a jointly deformed two-cloud family

Seed the new calculation with a saved primary and a small secondary. Continue
the secondary amplitude upward, solving for both fundamental profiles,
frequencies, and retained nonlinear components using the formulation from
milestone 1. Fix one amplitude normalization and phase convention per
independent fundamental in the nonresonant formulation.

Reuse the Newton solver and boundary machinery where applicable. Start with
selected configurations and short continuation paths rather than a dense
rectangular amplitude table. Include higher components according to an
explicit truncation whose effects can be measured.

**Checkpoint:** recover the primary-only limit and the corresponding
small-secondary truncation, then demonstrate a converged change in the primary
caused by the secondary. Compute total $E$ and $J$ from the joint field.
Agreement of Newton residuals alone does not establish radial, angular,
boundary, or Fourier convergence.

### 5. Compute actions, response derivatives, and slow evolution

Extend [branch_sensitivity.jl](../src/adiabatic/core/sensitivity.jl)
from a scalar amplitude derivative to multiple amplitude directions. If
$G(u,\mathbf A)=0$ is the discretized profile residual, with $u$ the vector of
profile unknowns, differentiate it as

$$
G_{,u}\,u_{,A_b}=-G_{,A_b}.
$$

Reuse one Jacobian factorization for all amplitude directions. Propagate the
response into the actions and other observables to obtain $D_{ab}$ and total
derivatives. Cross-check selected derivatives with finite differences between
independently corrected profile solutions.

Evaluate channel-resolved boundary fluxes and integrate the resulting slow
amplitude/action equations. Keep the profile solver separate from the time
integrator: the integrator should consume corrected profile data and its
derivatives, whether computed locally or later supplied by validated tables.
Begin with a short trajectory; defer a general high-dimensional interpolator
until sampling needs and errors are understood.

**Checkpoint:** on the same two-cloud family, compare action-derived rates
with an independently assembled $E,J$ balance solve. Check the differential
identities, each channel's angular/energy flux relation, and integrated total
balances. Do not obtain the actions by fitting the two balance equations and
then call their agreement independent validation. Track response-matrix
conditioning rather than hiding singularity with a pseudoinverse.

### 6. Add a third cloud and test truncation and averaging

First recover an established dilute three-cloud transfer system. Then allow
the third cloud to deform the joint configuration. Increase the number of
independent clouds separately from the number of Fourier and angular
components. Assess convergence of the target observable, not just populations;
a lightly populated component can still mediate substantial transfer.

Monitor relevant allowed resonances. With an integer interaction vector
$\mathbf k$, define the relative phase
$\chi=\mathbf k\cdot\boldsymbol\vartheta$ and detuning
$\dot\chi=\mathbf k\cdot\boldsymbol\omega$. An axisymmetric interaction
requires $\mathbf k\cdot\mathbf m=0$. Compare relevant detunings with coupling
and evolution rates, not just with the carrier frequencies. If a relative
phase remains slow enough to affect transfer, flag the amplitudes-only
approximation; adding its dynamics is a further extension. Likewise flag
emission thresholds requiring new boundary expansions, driven responses near
bound resonances, and loss of the timescale separation or branch stability
needed for adiabatic following.

**Checkpoint:** recover dilute individual transfer rates and total balances,
then demonstrate convergence of a short finite-amplitude three-cloud evolution.
Total $E,J$ balance alone cannot validate the extra amplitude direction; use
independent transfer benchmarks and, when feasible, selected short field
evolutions to test it.

## First development increment and completion criteria

The **behavior-preserving layout migration** above is complete; its validation
is recorded in the project handoff. The first new-physics increment is **component representation plus the two-phase
nonlinear projector**, alongside the written closure derivation. Its useful
output is tested infrastructure and explicit conventions, without committing
the joint solver to an unverified evolution law. Follow the incremental
development conventions in [CLAUDE.md](../CLAUDE.md).

The first physics result is a converged measurement of how a populated
secondary reshapes the primary. The first result establishing the new closure
is validated three-cloud evolution. Neither requires beginning with all modes
or a large amplitude grid.

Keep derivations and final conventions in design documentation, checkpoint
status in `HANDOFF.md`, and run commands, resolution studies, source revisions,
inputs, and saved evidence in campaign READMEs. Submit expensive solves through
Slurm under the existing project conventions. Merge each validated foundation
in a focused change, keeping experimental drivers opt-in until their physical
checks are complete.
