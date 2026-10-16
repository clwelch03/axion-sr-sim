# Multiphase adiabatic evolution: implementation plan

Updated 2026-10-15. The joint matching solver and finite-domain E/J observables
are implemented. The opt-in loader now includes observables, implicit
amplitude responses, and short E/J-balanced evolution. The user authorized
direct implementation of milestones 5's response, neighbor checks, and short
evolution, with explanatory comments and docstring skeletons.

The response uses a fixed reference coefficient scale and real-packed trial
coordinates. It factors one matching Jacobian for both amplitude directions,
then differentiates the observables along the induced profile changes.
Analytic vector-response and midpoint-evolution checks pass. The frozen-source
matched-cloud campaign completed both neighbor validations and 4/8-step cosine
trajectories. E/J derivatives agree with neighbors within 9.50e-9 relative;
charge-budget errors decrease fourfold when dt is halved. See
[the campaign record](../runs/adiabatic/multiphase_two_cloud/README.md).

The physical slow-growth approximation, profile/domain/support convergence,
and the later action-flux closure remain separate open questions.
[HANDOFF.md](../HANDOFF.md) remains the authority on current capabilities.

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

## Working approach

Work on `main` in the existing checkout, one small piece at a time.
The default teaching workflow uses small, already-run snippets. For the
October 15 response/evolution increment, the user explicitly authorized direct
library implementation, with readable comments and docstring skeletons.
Each piece still gets a decisive check, as in [CLAUDE.md](../CLAUDE.md).

Implement the new model under `src/adiabatic/multiphase/`, with matching tests
under `tests/adiabatic/multiphase/`. Add `components.jl`, `projection.jl`,
`profiles.jl`, `observables.jl`, and later action/evolution helpers as their
milestones are reached. Joint E/J and boundary fluxes currently share
`observables.jl`. Reuse `core/` and the existing Julia project; seed-generation and
comparison drivers may explicitly load the other models. Keep the multiphase
loader opt-in and existing models independent of it. The
[source guide](../src/adiabatic/README.md) owns the current file map and commands.

Create `runs/adiabatic/multiphase_two_cloud/` when physical runs are ready.
Follow the [campaign preservation rules](../runs/README.md), with distinct
checkpoint prefixes and recorded source revisions and uncommitted changes.
Use a stable source snapshot for long jobs.

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

### Independent clouds, driven components, and interaction channels

The length of every $\mathbf q$ is **the number $K$ of independent cloud
phases**, not the number of retained Fourier components or angular functions.
For two independent clouds, all four components in the table below have
two-entry vectors. Retaining another driven component adds another vector,
not another entry to each vector. Adding a third independent cloud adds a
third entry to every vector; existing components acquire a zero in that entry.

For example, adding a $411$ cloud with its own amplitude and phase to $211$
and $322$ makes $K=3$. Here the cloud labels denote the linear-limit quantum
numbers $(n,l,m)$, with $n$ the principal quantum number; $K$ avoids confusing
that quantum number with the number of clouds. Clouds can share $l,m$ and
still have independent phases and radial profiles.

A driven component is a sourced field response, not the interaction process
itself. For the two-cloud ordering $(211,322)$, a cubic source contains

$$
\bigl(R_1e^{-i\psi_1}\bigr)^2
\bigl(R_2^*e^{+i\psi_2}\bigr)
=R_1^2R_2^*e^{-i(2\psi_1-\psi_2)},
$$

where $R_a$ denotes the complex spatial coefficient of cloud $a$. This product
has $\mathbf q=(2,-1)$, frequency $2\omega_{211}-\omega_{322}$, and $m=0$.
It supplies the driven field associated with the familiar transfer channel
$211+211\rightarrow322+\mathrm{BH}$; BH denotes horizon absorption.
Its profile must be solved, but its phase is fixed by the two cloud phases.

If the $411$ population is also evolved independently, the proposed channel
$211+322\rightarrow411+\mathrm{BH}$ instead has a BH-directed component
with $\mathbf q=(1,1,-1)$ in cloud order $(211,322,411)$:

$$
\Omega_{\mathbf q}=\omega_{211}+\omega_{322}-\omega_{411},
\qquad m_{\mathbf q}=1+2-1=2.
$$

This is phase and azimuthal bookkeeping, not a calculation of that channel's
rate or proof that it matters. An interaction-generated bound population can
still require an independent cloud variable; its origin does not make it a
slaved response. Conversely, “driven” does not require propagation to infinity:
a below-threshold component can carry horizon flux. Near bound resonances,
whether a response can remain slaved needs the checks in milestone 6.

## Implementation milestones

### 1. Derive the closure and establish its approximation order

**Two-cloud sequencing decision (October 15):** retain total energy and
angular-momentum balance for the first mutually deformed two-cloud evolution
at fixed BH parameters. The two positive outer normalizations $(A_1,A_2)$
label the profile family; total phase-averaged $E$ and $J$ are computed from
the joint field, including interaction terms at the chosen approximation order.
They are outputs, not two additional independent evolution coordinates.
Use the full profile response in

$$
\begin{pmatrix}
\partial_{A_1}E & \partial_{A_2}E\\
\partial_{A_1}J & \partial_{A_2}J
\end{pmatrix}
\begin{pmatrix}\dot A_1\\\dot A_2\end{pmatrix}
=-\begin{pmatrix}F_E\\F_J\end{pmatrix}.
$$

This determines the rates only where the matrix is nonsingular and the
two-amplitude adiabatic description is valid. Establish how slow growth enters
the profile approximation before interpreting a matched family physically;
switching to balance equations does not resolve that question. Phase actions
remain a proposed later extension, to be derived before three-cloud evolution
and checked independently against the two-cloud balance result.

Keep the profile solver in `multiphase/profiles.jl`, with reconstruction and
projection in `projection.jl` and labels in `components.jl`. Separate joint
observables, profile responses, and the slow-time driver as those pieces are
implemented. Begin with nearby corrected profiles and a short continuation
path, rather than a dense two-dimensional table. Continuation can acquire its
own `multiphase/continuation.jl` when needed; campaign scripts and evidence
belong under `runs/adiabatic/multiphase_two_cloud/`.

**Observables interface (October 15):** reconstruct the joint real field and
its four Boyer-Lindquist coordinate derivatives, evaluate the scalar stress
tensor, then obtain energy, angular momentum, and boundary currents by
contraction and averaging. Separate analytic energy terms remain useful
independent checks. The first derivative prototype uses explicitly supplied
real carrier frequencies: its time derivative holds the radial profiles fixed
and omits slow amplitude/profile response. It rejects complex trial frequencies
rather than silently interpreting their imaginary parts. This is a frozen-field
kinematic evaluator, not a resolution of the slow-growth approximation.
Its 32 checks in `test_field_derivatives.jl` pass against the inserted helpers
in `projection.jl`, including internal polar-derivative construction from
`polar_values` and `grid`.
The potential-energy helper is already
inserted into `observables.jl` and passes its 10 checks when explicitly loaded;
the observable and derivative tests are now registered.

The next two candidates are `kerr_metric` and `scalar_stress_tensor`. Their
115 checks are in `core/test_kerr_metric.jl` (86) and
`multiphase/test_stress_tensor.jl` (29). Both helpers are inserted, and the
115 checks pass against the actual library.
Use signature $(-+++)$, coordinate order $(t,r,\theta,\varphi)$, and the
normalized tensor $T_{\mu\nu}/F_a^2$. The metric helper accepts a tortoise
position but returns Boyer-Lindquist components, matching the reconstructed
radial derivative. Raising the first index supplies $-T^t{}_t$ and
$T^t{}_\varphi$ for charge densities and $-T^r{}_t$ and $T^r{}_\varphi$ for
outward currents. Include $\rho^2$, the $2\pi$ spatial azimuthal integral,
and the independent-phase average in their angular integrals; multiply by
$dr/dr_*$ only for charges per unit tortoise radius. These local contractions
recover the existing single-cloud energy with real carriers to $3.1\times10^{-15}$
relative error for all three potentials, analytic joint channel fluxes, and
Schwarzschild angular momentum. The integrated evaluator is now implemented
as described below; physical endpoint/truncation convergence remains pending.

The inserted `two_phase_observables_at_radius` wraps these
contractions in one polar/phase quadrature loop and passed 44 checks in
`test_observables_at_radius.jl` against the actual library, with `observables.jl`
loaded explicitly.
Its four outputs are energy and angular momentum per unit tortoise radius,
`radial_energy_flux`, and `radial_angular_momentum_flux`. Both currents are
signed toward increasing radius; negative values describe inward transport.
The complete wrapper
recovers single-cloud energy to $2.0\times10^{-15}$ relative error and the
analytic joint-field currents, Schwarzschild angular momentum, and quartic
interaction energy. It also passes phase-origin and local differentiation
checks; the latter scales prescribed fields and does not compute the response
of a matched profile to either cloud amplitude. Loader inclusion and test
registration are complete; the test now uses the current `radial_*` field
names and passes against the library.

The inserted radial wrapper `two_phase_integrated_observables` reuses
`integrate_uniform` for both charges, permits either storage direction on each
uniform branch, checks that the branches meet, and negates the inner outward
currents to report horizon losses. All 49 checks in
`test_integrated_observables.jl` pass: exact cubic E/J integrals, odd/even
interval counts, reverse storage, single-cloud integrated energies for all
potentials, endpoint signs, component ordering, and saved-history/input checks.
The wrapper rejects nonfinite samples and nonfinite local observables.
The user-requested refactor extracting top-level
`ordered_radial_branch` and `integrate_two_phase_branch_observables` passes
the same 49 checks and is now inserted, together with explanatory comments
and a short radial-branch docstring. A branch here is a sequence of spatial
profile samples from one radial integration, not an amplitude family or time
series. It preserves the local
`radial_energy_flux` and `radial_angular_momentum_flux` names, positive toward
increasing radius. The tests are now registered. These are finite-domain
integrals; no boundary tail correction or physical-profile convergence is
established. The saved-history compatibility trial is deliberately dilute
after larger arbitrary boundary coefficients caused nonlinear growth and
observable overflow on its short test domain. The amplitude-response solve
is now implemented separately in `response.jl`, as described in milestone 5.

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

### 2. Component labels, radial equations, and boundary conditions

Represent each independent cloud by its linear-limit quantum numbers and
normalization convention. Represent each Fourier block by $\mathbf q$ and its
retained $l$ values; derive its frequency and azimuthal integer from the cloud
vectors. Start with a small representation containing `q` and `l_values`.
Keep the cloud frequency vector separate from these labels because frequencies
change along the nonlinear family. Require `length(q) == K` and distinct
retained integer $l\geq|m_{\mathbf q}|$; group blocks by `q`, not by $m$.

For the existing two-cloud example:

| Component | $\mathbf q$ | Frequency | Azimuthal integer for $(m_1,m_2)=(1,2)$ |
|---|---|---|---|
| Primary fundamental | $(1,0)$ | $\omega_1$ | 1 |
| Secondary fundamental | $(0,1)$ | $\omega_2$ | 2 |
| Driven dissipative component | $(2,-1)$ | $2\omega_1-\omega_2$ | 0 |
| Driven radiative component | $(-1,2)$ | $2\omega_2-\omega_1$ | 3 |

The single-cloud assembly (`linear_tortoise_rhs!` and `boundary_state` in
[single_cloud/profiles.jl](../src/adiabatic/single_cloud/profiles.jl)) derives
frequency from the harmonic multiple `m / fundamental_m`. The shared
`radial_operator_coefficients` and `tortoise_operator_coefficients` in
[core/radial_operator.jl](../src/adiabatic/core/radial_operator.jl), and
`horizon_boundary_data` and `outer_boundary_data` in
[core/boundaries.jl](../src/adiabatic/core/boundaries.jl), already accept the
actual component frequency. Reuse them in new multiphase block assembly;
the existing single-cloud assembly keeps its physical assumptions. Define and check
the conventions for zero/signed $m$, negative-frequency conjugates, and
bound/outgoing boundary selection. Do not silently merge distinct components
merely because they have the same $m$.

**Checkpoint:** reproduce existing radial operators and boundary states with
the new representation, including the driven $m=0$ block.

#### 2.1 The field represented by a component

Use $x=\cos\theta$. Let $Y_{lm}(x)$ denote the real polar function normalized by

$$
\int_{-1}^{1}Y_{lm}(x)Y_{l'm}(x)\,dx=\delta_{ll'}.
$$

For signed $m$, choose $Y_{lm}=Y_{l,|m|}$ from the existing polar table. This
fixes our basis convention; any alternative constant phase for negative $m$
would have to be absorbed into the radial coefficients. There is no extra
$1/\sqrt{2\pi}$ factor in these functions.

For each retained phase vector, define its spatial coefficient

$$
F_{\mathbf q}(r,x)=\sum_{l\in L_{\mathbf q}}
R_{l,\mathbf q}(r)Y_{l,m_{\mathbf q}}(x),
$$

where $L_{\mathbf q}$ is the chosen list of angular indices. The contribution
to the real scalar field is
$F_{\mathbf q}e^{-i\mathbf q\cdot\boldsymbol\psi}+\text{complex conjugate}$.
Each $R_{l,\mathbf q}$ is an unrescaled radial profile: no factor of $1/r$ is
included. Its value at one radius is distinct from a slow cloud coordinate
$A_a$, and from a leading boundary-series coefficient.

#### 2.2 The radial equation and the assembled matrix

Work in units $M=1$. Write $a_{\rm BH}$ for the dimensionless Kerr spin (the
code's `spin`) and $\alpha=M\mu$ for the scalar mass parameter. Define

$$
\Delta=r^2-2r+a_{\rm BH}^2,\qquad H=r^2+a_{\rm BH}^2,
\qquad \frac{dr_*}{dr}=\frac{H}{\Delta}.
$$

For a fixed component, abbreviate $\Omega=\Omega_{\mathbf q}$ and
$m=m_{\mathbf q}$. The spatial radial equation at frozen amplitudes is

$$
\frac{d}{dr}\left(\Delta\frac{dR_{l,\mathbf q}}{dr}\right)
+B_{l,\mathbf q}R_{l,\mathbf q}
+a_{\rm BH}^2(\Omega^2-\alpha^2)
\sum_{l'}C^m_{ll'}R_{l',\mathbf q}
+N_{l,\mathbf q}=0,
$$

$$
B_{l,\mathbf q}=
\frac{[H\Omega-a_{\rm BH}m]^2}{\Delta}
+2a_{\rm BH}m\Omega-a_{\rm BH}^2\Omega^2
-\alpha^2r^2-l(l+1),
\qquad
C^m_{ll'}=\int_{-1}^{1}Y_{lm}(x)x^2Y_{l'm}(x)\,dx.
$$

$N_{l,\mathbf q}$ is the nonlinear source defined in milestone 3. The
linear terms mix angular indices within a fixed $\mathbf q$; the nonlinear
source couples different $\mathbf q$ values. Assembling this spatial equation
does not yet supply the slow-growth terms required by milestone 1.

Collect the retained radial profiles into $\mathbf R_{\mathbf q}$ and define
$\mathbf Q_{\mathbf q}=d\mathbf R_{\mathbf q}/dr_*$. Then

$$
\frac{d\mathbf R_{\mathbf q}}{dr_*}=\mathbf Q_{\mathbf q},\qquad
\frac{d\mathbf Q_{\mathbf q}}{dr_*}
=-p\mathbf Q_{\mathbf q}-V_{\mathbf q}\mathbf R_{\mathbf q}
-s\mathbf N_{\mathbf q},
$$

$$
s=\frac{\Delta}{H^2},\qquad p=2rs,\qquad
(V_{\mathbf q})_{ll'}=
s\left[B_{l,\mathbf q}\delta_{ll'}+
a_{\rm BH}^2(\Omega^2-\alpha^2)C^m_{ll'}\right].
$$

An **operator** here is the rule taking profiles and their derivatives to
the left-hand side of the radial equation. The implemented linear assembly
returns $p$ and $V_{\mathbf q}$. Near the horizon, evaluate these through the
existing regularized core helper, which cancels the apparent $1/\Delta$
singularity analytically. Do not numerically form divergent $B$ and then
multiply by a vanishing $s$.

#### 2.3 Boundary series and their coefficients

Let $r_\pm=1\pm\sqrt{1-a_{\rm BH}^2}$. At the horizon, use the ingoing series

$$
\mathbf R_{\mathbf q}(r)\simeq
(r-r_+)^{\nu_H}\sum_{j=0}^{J_H}\mathbf b_{H,j}(r-r_+)^j,
\qquad
\nu_H=-i\frac{2r_+\Omega-a_{\rm BH}m}{r_+-r_-}.
$$

At large radius, use

$$
\mathbf R_{\mathbf q}(r)\simeq e^{ikr}r^{\nu_\infty}
\sum_{j=0}^{J_\infty}\frac{\mathbf b_{\infty,j}}{r^j},
\qquad
\nu_\infty=-1+i\frac{2\Omega^2-\alpha^2}{k}.
$$

Here $J_H,J_\infty$ are series truncation orders. For real positive $\Omega$,
choose $k=i\sqrt{\alpha^2-\Omega^2}$ below the mass threshold (decaying),
or $k=\sqrt{\Omega^2-\alpha^2}$ above it (outgoing). The existing core helper
defines the continuation for complex frequency. Exact $k=0$ requires a
different expansion.

The input `leading_coefficients` is $\mathbf b_{H,0}$ or
$\mathbf b_{\infty,0}$, with one entry per retained $l$. Higher coefficients
follow from the radial equation. Evaluate both the series and its tortoise
derivative to obtain the boundary state $[\mathbf R;\mathbf Q]$.
The current wrapper supports signed/zero $m$ and finite
$\operatorname{Re}\Omega>0$; static and negative-frequency representatives
remain outside its supported boundary convention.

**Checks already used in prototypes:** compare linear and boundary assembly
to the single-cloud implementation; compare $m=0$ to analytic expressions;
permute angular indices and require the same permutation of states and both
matrix axes. Preserve separate blocks when two components share $m$.

### 3. Generalize nonlinear projection

**Purpose:** compute $N_{l,\mathbf q}$ from the joint real field. Projection
means extracting the coefficient of a chosen phase and angular pattern.
The real-field reconstruction convention, angular basis, and potential
functions already have tests in the existing model. Reuse them. The new
mathematics is the independent-phase average, which distinguishes components
that a fixed-time azimuthal average combines.

Start at one radius, with prescribed complex radial coefficients. No profile
solve, boundary matching, or time evolution is needed to check this operation.

#### 3.1 Retained components and real-field reconstruction

Let $\mathcal Q$ be a finite set of nonzero integer phase vectors, containing
at most one member of each conjugate pair $\{\mathbf q,-\mathbf q\}$.
Together with the angular lists $L_{\mathbf q}$, it defines the truncation.
At a fixed radius, reconstruct

$$
\Phi(x,\boldsymbol\psi)
=\sum_{\mathbf q\in\mathcal Q}
\left[F_{\mathbf q}(x)e^{-i\mathbf q\cdot\boldsymbol\psi}
+F_{\mathbf q}(x)^*e^{+i\mathbf q\cdot\boldsymbol\psi}\right]
=2\operatorname{Re}\sum_{\mathbf q\in\mathcal Q}
F_{\mathbf q}(x)e^{-i\mathbf q\cdot\boldsymbol\psi}.
$$

The factor of two supplies the conjugate partner; it must not appear a second
time in the projector. In the signed polar convention of section 2.1,
$F_{-\mathbf q}=F_{\mathbf q}^*$. A zero vector would be its own conjugate
and needs a separate real coefficient; exclude it from this first prototype.
A component with $m=0$ and nonzero $\mathbf q$, such as $(2,-1)$, uses the
same formula and can oscillate in time: $m=0$ does not imply $\Omega=0$.

The projector itself does not need numerical carrier frequencies. It uses
the integer phase labels, the cloud azimuthal integers, and the radial
coefficients. When using these components in the boundary problem, choose
representatives compatible with its positive-frequency convention and detect
any approach to zero frequency separately.

#### 3.2 The continuous projector and its normalization

Define the normalized average over independent phases by

$$
\langle G\rangle_\psi=
\frac{1}{(2\pi)^K}\int_0^{2\pi}\!d\psi_1\cdots
\int_0^{2\pi}\!d\psi_K\,G(\boldsymbol\psi).
$$

Integer Fourier patterns obey

$$
\left\langle e^{i(\mathbf q-\mathbf p)\cdot\boldsymbol\psi}
\right\rangle_\psi=\delta_{\mathbf q,\mathbf p}.
$$

For a real polar/phase field $G$, define

$$
\mathcal P_{l,\mathbf q}[G]=
\int_{-1}^{1}dx\,Y_{l,m_{\mathbf q}}(x)
\left\langle e^{+i\mathbf q\cdot\boldsymbol\psi}
G(x,\boldsymbol\psi)\right\rangle_\psi.
$$

Combining phase orthogonality with polar orthonormality gives
$\mathcal P_{l,\mathbf q}[\Phi]=R_{l,\mathbf q}$ for retained entries, and
zero for unpopulated entries. This fixes the sign and normalization of the
new projector. There is no additional azimuthal integral or factor of $2\pi$
in this coefficient extraction.

This is an algebraic identity on the phase torus (the product of the phase
circles). Interpreting a torus average as the appropriate dynamical average
still requires the nonresonance and timescale checks in milestones 1 and 6.
Distinct labels must not be assumed dynamically incoherent merely because
this mathematical projection separates them.

#### 3.3 The first decisive test: equal angular pattern, different phases

Use cloud azimuthal integers $(m_1,m_2)=(1,2)$ and retain just

$$
\mathbf q_A=(3,0),\qquad \mathbf q_B=(-1,2),\qquad
m_{\mathbf q_A}=m_{\mathbf q_B}=3,\qquad l_A=l_B=3.
$$

Choose complex radial coefficients $u=0.2+0.1i$ and $v=-0.07+0.04i$.
Then

$$
\Phi=2Y_{33}(x)\operatorname{Re}
\left[u e^{-3i\psi_1}+v e^{-i(-\psi_1+2\psi_2)}\right].
$$

The required results are

$$
\mathcal P_{3,\mathbf q_A}[\Phi]=u,\qquad
\mathcal P_{3,\mathbf q_B}[\Phi]=v.
$$

If instead we set both temporal phase origins to zero and sample only the
physical azimuth, then $\psi_1=-\varphi$, $\psi_2=-2\varphi$. Both factors
become $e^{3i\varphi}$: azimuthal projection recovers $u+v$. This is precisely
the distinction the existing one-phase reconstruction cannot supply.

Use the existing polar tables and a modest, adequately resolved phase grid;
16 polar nodes and $9\times9$ phase nodes suffice for this particular linear
round trip. Require recovery to floating-point quadrature accuracy, including
the imaginary parts. These are bookkeeping test coefficients, not a solved
physical profile.

**Library reconstruction and projector checked (October 13):**
`tests/adiabatic/multiphase/test_phase_separation.jl` reconstructs the field
using `reconstruct_two_phase_field_at_radii` (one radius per call), and
evaluates the projection sums explicitly. All three separation checks pass:
the two coefficients are recovered
separately with maximum absolute error $3.11\times10^{-16}$, and the
azimuth-only projection recovers their sum. Five additional checks cover
unequal phase-grid sizes, the phase sign, an empty input field, and dimension
mismatches. Five further checks exercise `project_two_phase_field`: recovery
of both known coefficients, zero for an absent component, conjugate recovery
at negative $\mathbf q$, and projection of an independently constructed sine
field on unequal phase grids. Sections 3.1 and 3.2 now have tested library
implementations using the discrete quadrature in section 3.6.

#### 3.4 Evaluate the nonlinear remainder, then project it

Let $U$ be the dimensionless potential, excluding its overall $\alpha^2$.
Define $h(\Phi)=\Phi-U'(\Phi)$, so the existing potentials give

$$
h(\Phi)=
\begin{cases}
\Phi-\sin\Phi,&\text{cosine},\\
\Phi^3/6,&\text{quartic},\\
\Phi^3/6-\Phi^5/120,&\text{sextic}.
\end{cases}
$$

The source in the radial equation is

$$
\boxed{N_{l,\mathbf q}(r)=
\alpha^2\mathcal P_{l,\mathbf q}
\left[(r^2+a_{\rm BH}^2x^2)h(\Phi)\right].}
$$

Evaluate $h$ on the **joint** field; applying it to each cloud separately
loses the cross terms. Include the geometric factor inside the polar
projection: it can mix angular indices. This source enters the radial
equation with a plus sign, hence contributes $-sN_{l,\mathbf q}$ to
$dQ_{l,\mathbf q}/dr_*$. Apply that tortoise factor once, after projection.

The output component set may include selected labels absent from the input
field. Their projected sources measure which omitted components are driven.
This does not mean they have already acquired solved radial profiles.

**Source checkpoint (October 14):** `nonlinear_two_phase_source` implements
this expression. All five checks in
`tests/adiabatic/multiphase/test_nonlinear_source.jl` pass against the library:
the two integrated quartic formulas below, a parity-forbidden angular output,
a phase label absent from the cubic convolution, and a zero field. The source
is now connected to the tested radial RHS (section 3.8). Twelve additional checks compare the
single-cloud limit for quartic, sextic, and cosine potentials: with only
$R_{1,(1,0)}$ populated, fundamental and generated third-harmonic sources agree
with the single-cloud implementation, the third harmonic is nonzero, and the
$(l,\mathbf q)=(3,(-1,2))$ source vanishes. These comparisons use a fixed
24-by-64-by-11 polar/phase grid, with $\psi_1=-\varphi$ to match the old
azimuthal convention. Subsequent checks in `test_source_validation.jl` cover the
complete quartic table, phase origins, and cosine quadrature as recorded below.

#### 3.5 Exact quartic benchmarks before the cosine potential

For the quartic potential, the nonlinear field is a finite Fourier
convolution. Extend the input set by its conjugates to
$\mathcal S=\mathcal Q\cup(-\mathcal Q)$. At each polar coordinate,

$$
\widehat h_{\mathbf q}(x)
\equiv\left\langle e^{i\mathbf q\cdot\boldsymbol\psi}
\frac{\Phi^3}{6}\right\rangle_\psi
=\frac16\sum_{\substack{\mathbf p,\mathbf s,\mathbf t\in\mathcal S\\
\mathbf p+\mathbf s+\mathbf t=\mathbf q}}
F_{\mathbf p}(x)F_{\mathbf s}(x)F_{\mathbf t}(x).
$$

The sum is over ordered triples; their multiplicities determine the
coefficients. It provides a reference independent of phase quadrature.
In particular, let the input contain only the two fundamentals,
$F_{(1,0)}=f_1(x)$ and $F_{(0,1)}=f_2(x)$. Selected exact outputs are:

| Output label | Exact phase coefficient $\widehat h_{\mathbf q}$ |
|---|---|
| $(1,0)$ | $\tfrac12\lvert f_1\rvert^2f_1+\lvert f_2\rvert^2f_1$ |
| $(0,1)$ | $\tfrac12\lvert f_2\rvert^2f_2+\lvert f_1\rvert^2f_2$ |
| $(2,-1)$ | $\tfrac12 f_1^2f_2^*$ |
| $(-1,2)$ | $\tfrac12 f_2^2f_1^*$ |
| $(3,0)$ | $f_1^3/6$ |
| $(0,3)$ | $f_2^3/6$ |
| $(1,1)$ | $0$ |

First check these phase coefficients, then integrate them against
$\alpha^2Y_{l,m_{\mathbf q}}(r^2+a_{\rm BH}^2x^2)$ to check the complete
radial source. For pure $211$/$322$ angular shapes, set
$f_1=uY_{11}$ and $f_2=vY_{22}$; $u,v$ are local complex radial coefficients.
Two fully integrated benchmarks in our polar normalization are

$$
N_{0,(2,-1)}=
\alpha^2u^2v^*\frac{\sqrt{30}}{20}
\left(r^2+\frac{a_{\rm BH}^2}{7}\right),
$$

$$
N_{3,(-1,2)}=
\alpha^2v^2u^*\frac{\sqrt{210}}{42}
\left(r^2+\frac{a_{\rm BH}^2}{11}\right).
$$

These follow by integrating the explicit polar polynomials and square-root
factors, and test the driven-source normalization. They test source terms,
not the driven profiles or their fluxes. Also check an angular projection
forbidden by polar parity, and a phase label absent from the convolution.

**Checked (October 14):** all seven rows of the phase-coefficient table pass
at every sampled polar node, and their integrated sources agree with the
library (14 checks). The two fully analytic integrated benchmarks and forbidden
outputs also pass in `test_nonlinear_source.jl`.

#### 3.6 Discrete quadrature and aliasing

For $K=2$, use Gauss–Legendre polar nodes $x_i$ with weights $w_i$ and uniform
phase nodes $\psi_{a,j}=2\pi j/n_a$, $j=0,\ldots,n_a-1$. Here $n_a$ is a
phase-grid size, not a cloud principal quantum number. Approximate

$$
\mathcal P_{l,\mathbf q}[G]\simeq
\sum_i w_iY_{l,m_{\mathbf q}}(x_i)
\frac{1}{n_1n_2}\sum_{j=0}^{n_1-1}\sum_{k=0}^{n_2-1}
e^{i(q_1\psi_{1,j}+q_2\psi_{2,k})}G(x_i,\psi_{1,j},\psi_{2,k}).
$$

Discrete phase orthogonality distinguishes integer vectors only modulo the
grid sizes. Frequencies separated by an integer multiple of $n_a$ in every
phase direction are aliased. Derive the grid requirement from input and
requested output supports, including the implicit conjugates.

For a polynomial term of degree $d$, let
$Q_a=\max_{\mathbf q\in\mathcal S}|q_a|$ and
$T_a=\max_{\mathbf q\in\mathcal Q_{\rm out}}|q_a|$, where
$\mathcal Q_{\rm out}$ is the set of requested output phase labels. A simple
sufficient condition for alias-free extraction of those outputs is

$$
n_a>dQ_a+T_a\quad\text{for every phase direction }a.
$$

For reconstruction/projection on the same support this gives $n_a>2Q_a$.
For a cubic remainder projected onto that support it gives $n_a>4Q_a$; for
the quintic term, $n_a>6Q_a$. Requesting higher omitted components can require
a larger grid. This bound concerns phase aliasing only; polar quadrature
must be checked separately.

The cosine remainder has unbounded Fourier support. Reuse its stable
small-field evaluation and refine phase resolution until target sources
converge. Vary $n_1$, $n_2$, and the polar node count independently while
holding the retained component set fixed. Then vary the retained Fourier
and angular sets separately. Quadrature convergence alone does not establish
truncation convergence of a solved field.

**Quadrature checkpoint (October 14):** `test_source_validation.jl` uses
$R_{1,(1,0)}=2.3+0.7i$, $R_{2,(0,1)}=1.6-0.5i$, $r=8$, $a=0.7$, and
$\alpha=0.2$. Seven requested outputs include both fundamentals, both driven
components, both third harmonics, and $(l,\mathbf q)=(5,(1,0))$. Against a
$64^3$ reference, varying one grid count while holding the others at 64 gives:

| Refined direction | 8 nodes | 16 nodes | 32 nodes |
|---|---:|---:|---:|
| Polar | $6.33\times10^{-3}$ | $6.98\times10^{-12}$ | $8.28\times10^{-14}$ |
| First phase | $7.39\times10^{-2}$ | $7.25\times10^{-7}$ | $9.62\times10^{-14}$ |
| Second phase | $4.51\times10^{-2}$ | $3.79\times10^{-8}$ | $1.24\times10^{-13}$ |

Entries are maximum absolute source errors over the seven outputs. The joint
$32^3$ grid also agrees within $10^{-10}$. On this host a warmed evaluation
including table setup took about 0.020 seconds and allocated 4.98 MB; one real
field occupies 256 KiB. Allocated bytes are not peak process memory. These
13 regression checks establish convergence for this fixed input/output support;
grid adequacy and retained-support convergence must be revisited for solved fields.

#### 3.7 Phase-origin and single-phase checks

For a vector of phase shifts $\boldsymbol\delta$, transform the input by

$$
\widetilde R_{l,\mathbf q}=e^{-i\mathbf q\cdot\boldsymbol\delta}R_{l,\mathbf q}.
$$

This gives $\widetilde\Phi(\boldsymbol\psi)=
\Phi(\boldsymbol\psi+\boldsymbol\delta)$, so the projected source must obey

$$
\widetilde N_{l,\mathbf q}=
e^{-i\mathbf q\cdot\boldsymbol\delta}N_{l,\mathbf q}.
$$

Phase-averaged quantities such as $\langle U(\Phi)\rangle_\psi$ remain
unchanged, up to quadrature error. Check shifts in both phase directions.

For comparison to the established single-cloud projector, retain only
$\mathbf q=(j,0)$ with positive integer $j$, and take $m_1=1$.
At zero temporal phase, $\psi_1=-\varphi$, giving
$e^{-ij\psi_1}=e^{+ij\varphi}$ as in the existing implementation. Map
$R_{l,(j,0)}$ to its existing $(l,m=j)$ input and compare projected sources
for all three potentials. The second phase average integrates a constant.
This verifies compatibility with the old convention rather than requiring
the old angular machinery to be rebuilt.

**Checked (October 14):** for all three potentials, source covariance and
phase-averaged potential invariance pass for shifts $(0.37,0)$, $(0,-0.53)$,
and $(0.37,-0.53)$ (18 checks). The shifts are not integer grid translations.
Potential invariance is tested at each polar node before polar integration.
The single-cloud comparisons above also pass.

#### 3.8 Order of implementation and acceptance

1. Specify the retained input pairs and requested output labels; establish
   polar and phase quadrature with the normalization above.
2. Recover the two coefficients in the equal-$l,m$ test in section 3.3.
3. Recover the quartic convolution coefficients and integrated driven-source
   benchmarks in section 3.5, including forbidden outputs.
4. Recover the single-phase sources and phase-origin covariance.
5. Demonstrate separate polar/phase convergence for the cosine potential,
   and measure selected omitted-component sources.
6. Feed the source into $d\mathbf Q/dr_*=-p\mathbf Q-V\mathbf R-s\mathbf N$
   with the same ordering as the linear block. Check the single-cloud limit
   of the complete radial right-hand side.

Each item is a small development step with an explicit expected result.
Measure runtime and memory for the initial two-phase grid; its point count
is $n_{\rm polar}n_1n_2$. Tensor phase grids grow rapidly with $K$, so defer
a scalable many-cloud projector until this prototype's cost is known.
Passing these checks establishes a nonlinear source evaluator, not a joint
stationary solution, an action closure, or physical evolution.

**RHS checkpoint (October 14):** the implemented `nonlinear_two_phase_rhs!` computes
one joint source and assembles each linear component block in state order
$[R_1,\ldots,R_n,Q_1,\ldots,Q_n]$. It passes 21 library checks in
`test_nonlinear_rhs.jl`: full single-cloud RHS comparisons for all three
potentials, component/angular permutations, state preservation, dimension
checks, analytic driven-source signs and tortoise scaling, and a zero-field
test with nonzero slopes. The function is in `profiles.jl`, and the tests are
registered in the multiphase runner; all 115 checks pass. This completes the
local source/RHS checkpoint. The implemented matching residual below connects
the RHS to the shared radial integrator and passes short integration and
step-refinement checks. The implemented Newton wrapper now produces a matched trial.

### 4. Construct a jointly deformed two-cloud family

**Numerical prototype (October 14):** follow the existing single-cloud matching
structure, prescribing real positive outer coefficients for the two fundamental
components. If $d$ is the number of retained $(l,\mathbf q)$ pairs, the trial
unknowns are $d$ complex horizon coefficients, $d-2$ complex outer coefficients,
and two complex trial frequencies: $4d$ real unknowns for $4d$ real profile/slope
matching equations. For the minimal $211$, $322$, driven $m=0$, and driven $m=3$
components with one angular term each, this is a 16-real-variable problem.
The imaginary trial frequencies extend the existing numerical formulation;
they do not by themselves establish nonlinear amplitude rates or the
slow-growth approximation required by milestone 1. Start with nonzero secondary
normalization so its frequency is identifiable.

**Matching-residual checkpoint:** `two_phase_matching_residual` assembles each
boundary in global [R; Q] order, integrates both sides, and returns their
mismatch. The 24 registered checks in `test_matching_residual.jl` pass against the library:
single-cloud histories for all potentials, saved/plain returns, joint component
ordering, input validation, and short radial integration in both directions.
For radial steps 0.4/0.2/0.1, successive endpoint differences shrink by factors
17.60 outward and 14.28 inward, consistent with fourth-order RK4. These checks
do not establish asymptotic-boundary accuracy or a physical matched profile.
The function is in `profiles.jl`; with the Newton checks below, the full
registered suite has 159 passing checks.

**Newton-wrapper library checkpoint (October 15):** `solve_two_phase_mode` fixes the
first retained l coefficient in each fundamental block to the corresponding
positive real normalization. Its trial vector contains two complex frequencies,
all horizon coefficients, and free outer coefficients. Coefficients and residuals
are scaled by the smaller normalization. The two fundamentals remain bound;
all retained component frequencies remain positive and retain their initial
bound/outgoing classifications during Newton corrections.

All 20 registered checks in `test_nonlinear_matching.jl` pass against the
inserted library implementation, reproducing the October 14 candidate. At spin=0.99,
alpha=0.42, normalizations $(10^{-4},10^{-5})$, the four-component quartic trial
with l=(1,2,0,3), radial endpoints (-30,80), match=5, step=0.4, and quadrature
8x9x9 converges in three corrections from scaled mismatch 13.3874 to
$1.2248\times10^{-10}$. The resulting trial frequencies are approximately
$0.408805528742+1.50386627\times10^{-7}i$ and
$0.415628386519+6.37117738\times10^{-11}i$. Both driven outer coefficients
are nonzero. The checks verify fixed normalizations, a freshly evaluated
matching residual, input preservation, component reordering, and invalid inputs.
This is a numerical matching test at modest resolution, not a validated physical
family or nonlinear growth prescription. The wrapper is in `profiles.jl`, and
its tests are registered in the multiphase runner.

Seed the new calculation with a saved primary and a small secondary. Continue
the secondary amplitude upward, solving for both fundamental profiles,
frequencies, and retained nonlinear components using the formulation from
milestone 1. Fix one amplitude normalization and phase convention per
independent fundamental in the nonresonant formulation.

Reuse the Newton solver and boundary machinery where applicable. Start with
selected configurations and short continuation paths rather than a dense
rectangular amplitude table. Include higher components according to an
explicit truncation whose effects can be measured.

Once milestone 1 has specified a consistent profile problem, its numerical
matching part will compare inner and outer integrations at a common radius:

$$
\mathbf G_{\rm match}=
\begin{pmatrix}
\mathbf R_H(r_{*,\rm match})-\mathbf R_\infty(r_{*,\rm match})\\
\mathbf Q_H(r_{*,\rm match})-\mathbf Q_\infty(r_{*,\rm match})
\end{pmatrix}=0.
$$

Here the vectors concatenate all retained $(l,\mathbf q)$ entries in a fixed
order. The inner and outer integrations use the joint nonlinear source,
not separate per-cloud nonlinearities. The unknowns include free leading
boundary coefficients and the carrier frequencies; any slow-growth or
solvability variables must come from milestone 1. Convert complex residuals
to real and imaginary parts for Newton differentiation, because reconstruction
depends on both coefficients and their conjugates.

Define the cloud normalization explicitly before counting unknowns. For
example, if a selected fundamental outer coefficient defines $A_a$, choose
that coefficient real and equal to $A_a$ and eliminate it from the free
unknowns. This is a possible coordinate choice, not an identification of
$A_a$ with cloud energy. Other normalization choices require their own
constraints and derivatives. Matching the frozen spatial equation alone
does not resolve the open-boundary solvability issue or prove that a real-
frequency periodic family exists.

**Checkpoint:** recover the primary-only limit and the corresponding
small-secondary truncation, then demonstrate a converged change in the primary
caused by the secondary. Compute total $E$ and $J$ from the joint field.
Agreement of Newton residuals alone does not establish radial, angular,
boundary, or Fourier convergence.

### 5. Compute amplitude responses and E/J evolution, then extend to actions

**Current implementation (October 15):** the vector-amplitude extension of
`core/sensitivity.jl` and the multiphase adapters in `response.jl` are present.
Production derivatives use implicit differentiation; independently corrected
neighbors provide validation. `evolution.jl` supplies the scaled E/J balance
solve and midpoint evolution with corrected profiles at every stage. The
matched-cloud validation campaign completed the neighbor checks and 4/8-step
cosine trajectories; its results and limits are recorded separately. The equations below describe the implemented approach
and the later extension to actions.

First expose reusable unpacking and the real-packed matching residual
G(u,A1,A2), where u contains real and imaginary parts of the two trial
frequencies and free boundary coefficients. Reconstruction uses conjugates,
so this is a real derivative problem. Handle the existing normalization
scaling consistently: the solver scales coefficients and residuals by
`minimum(normalizations)`. Decide on a fixed reference scale or include all
scale dependence before differentiating. Check ForwardDiff support through
the boundary series and joint radial integration; do not differentiate the
Newton iterations.

For the immediate two-cloud implementation, propagate the state response into
O=(E,J) with total derivative O_A + O_u u_A, including frequencies, fixed outer
normalizations, and induced radial-profile changes. Keep the background,
quadrature, radial grid, and retained components fixed. The observable adapter
must explicitly choose real carriers from the complex matching frequencies.
This computes a derivative of the discretized numerical family; it does not
validate the slow-growth approximation or spatial convergence.

After neighbor checks pass, assemble the 2x2 balance matrix in milestone 1,
using horizon plus outer losses for each right-hand side, monitor conditioning,
and develop a short slow-time trajectory. The action formulas and independent
checks below describe the later extension, not a prerequisite for this first
two-cloud E/J response implementation.

Extend [core/sensitivity.jl](../src/adiabatic/core/sensitivity.jl)
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

In particular, the action response must include both explicit and induced
profile dependence:

$$
D_{ab}=\left.\frac{\partial I_a}{\partial A_b}\right|_u
+\sum_j\frac{\partial I_a}{\partial u_j}\frac{\partial u_j}{\partial A_b}.
$$

Apply the same chain rule to $E$ and $J$. On a family for which the milestone 1
identities have been established, the independently evaluated derivatives
should satisfy

$$
\frac{dE}{dA_b}=\sum_a\omega_aD_{ab},\qquad
\frac{dJ}{dA_b}=\sum_a m_aD_{ab}.
$$

This explains the proposed closure's balance check: combining $D\dot{\mathbf A}
=\mathbf S$ with these identities yields

$$
\dot E=-\sum_{\mathbf q}F_{E,\mathbf q},\qquad
\dot J=-\sum_{\mathbf q}\frac{m_{\mathbf q}}{\Omega_{\mathbf q}}
F_{E,\mathbf q}.
$$

These are consequences conditional on the action identities and proposed
sources, not a derivation of them. Compute each channel's angular flux
independently to check its ratio to energy flux. In the single-cloud limit,
the same identities give the familiar rate
$\dot A=-F_E/(dE/dA)$, with $F_E$ the total signed loss. For the dilute
two-cloud channels, the proposed action source assigns changes proportional
to $(-2,+1)$ for $\mathbf q=(2,-1)$ and $(+1,-2)$ for $\mathbf q=(-1,2)$,
with proportionality factor $F_{E,\mathbf q}/\Omega_{\mathbf q}$ for each.
These supply explicit comparison targets for the closure derivation.

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

The first increment is **component representation plus the two-phase nonlinear
projector**, alongside the written closure derivation. Its useful
output is tested infrastructure and explicit conventions, without committing
the joint solver to an unverified evolution law. Follow the incremental
development conventions in [CLAUDE.md](../CLAUDE.md).

The labels and linear radial assembly are present, with repository unit tests.
The positive-frequency boundary wrapper is present; transfer its prototype
checks into the test suite and run them against the actual implementation.
The two-phase projector now passes the same-$l,m$ benchmark in section 3.3.
Nonlinear-source evaluation now passes the quartic table, single-cloud limit,
phase-origin, and fixed-support cosine convergence checks recorded in section 3.
The full RHS, matching residual, and Newton wrapper pass their registered
library checks. Next register the boundary-wrapper checks and establish
matched-profile convergence before physical continuation, as described in section 4.
The closure derivation remains open; the
milestone numbering does not require finishing it before these projection
checks, but the joint solver and production evolution depend on it.

The first physics result is a converged measurement of how a populated
secondary reshapes the primary. The first result establishing the new closure
is validated three-cloud evolution. Neither requires beginning with all modes
or a large amplitude grid.

Keep derivations and final conventions in design documentation, checkpoint
status in `HANDOFF.md`, and run commands, resolution studies, source revisions,
inputs, and saved evidence in campaign READMEs. Submit expensive solves through
Slurm under the existing project conventions. Keep each validated foundation
in a focused change on the agreed working branch, keeping experimental drivers
opt-in until their physical checks are complete.
