**Gravitational backreaction for the coupled scalar towers — research note, 2026-09-15**

The current angular spectral reduction can also represent gravitational
backreaction. The useful first extension is a weak, slowly varying gravitational
potential represented by radial multipoles, together with a separate treatment
of black-hole mass and angular-momentum exchange. Evolving a general metric is
possible, but requires a new gravitational solver and a rederived scalar
operator. A small number of purely radial metric functions cannot exactly
describe the rotating real axion cloud.

This note distinguishes established results in the cited papers from proposed
implementation choices for this repository. No gravitational evolution has
been implemented or tested here.

**What the repository currently evolves.**

[solver.py](../src/solver.py) implements Dolan's fixed-Kerr scalar equations as
radial towers, with a tridiagonal solve over angular channels at each radial
point. `Tower.__init__` fixes `mass = 1`, constructs `r(x)`, and caches the metric
coefficients, potential and horizon angular velocity.
[nonlinear_source.py](../src/nonlinear_source.py) reconstructs

\[
\Phi(t,r,\theta,\varphi)=\frac{2}{r}\operatorname{Re}
\sum_{m>0,\ell}\psi_{\ell m}(t,x)Y_{\ell m}(\theta,\varphi),
\qquad \Phi=\phi/f_a,
\]

and projects the sine remainder into the scalar RHS before the tridiagonal
solve. [nonlinear_tower.py](../src/nonlinear_tower.py) advances positive odd
\(m\) towers together, recomputing their common source at each RK4 stage.
[sine_remainder.jl](../src/sine_remainder.jl) is an arithmetic kernel for that
potential; gravitational coupling belongs in the common RHS and geometry
machinery. The nonlinear progress note identifies physical energy,
angular-momentum and flux diagnostics as outstanding work.

The reduction retains angular information; it is not spherical symmetry or
intrinsically two-dimensional gravity. Dolan's azimuth is horizon adapted:
\(d\varphi=d\phi_{\rm BL}+a\,dr/\Delta\). Stress-energy calculations must use
the metric and derivatives in consistent coordinates. In particular,

\[
\left.\partial_r\right|_{\phi_{\rm BL}}=
\frac{r^2+a^2}{\Delta}\left.\partial_x\right|_{\varphi}
+\frac{a}{\Delta}\partial_\varphi.
\]

The radial derivative of the reconstructed field also differentiates its
\(1/r\) factor. See [Dolan (2013), Sec. III, especially Eqs. (18)–(36)](https://arxiv.org/html/1212.1477#S3).

**The physical coupling that is currently absent.**

With signature \((-+++)\), \(c=\hbar=1\), and
\(V(\phi)=\mu^2 f_a^2[1-\cos(\phi/f_a)]\), the system is

\[
\Box_g\Phi-\mu^2\sin\Phi=0,
\qquad G_{ab}=8\pi G T_{ab},
\]

\[
T_{ab}=f_a^2\left[
\nabla_a\Phi\nabla_b\Phi
-g_{ab}\left(\frac12\nabla_c\Phi\nabla^c\Phi
+\mu^2[1-\cos\Phi]\right)\right].
\]

Thus the additional dimensionless parameter can be
\(\kappa_g=8\pi Gf_a^2=(f_a/\bar M_{\rm Pl})^2\), where
\(\bar M_{\rm Pl}=(8\pi G)^{-1/2}\) is the reduced Planck mass. Papers often
instead use the unreduced \(M_{\rm Pl}=G^{-1/2}\); record the convention.
The existing dimensionless field history does not determine its physical
cloud mass until this normalization is specified. Large \(\Phi\) and strong
sine nonlinearities can coexist with a small metric perturbation.
Weak metric deformation alone does not justify Poisson gravity: relativistic
pressure, currents or rapid collapse can require the linearized Einstein
equations even while \(|h|\ll1\).

Construct stress-energy from the real field and all its first derivatives,
including angular derivatives, using the full cosine potential. The stored
`pi` is \(\partial_t\psi\), not the normal-observer scalar momentum. The PML
auxiliary `chi` is not matter. A complex scalar's stress tensor cannot simply
be substituted for the reconstructed real field's stress tensor; doing so can
change normalization and discard oscillating terms.

**Three different effects should be distinguished.**

| Effect | Additional unknowns | What it captures | Principal limitation |
|---|---|---|---|
| Black-hole response | \(M_{\rm BH}(t),J_{\rm BH}(t)\) | Spin-down, mass loss, changing superradiant threshold | Does not describe the cloud's gravitational field |
| Conservative cloud gravity | Potential multipoles \(U_{Lq}(t,r)\), or metric perturbations | Frequency shifts, changes in shape, gravitational mode coupling | Poisson gravity omits gravitational radiation and relativistic stress/current effects |
| Dynamical spacetime | Tensor metric variables and their momenta | Retardation, gravitational waves, strong deformation and collapse | Requires gauge, constraints and appropriate boundaries |

These can be combined, with a consistent division between the black-hole
parameters and the cloud-sourced perturbation. The gravitational monopole of
an extended cloud is not equivalent to adding its entire mass to the central
black hole.

**Black-hole evolution through horizon fluxes.**

Define \(F_E^H,F_J^H\) positive into the hole. In geometrized units,

\[
\dot M_{\rm BH}=F_E^H,\qquad
\dot J_{\rm BH}=F_J^H,\qquad
\chi_{\rm BH}=J_{\rm BH}/M_{\rm BH}^2.
\]

Consequently,

\[
\dot\chi_{\rm BH}=\frac{F_J^H}{M_{\rm BH}^2}
-2\chi_{\rm BH}\frac{F_E^H}{M_{\rm BH}}.
\]

Both horizon fluxes are negative for a superradiant positive-frequency mode.
For one nearly monochromatic mode,
\(F_J^H=(m/\omega_R)F_E^H\); a general nonlinear field requires the full
fluxes. Their stationary-background currents are
\(j_E^a=-T^a{}_b(\partial_t)^b\) and
\(j_J^a=T^a{}_b(\partial_\varphi)^b\).
This is the quasiadiabatic approach of
[Brito, Cardoso & Pani (2015), Sec. 3](https://arxiv.org/html/1411.0686#S3).

In this code, an evolving Kerr approximation requires more than assigning
`tower.spin`. Separate \(a=J/M\) from \(\chi=J/M^2\); they coincide
numerically only when \(M=1\). Use fixed reference units based on the initial
mass \(M_0\). The particle mass \(\mu\) stays constant while
\(\alpha(t)=M_{\rm BH}(t)\mu\) changes. The existing use of `alpha**2` as
the mass term is valid because the present code fixes \(M=1\).

Also choose a coordinate prescription before rebuilding geometry. Keeping the
same array indices while recomputing a new tortoise map relocates the field
without transporting it. A time-dependent tortoise/azimuthal map introduces
coordinate-advection terms. A fixed coordinate chart with a rederived
time-dependent operator, or an explicitly controlled adiabatic remapping
scheme, is needed. Simply substituting time-dependent Kerr parameters is an
approximation, not an exact solution of Einstein's equations.

Measure fluxes on a physical extraction surface outside the PML. Relating that
flux to the horizon requires accounting for propagation and energy stored
between the surface and the hole, or a justified slow-evolution limit. PML
dissipation itself is not the horizon flux. Verify the fixed-background
balance \(\dot E_{\rm ext}=-F_E^{\rm inner}-F_E^{\rm outer}\), and its angular
momentum counterpart, before using the fluxes to move the black hole.

**A radial multipole Poisson solve is the simplest gravitational field solve.**

In the dilute, weak-field cloud region, let \(\rho_c\) be the slowly varying
mass density. For a nonrelativistic envelope defined by
\(\phi=(\Psi e^{-i\mu t}+\Psi^*e^{i\mu t})/\sqrt{2\mu}\), its leading value
is \(\rho_c=\mu|\Psi|^2\). Expand the cloud's density and its additional
potential in ordinary spherical harmonics:

\[
\rho_c=\sum_{Lq}\rho_{Lq}(t,r)Y_{Lq},\qquad
U=\sum_{Lq}U_{Lq}(t,r)Y_{Lq}.
\]

Projecting \(\nabla^2U=4\pi G\rho_c\) gives, independently for each multipole,

\[
\frac{1}{r^2}\partial_r(r^2\partial_r U_{Lq})
-\frac{L(L+1)}{r^2}U_{Lq}=4\pi G\rho_{Lq}.
\]

The isolated Newtonian solution is

\[
U_{Lq}(r)=-\frac{4\pi G}{2L+1}\left[
r^{-L-1}\int_0^r s^{L+2}\rho_{Lq}(s)\,ds
+r^L\int_r^\infty s^{1-L}\rho_{Lq}(s)\,ds\right].
\]

This needs two cumulative radial integrals per multipole, so the radial
solve scales as \(O(N_rN_{\rm multipoles})\). Angular reconstruction and
projection remain separate costs. On the existing tortoise grid use
\(dr=\Delta\,dx/(r^2+a^2)\); integrating with uniform `dx` alone is wrong.
A separate gravity grid in the cloud region is another option.

The limits above describe a flat-space isolated-source problem. On a truncated
exterior domain, state what happens to omitted matter. For an empty inner
cavity, the regular cloud-generated solution behaves as \(r^L\); for no
exterior sources the outer solution behaves as \(r^{-L-1}\). Do not set both
potential boundaries to zero. These are Newtonian matching assumptions, not
relativistic black-hole horizon conditions. Include the potential of exterior
shells: even their spatially constant interior contribution changes frequencies
relative to time at infinity.

For the dilute \(211\) mode, \(\rho_c\propto|Y_{11}|^2\), and

\[
|Y_{11}|^2=\frac{1}{4\pi}[1-P_2(\cos\theta)].
\]

Therefore \(L=0,q=0\) and \(L=2,q=0\) are the exact density multipoles in
that approximation. Starting with the monopole provides a useful checkpoint;
adding the quadrupole restores this angular density. On Kerr, with multiple
modes or finite sine amplitudes, converge the retained gravitational
multipoles independently of the scalar truncation and angular quadrature.

An analytic benchmark follows directly from the isolated Green function. With
normalized hydrogenic \(R_{21}=r e^{-r/(2a_0)}/(2\sqrt6\,a_0^{5/2})\),
\(a_0=(\mu\alpha)^{-1}\), and \(G=c=\hbar=1\),

\[
\delta\omega_{211}^{L=0}=-\frac{93}{512}\mu\alpha^2\frac{M_c}{M},
\qquad
\delta\omega_{211}^{L=0,2}=-\frac{237}{1280}\mu\alpha^2\frac{M_c}{M}.
\]

I checked the dimensionless radial integrals using exact rational arithmetic:
the contributions are \(93/512\) and \(9/2560\). This checks the proposed
benchmark, not a numerical gravity implementation. The frequency correction
is the expectation of \(\mu U\); the total self-binding energy contains a
factor of one half and must not be confused with the single-mode frequency.
The multipole treatment and rounded \(211\) coefficient are consistent with
[Baryakhtar et al. (2021), Appendix G](https://arxiv.org/html/2011.11646#A7).

**How the potential couples back to the present relativistic PDE.**

For a nonrelativistic envelope the coupling is simply
\(i\partial_t\Psi=[-\nabla^2/(2\mu)-M\mu/r+\mu U+\cdots]\Psi\).
That envelope equation by itself has no Kerr horizon superradiance; a separate
horizon matching or growth/absorption prescription is required. It is useful
for slow cloud evolution but does not replace the current relativistic solver
without loss of physics.

For the current real-field solver, instead choose a weak metric perturbation
\(g=g_K+h\) and derive its effect on the wave operator:

\[
\Box_{g_K}\Phi-\mu^2\sin\Phi=-\delta\Box[h]\Phi+O(h^2),
\]

\[
\delta\Box[h]\Phi=-h^{ab}\nabla_a\nabla_b\Phi
-\left(\nabla_a h^{ab}-\tfrac12\nabla^b h\right)\nabla_b\Phi,
\qquad h=g_K^{ab}h_{ab}.
\]

Indices and covariant derivatives here use the background Kerr metric. In a
weak-field chart a Newtonian potential determines both temporal and spatial
metric perturbations. Their extension and matching to the near-hole region
must be specified; adding a flat-space potential everywhere is not a
controlled strong-field metric solution.

Unlike the sine remainder, this correction contains field derivatives,
including second time derivatives. Schematically the projected equation is

\[
(A_K+\delta A)\dot\pi=b_K+b_{\rm NL}+\delta b.
\]

At first order in \(h\), one implementation option preserves the existing
tridiagonal solve:

\[
a_0=A_K^{-1}(b_K+b_{\rm NL}),\qquad
\dot\pi=a_0+A_K^{-1}(\delta b-\delta A\,a_0)+O(h^2).
\]

This is an order-reduced approximation in the physical interior. Otherwise
solve the perturbed acceleration matrix, which need not remain tridiagonal.
`additional_b` can carry a derived first-order correction, but is not by
itself a generic metric-coupling interface. The geometric projection factors
must be derived as for the existing nonlinear-source normalization test.

At each RK stage the proposed sequence is: reconstruct the field and its
derivatives; construct the chosen slow gravitational source; project density;
solve for \(U\); evaluate and project the gravitational correction; advance
all scalar towers together. A slow potential updated less often needs its own
time-interpolation or multirate accuracy check. If a cycle average is used,
define it explicitly: averaging over fast \(2\mu\) oscillations does not
automatically remove slow intermode beats, which can be nonaxisymmetric.

Two particularly relevant perturbative references are
[Cannizzaro, Sberna, Green & Hollands (2024)](https://arxiv.org/abs/2309.10021),
which develops relativistic mode perturbation theory, and
[Takahashi (2026), Sec. III.3](https://arxiv.org/html/2604.21239v2#S3.SS3),
which includes spatial metric perturbations, a relativistic source and higher
potential multipoles. The latter still uses a Newtonian gravitational solve;
its comparison retains about a 10% discrepancy in the self-gravity frequency
shift at \(M\mu=0.4\). This is a relevant accuracy warning for the repo's
\(\alpha\simeq0.42\) regime, not a measured error for this proposed solver.

**Evolving gravitational perturbations or the full metric.**

For small metric deformation, solve the sourced linearized Einstein equations,

\[
\delta G_{ab}[h]=8\pi G T_{ab}[\Phi;g_K],
\]

with a gauge choice, constraint-compatible data and physical boundaries.
Tensor or spin-weighted angular harmonics again leave coefficients depending
only on \(t,r\). Kerr produces angular coupling, and nonlinear feedback
produces further couplings. This retains the radial-tower architecture in
principle, but requires a gravitational evolution system; it is not another
instance of the existing massive-scalar `Tower`.

The gravity sector needs different azimuthal modes. Already at quadratic order,
a real \(m=1\) scalar sources stationary \(q=0\) and oscillatory \(q=\pm2\)
stress-energy. With the current odd-m scalar symmetry, stress-energy has even
azimuthal numbers including zero. The scalar-only positive-odd-m validator
therefore cannot be reused for gravitational fields. The stationary rotating
cloud also has a momentum current, whose frame-dragging response is absent
from a scalar Poisson potential. At finite amplitude the full cosine energy
has arbitrarily high even harmonics; a quadratic \(2\ell_{\max}\) rule is not
an exact truncation prescription for the full problem.

The sourced spin-2 Teukolsky equation is a useful intermediate project for GW
generation. Computing outgoing curvature/waves from the scalar stress tensor
does not alone return the cloud's complete gravitational field or feed a
radiation reaction into the scalar evolution. Metric reconstruction and the
nonradiative mass/angular-momentum pieces are also needed for full feedback.
For sources distributed through the cloud, vacuum reconstruction alone is
insufficient. [Hollands & Toomani, *Metric Reconstruction in Kerr Spacetime*](https://arxiv.org/abs/2405.18604)
describes a sourced initial-value construction with reconstruction, transport
corrections and an ADM-determined contribution. It supplies a formal route,
but is a substantial implementation project.

For strong gravitational deformation, evolve the Einstein-scalar system with
a well-posed formulation such as generalized harmonic or BSSN/CCZ4, including
initial constraints, lapse/shift gauge, horizon treatment and outgoing
boundaries. An angular spectral representation remains possible, with
convergence in the number of modes. Exact spherical symmetry admits radial
Einstein-scalar systems but removes the rotating \(m=1\) cloud and Kerr
superradiance. Axisymmetry of a real rotating field's stress tensor requires
an additional approximation; a complex single-frequency field can instead
have an axisymmetric stress tensor exactly.

[May, East & Siemonsen (2025)](https://arxiv.org/abs/2410.21442) takes another
useful route: stationary axisymmetric Einstein-matter solutions with a
synchronized complex cloud, interpreted as an oscillation-averaged model of
real clouds. It provides relativistic self-gravity frequency benchmarks and
updates SuperRad. This is not a full dynamical, multimode, cosine-axion
evolution; it is particularly useful for validation near saturation.
[East & Pretorius (2017)](https://arxiv.org/abs/1704.04791) demonstrates
self-consistent dynamical growth and saturation in Einstein-Proca evolution,
but uses a complex vector field and axisymmetric stress-energy, not the real
scalar system here.

**Other relevant literature and a staged recommendation.**

| Paper | Relevance and boundary of applicability |
|---|---|
| [Yoshino & Kodama (2012), *Bosenova collapse of axion cloud around a rotating black hole*](https://arxiv.org/abs/1203.5070) | Nonlinear sine-Gordon cloud dynamics; distinguishes self-interaction collapse from gravitational metric collapse. |
| [Yoshino & Kodama (2014), *Gravitational radiation from an axion cloud around a black hole*](https://arxiv.org/abs/1312.2326) | Sourced gravitational-wave calculation during superradiant growth; a benchmark for a future radiation sector. |
| [Omiya, Takahashi & Tanaka (2022), *Adiabatic evolution of the self-interacting axion field around rotating black holes*](https://arxiv.org/abs/2201.04382) | Nonlinear stationary cloud sequences and flux-driven slow evolution. A way to avoid resolving every oscillation over the growth time; its original treatment omits BH spin-down and additional independent superradiant modes. |

My proposed order for this repository is:

1. Add physical normalization and stress-energy diagnostics. Check real-field
   factors, the cosine potential, derivative coordinates and physical
   inner/outer flux balance on fixed Kerr. Report \(M_c/M\), angular momentum
   transfer and an estimate of \(|U|\) before deciding the required gravity
   approximation.
2. Build a gravity diagnostic on frozen, dilute field profiles: monopole then
   quadrupole Poisson integrals. Recover the analytic \(211\) shifts above,
   check vacuum tails and radial/angular convergence, and compare a
   relativistic mode calculation with the published benchmarks.
3. Derive the first-order metric correction to the scalar RHS and feed the
   potential back at the RK stages. Recover the current nonlinear evolution
   when gravity is disabled. At small amplitude verify the gravity-induced
   frequency shift scales with \(f_a^2\) and field amplitude squared.
4. Couple horizon fluxes to \(M_{\rm BH},J_{\rm BH}\), after settling the fixed
   units and moving-geometry treatment. Check flux signs, spin-down and
   approach to the superradiant threshold. Avoid double-counting BH mass/spin
   changes in both Kerr parameters and the perturbation.
5. Add sourced gravitational radiation or fully relativistic metric
   perturbations if the science needs them. A weak self-gravity calculation
   can accumulate a significant phase shift while the metric remains small;
   small cloud mass alone is not a reason to ignore gravity in phase studies.

For long nonlinear evolutions, the current frozen outer boundary also needs
attention: self-interactions generate escaping scalar waves, and a gravity
sector would generate escaping GWs. The bound-state falloff argument used in
the linear solver does not absorb either radiation channel. Relativistic
gravity sources must be conserved to the retained perturbative order; sourcing
them with PML-damped matter throughout the computational domain would violate
that requirement. These boundary and diagnostic questions should be settled
before interpreting a long evolution as physical saturation.
