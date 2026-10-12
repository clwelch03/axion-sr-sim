**Cloud self-gravity in the adiabatic framework — literature review, 2026-09-28**

The best fit to this repository is a perturbation in gravitational coupling
around each already nonlinear scalar branch point. Solve for the gravitational
field sourced by that configuration, then differentiate the scalar matching
problem to obtain its response. The cosine potential need not be expanded in
field amplitude. A multipole Poisson calculation is a useful first diagnostic;
a sourced metric perturbation on Kerr is the systematic relativistic extension.
These are proposals, not implemented or validated gravity calculations.

**October 7 update:** the standalone dilute Newtonian $211$ benchmark below
is now checked. Gravity on the saved nonlinear profiles and gravitational
feedback remain unimplemented.

This supplements the earlier [time-domain backreaction review](gravitational-backreaction-review.md).
The present implementation already has nonlinear matching, energy/flux
observables, implicit branch derivatives, and an amplitude integrator in
`src/adiabatic/`. Some older design prose predates the last two capabilities.

**Literature most relevant to the choice of method.**

| Reference | What it contributes and how I would use it |
|---|---|
| [Omiya, Takahashi & Tanaka, 2022, *Adiabatic evolution of the self-interacting axion field around rotating black holes*](https://arxiv.org/abs/2201.04382), Sec. 3 | The existing method: nonlinear scalar configurations with an approximate helical symmetry, followed through slow energy balance. This is the scalar foundation to retain. |
| [Baryakhtar et al., 2021, *Black hole superradiance of self-interacting scalar fields*](https://arxiv.org/html/2011.11646#A7), Appendix G | Poisson gravity, a multipole Green function, and hydrogenic frequency shifts. The dilute-limit normalization and sign benchmark. |
| [Siemonsen, May & East, 2023, *SuperRad*](https://arxiv.org/html/2211.03845#S4), Sec. IV | Compares Newtonian, quasi-relativistic, and relativistic frequency shifts. Useful for benchmarking; its treatment neglects nongravitational self-interactions. |
| [Cannizzaro, Sberna, Green & Hollands, 2024, *Relativistic perturbation theory for black-hole boson clouds*](https://arxiv.org/abs/2309.10021) | A relativistic bilinear product for dissipative Kerr modes and frequency perturbation theory. The self-gravity example still uses a semi-Newtonian, monopole treatment; a relativistic product alone does not solve the Einstein equations. |
| [May, East & Siemonsen, 2025, *Self-gravity effects of ultralight boson clouds formed by black hole superradiance*](https://arxiv.org/html/2410.21442v2), Sec. III and Appendix B | Stationary, axisymmetric Einstein–complex-scalar/vector solutions, with synchronization and no nongravitational self-interaction. Provides relativistic frequency benchmarks and explicit elliptic equations. Newtonian shifts can underestimate the relativistic result by an order-one factor. |
| [Takahashi, 2026, *Relativistic frequency shifts in gravitational waves from axion clouds*](https://arxiv.org/html/2604.21239v2#S3.SS3), Secs. III.3–IV.1 | Especially relevant: a general metric-induced Hamiltonian perturbation, then a practical semi-Newtonian implementation with spatial metric terms, relativistic source, and higher multipoles. Reports about 10% discrepancy in the self-gravity shift at alpha = 0.4 compared with the relativistic benchmark. |
| [Wardell, Kavanagh & Dolan, 2025, *Sourced metric perturbations of Kerr spacetime in Lorenz gauge*](https://arxiv.org/abs/2406.12510) | Reconstruction from six sourced scalar equations. Crucially, static contributions require separate treatment: the formula is not a ready-made solution for the stationary cloud potential. |
| [Hollands & Toomani, *Metric Reconstruction in Kerr Spacetime*](https://arxiv.org/abs/2405.18604) | Sourced reconstruction with transport correctors and an ADM-determined algebraically special contribution. Relevant because matter fills the cloud region; vacuum reconstruction alone is insufficient. |
| [Kim & Lenoci, 2025, *Self-Gravity in Superradiance Clouds*](https://arxiv.org/abs/2508.08367) | Density-dependent level shifts can reorder hyperfine levels and alter binary resonances. More relevant to a future multimode/binary extension than the present isolated single-pattern branch. |

The first three papers to read for implementation are Takahashi Sec. III.3,
May–East–Siemonsen Appendix B, and Cannizzaro et al. The distinction between
relativistic scalar perturbation theory and relativistic gravitational sourcing
is essential at our alpha = 0.42. The cited 10% is an error in the *shift* in a
particular near-threshold comparison, not in the total frequency and not an
error estimate for our nonlinear branch.

**Proposed expansion and physical normalization.**

Let Phi = phi/f_a denote the dimensionless real field used by the solver. For
the cosine potential,

\[
\Box_g\Phi-\mu^2\sin\Phi=0,\qquad
T_{ab}=f_a^2\left[\partial_a\Phi\partial_b\Phi
-g_{ab}\left(\tfrac12(\nabla\Phi)^2+\mu^2(1-\cos\Phi)\right)\right].
\]

The additional parameter is kappa = 8 pi G f_a^2. It equals
(f_a / reduced-Planck-mass)^2. Existing energies are normalized by f_a^2;
specifying A0 alone therefore does not specify the physical cloud mass.
At fixed dimensionless profile, changing f_a changes self-gravity without
changing the existing fixed-Kerr scalar solution.

Write g = g_K + h and Phi = Phi_0 + delta-Phi, retaining first order in
the gravitational perturbation. The proposed sequential equations are

\[
\delta G_{ab}[h]=8\pi G T_{ab}[\Phi_0;g_K],
\]

\[
\left(\Box_K-\mu^2\cos\Phi_0\right)\delta\Phi
=-\delta\Box[h]\Phi_0.
\]

The second equation is shorthand for the linearized periodic boundary-value
problem: its unknowns also include the frequency correction, and its amplitude
and phase conventions must be fixed. All covariant derivatives in

\[
\delta\Box[h]\Phi=-h^{ab}\nabla_a\nabla_b\Phi
-\left(\nabla_a h^{ab}-\tfrac12\nabla^b h\right)\nabla_b\Phi
\]

use Kerr. Thus metric feedback changes derivative terms as well as any
effective potential. A shift of the scalar mass term alone is insufficient.

This expansion treats self-interaction nonperturbatively and gravity
perturbatively. Corrections to T from delta-Phi enter the metric only at the
next order in kappa. Consequently, away from singular scalar response, one
metric solve followed by one scalar response solve gives the first-order
answer; iterative self-consistency is optional and resums selected higher
orders. It does not make linearized gravity exact at higher orders.

**The implementation opportunity: reuse the branch Jacobian.**

Let F(z,A0,kappa) = 0 be the normalized, real-packed matching residual, with z
containing the frequency and shooting amplitudes. For each saved gravity-free
solution z0, obtain h = kappa h1 from its retained profile. Then

\[
J z_1=-F_{\kappa},\qquad
J=\left.\frac{\partial F}{\partial z}\right|_{z_0,\kappa=0},
\qquad z=z_0+\kappa z_1+O(\kappa^2).
\]

This is the same implicit differentiation pattern already used by
[branch_sensitivity.jl](../src/adiabatic/core/sensitivity.jl).
The new work is the metric source/solve and its contribution to the radial
operator and boundary data. It is not just an extra invocation of the existing
derivative routine. Frequency dependence of the boundary series must remain
inside F, together with explicit gravitational boundary corrections.

For an observable O, compute

\[
\frac{dO}{d\kappa}=O_z z_1+O_{\kappa}.
\]

This produces changes in omega, the profile, energy, and fluxes. A linear-mode
matrix element is a useful dilute-limit cross-check, but inserting a strongly
nonlinear profile into a linear eigenmode formula does not capture its full
response. In particular, the cosine term couples variations of all retained
harmonics. Near a poorly conditioned J, a small metric can cause a large scalar
response; check the response norm and continuation, not just the size of h.

**Which gravitational equations to solve first?**

For a first diagnostic, reconstruct the real field and derivatives, average
over its fast phase, and project the chosen density onto Y_L0. In the
nonrelativistic limit rho = mu |Psi|^2. Then

\[
U_L(r)=-\frac{4\pi G}{2L+1}
\left[r^{-L-1}\int_0^r s^{L+2}\rho_L(s)ds
+r^L\int_r^\infty s^{1-L}\rho_L(s)ds\right].
\]

This requires two radial cumulative integrals per multipole. For a hydrogenic
211 density, L = 0 and 2 suffice; nonlinear distorted profiles require an
independent multipole convergence check. The isolated Newtonian expression
assumes an inner continuation and vanishing external sources. On our exterior
grid, these must become explicit matching conditions. Include the radial
Jacobian dr/dr_* and the potential of shells outside the evaluation radius.

The frequency perturbation in this limit is the expectation of mu U, while
the gravitational binding energy is one half of integral rho U. The earlier
review records the analytic hydrogenic benchmark

\[
\frac{\delta\omega_{211}}{\mu}
=-\frac{237}{1280}\alpha^2\frac{M_c}{M_{\rm BH}}
\]

including L = 0,2. It agrees with the rounded coefficient in Baryakhtar et al.
This is a normalization benchmark, not a reliable finite-alpha prediction.

To reproduce Takahashi's diagnostic, use that paper's specified T_tt source
and temporal *and spatial* metric perturbations. T_tt, normal-observer energy
density, and the Killing energy integrand are different on Kerr; choosing
among them is part of the approximation, not interchangeable notation.
Extending that prescription to the full cosine source is our proposed model,
and needs its own validation.

For controlled relativistic conservative gravity, my preferred next step is
a stationary axisymmetric *linearized Einstein* boundary-value solve sourced
by the phase-averaged stress tensor. It must include stresses and angular
momentum currents, hence cloud-induced frame dragging. A direct gauge-fixed
solve in (r,theta), or a tensor-harmonic radial system, seems more direct than
beginning with general radiative metric reconstruction. It needs regular
horizon data, asymptotic time normalization, constraints, and a convention
separating BH mass/spin perturbations from cloud gravity. Averaging and dropping
slow evolution must respect source conservation to the retained adiabatic order.

An extended cloud cannot be modeled by replacing M_BH with M_BH + M_c in
the Kerr operator. Its potential inside the cloud and its horizon response
are different from a central point mass.

**A real cloud also sources oscillating gravity.**

With chi = omega_R t - m0 varphi, the current odd scalar harmonics give

\[
T_{ab}=T^{(0)}_{ab}+
\sum_{p=2,4,\ldots}\left[T^{(p)}_{ab}e^{-ip\chi}+\mathrm{c.c.}\right].
\]

Thus the metric needs zero and even harmonics, with gravitational frequency
p omega_R and azimuthal number p m0. At leading quadratic order a single
fundamental produces p = 0 and 2; the full cosine source produces further
even harmonics. A phase average is axisymmetric for this single locked pattern.
Independent modes can retain slow, nonaxisymmetric beat terms.

The nonzero sectors can also be solved in the frequency domain, with ingoing
horizon and outgoing infinity conditions. No evolution over every scalar
oscillation is necessary. Their scalar feedback follows by projecting
delta-Box[h] Phi_0. In particular, h^(2) times the conjugate fundamental can
feed the fundamental: discarding oscillating gravity is an additional
approximation, not a consequence of adiabaticity alone. Evaluate its size
before claiming a complete relativistic first-order real-field result.

Gravitational wave fluxes can be obtained from a sourced Teukolsky calculation,
but curvature/wave extraction alone does not supply the full conservative
metric. This is where the sourced reconstruction references become relevant.
An exactly periodic outgoing solution is only a local-in-slow-time model:
perpetual scalar radiation tails cannot be integrated to infinity as a finite,
stationary self-gravitating cloud. Separate the bound near zone and outgoing
radiation with a consistent matching/extraction prescription.

**Energy balance and the useful scientific outputs.**

Once branch corrections are available, preserve the existing amplitude
evolution architecture. For fixed BH parameters and a consistently defined
cloud-plus-binding energy E_cons, write schematically

\[
E_{{\rm cons},A}\dot A
=-F_H^{\phi}-F_\infty^{\phi}
-F_H^{\rm GW}-F_\infty^{\rm GW}.
\]

Horizon fluxes are positive into the BH and infinity fluxes positive outward.
Use a near-zone energy/boundary-charge definition consistent with the metric
approximation; a scalar stress-energy integral on fixed Kerr is not the full
self-gravitating energy. In a Newtonian implementation the added binding term
has the one-half factor above. Avoid counting it twice through matter energy
and a separate potential term. GW flux, quadratic in h, can still contribute
at the leading relative gravitational order to cloud loss.

If BH parameters evolve, the left side also includes E_M dot-M and E_J dot-J,
and BH evolution uses the corresponding total horizon fluxes. Keep mu fixed
and use fixed reference units when M_BH changes. Fix how the background Kerr
parameters and homogeneous metric mass/spin pieces divide this evolution.

The existing preliminary cosine branch has a scalar flux zero near A0 = 12.96,
with spatial convergence still outstanding in the
[September 28 note](adiabatic-branch-work-2026-09-28.md). A particularly useful
new measurement is the displacement of that zero. If L(A,kappa) denotes net
loss and L_A is nonzero at the original zero,

\[
\left.\frac{dA_{\rm sat}}{d\kappa}\right|_0=-\frac{L_\kappa}{L_A}.
\]

Also measure the change in the branch's energy slope and any turning point.
A turning point alone is not a complete dynamical stability calculation.
These quantities connect self-gravity directly to the present project rather
than only to gravitational-wave phase.

**Recommended sequence.**

1. Specify f_a/Planck-mass, normalize the source, and obtain M_c/M_BH and a
   weak-field potential estimate from a few existing profiles.
2. Implement the frozen-profile monopole/quadrupole diagnostic, recover the
   hydrogenic result, and compare the dilute relativistic shift with published
   benchmarks at matching spin and normalization.
3. Add the metric-induced matching residual and solve J z1 = -F_kappa.
   Validate the linear response against small finite changes in kappa within
   the same gravitational approximation. Track energy, flux and saturation
   changes as well as frequency.
4. Replace the semi-Newtonian gravity solve with a stationary Kerr metric
   perturbation solve if precision at alpha = 0.42 is required. Quantify the
   omitted oscillating sectors; add sourced waves and their feedback as needed.

The distinctive extension is a gravitational response of the *nonlinear*
adiabatic branch. The papers reviewed supply the ingredients and benchmarks,
but I did not find a ready-made calculation combining our full cosine branch,
its radiative harmonics, and self-consistent relativistic metric feedback.

## 211 benchmark checkpoint — October 7, 2026

This first piece solves the Newtonian potential of an analytic, dilute,
single-mode cloud. It uses no Kerr profiles or horizon boundary conditions.
The derivation follows the Poisson and frequency-perturbation prescription in
[Baryakhtar et al., Appendix G](https://arxiv.org/html/2011.11646#A7);
the rational coefficients below are independently evaluated here.

Use $c=\hbar=1$ with $G$ explicit, $\alpha=GM_{\rm BH}\mu$, and define the
Bohr radius $a_B=1/(\mu\alpha)$. This notation avoids confusing the cloud
length scale with Kerr spin. Normalize the single-particle wavefunction by
$\int|\psi_{211}|^2d^3r=1$:

$$
\psi_{211}=R_{21}(r)Y_{11}(\theta,\varphi),\qquad
R_{21}=\frac{1}{2\sqrt6\,a_B^{3/2}}
        \frac{r}{a_B}e^{-r/(2a_B)},\qquad
|Y_{11}|^2=\frac{3}{8\pi}\sin^2\theta.
$$

For occupation $N$, the leading nonrelativistic mass is $M_c=N\mu$ and the
fast-oscillation-averaged density is $\rho=M_c|\psi_{211}|^2$.
Define $x=r/a_B$ and the dimensionless density $f=a_B^3\rho/M_c$.
With $P_L$ denoting Legendre polynomials, rather than normalized spherical
harmonics,

$$
f(x,\theta)=\frac{x^2e^{-x}}{96\pi}
             [1-P_2(\cos\theta)]
           =f_0(x)+f_2(x)P_2(\cos\theta),
\qquad f_2=-f_0.
$$

The normalization is $4\pi\int_0^\infty x^2f_0dx=1$.
Write the cloud potential as
$U=(GM_c/a_B)\sum_{L=0,2}v_L(x)P_L(\cos\theta)$.
Regularity at the origin and $U\to0$ at infinity give

$$
v_L(x)=-\frac{4\pi}{2L+1}\left[
 x^{-L-1}\int_0^x y^{L+2}f_L(y)dy
 +x^L\int_x^\infty y^{1-L}f_L(y)dy\right].
$$

The two integrals include matter inside and outside the evaluation radius.
In particular, $v_0(0)=-1/4$, $v_2(0)=0$, $v_0\sim-1/x$, and
$v_2\sim6/x^3$ at large $x$. The BH potential is already responsible for the
hydrogenic wavefunction; it is not included again in $U$.

Angular orthogonality gives the frequency shift

$$
\frac{\delta\omega}{\mu}
=\frac{GM_c}{a_B}\,4\pi\sum_{L=0,2}\frac{1}{2L+1}
  \int_0^\infty x^2 f_L(x)v_L(x)dx
=-\frac{GM_c}{a_B}\left(\frac{93}{512}+\frac{9}{2560}\right).
$$

For an independent exact evaluation, use the radial probability density
$w(x)=x^4e^{-x}/24$. The symmetric radial integrals are

$$
J_L=\frac{2}{24^2}\int_0^\infty x^{3-L}e^{-x}
       \int_0^x y^{4+L}e^{-y}dy\,dx.
$$

The monopole coefficient is $J_0=93/512$; the quadrupole coefficient is
$J_2/25=9/2560$, with the factor $1/25$ supplied by the angular projection
and Green function. The inner integral is a finite incomplete-gamma
polynomial, so these values can be obtained with exact rational arithmetic.
Consequently,

$$
\boxed{\frac{\delta\omega_{211}}{\mu}
=-\frac{237}{1280}\alpha^2\frac{M_c}{M_{\rm BH}}.}
$$

There is no factor $1/2$ in the frequency shift. The total gravitational
binding energy instead is $W=\tfrac12\int\rho U\,d^3r
=\tfrac12N\delta\omega$ for this fixed profile. Differentiating
$W\propto N^2$ with respect to $N$ recovers $\delta\omega$.
For the fixed-background code normalization $M_c/M_{\rm BH}=\epsilon_f E_1$,
this becomes $GM_{\rm BH}\delta\omega=-(237/1280)\alpha^3\epsilon_f E_1$,
at the same leading dilute order, with the unreduced-Planck convention
$\epsilon_f=(f_a/M_{\rm Pl})^2$.

Run the standalone check from the repository root:

```bash
julia --project=. tests/adiabatic/core/check_self_gravity_211.jl
```

The script projects the density using the repo's normalized polar harmonics,
evaluates the potential with cumulative radial trapezoids, and integrates its
expectation value. It independently evaluates $J_L$ with rational arithmetic.
All 16 checks passed. At $x_{\max}=40$ and 8 polar nodes:

| Radial spacing | Monopole coefficient | Quadrupole coefficient | Absolute error in their sum |
|---|---:|---:|---:|
| 0.08 | 0.1816458337 | 0.0035166666 | $6.250\times10^{-6}$ |
| 0.04 | 0.1816419271 | 0.0035158854 | $1.563\times10^{-6}$ |
| 0.02 | 0.1816409505 | 0.0035156901 | $3.906\times10^{-7}$ |
| 0.01 | 0.1816407064 | 0.0035156413 | $9.766\times10^{-8}$ |

The exact sum is $0.18515625$. Refinement gives second-order radial
convergence. Separate checks increase the polar nodes from 8 to 16 and
$x_{\max}$ from 40 to 60; coefficient changes satisfy relative tolerances
of $10^{-13}$ and $10^{-11}$, respectively. Normalization, absence of an
$L=4$ source, origin limits, and far-field moments also pass.
This validates the analytic-source Newtonian calculation, not the density
prescription, boundaries, or frequency response for nonlinear Kerr clouds.
