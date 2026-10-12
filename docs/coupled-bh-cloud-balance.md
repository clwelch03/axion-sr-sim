# Coupled BH–cloud balance: first scalar-only model

This is the October 8, 2026 balance checkpoint following the local stationary
[background-response checks](../runs/adiabatic/background_continuation/README.md).
It supplies a leading small-cloud-mass evolution law and separates the
additional terms needed for a cloud on a genuinely time-dependent metric.
That checkpoint did not supply a time-integrated trajectory or calculate those
additional response terms. Its equations impose no QCD relation between mass
and decay constant; the later campaign below selects that relation explicitly.

The instantaneous balance functions now live in
[`src/adiabatic/multimode/coupled_balance.jl`](../src/adiabatic/multimode/coupled_balance.jl),
loaded through `src/adiabatic/load.jl`. The campaign's
`check_coupled_balance.jl` supplies saved inputs and validation only.
This was the first implementation piece; the local alpha=0.42 integrator draft
remains paused for incremental review. Later on October 8, the user explicitly
authorized completing a separate first graph directly. The
[QCD spin-down campaign](../runs/adiabatic/qcd_spindown/README.md) now integrates
the leading law for a 20-solar-mass BH, using fresh stationary inputs and
explicit classical seeds. It does not supply the higher-order response terms
derived below.

## 1. Fixed reference units and state

Use $G=c=\hbar=1$. Let $M_0$ be the initial BH mass, $\mu$ the fixed particle
mass, and $\epsilon_f=Gf_a^2=(f_a/M_{\rm Pl})^2$, with the **unreduced** Planck
mass. Define the following notation for this implementation:

$$
\tau=t/M_0,\quad b=M_{\rm BH}/M_0,\quad
j=J_{\rm BH}/M_0^2,\quad
q_i=E_i^{\rm physical}/(\epsilon_f M_0).
$$

The four independent state variables are $(b,j,q_1,q_2)$. Initially $b=1$ and
$j=\chi_0$. The quantities passed to the instantaneous $M_{\rm BH}=1$ solver
are

$$
\boxed{\chi=j/b^2,\qquad\alpha=\alpha_0 b,\qquad e_i=q_i/b,}
\qquad \alpha_0=M_0\mu.
$$

Here $e_i$ are the existing code cloud energies, called `primary_energy` and
`secondary_energy` in the fixed-background routines. A prime below means
$d/d\tau$. Holding $q_i$ constant as $b$ changes does not hold $e_i$ constant.
The dimensionless cloud mass fraction is $\eta=\epsilon_f(e_1+e_2)$.

## 2. Boundary fluxes

At the current $(e_1,e_2,\alpha,\chi)$, let $F_H,F_o$ be the summed **code
energy fluxes**, positive into the BH and outward, respectively. Let $T_H,T_o$
be the corresponding **code angular-momentum fluxes**. The subscript $o$
uses the existing finite outer boundary for the primary, secondary, and $m=0$
fields, and the infinity current for $m=3$; it is not uniformly an exact
infinity extraction.

For each boundary separately, sum the retained channels:

| Channel | Frequency in instantaneous units | Angular/energy flux ratio |
|---|---|---|
| Primary and its locked harmonics | $\hat\omega_1=M\omega_{1,R}$ | $1/\hat\omega_1$ |
| Secondary | $\hat\omega_2=M\omega_{2,R}$ | $2/\hat\omega_2$ |
| Driven $m=0$ | $2\hat\omega_1-\hat\omega_2$ | $0$ |
| Driven $m=3$ | $2\hat\omega_2-\hat\omega_1$ | $3/(2\hat\omega_2-\hat\omega_1)$ |

Thus $m=0$ can return mass to the BH without returning angular momentum.
The physical energy flux is $\epsilon_f F$, whereas the physical torque is
$\epsilon_f M_{\rm BH}T$. The extra mass factor in torque matters.

At leading order in the cloud's gravitational backreaction, with no accretion,

$$
\boxed{b'=\epsilon_f F_H,\qquad j'=\epsilon_f bT_H,}
$$

$$
\boxed{\chi'=\frac{\epsilon_f}{b}(T_H-2\chi F_H),\qquad
\alpha'=\alpha_0\epsilon_f F_H.}
$$

Only horizon flux enters these BH equations. Scalar emission outward changes
the cloud, not the BH directly. GW contributions are absent in this model.

## 3. Leading instantaneous-background cloud law

Our fixed-Kerr two-cloud model uses

$$
g(e_1,e_2;\alpha,\chi)
=g_1(e_1;\alpha,\chi)+\frac{2e_2}{\hat\omega_2(e_1;\alpha,\chi)},
\qquad \partial_{e_1}g_1=\frac1{\hat\omega_1}.
$$

Here $g$ is cloud angular momentum in units $\epsilon_f M_{\rm BH}^2$.
This is the existing Sec. III approximation from
[Omiya et al.](https://arxiv.org/abs/2211.01949), PDF Eqs. (62)–(63), (72),
as implemented in `two_cloud_rates`; it neglects secondary deformation of
the primary. In particular, $g_1\ne e_1/\hat\omega_1$ on a nonlinear branch.

Define

$$
a_1=\frac1{\hat\omega_1}
-\frac{2e_2}{\hat\omega_2^2}\partial_{e_1}\hat\omega_2,
\qquad a_2=\frac2{\hat\omega_2}.
$$

The leading instantaneous-background extension is

$$
\boxed{
\begin{pmatrix}1&1\\a_1&a_2\end{pmatrix}
\begin{pmatrix}q_1'\\q_2'\end{pmatrix}
=-\begin{pmatrix}F_H+F_o\\T_H+T_o\end{pmatrix}.}
$$

All coefficients are evaluated on the **current** background. There is no
extra $b$, $\epsilon_f$, or growth-rate factor multiplying the solution of this
system. For fixed $b$, $E_i^{\rm physical}=\epsilon_f M_{\rm BH}e_i$ and
$dt=M_{\rm BH}\,d(t/M_{\rm BH})$ cancel the mass factors in the energy rate.
The code energies, if needed as dependent variables, obey

$$e_i'=\frac{q_i'-e_i b'}b.$$

Changing background coefficients without this energy-unit conversion would
introduce spurious cloud mass changes. The conserved-energy diagnostic of
this **truncated model** is

$$[b+\epsilon_f(q_1+q_2)]'=-\epsilon_f F_o.$$

This identity is exact algebraically for this model; it is not a claim of
exact gravitational energy accounting. Monitor the balance matrix's
conditioning and avoid continuing through its singularities.

## 4. What has been dropped, and why it matters

An instantaneous stationary family and a scalar field on a time-dependent Kerr
metric are not the same approximation. From $\nabla_a T^a{}_b=0$,

$$
\partial_a(\sqrt{-g}T^a{}_t)
=\frac12\sqrt{-g}T^{ab}\partial_tg_{ab}.
$$

For energy density $-\sqrt{-g}T^t{}_t$, the cloud balance therefore contains

$$
\dot E_{\rm cl}^{\rm physical}=-F_E^{H,\rm physical}
-F_E^{o,\rm physical}+W,
\qquad W=-\frac12\int d^3x\sqrt{-g}\,\langle T^{ab}\rangle\partial_tg_{ab}.
$$

The brackets denote the same fast-phase averaging as the cloud model.
This formula assumes a specified common coordinate system, time normalization,
and fixed coordinate integration boundaries. Moving horizon/extraction
surfaces add transport terms; the instantaneous tortoise grids are **not**
already that common coordinate construction. A physical implementation must
make these choices consistently, using regular horizon coordinates as needed.
Axisymmetry gives no analogous bulk angular-momentum work term.

Write $w=W/\epsilon_f$ and

$$K(q_1,q_2,b,j)=b^2g(q_1/b,q_2/b;\alpha_0b,j/b^2).$$

The cloud angular momentum in fixed units is $\epsilon_f K$. Within this
cloud-charge approximation, the fuller projected balance is

$$
\boxed{
\begin{pmatrix}1&1\\a_1&a_2\end{pmatrix}
\begin{pmatrix}q_1'\\q_2'\end{pmatrix}
=\begin{pmatrix}
-(F_H+F_o)+w\\
-(T_H+T_o)-(K_b b'+K_j j')/b
\end{pmatrix}.}
$$

The partial derivatives $K_b,K_j$ hold $q_1,q_2$ fixed. They include the
explicit mass scaling and the profile/frequency response, not just a change
in $\alpha$. Writing $g_\alpha,g_\chi$ at fixed $e_i$,

$$K_b=b[2g-e_1a_1-e_2a_2+\alpha g_\alpha-2\chi g_\chi],
\qquad K_j=g_\chi.$$

Neither $g_1$'s background derivatives nor $w$ are supplied by the current
energy/flux tables. In particular, differentiating energies at fixed $A_0$
does not measure work at fixed adiabatic action. Setting
$W=E_M\dot M+E_J\dot J$ without specifying what cloud state is held fixed
would not resolve this issue.

For a light cloud with regular parameter response and self-driven BH changes,
$b',j'$ are $O(\epsilon_f)$ at fixed code configuration. Both omitted cloud
terms in the boxed equation are then $O(\epsilon_f)$ relative to its ordinary
code-rate scale; their contributions to physical charge rates are
$O(\epsilon_f^2)$. The relevant practical control is the cloud mass fraction
times its parameter sensitivity, not $\epsilon_f$ alone. Near a resonance,
ill-conditioned balance, or cancellation of large fluxes this counting need
not control the small net rate. It requires an error estimate there.

The leading law above drops **both** terms at that order. Its actual family
angular-momentum audit has remainder
$[j+\epsilon_f K]'+\epsilon_f bT_o
=\epsilon_f(K_b b'+K_j j')$. Including the response term cancels it.
With metric work retained but BH mass still driven only by the leading
horizon flux, the energy audit instead has remainder $\epsilon_f w$.
Completing total gravitational energy at that order also needs consistent
metric/backreaction and boundary-charge accounting; arbitrarily subtracting
$W$ from the BH horizon mass equation is not a derived solution.

## 5. Checks and next implementation step

The runnable checkpoint is
`runs/adiabatic/background_continuation/scripts/check_coupled_balance.jl`.
It checks recovery of the existing nonlinear fixed-BH RHS, all four dilute
channel transfer signs and multiplicities, nonunit-mass unit conversion,
the BH spin derivative, and the energy/angular audits above. An independent
two-mode model $q_i=\Omega_i(b,j)N_i$ with conserved occupations tests that
the work and parameter-response terms recover the correct energy changes
under an externally prescribed slow background variation. It also explicitly
shows that setting work to zero fails that test.

The script evaluates a local leading-order RHS at the saved baseline, with
an **illustrative** $e_2=0.1$, $f_a=10^{16}$ GeV, $b=1$, and $\chi=0.99$.
These are test inputs, not selected initial conditions for a physical run.
Results and check counts belong in the campaign README.

The first trajectory can use the leading instantaneous-background system,
provided it is labeled as such and remains in the light-cloud regime with
regular response. It still needs an adaptive sequence of stationary solves
or interpolation in amplitude, alpha, and spin, and physical seeds and $f_a$.
Going beyond that approximation requires the cloud angular-momentum response
and metric-work/boundary calculation above. The algebraic checkpoint does
not validate those unimplemented corrections.
