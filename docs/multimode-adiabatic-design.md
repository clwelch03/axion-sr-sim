# Multimode adiabatic evolution: design and Sec. III reproduction

Equation numbers refer to the PDF
[2211.01949v1](papers/2211.01949v1.pdf) (Omiya, Takahashi, Tanaka & Yoshino),
not the arXiv HTML, whose numbering differs. Fixed Kerr background throughout.
Current capabilities, caveats, and next steps: [project handoff](../HANDOFF.md).
Results and commands: [secondary campaign](../runs/adiabatic/multimode_secondary/README.md).
This document owns equations, implementation decisions, and the component map.

See the [September 28 literature review](multimode-literature-review-2026-09-28.md)
for prior work and the narrower novelty claim: mutual nonlinear profile
response and its effect on evolution.

## Reproduction target: Sec. III

The field is split as ϕ = ϕ₁ + ϕ₂ + ϕ_r (Eq. 46) with, by Eqs. (47)–(48),

\[
\phi_1=\phi_1(A_1),\qquad \phi_2=A_2\,\phi_2(A_1),
\]

neglecting the secondary's deformation of the primary. Then

- Eq. (49): □ϕ₁ − μ² sin ϕ₁ = 0, the single-cloud problem already solved by
  `continue_branch.jl` (ansatz Eq. 53, harmonics of ω₁t − m₁φ);
- Eq. (50): □ϕ₂ − μ² (cos ϕ₁) ϕ₂ = 0, projected as Eq. (61) with ansatz (60);
- Eqs. (51)–(52): driven radiative fields with source μ²(cos ϕ₁)ϕ₂ − ½μ²(sin ϕ₁)ϕ₂² + …,
  projected onto (2ω₁ − ω₂, m = 0) and (2ω₂ − ω₁, m = 3);
- Eqs. (62)–(63): Ȧ₁, Ȧ₂ from total energy and angular-momentum balance, with
  E = E₁(A₁) + Ê₂(A₁)|A₂|² (Eqs. 68–69) and fluxes from Eqs. (64)–(67).

The paper's numerics are otherwise unstated beyond r* = ±100 for the primary at
Mμ = 0.42 and outer series order a₇ (p. 11).

## Decisions and their reasons

**Averaged reading of Eq. (61).** The integral in (61) is written at fixed t.
With two incommensurate clouds, (1 − cos ϕ₁)ϕ₂ contains two pieces with m = 2:
⟨1 − cos ϕ₁⟩R at ω₂ (the period average, ⟨·⟩ over one primary period, equal
to the azimuthal average at t = 0) and a piece at 4ω₁ − ω₂ proportional to R*.
We keep only the first. Reasons:

1. *Phase test.* With the fixed-time term, fixing the (2,2) outer amplitude to
   e^{iθ} instead of 1 changes ω₂: at A₁ = 2.05, Im ω₂ = +5.8e-8, −1.4e-4,
   −2.7e-8, +1.4e-4 for θ = 0, π/4, π/2, 3π/4. A physical frequency cannot
   depend on a normalization phase. The averaged equation is complex-linear and
   phase independent.
2. *The paper's bookkeeping.* Sec. II B (after Eq. 22) classifies every term by
   its (ω, m) pair; the components acting back on a cloud are those at its own
   (ω, m), e.g. (ω₂ + ω₁ − ω₁*, 2). The R* piece is a separate (4ω₁ − ω₂, 2)
   component, which would radiate (4ω₁ − ω₂ > μ).

Both readings agree at O(A₁²); the R* piece starts at O(A₁⁴). Code:
`secondary_background_source(...; reading=:averaged)` (default), with
`:fixed_time` kept only to reproduce the phase test.

**Missing μ².** Linearizing our radial equation (adiabatic-design.md §3) about
ϕ₁ gives the source α²(r² + a²x²)⟨1 − cos ϕ₁⟩ϕ₂. The printed (61) omits μ².

**Projection sign.** The printed weight e^{+im₂φ} selects the conjugate
component; the code projects with e^{−im₂φ}, consistent with the rest of (61).

**Secondary truncation and domain.** (l, m₂) = (2,2), (4,2), (6,2), the
analogue of the primary's l = 1, 3, 5 (the paper says only that l = 2
dominates). The secondary uses the primary's radial domain at twice its step,
so every RK4 stage reads a saved primary sample; no interpolation and no
extrapolation beyond the primary's outer boundary.

**Normalization of ϕ₂.** Fixed by unit (2,2) outer amplitude. The paper scales
ϕ₂ to unit energy instead. Eqs. (62)–(63) are balances of physical E and J, so
redefining A₂ → f(A₁)A₂ changes variables but not the trajectory in (E₁, E₂).

**Secondary energy.** The second-order term in ϕ₂ of μ²(1 − cos(ϕ₁ + ϕ₂)),
period-averaged, is ½μ²⟨cos ϕ₁⟩ϕ₂²: a mass term with effective mass squared
μ²⟨cos ϕ₁⟩. Kinetic and gradient terms reuse `energy_per_tortoise_radius` with
α = 0. On the time-independent averaged background this energy is conserved up
to boundary fluxes, so 2Γ₂E₂ = −(F_H + F_out) exactly. The campaign records
the numerical check of this identity for each saved branch.

## Reproduction pieces

The pieces identify components of the calculation; completion and validation
are tracked in the handoff and campaign. Paths below are relative to
`src/adiabatic/`.

| Piece | Content | Implementation |
|---|---|---|
| 1 | Period-averaged primary background | `multimode/background.jl` |
| 2 | Secondary eigenproblem (61) along the primary branch | `multimode/secondary_mode.jl` |
| 3 | Secondary energy, fluxes, and energy balance | `multimode/secondary_energy.jl` |
| 3b | Repeat pieces 1–3 at the Fig. 7 parameters | Same routines, separate campaign inputs |
| 4 | Driven $m = 0$ field and horizon flux | `multimode/dissipative_mode.jl`; [piece 4 guide](multimode-piece4-guide.md) |
| 5 | Driven $m = 3$ field and flux to infinity | `multimode/radiative_mode.jl` |
| 6 | Energy/angular-momentum balance (62)–(63) | `multimode/two_cloud_evolution.jl` |

**Piece 6 design (confirmed 2026-10-05).** Variables are the cloud energies
$E_1$ (labels the primary branch; monotonic on both saved branches) and
$E_2 = \hat E_2(A_1)|A_2|^2$ (Eq. 69). $E_2 = |A_2|^2$ in the paper's
normalization (unit-energy secondary), so this is their (62)–(63); it avoids
differentiating our normalization-dependent $\hat E_2$. With
$J = J_1(E_1) + 2E_2/\omega_2(E_1)$ and $dJ_1 = dE_1/\omega_1$ (Eq. 72):

$$
\dot E_1 + \dot E_2 = -F_E, \qquad
\Big(\frac{1}{\omega_1} - \frac{2E_2}{\omega_2^2}\frac{d\omega_2}{dE_1}\Big)\dot E_1
+ \frac{2}{\omega_2}\dot E_2 = -F_J .
$$

$F_E$ sums the primary's own flux $F_1(E_1)$, the secondary's own flux
$\propto E_2$, the $m = 0$ flux $\propto E_2$, and the $m = 3$ flux
($F_\infty + F_H$, Eq. (65)) $\propto E_2^2$. In $F_J$ these carry
$1/\omega_1$, $2/\omega_2$, $0$, $3/\omega_3$. All frequencies are real parts.
Coefficients are not-a-knot cubic splines in $E_1$ of each flux divided by
its leading small-field power; $d\omega_2/dE_1$ comes from the $\omega_2$
spline. (Splines in $\ln E_1$, first proposed, were 5–100 times worse when
built from every other point, and their $d\omega_2/dE_1$ was off by 120% at
small $E_1$: the coefficients are smooth in $E_1 \propto A_1^2$, not in
$\ln E_1$.) Initial data and time unit as in Fig. 6.

Validation identities for pieces 4–5: field amplitude ∝ A₂ (m = 0) and ∝ A₂² (m = 3);
small-A₁ powers against the perturbative F₀, F₃ (Eqs. 38–39, (A11), (A14)); F_J = (m/Ω)F_E
per channel. The campaign records which checks have been performed;
open validation work is listed in the handoff.

### What a configuration and its fluxes mean

In the current model, a *cloud configuration* is the instantaneous pair
$(E_1,E_2)$ at fixed BH parameters. On the selected monotonic primary branch,
$E_1$ selects the nonlinear primary profile, the secondary eigenmode on that
background, and the associated driven responses. $E_2$ sets the secondary's
population through $|A_2|^2=E_2/\hat E_2(E_1)$. The spatial profiles follow
these stationary solutions adiabatically; the secondary does not yet reshape
the primary, and the driven fields have no independently evolved populations.

The boundary currents were calculated when constructing the stationary
tables. `multimode/scalar_energy_accounting.jl` reconstructs them along a saved
trajectory by interpolating coefficients in $E_1$ and restoring the energy
factors:

| Channel | Energy flux at each boundary |
|---|---|
| Primary | $C_1(E_1)E_1$ |
| Secondary | $C_2(E_1)E_2$ |
| Driven $m=0$ | $C_0(E_1)E_1^2E_2$ |
| Driven $m=3$ | $C_3(E_1)E_1E_2^2$ |

Each boundary has its own coefficient. These coefficients retain the
nonlinear dependence on the primary profile; factoring out the leading
small-field powers does not replace that dependence with a perturbative
approximation. Dividing the tabulated fluxes by the corresponding reference
energies removes the arbitrary secondary amplitude normalization.

Summing the channels gives $F_H$, positive into the BH, and $F_{\rm out}$,
positive outward. Thus negative $F_H$ means net extraction, including any
energy returned through the dissipative channel. Outward accounting uses
finite-radius currents for the primary, secondary, and $m=0$ fields, and
the infinity flux for $m=3$. This radiative channel is scalar radiation;
GW losses are not included. Multiplying each channel's energy flux by
$m/\omega$ gives its angular-momentum flux, so $m=0$ transports energy but
no angular momentum.

Fluxes are rates per $x=t/M$. Convert saved times $\tau=\gamma_1 x$ back to
$x$ before integrating. From the initial saved time $x_0$, define

$$
E_{\rm extracted}(x)=-\int_{x_0}^{x}F_H\,dx',\qquad
E_{\rm outward}(x)=\int_{x_0}^{x}F_{\rm out}\,dx'.
$$

The energy check is
$E_{\rm cloud}(x)-E_{\rm cloud}(x_0)=E_{\rm extracted}-E_{\rm outward}$,
with $E_{\rm cloud}=E_1+E_2$. Angular-momentum fluxes are integrated with
the same signs. This checks consistency with the evolution's balance
equations, not the validity of the underlying adiabatic approximation.
`observables/physical_normalization.jl` supplies the conversion to physical
units; the normalization script reports powers, torques, and cumulative
fractions of the initial BH reservoirs. This remains postprocessing on a
fixed background: a cloud-energy plateau can coexist with continuing BH
extraction balanced by outward radiation. Evolving the BH mass and spin
requires feeding those horizon fluxes back into the dynamics.

## Code identities to keep separate

The assembled stationary routines assign frequencies as
`frequency * m / fundamental_m`; their `(l, m)` entries are harmonics of one
pattern, not independent clouds. The secondary therefore runs as its own
block with `fundamental_m = m₂`. Keep three identities distinct:

- a bound cloud, labeled by its linear-limit (n, l, m) and continued eigenbranch;
- a Fourier component, labeled by an integer vector **q** with
  Ω_q = **q**·**ω** and m_q = **q**·**m**;
- the angular index l within each component.

The $m = 0$ channel is a zero-$m$ block, which `linear_tortoise_rhs!` and
`boundary_state` reject (they require $m > 0$). `multimode/dissipative_mode.jl`
therefore assembles it directly from `tortoise_operator_coefficients` and the
series helpers `horizon_boundary_data` / `outer_boundary_data`, which accept
$m = 0$.

## Beyond Sec. III

These are extension ideas. See [the handoff](../HANDOFF.md#next-steps) for
project priorities and the [October 6 proposal](further_work_1006.md) for the
comparison of possible research directions.

**Two mutually nonlinear clouds.** A quasi-periodic field Φ(r, θ; ψ₁, ψ₂) with
ψ_a = ω_a t − m_a φ, expanded in Fourier vectors **q**, generalizes the
single-cloud harmonic expansion. Sec. III is a truncation of it (|q₂| ≤ 1 plus
the O(A₂²) source, no feedback onto the (k, 0) components). For two clouds the
two relative phases are removable by time translation and rotation when
m₁ω₂ ≠ m₂ω₁, so the family is labeled by (A₁, A₂) and Eqs. (62)–(63) still close.

**Three or more clouds.** E and J balances give at most two equations:

\[
\begin{pmatrix}
\partial_{A_1}E&\cdots&\partial_{A_K}E\\
\partial_{A_1}J&\cdots&\partial_{A_K}J
\end{pmatrix}\dot{\boldsymbol A}
=-\begin{pmatrix}F_E\\F_J\end{pmatrix}.
\]

A candidate closure is Whitham's averaged-Lagrangian theory: one action
I_a = ∂⟨L⟩/∂ω_a per phase, with dI_a/dt = −Σ_q q_a 𝒩_q from the number fluxes
of driven components, J = Σ m_a I_a and dE = Σ ω_a dI_a on shell. In the weak
field this reduces to the stoichiometric rate network
ν_ac = δ_ak − δ_ai − δ_aj for processes i + j → k + emitted field, with
Σ_a ω_a ν_ac = −Ω_c and Σ_a m_a ν_ac = −m_c as cheap checks. It needs a written
derivation for open Kerr boundaries and complex ω before use. Near-resonant
**q**·**ω** (small divisors) would make a relative phase physical.

[Witte and Mummery, Secs. II.3–II.4](https://arxiv.org/html/2412.03655v3#S2.SS4)
give a quartic rate-network comparison with higher bound levels.
