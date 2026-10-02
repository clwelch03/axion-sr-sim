# Multimode adiabatic evolution: design and Sec. III reproduction

Revised 2026-10-02. Equation numbers refer to the PDF
[2211.01949v1](papers/2211.01949v1.pdf) (Omiya, Takahashi, Tanaka & Yoshino),
not the arXiv HTML, whose numbering differs. Fixed Kerr background throughout.
Results and commands: [secondary campaign](../runs/adiabatic/multimode_secondary/README.md).

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
to boundary fluxes, so 2Γ₂E₂ = −(F_H + F_out) exactly; the code satisfies it to
3e-6 at every branch point.

## Reproduction pieces

| Piece | Content | Status |
|---|---|---|
| 1 | ⟨1 − cos ϕ₁⟩ from a saved primary | done; exact quadrature and small-field checks |
| 2 | Secondary eigenproblem (61) along the primary branch | done at Mμ = 0.42; Leaver linear limit, A₁² scaling, step convergence at A₁ = 2.05 |
| 3 | Ê₂(A₁), fluxes, energy balance | done at Mμ = 0.42 |
| 3b | Pieces 1–3 at Mμ = 0.3 for Fig. 7 | done; Fig. 7 peaks reproduced to reading precision |
| 4 | Driven m = 0 field at 2ω₁ − ω₂, horizon flux | next |
| 5 | Driven m = 3 field at 2ω₂ − ω₁, flux to infinity | 2ω₂ − ω₁ ≈ 0.423 vs μ = 0.42: radial wavelength ~120M, outer-boundary care needed |
| 6 | J, Eqs. (62)–(63), Fig. 6 | J₂ = m₂E₂/Re ω₂ for the single-frequency secondary; J₁ from δE₁/ω₁ = δJ₁/m₁ (Eq. 72) |

Planned checks for pieces 4–5: field amplitude ∝ A₂ (m = 0) and ∝ A₂² (m = 3);
small-A₁ powers against the perturbative F₀, F₃ (Eqs. 41–42); F_J = (m/Ω)F_E
per channel.

## Code identities to keep separate

The assembled stationary routines assign frequencies as
`frequency * m / fundamental_m`; their `(l, m)` entries are harmonics of one
pattern, not independent clouds. The secondary therefore runs as its own
block with `fundamental_m = m₂`. Keep three identities distinct:

- a bound cloud, labeled by its linear-limit (n, l, m) and continued eigenbranch;
- a Fourier component, labeled by an integer vector **q** with
  Ω_q = **q**·**ω** and m_q = **q**·**m**;
- the angular index l within each component.

The m = 0 channel needs a zero-m block, which the existing boundary routines
reject (they require m > 0).

## Beyond Sec. III

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
