# Project handoff

Updated 2026-10-02. Active work: reproducing Omiya, Takahashi, Tanaka & Yoshino,
*Impact of multiple modes…* ([2211.01949v1](docs/papers/2211.01949v1.pdf)),
Sec. III — a nonlinear l = m = 1 primary with a linear l = m = 2 secondary on it,
evolved by energy and angular-momentum balance (Eqs. 62–63). The single-cloud
adiabatic work (OTT 2201.04382) is the primary's input; see
[adiabatic-repro.md](docs/adiabatic-repro.md).

Read first: [multimode-adiabatic-design.md](docs/multimode-adiabatic-design.md)
(equations, decisions, piece plan) and the
[secondary campaign README](runs/adiabatic/multimode_secondary/README.md)
(commands, results, checks). Cite equation numbers from the PDF in
`docs/papers/`, not the arXiv HTML.

## Where things stand

Pieces 1–3 of six are done at Mμ = 0.42, a/M = 0.99:

1. Period average ⟨1 − cos ϕ₁⟩ of the primary (`src/adiabatic/multimode/background.jl`).
2. Secondary eigenproblem, Eq. (61), continued along all 232 primary points
   (`multimode/secondary_mode.jl`, driver `drivers/continue_secondary.jl`).
3. Secondary energy Ê₂(A₁) and boundary fluxes (`multimode/secondary_energy.jl`,
   driver `drivers/secondary_energies.jl`).

Validated: linear limits against spheroidal Leaver (both clouds, both α);
small-A₁ coefficients (Re ω₂ shift −0.0101 A₁², Γ₂/Γ₂^lin − 1 = 42 A₁²);
exact energy balance 2Γ₂E₂ = −F to 3e-6 at every point; radial-step
convergence of the secondary at A₁ = 2.05 (RK4 ratio 16).

Pieces 1–3 also run at Mμ = 0.3 and **reproduce the paper's Fig. 7
quantitatively**: at E₁ = E₂ ≈ 660, nonlinear/linear peak |R| are
0.125/0.116 (l = m = 1) and 0.034/0.030 (l = m = 2), matching the figure to its
reading precision, with the same peak positions. Commands and the table are in
the campaign README.

## Decisions made (with reasons in the design doc)

- **Averaged reading of Eq. (61).** The background is ⟨1 − cos ϕ₁⟩, not the
  fixed-time product. The fixed-time version couples R to R* and made ω₂
  depend on the secondary's normalization phase (Γ₂ swung ±1e-4 at A₁ = 2);
  that run is kept as `results/cosine_secondary_fixed_time.*`. The paper's
  (ω, m) classification (Sec. II B) implies the average.
- **μ² restored** in the projection term that the printed Eq. (61) omits.
- **Secondary on the primary's domain** at twice its step; normalization by
  unit (2,2) outer amplitude. A₂'s normalization does not affect Eqs. (62)–(63).
- **Mμ = 0.3 Newton tolerance** follows the residual's change from one ulp of
  Re ω₂ (`AXION_SECONDARY_ULP_TOLERANCE_FACTOR=2`); a fixed tolerance fails
  partway along the branch because that floor grows as the secondary contracts.

## Known caveats

- The Mμ = 0.42 primary branch is **not step-converged above A₁ ≈ 9**: its
  step-0.1 solutions re-evaluated at step 0.05 have scaled residual 4e-5 at
  A₁ = 9.1 rising to 4e-2 at 13. This bears on the cosine flux zero near
  A₁ ≈ 12.96. It does not affect Fig. 6 at Mμ = 0.42, which stays below
  E₁ ≈ 1000 (A₁ ≲ 0.8, step-converged to ~1e-7).
- Unchecked: secondary l-truncation (l ≤ 6) and outer boundary at large A₁.
  The paper reports neither.

## Next step: piece 4

The m = 0 dissipative field at ω = 2ω₁ − ω₂ (Eqs. 51–52): the part of
μ²(1 − cos ϕ₁)ϕ₂ varying as e^{−i(2ω₁−ω₂)t} with m = 0, i.e. the e^{−2i(ω₁t−φ)}
part of 1 − cos ϕ₁ times ϕ₂*. Ingoing at the horizon, decaying at infinity
(2ω₁ − ω₂ < μ). Checks: horizon flux ∝ A₂²; small-A₁ power matches the
perturbative F₀ of Eq. (41). Then piece 5 (m = 3 at 2ω₂ − ω₁, just above μ, so
an outer-boundary issue) and piece 6 (Eqs. 62–63, Fig. 6).

## Working with the user

Build incrementally with a checkpoint per piece; the user types library code
from pasted snippets (paste only code that has been run), Claude may write
drivers, scripts, and runs. Define new notation in place, and prefer naming
physical quantities over bookkeeping symbols. Size checks to what the next
result needs; keep checks that can change a conclusion. Heavy runs go through
Slurm; check `git status` first — this is a shared workspace.

The previous handoff (2026-09-24, amplitude-evolution interpolation checks) is
in git history; its state is superseded by
[adiabatic-branch-work-2026-09-28.md](docs/adiabatic-branch-work-2026-09-28.md).
