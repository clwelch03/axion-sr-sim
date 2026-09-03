# Reproducing Dolan 1212.1477 Sec. III from scratch

Branch: `repro/dolan-coupled-tower`. Living handoff doc — kept current so this
can be picked up cold at any point. See `CLAUDE.md` ("Build from-scratch
reproductions incrementally...") for the working rule this branch follows:
small pieces, built and checked together, not autonomous whole-module drops.

## Why this branch exists

The prior `radial1d/` build (single-mode ℓ=m=1 truncation, on branch
`feat/eigenmode-seed-and-backreaction`) got the right answer but drifted from
actually reproducing Dolan's method: it approximated away the coupled-ℓ tower
(Sec. III, Eq. 25/38) rather than implementing it, and picked up backreaction/
nonlinear/multilevel extensions before the foundation was faithfully rebuilt.
Decision: start over, reproduce Dolan's Sec. III (equations, discretization,
PML, boundaries) exactly, validate against his own numbers, and only then
revisit backreaction/detuning/multilevel on top of the faithful version.

Old code is deleted from this branch (`git rm -r radial1d/`, not yet
committed) but still lives on `feat/eigenmode-seed-and-backreaction` and in
history if anything needs to be checked against it.

## What we have from the paper so far

- **Sec. III.D (finite-difference method):** grid spacing `Δx = M/n`, `n =
  4..16`; `Δt = κΔx`, `κ=0.8`; Method of Lines, RK4 in time, standard 4th-order
  5-point stencils (his Eq. 37, matches what we had before). Grid indexed
  `ψ_(i) = ψ_j(x = x_+ + iΔx)`, integer `i` (`x_+` = chosen left edge of the
  truncated domain, not the literal horizon — tortoise `x` diverges there).
- **Eq. 38, `Aẋ = b`:** the coupled-ℓ system. `x^T = [π_k, π_{k+2}, ...]`
  (same-parity ℓ only), `A` real symmetric diagonally-dominant tridiagonal (in
  ℓ), solved every grid point / every RK stage via the Thomas algorithm. This
  is the piece the old single-mode code never built.
- **Sec. III.E (PML):** near-horizon region `x/M ≪ 0` decouples in ℓ (the
  `~Δ` coupling terms vanish there), reducing to the per-ℓ advective equation
  `ψ̈ - ψ'' + 2imΩψ̇ - 2imΩψ' = 0` (his Eq. 43). PML via complex-coordinate
  stretching gives (Eq. 44, subscript ℓ dropped):
  ```
  ψ̇ = π
  π̇ = ψ'' + 2imΩ(ψ' - π) - χ' - imΩχ - γ(π + imΩψ)
  χ̇ = γ(ψ' - χ + imΩψ) - imΩχ
  ```
  matches what the old build already had.
- **Sec. III.F (boundaries):** hard `ψ=0=π` at *both* edges, no separate
  outer sponge — relies on PML attenuation (left) and bound-state exponential
  falloff at `r_max ~ 500M` (right, `Mν ~ (M/r_max)^4`). Simpler than the old
  build's extra `gout` knob; drop that.

Still needed from the paper: Sec. III.A-C (Eq. 20/25 derivation, the `c^m_jl`
coupling coefficients — we have these from the old build's README, should
re-derive/confirm from the paper directly), and Sec. IV for his own benchmark
numbers to validate against (in addition to the Leaver cross-check we already
trust).

## Plan (small steps, in order)

1. Grid + single decoupled channel (ℓ=1, no tridiagonal solve, no PML yet) —
   just the interior equation and stencils on Dolan's exact grid convention.
2. Add the PML (Eq. 44) and hard boundaries (III.F). Validate against the
   known Leaver damped benchmark (a=0.9, α=0.5, Γ=-5.205e-4) as a first
   checkpoint — this re-derives the old single-mode result but under Dolan's
   exact discretization.
3. Add the ℓ=3 channel and the tridiagonal `Aẋ=b` solve (Thomas algorithm).
   Check convergence of Γ as the tower is truncated at different ℓ_max.
4. Add ℓ=5, repeat the convergence check.
5. Compare against any explicit benchmark Dolan gives in Sec. IV, not just
   Leaver.

## Status

Nothing built yet. Old `radial1d/` staged for deletion (uncommitted). This
README is step 0.
