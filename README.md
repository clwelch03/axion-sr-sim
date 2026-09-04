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

**Design decision (changed from the original plan):** the code is *tower-native
from day one*, not "single channel now, tower bolted on later." One `Tower`
object holds all same-parity channels `l ∈ {|m|, |m|+2, …, l_max}`; the
tridiagonal `Aẋ=b` solve runs for every `l_max`. `l_max=|m|` is just the
one-rung case and goes through the *identical* `rhs` — so there is no
single-channel version to replace, and the ℓ-tower cannot get orphaned the way
it did last time. The build order below still validates at `l_max=1` first,
because that isolates the grid/stencil/PML machinery from the coupling with a
cheap decisive checkpoint — but that is a *validation* order, not a "defer the
hard part" order. The coupling code exists and is unit-tested from the start.

1. **[DONE]** `solver.py`: grid, Kerr tortoise `r(x)`, coupling matrix
   `c^m_{jl}` (via the `cosθ` ladder, no Clebsch–Gordan needed), Eq. (25)/(26)
   coefficients, 4th-order stencils (Eq. 37), Thomas solve, RK4/MoL. All
   unit-tested: stencils show exact 4th-order (16× per halving); Thomas matches
   `np.linalg.solve` to machine precision at `n_ch=1,3,5`; single-channel `rhs`
   exactly equals the plain `b/D(r)` divide; tower run leaks ℓ=1→ℓ=3 through
   the coupling as expected.
2. **[NEXT]** Add the PML (Eq. 44) near the left edge + hard boundaries
   (III.F). Validate at `l_max=1` against the Leaver damped benchmark (a=0.9,
   α=0.5, Γ=-5.205e-4). NOTE: single-channel will *not* hit Leaver exactly —
   Leaver's continued fraction already resums all ℓ-coupling, so single-channel
   lands close and the *tower* closes the residual gap. That gap closing is the
   real coupling validation (step 3).
3. Run `l_max=3`, then `5` — zero new code, pure runtime — and watch Γ converge
   toward Leaver. This convergence *is* the proof the coupling is right.
4. Compare against any explicit benchmark Dolan gives in Sec. IV, not just
   Leaver.

## Status

**Step 1 DONE, validated.** `solver.py` holds the full tower solver (see step 1
above for the test list). Currently hard Dirichlet at *both* edges, no PML yet,
so reflections contaminate any long run — the machinery is proven to propagate,
but no physics until the PML (step 2) is in. Next: PML (Eq. 44).
