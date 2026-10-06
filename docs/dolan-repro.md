# Reproducing Dolan 1212.1477 Sec. III from scratch

## Current entry points — 2026-09-24

The solver is packaged under `src/dolan/`, with independent benchmarks in
`src/common/`. The original reproduction branch and extraction tasks below
are historical. Working conventions are in [CLAUDE.md](../CLAUDE.md);
current project status and priorities are in [HANDOFF.md](../HANDOFF.md).

- [Dolan Fig. 13 campaign](../runs/dolan/fig13/README.md): peak runs and
  alpha/resolution sweeps, including the `t=10^6 M` array recipe.
- `PYTHONPATH=src python -m dolan.fig13_growth_rates`: Leaver/eigenvalue plot
  with available time-domain points (default resolution `M/6`).
- `PYTHONPATH=src python -m dolan.analyze_timedomain_run <run.npz>`:
  frequency-filter analysis of a saved record.
- [Nonlinear progress](nonlinear-progress.md) and
  [performance experiments](nonlinear-performance.md): positive odd-m source
  and synchronized RK4 checks, including the optional Julia benchmark bridge.
- [Adiabatic validation record](adiabatic-repro.md): the subsequent stationary-profile
  approach to slow cloud evolution.

## Original linear reproduction record

The following design notes and dated measurements preserve the original
validation history. The September 8 next-step list is not the current project
backlog; use the campaign recipes above for current paths and run settings.

## Why this branch exists

The prior `radial1d/` build (single-mode ℓ=m=1 truncation, on branch
`feat/eigenmode-seed-and-backreaction`) got the right answer but drifted from
actually reproducing Dolan's method: it approximated away the coupled-ℓ tower
(Sec. III, Eq. 25/38) rather than implementing it, and picked up backreaction/
nonlinear/multilevel extensions before the foundation was faithfully rebuilt.
Decision: start over, reproduce Dolan's Sec. III (equations, discretization,
PML, boundaries) exactly, validate against his own numbers, and only then
revisit backreaction/detuning/multilevel on top of the faithful version.

The predecessor lived on `feat/eigenmode-seed-and-backreaction` in the
GRChombo fork. This repository was subsequently extracted as a standalone
project; it does not require GRChombo.

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

## How the PML fits the pipeline (design notes)

**State = three fields**, each shaped `(n_channels, n_points)`:
- `ψ` — the scalar field per channel. The only observable; everything else is
  machinery to evolve it.
- `π ≡ ∂_t ψ` — its "velocity." Exists because the wave equation is 2nd-order
  in time; we split into two 1st-order eqs (`ψ̇ = π` and an eq for `π̇`).
- `χ` — the PML auxiliary field. Pure bookkeeping, nonzero only in the PML
  strip, never an observable. It's the local-in-time "memory" variable that
  makes complex-coordinate stretching work.

**Everything is Method of Lines:** one RHS function maps `(ψ,π,χ) → (ψ̇,π̇,χ̇)`,
and RK4 marches all three together. Inside one RHS call:
1. stencils compute `ψ'`, `ψ''`, `χ'` (globally, 4th-order Eq. 37);
2. `ψ̇ = π` — free, no solve; *this* is the equation that advances `ψ`;
3. `π̇` = `thomas_solve(A, b)` (full Eq. 25, channel-coupled) **plus** the PML
   source `−χ' − imΩχ − γ(π + imΩψ)` added on top *only where `γ>0`*;
4. `χ̇ = γ(ψ' − χ + imΩψ) − imΩχ` — its own algebraic formula, auto-zero where
   `γ=0`.

So the Thomas solve and the PML don't "combine" — they live in **different
parts of the same `π̇`**: Thomas gives `π̇` in the coupled interior, the PML
terms dominate `π̇` in the decoupled absorbing strip. `χ` only feeds `π̇`
through `−χ'−imΩχ`. RK4 steps `(ψ,π,χ)` forward at once; we read `ψ`, discard `χ`.

**Two design decisions that make this clean:**
- **The PML does NOT enter the Thomas solve.** The channel coupling *is* a `~Δ`
  term (`A` off-diagonals `= a²Δ c_{l,l±2}`; `b` coupling `= a²μ²Δ c_{l,l±2}`),
  and the PML is placed where these terms are negligible. The implementation
  still applies the Thomas solve everywhere; the off-diagonal terms approach
  zero near the horizon.
  The PML is then added to that result, without an explicit domain split.
- **No two-region solver / no seam matching.** As `Δ→0` the interior `π̇ = b/A`
  limits *exactly* onto Eq. (44) with `γ=χ=0` (since `A→(r²+a²)²`, the `ψ''`
  coeff `→1`, `2iam(r²+a²)/A→2imΩ`, `V_l→0`). So we run the full Eq. (25)
  everywhere and just *add* the PML source where `γ>0`. `γ(x)` itself is the
  switch (PML region ≡ `γ>0`); it must ramp **smoothly from exactly 0** at the
  interior edge `x_pml_start` (a hard jump reflects) — polynomial ramp
  `γ_max·((x_s−x)/(x_s−x_left))^p`, `p=2–3` — and `x_s` sits deep enough that
  `Δ` is already negligible, making the unit-`ψ''` normalization self-consistent
  (this *is* Dolan's "neglect `~Δ`" statement). `Ω_H = a/(r₊²+a²)`.

**Only change from the step-1 solver:** state grows `(ψ,π) → (ψ,π,χ)`; `rhs`
gains the `χ̇` equation and the PML terms on `π̇`. Stencils, Thomas solve, and
`ψ̇ = π` are untouched.

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
2. **[DONE]** PML (Eq. 44) as an additive source on the Eq. 25 solver
   (`x_pml_center`, `pml_width`, `pml_gamma`, `Omega_H`, three-field `(ψ,π,χ)`
   rhs). Dolan's geometry (Sec. III.E/F): a *symmetric bump* centered at x=−100M
   (width 10M) with a hard ψ=π=0 wall at x=−200M and a dead buffer between.
   Reduces exactly to the two-field solver in the interior. Deep placement forced
   the tortoise inversion to be recast in the horizon displacement `s=r−r₊`
   (log-space Newton), so it converges at x≲−100 instead of crashing.
3. **[DONE] Γ measurement.** `qbs_eigenmode.py` (`solve_qbs` — the QBS as an
   eigenvector of the linear generator) and `measure_gamma.py` (Dolan's
   *time-domain* method: evolve broadband Gaussian initial data, then extract Γ
   with the frequency filter — see the Status section and `frequency_filter.py`;
   this replaced the earlier eigenmode-seed/log-slope driver). Benchmark is
   `leaver.py` (Dolan-2007 continued fraction), evaluated at **runtime** via
   `leaver_benchmark(tower)` — no hardcoded constant.

## What we learned (corrects the earlier premise)

The organizing assumption above — "single channel lands close, the *tower*
closes the gap to Leaver" — **is wrong at a=0.9, α=0.5, and we verified it
directly.** Findings:

- **The ℓ-tower is negligible here.** The ℓ=3 channel populates at only ~4e-4
  of ℓ=1 (ℓ=5 ~0), so its effect on Γ is ~1e-7. `l_max`=1, 3, 5 give *identical*
  Γ to 5 sig figs. The coupling is correctly implemented — it is just physically
  tiny at this spin/coupling (weak spheroidal mixing, `aω≈0.43`).
- **There was no gap: the eigenvalue Γ matches Leaver to ~0.01%.** The earlier
  "~1% short of Leaver" compared against the *flat* Λ=l(l+1) continued fraction.
  Our operator carries the spheroidal angular correction (the `a²μ²cos²θ` term),
  so the apples-to-apples benchmark is `leaver.py` with `spheroidal=True`.
  Against that, the **eigenvalue** Γ=Im(ω) from `solve_qbs` reproduces *both* ω_R
  and Γ to ~0.01% (a=0.9, α=0.5), and is *flat* against PML width, PML reflection,
  x_right, and dx — a converged, robust operator property (verified 2026-09-08
  with the deep-bump PML). The old "push `x_left` deeper → 0.95→0.99" trend was
  the shallow PML eating the near-horizon mode; the deep bump removes it.
- **The time-domain *fit* is separately PML-contaminated.** The fitted slope Γ
  runs ~2.5% high vs. spheroidal Leaver and *drifts with PML width* — a limitation
  of the slope-fit extraction, not the operator. For rates, prefer the eigenvalue
  method; a faithful time-domain extraction needs Dolan's frequency-filter method
  (Sec. IV A 2), not a plain log-slope.
- **We still need the tower** — just not as the Leaver gap-closer here. It is
  the faithful equation, and it is expected to matter at higher spin/coupling
  (stronger spheroidal mixing) and for multi-mode/nonlinear work. Confirming
  where it becomes non-negligible is a live thread.

## Status

**Dolan Fig. 13 reproduced by BOTH methods** (2026-09-08). The linear
reproduction is complete; what remains is polish and the eventual nonlinear step.

**Tier A — eigenvalue Fig. 13: DONE.** `fig13_growth_rates.py` sweeps
α∈[0.30, 0.44] at a=0.99 and overlays our `solve_qbs` eigenvalue points on
spheroidal Leaver for *both* the n=0 fundamental and n=1 overtone. The points
sit on the Leaver curves through the peak (α≈0.42, Mν≈1.5e-7), the rollover, and
the n=1 sign change past α≈0.44 (`plots/fig13_growth_rates.png`). `solve_qbs`
gained a complex `omega_guess` (to select the right overtone) and a `maxiter`
cap (bad points raise `ArpackNoConvergence` fast instead of spinning). The
eigenvalue operator also reproduces full-spheroidal Leaver (ω_R *and* Γ) to
~0.01% at a=0.9, α=0.5, flat against every numerical knob (PML width/reflection,
x_right, dx).

**Tier B — time-domain frequency-filter: DONE (proof point).** This is the
method that, unlike the eigenvalue solve, *extends to the nonlinear problem*.
A checkpointed 10⁵ M evolution of Dolan's Gaussian initial data (Eq. 45) ran at
the a=0.99, α=0.42 peak (job 25713497, lr6, ~1h07m; recipe in `runs/README.md`).
- The raw probe spectrum is a textbook Dolan result: n=0 and n=1 peaks land
  *exactly* on Leaver, with a resolved tower of higher overtones.
- Dolan's frequency filter (`frequency_filter.py`, Eq. 46/49) recovers the
  growth rate: **n=0 Γ = 1.49e-7 vs Leaver 1.50e-7 (ratio 0.992, 0.8%)**;
  n=1 Γ = 4.68e-8 vs 5.36e-8 (ratio 0.873, 13%).
- Reproduce the spectrum + Γ table from the saved probe record with
  `analyze_timedomain_run.py <run.npz>`.

This supersedes the old plain log-slope extraction (which was PML-contaminated,
~2.5% and width-dependent) — the frequency filter is Dolan's actual Sec. IV A 2
method and is robust.

**The ℓ-tower stays negligible here** (ℓ=3 populates at ~4e-4 of ℓ=1), as at
a=0.9; `l_max`=1,3,5 give identical Γ. It is the faithful equation and is
expected to matter at higher spin/coupling and for multi-mode/nonlinear work.

**Immediate next steps (agreed order, 2026-09-08):**
1. **[DONE] Consolidate & document** — this Status rewrite, plus `runs/README.md`
   (run recipe) and `analyze_timedomain_run.py` (rerunnable spectrum + Γ table).
2. **Diagnose the n=1 13% deficit** — likely frequency-filter bleed (the n=0/n=1
   peaks are only ~1.5 FFT bins apart at Δω≈6.3e-5) or the fit window. Cheap
   knobs: `filter_bin_count`, `central_fraction`, probe position, longer record.
3. **Full time-domain Fig. 13 curve** — a SLURM array of ~15 α points at 10⁵ M,
   extracting Γ(α) time-domain for both overtones, to confirm the filter method
   tracks Leaver across α (not just at the peak). ~15×1h.

After that, the time-domain extractor is validated and the path opens to the
actual science: turn on the `|ψ|²ψ` self-interaction and measure saturation.

**Future upgrade path (parked):** hyperboloidal slicing replaces the PML+wall
entirely — a real coordinate change to slices reaching the horizon and scri⁺,
exactly reflectionless (no γ tuning, no reflection floor) and QBS-regular; now
the community default for QNM/QBS/superradiance time-domain work (Zenginoğlu
2008/2011; Panosso Macedo 2020; review arXiv:2409.11478). Worth it only if the
PML placement gets fiddly — not a mid-reproduction switch.
