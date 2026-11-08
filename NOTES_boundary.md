# Outer boundary: integrated data for propagating components (branch `boundary-integration`)

To be merged into `NOTES.md` (Day 2 / Day 3).

### Done and verified
- **Design** (as agreed after `validation/outer_integration_check.jl`): for every Fourier component
  with Re ω_q ≥ α the outer endpoint data are no longer the truncated 1/r series at the mesh endpoint.
  The series is evaluated at a far radius where it converges, 2 k r_far ≥ |η|² + 40 with
  η = −(2ω² − α²)/k (r_far ≈ 1370 for the (−1, 2) wave at a = 0.99, α = 0.3; last series term
  4e-20 there at order 40), and the **linear** outer system (full Kerr B, a²(ω²−α²) C mixing, no
  nonlinear source) is integrated inward to the endpoint. The leading outer coefficients stay the
  unknowns; layout, amplitude definition, fluxes and `charge_jacobian` are untouched.
- **Envelope integration**: R = e^{ikr} r^ν T(r) turns the 300 M-wavelength oscillation into a
  slowly varying T, integrated with 3-stage Gauss–Legendre steps (order 6, A-stable, one dense
  6×6 solve per step for one channel). Plain RK4 on R needed 15000+ steps; RK4 on T needs ≈ 800 for
  1e-10, Gauss order 6 needs ≈ 100 for 1e-11 (scratch comparison against a 40000-step reference
  integration converged to 5e-12). Default 48 steps per wavelength (158 steps from r = 1371 to 400),
  ≈ 4 ms per evaluation, ForwardDiff-transparent (AD gradient vs finite difference 1e-6–1e-5, the
  finite-difference error). Step count and r_far come from `ForwardDiff.value`s only.
- **Settings**: `OuterBoundary(; data=:integrated | :series, include_bound=false, convergence_margin=40,
  far_order=40, steps_per_wavelength=48)`, stored in `RadialProblem.outer_boundary` and carried by
  `with_background` / `with_potential`; `RadialProblem(...; outer_data=:series)` is the shorthand
  that reproduces the old behaviour exactly (the series code path is unchanged).
  `include_bound=true` also integrates bound components (stable inward; for them r_far is usually
  inside the domain, so this amounts to the order-40 series at the endpoint).
- **Default is `:integrated`.** It does shift existing numbers by more than 1e-10, because it
  corrects the series error: at a = 0.9, r* ∈ [−80, 180] (default mesh) the linear-point A′₂₁₁ moves
  by 8.5e-7 relative (1.4778702e-14 → 1.4778689e-14), and at A = (1e-2, 1e-4) A′₂₁₁ changes sign,
  −1.13e-9 → +3.51e-9 (`validation/logs/background_derivatives_integrated.log` vs
  `background_derivatives.log`). **The day-2 long fixed-background evolution (a = 0.9, r*_max = 180)
  was run with the series boundary and should be rerun** before its saturation numbers are quoted.
- `domain_estimate(...; integrated_outer=true)` (default): `rstar_max` from containment alone,
  `capped = false`; the new field `far_radius` reports where the integration starts, and a warning
  fires when it exceeds 20 `max_radius` (near threshold; cost ∝ 1/k², k → 0 still needs a different
  expansion). `integrated_outer=false` keeps the old sizing (and its tests).
- **Tests**: 107/107 pass (84 old + 23 new; `validation/logs/runtests_boundary.log`). New
  "Integrated outer data" set: (a) far-start integration over a dozen wavelengths reproduces a
  converged series at the endpoint to 1e-11, and r_far ≤ endpoint short-circuits to the series;
  (b) the Gauss stepper reproduces the Tricomi-U Coulomb function H⁺₃(η = −4.28, ρ) from ρ = 60 to
  8.5 to 1e-9 and shows order-6 error ratios; (c) r_far / step / far-order independence to 1e-10 for
  the (−1, 2) component; (d) frequency AD vs finite differences; a bound-component integration
  is step-independent and within 1e-3 of the J_B = 7 series at κr ≈ 12. Domain-estimate tests extended.
- **Dual background**: `validation/background_derivatives.jl` rerun with the integrated default:
  IFT vs four-re-solve finite differences still agree to 1e-8 at A = (1e-2, 1e-4)
  (`validation/logs/background_derivatives_integrated.log`).

### Acceptance: endpoint table at a = 0.99, α = 0.3, E₁ = E₂ = 650 (cosine)
`validation/boundary_endpoint_table.jl`, log in `validation/logs/`; stretched meshes
(r*_min = −100, step 0.4 to 150, ×1.03, max 3.7), amplitudes recalibrated per row.

| outer data | r*_max | Im ω₂₁₁ | Im ω₃₂₂ | L_E | L_J |
| --- | --- | --- | --- | --- | --- |
| series (old) | 200 | −4.638e-6 | 1.537e-6 | 4.029e-3 | 6.710e-3 |
| series (old) | 300 | −4.069e-6 | 0.970e-6 | 4.021e-3 | 9.152e-3 |
| series (old) | 400 | −3.902e-6 | 0.625e-6 | 4.237e-3 | 1.129e-2 |
| integrated | 200 | −3.998e-6 | 0.5308e-6 | 4.507e-3 | 1.275e-2 |
| integrated | 300 | −3.8905e-6 | 0.5952e-6 | 4.272e-3 | 1.162e-2 |
| integrated | 400 | −3.8887e-6 | 0.5962e-6 | 4.257e-3 | 1.148e-2 |

- 300 → 400: Im ω₃₂₂ changes by 0.17 % (series: 36 %), Im ω₂₁₁ by 0.05 %, L_E 0.4 %, L_J 1.2 %.
  The best converged secondary rate is **Im ω₃₂₂ ≈ 5.96e-7**, i.e. even the old r*_max = 400 value
  was 5 % high and the r*_max = 200 value 2.6× too high.
- (J_H, J_B) = (3, 7) → (5, 9) with the integrated boundary: identical to all printed digits
  (series boundary: −3.0 %). Doubling the step density and margin 40 → 100: identical.
- `include_bound=true` gives identical rows, so the remaining 200 → 300 difference (11 % in Im ω₃₂₂)
  is not series truncation of the bound components. It is the |322⟩ cloud not being contained at
  r* = 200 (12 % of its peak), where the boundary data neglect the nonlinear source: the
  domain must still reach the containment radius (≈ 500 at α = 0.3), which `domain_estimate` now
  asks for.
- Cost: residual 0.07–0.13 s for both boundary types (no measurable change); full Jacobian
  0.9 s vs 0.6 s, the difference being ForwardDiff through the integration in the outer
  boundary block (see open issues).
- The r*_max = 400 calibrations log several "did not converge" warnings from intermediate secant
  steps (|F| stalling at 1e-7 – 1e-4: the normalization floor of day 2, |211⟩ ≈ 1e-8 of its peak
  near r ≈ 450), with both boundary types. The **final solutions are converged**
  (`validation/boundary_convergence_check.jl`): |F| = 1.2e-15 (integrated) and 2.5e-15 (series) after
  3 iterations, a re-solve with a fresh factorisation moves ω by 1e-20, and E = 650.000 ± 5e-4.

### Open issues
1. **Jacobian cost of the boundary block**: ForwardDiff carries all node-block and parameter
   partials through the integration: outer block 0.26 s vs 0.08 s with the series (r*_max = 400 mesh). The integration is linear in the leading coefficients, so a
   transfer matrix with only the two frequency partials (plus the background's, when dual) would
   remove most of the 0.3 s. Not done.
2. Threshold k → 0 along an evolution: r_far ∝ 1/k³, steps ∝ 1/k²; fine to k ≈ 3e-3, beyond that a
   different expansion (outline footnote) is needed. The warning in `domain_estimate` flags it.
3. The normalization floor near r ≈ 450 still limits how far the mesh can go; with integrated data
   it no longer needs to go further than containment, so this is now benign at α = 0.3.
4. Long a = 0.9 evolution and any rates quoted from r*_max ≤ 200 meshes should be recomputed.
