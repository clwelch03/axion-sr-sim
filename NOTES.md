# Progress notes

## Day 3 (2026-11-07, evening session)

Repository initialised at the project root (`main`); the outer-boundary work runs in parallel on the
`boundary-integration` branch (worktree `../can-codex-boundary-worktree`), see `NOTES_boundary.md` once merged.

### Done and verified (main)
- **Faster source Jacobian** (`field.jl`): the phase averages ⟨e^{iq·ξ} h′(Φ) 2cos/sin(q′·ξ)⟩ are
  ĥ′(q+q′) ± ĥ′(q−q′), so all of them come from one real BLAS product of a precomputed averaging table
  (one row per distinct ±k; 25 rows for the two-mode set) with h′(Φ)[point, polar node]. Matches AD to 2e-14
  for cosine, quartic (n_angular = 2) and chiral. Per node 275 μs → 61 μs; production Jacobian (751 nodes)
  0.58 s → 0.21 s; a Newton step is now dominated by the sparse LU (≈ 0.35 s).
- **Integrator tuning** (`evolution.jl`, `validation/evolution_step_control.jl`): the day-2 explanation was
  wrong: the Newton-tolerance noise in d ln A/dt is 1e-12 of Im ω₂₁₁, negligible. The waste came from Tsit5's
  initial-step heuristic (it picked dt ≈ 1, then grew ×10 per step: ≈ 40 evaluations to reach 1e6), and from
  integrating an exponential in A (steps limited to γ dt ≈ 0.1). Now the state is ln A by default
  (`log_amplitudes=true`, reltol = abstol = 1e-7) and the first step is 0.05 / max|A′/A| (`initial_step_fraction`).
  Three linear e-folds: 141 → 103 (first step only) → **31** evaluations, identical accuracy. The 5.5e-5
  remaining error in all three runs is not integration error but the energy-flux bias below.
  `amplitude_trajectory(result)` reads amplitudes back in either convention (scripts updated).
- **Nonlinear boundary terms** (`validation/nonlinear_boundary_terms.jl`): local ratio |N|/|B R| per channel.
  Horizon endpoints: ≤ 4e-7 (a = 0.99 joint E = 650 point), ≤ 4e-9 (a = 0.9 peak and saturated states).
  Outer endpoint (r* = 200 resp. 180): fundamentals 3e-5–6e-4 (211) and 3e-4–4e-3 (322); second-order
  components (2,−1) 4e-3–6e-3 and (−1,2) 1e-3–4.5e-3; the third-order components (3,0), (2,1), (1,2), (0,3)
  are O(1) — they are entirely source-driven there (and carry ≤ 1e-12 of the flux). So the outline's
  "sufficiently small" holds at the few-1e-3 level for every component that matters, not uniformly.
- **Wave actions and the action balance** (`charges.jl`, `validation/action_balance.jl`): N_i = ∂⟨L⟩/∂ω_i =
  Σ_q q_i n_q with n_q the Klein–Gordon charge of component q. Verified to round-off at linear and nonlinear
  points: J = Σ m_i N_i, L_E = Σ ω_i L_{N_i}, L_J = Σ m_i L_{N_i}. E − Σ ω_i N_i is the outer surface term
  (plus potential energy): −3.6e-4 E at r*_max = 180, −4e-5 at 250, −3e-6 at 320 (linear point).
  **Energy-flux bias found**: with a complex trial frequency the energy flux is
  4π(r²+a²) Σ [ω Im(R*Q) − Im ω̃ Re(R*Q)]; `outward_flux` drops the second term, which is the same order as
  the losses wherever the cloud is not contained. Linear point, default mesh: A′/(γA) − 1 =
  (7.4e-5, −4.2e-4) with the energy balance as coded, (−1.5e-6, −1.2e-5) with the complete flux, and
  identically (−1.5e-6, −1.2e-5) with the action balance, whose flux contains no time derivative. At the
  nonlinear point (1e-2, 1e-4) the three balances agree to 1e-3. Not changed in `outward_flux` (the γ that
  belongs there in the nonlinear adiabatic setting is ambiguous); the action balance sidesteps it.
  `amplitude_rates(...; balance=:action)` and `evolve_*(...; balance=:action)` implement it for any
  number of modes; the default stays `:energy` so all existing results reproduce.
- **Metric work, settled** (`evolution.jl`, `validation/metric_work_containment.jl`): actions are adiabatic
  invariants, so dE = Σ ω_i dN_i + w with
  w = (∂_b E − Σ ω_i ∂_b N_i)|_A b′ + (∂_j E − Σ ω_i ∂_j N_i)|_A j′, using the same fixed-A, fixed-r*/M
  derivatives as the drift terms (one extra dual pass); the J-counterpart vanishes identically (J = Σ m N),
  matching the outline's (w, 0). Independent of r*_min to 7 digits (−4.011009e-3 E b′ for r*_min = −45 … −120,
  where the Boyer–Lindquist version runs −1e-2 … +5e5). Equals the linear-regime prediction Σ N_i dω_i/dτ once
  the cloud is contained: w/(E b′) = −4.01e-3, −9.52e-3, −1.0561e-2, −1.0660e-2 vs −1.0664e-2 for
  r*_max = 180, 250, 320, 400; and there energy balance + w agrees with the action balance to 1e-7
  (without w: 1e-4). `include_metric_work=true` now uses this w; still off by default.
  `boyer_lindquist_metric_work` kept for comparison. `BackgroundDerivatives` (struct) replaces the pair
  returned by `background_charge_derivatives`, adding action and frequency derivatives.
- `Pkg.test()`: 95/95 (≈ 3.5 min; new sets for actions/balances and a short log-amplitude evolution).

### Open (Catherine)
1. **Switch the evolution to the action balance?** It is exact for the region, needs no metric work, removes
   the two-mode restriction (outline "Future work"), and avoids the energy-flux bias. Energy and angular
   momentum remain available as diagnostics. My recommendation: yes, once the boundary branch is merged
   (contained domains make the two formulations agree to 1e-7 anyway).
2. The drift terms ∂_b Q|_A are large on non-contained domains (the fixed outer coefficient rescales the whole
   tail as κ changes); with contained domains this is benign. Same remedy as above.

## Day 2 (2026-11-07, afternoon session)

### Done and verified
- **Analytic Jacobian** replaces per-node ForwardDiff (`rhs_derivatives` in `bvp.jl`,
  `nonlinear_source_jacobian` in `field.jl`, `nonlinearity_derivative` in `potential.jl`).
  Matches AD to 1e-16 for all four potentials, with `n_angular = 2` and the frequency columns
  (new test set; `rhs_derivatives_autodiff` kept as the reference). Node right-hand sides are
  now evaluated once per residual/Jacobian instead of twice.
  Production two-mode problem (751 nodes, 24064 unknowns): Jacobian 7.3 s → 0.7 s,
  residual 0.12 s → 0.07 s, sparse LU 0.33 s. Newton iteration ≈ 1.1 s (was ≈ 8 s).
- `Pkg.test()`: 84/84 pass (≈ 2 min) after all of today's changes (the new non-uniform Simpson test caught a wrong trailing-interval weight, fixed; the domain-estimate and singularity-guard changes are covered).
- **Yesterday's "first solve did not converge" was a misreading**: the warning came from the
  A = 1e-3 step (stderr is unbuffered). The hydrogenic seed converges in 6–8 Newton steps.
- **Two-mode cosine solve at a = 0.99, α = 0.3, E₁ = E₂ = 650, r* ∈ [−100, 200], step 0.4**
  converges: A = (4.209e-2, 1.305e-4), ω₂₁₁ = 0.295948 − 4.64e-6 i, ω₃₂₂ = 0.298350 + 1.54e-6 i,
  peak |R| = (0.123, 0.035), E_total = 1302 (fundamental free-field energies 650 each).
  → the primary *decays* and the secondary grows in the joint solution.
  Linear (free) two-mode solve at the same point: ω₂₁₁ = 0.296304 + 2.68e-8 i, ω₃₂₂ = 0.298452 + 1.80e-12 i
  (both positive; yesterday's negative Im ω₃₂₂ was from the ramped-amplitude run, not the linear problem).
- Single |211⟩ cosine at E = 650: A = 0.03629 (linear estimate 0.0264 → E/target 0.66: strong deformation).
- **Fixed-background evolution, linear regime** (a = 0.9, default mesh, A = (1e-6, 1e-7)):
  `amplitude_rates` = γA to 7e-5 (211) and 4e-4 (322). Tsit5 over three e-folds of |211⟩
  (t ≤ 2.03e8) reproduces A = A₀e^{γt} to 5e-5 relative: 141 rhs evaluations in 744 s.
  Warm timings: `charge_jacobian` 0.13 s, `cloud_charges` 0.16 s, warm solve ≈ 4 s
  (the 14 s seen earlier was compilation).
- **Chord iteration** in `solve_radial` (`reuse_factorization=true`, `chord_threshold=1e-3`):
  warm-started solves reuse the previous factorised Jacobian while each step halves the residual,
  and refactorise once at the end for the sensitivities. Matches full Newton to 1e-12 in dA/dt
  along an amplitude sequence; saves 1–2 s per solve (`scratchpad/chord_check*.log`).
- `amplitudes_for_energy` now starts from a linear estimate (one free-field solve, rescaled per mode),
  iterates with a per-mode secant on d ln E/d ln A (≈ 1.3–1.6 for deformed clouds, 2 for linear),
  backs off when a step fails to converge, and warm-starts every step. Default `tol = 1e-6`
  (the diagonal secant oscillates at the 1e-8 level for the coupled two-mode problem).
- `with_potential(problem, potential)`, `fundamental_energies(solution)` added and exported;
  `verbose=true` on the evolution drivers prints (and flushes) one line per rhs evaluation.
- **Comparison plot reproduced**: `scripts/profile_comparison.png` (a = 0.99, α = 0.3, E₁ = E₂ = 650,
  r* ∈ [−100, 200]). Primary single-cloud and joint profiles nearly coincide (taller, narrower
  than linear); the joint secondary is pulled inward and peaks ≈ 10 % higher than alone.
  Amplitudes: linear (2.639e-2, 8.43e-5), single (3.629e-2 | 9.62e-5), joint (4.209e-2, 1.305e-4).
- **Endpoint sensitivity of the joint rates** (same point, step 0.4):

  | r* domain | Im ω₂₁₁ | Im ω₃₂₂ | L_E | L_J |
  | --- | --- | --- | --- | --- |
  | [−100, 200] | −4.64e-6 | 1.54e-6 | 0.00403 | 0.00671 |
  | [−100, 300] | −4.07e-6 | 0.97e-6 | 0.00402 | 0.00915 |
  | [−150, 300] | −4.07e-6 | 0.97e-6 | 0.00402 | 0.00915 |
  | [−100, 400] | −3.90e-6 | 0.62e-6 | 0.00424 | 0.01129 |

  Horizon endpoint converged (−100 vs −150: 1e-8 level). Outer endpoint **not** converged.
  **Confirmed by the per-component diagnostic** (`scratchpad/components.log`, r* ∈ [−100, 200]):
  - q = (−1, 2), m = 3: ω_q = 2ω₃₂₂ − ω₂₁₁ = 0.30075 > α, k = 0.0213 (wavelength ≈ 300 M),
    ν = −0.98 + 4.28i. It carries the dominant *outgoing* energy flux, F_E = 6.4e-4. Its outgoing
    series needs k r ≫ |ν| (r ≫ 200), so none of the uniform domains is in its asymptotic region;
    L_J grows 0.0067 → 0.0092 → 0.0113 with r*_max. Hydrogenic estimate 2ω₃ − ω₂ − α = α³/72 > 0:
    this component is *always* just above threshold, and amplitude-driven frequency shifts (≈ 1e-3)
    can push it through k = 0 during an evolution (outline footnote on the k = 0 expansion).
  - q = (2, −1), m = 0: bound (κ = 0.062), peak |R| = 0.037, absorbed at the horizon with
    F_E(in) = −3.4e-3: the **dominant loss channel** (211 × 211 → 322 × BH).
  - |322⟩ itself: |R| at r* = 200 is 12 % of its peak (κ = 0.031 with the r² prefactor), so the
    cloud is not contained either; F_E(out) = 6.7e-5 through the boundary.
  - The order-3 propagating components (3ω) are negligible (fluxes ≤ 1e-12).
  → Implemented: stretched meshes (`RadialMesh(...; stretch_from, stretch_factor, max_step)`) with
  non-uniform composite Simpson weights in `charges.jl` (tested). Stretched r* ∈ [−100, 400]
  (step 0.4 to 150, ×1.03 per interval, max 4; 736 nodes) reproduces the uniform 1251-node
  result to all printed digits (ω, L, A). Domain convergence r*_max = 400 / 1000 / 2000 / 4000:
  `scratchpad/stretched.log`.
  Alternative if the series still misbehaves at small k: exact Coulomb-wave outer data (H⁺_l(η, kr))
  for propagating components instead of the 1/r series.
  **It does misbehave** (`validation/truncation_convergence.jl`, stretched r* ∈ [−100, 400]):
  (J_H, J_B) = (3, 7) → (5, 9) → (7, 11) shifts Im ω₃₂₂ by −3.0 % then −5.2 % (Im ω₂₁₁ by −0.2 %,
  −0.4 %), i.e. the outer expansion is asymptotic, not convergent, at k r ≈ 8.5 for the (−1, 2) wave.
  Pushing r*_max out (domain_convergence.log) or switching to Coulomb-wave data are the only fixes.
  **Pushing r*_max out does not work with the present normalization** (`validation/domain_outer_debug.jl`):
  on the stretched r* ∈ [−100, 1000] mesh the |211⟩ profile at the boundary is e^{−κr} ≈ 1e-21 of its
  peak, so the row imposing A₂₁₁ through the leading outer coefficient constrains nothing at the
  Newton tolerance; the system is singular, ω₂₁₁ comes out as 0.29487 + 1.7e-7 i (wrong) and the
  calibration estimate explodes (A₁ ≈ 4e7). The usable outer radius is where the fastest-decaying
  fundamental is still ≈ 1e-8 of its peak (r ≈ 450 at α = 0.3), now built into `domain_estimate`
  (`normalization_radius`, `normalization_floor=1e-8`) and `solve_radial` stops with a warning when
  the residual sits at the floor while the step does not vanish. Hence the (−1, 2) wave *cannot* be
  made asymptotic with the 1/r series and outer-coefficient amplitudes: the remedy must be
  **Coulomb-wave outer data for propagating components** (boundary at r ≈ 400–500), possibly with
  exponentially rescaled tail unknowns if even larger domains are ever needed. Alternatively define
  the amplitudes by the fundamental free-field energies (dense Jacobian rows) instead of the outer
  coefficients. Catherine's decision.
  **Coulomb-wave prototype works** (`validation/coulomb_prototype.jl`, log in `validation/logs/`):
  H⁺_l(η, ρ) = e^{iθ_l} (−2iρ)^{l+1+iη} U(l+1+iη, 2l+2, −2iρ) with Tricomi's U from the integral
  representation on a rotated contour and the substitution x = e^v (trapezoidal rule, exponentially
  convergent, h = 0.1). Agrees with the present asymptotic series to 1e-14 wherever the series
  converges (ρ ≥ 21), costs 21 μs, works under ForwardDiff. At ρ = k r = 8.5 (r* = 400 for the
  (−1, 2) wave) the series is off by 1.3e-3, which is the error behind the endpoint drift.
  Remaining work for production use: the Kerr O(1/r²) corrections (a², 2M, angular mixing) on
  top of the Coulomb leading factor, i.e. R = H⁺_l(η, kr)/r × Σ c_j r^{−j} with a recursion for c_j,
  then the boundary residual and its AD. Note on the outline's "converged by r* = 180": true at a
  point with a weak secondary, where the (−1, 2) flux is negligible; at E₁ = E₂ = 650 it is the
  dominant outer flux and nothing below r* ≈ 400 is converged.
- **Spindown balance terms** (a = 0.9, A = (1e-6, 1e-7), `scratchpad/spindown_check.log`):
  `balance_terms` costs ≈ 45 s (four re-solves); the background drift ∂_b Q b′ + ∂_j Q j′ is ≈ 1 % of
  the losses and shifts A′₂₁₁ by 0.5 %. `metric_work` was finite at r*_min = −80 (w ≈ 2.5 E b′) but
  NaN at −120: Δ was formed from r (catastrophic cancellation at s ≈ 1e-17). Rewritten in terms of
  the horizon offset, which shows the real problem: **w diverges like 1/s towards the horizon**
  (`scratchpad/metric_work_check.log`: w/(E b′) = −0.006, 2.5, 1.1e3, 4.6e5, 8.7e10 for
  r*_min = −60, −80, −100, −120, −160). ∂_τ g^{tt} at fixed Boyer–Lindquist r goes like
  ∂_τ(1/Δ) ∝ 1/Δ² and the measure supplies only one Δ. The BL-fixed convention is simply the wrong
  chart at a moving horizon: the derivative should be taken at fixed horizon-comoving coordinate
  (r*/M, the same convention as ∂_b E), which introduces ∂_τ r|_{r*} cross terms. Not done; the term
  is now **optional and off by default** (`include_metric_work=false`), leaving the leading-order
  balance [∂_A E] A′ = −L − ∂_b E b′ − ∂_j E j′ in which the cloud energy changes only by fluxes.
  **Physics gap fixed**: the hole's mass/spin rates were the horizon fluxes in f_a² units, i.e.
  f_a = M_Pl. `evolve_with_spindown` now takes `decay_constant_ratio = f_a/M_Pl` (non-reduced) and
  multiplies the horizon fluxes by its square (`EvolutionCache.coupling_to_hole`). Consequently the
  background-drift and metric-work corrections are O((f_a/M_Pl)² E) relative to the losses.
  **Open formalism point**: ∂_b E is taken at fixed r*/M mesh while w is ∂_τ g at fixed BL r;
  the conventions must match for the balance to be exact (outline calls this part schematic).
- **Spindown run end to end** (`scratchpad/spindown_short.log`: a = 0.9, A = (1e-6, 1e-7),
  f_a/M_Pl = 1e-2, t ≤ 2e7): 57 rhs evaluations in 1092 s (19 s each; `balance_terms` re-solves four
  times). A₂₁₁(2e7) = 1.343891e-6 vs fixed-background 1.343868e-6 (2e-5 apart, the background drift);
  1 − b = 7.9e-11 = (f_a/M_Pl)² × ΔE_cloud, j drops by 2.9e-10.
- **Exact background derivatives** (`background_charge_derivatives`, default in `balance_terms`;
  `resolve_background=true` keeps the four-re-solve finite differences for validation):
  dQ_ref/db = ∂Q/∂b|_u − (∂Q/∂u) J⁻¹ ∂F/∂b|_u with ∂F/∂b by ForwardDiff through the background,
  mass coupling and mesh offsets (`RadialProblem`/`RadialMesh` now accept dual-number backgrounds;
  `with_background` builds them). A central difference for ∂F/∂b was 1 % off at relative step 1e-5
  (clean step² error): the |322⟩ outer-boundary rows vary as e^{−κr} with d(κr)/dα ≈ 1800, and J⁻¹
  maps their error straight onto the profile near the |322⟩ peak. With the exact partial, IFT and
  re-solve agree to 1e-7 (linear point) and 2e-8 (A = (1e-2, 1e-4)); `balance_terms` 0.6 s vs 30 s
  (`validation/background_derivatives.jl`, log in `validation/logs/`). The short spindown run now
  takes 222 s for 57 rhs evaluations (3.9 s each, was 19 s) with identical output.
- **Radial step and quadrature converged**: halving the step from 0.4 to 0.2 on the stretched a = 0.99
  mesh changes the joint rates by 3e-8 relative; n_polar 12 → 16 and n_phase 16 → 24 by 1e-12
  (`validation/logs/truncation_convergence.log`). Only the outer series/endpoint issue remains.

- **First long fixed-background evolution** (a = 0.9, α = 0.3, cosine, default mesh r* ∈ [−80, 180],
  A₀ = (1e-6, 1e-7), t ≤ 1.5e9): 711 rhs evaluations in 3689 s (5.2 s each, pre-chord code).
  Data: `results/fixed_background_a0.9_alpha0.3_{diagnostics,solution}.csv`, plot
  `results/fixed_background_a0.9_alpha0.3.png`, script `scripts/evolution_fixed_background.jl`,
  plotting `scripts/plot_evolution.jl`.

  | t/M | A₂₁₁ | A₃₂₂ | E/f_a²M | Im ω₂₁₁ | Im ω₃₂₂ |
  | --- | --- | --- | --- | --- | --- |
  | 0 | 1.0e-6 | 1.0e-7 | 8.6e-4 | 1.48e-8 | 1.4e-12 |
  | 6.58e8 (peak of |211⟩) | 1.905e-2 | 9.6e-6 | 266 | 6.5e-10 | 2.8e-7 |
  | 1.5e9 (saturated) | 5.52e-3 | 3.44e-5 | 122 | 1.2e-12 | −2.0e-12 |

  |211⟩ grows linearly until E ≈ 100, saturates at E ≈ 266 (its eigen-rate suppressed 20×, the m = 0
  component being absorbed), then |322⟩ grows explosively (1e-7 → 3.8e-5 within ≈ 1.5e8) and takes
  half of the energy; afterwards the pair sits in a quasi-equilibrium with net L_E ≈ 3e-10 and
  frequency shifts of −2e-5 (211) and −7e-6 (322). This is the 211 × 211 → 322 × BH saturation picture.
  Caveats: default domain (r*_max = 180, the |322⟩ profile is ≈ 10 % of its peak there, and the
  near-threshold (−1, 2) wave is not resolved); Tsit5 with reltol 1e-6 took tiny steps at t ≲ 1e4
  (the |F| floor of the Newton solve injects 1e-13 noise into A′) — a larger `abstol` or `dtmin`
  would save ≈ 100 evaluations.

- **Physically motivated domains** (Catherine's suggestion): `domain_estimate(background, α, mode_set)`
  and `RadialMesh(background, α, mode_set; ...)` set the endpoints from hydrogenic estimates:
  r*_min where s = 1e-8 (r₊ − r₋) (Frobenius convergence radius), r*_max from the slowest bound
  component's containment (`r^ν e^{−κr}` below 1e-4 of its peak) and from the propagating
  components' asymptotic requirement `2 k r ≥ 10 |ν|²`, capped at r = 5000; stretching starts at
  twice the largest bound peak radius and the step is capped at a quarter wavelength of the fastest
  propagating component (`validation/domain_estimates.jl`):

  | a | α | r*_min | containment r | asymptotic r | nodes |
  | --- | --- | --- | --- | --- | --- |
  | 0.9 | 0.3 | −60 (ratio 1e-8; −45 at 1e-6) | 506 | 5000 (capped; needs ≈ 12000) | 3144 |
  | 0.99 | 0.3 | −150 (−110 at 1e-6) | 506 | capped | 3372 |
  | 0.99 | 0.2 | −150 | 1142 | capped | 2874 |
  | 0.999 | 0.3 | −434 | 506 | capped | 4081 |

  **Inner endpoint converged far earlier than assumed**: for the linear two-mode problem at a = 0.9
  the rates at r*_min = −30 (s/gap = 8e-5) agree with −120 to 1e-10 relative, so the default
  `horizon_offset_ratio` is 1e-6 (r*_min ≈ −45 at a = 0.9, −110 at a = 0.99) and the old −80/−100
  were overkill (the day-1 "convergence issues at the horizon" were presumably something else).
  The near-threshold (−1, 2) wave (hydrogenic k = α²/6 ≈ 0.015) dominates the outer requirement at
  every point and costs ≈ 4× the nodes of the default mesh. **Recommendation**: replace the 1/r series
  by exact Coulomb-wave outer data `H⁺_l(η, kr)` for propagating components; then the outer boundary
  only needs the containment radius (≈ 500 here) and the small-k problem disappears (the k → 0
  threshold still needs thought). Inner-endpoint convergence with the gap criterion and outer
  convergence on stretched meshes: `validation/domain_convergence.jl` → `validation/logs/domain_convergence.log`.

- **Outer-boundary options quantified on the Kerr equation** (`validation/outer_integration_check.jl`,
  log in `validation/logs/`; (−1, 2) component at the a = 0.99, α = 0.3 joint point, k = 0.0213, η = −4.28).
  Reference: the linear Kerr radial equation integrated inward (RK4 in r) from r = 1000–1500, where the
  1/r series is converged to 1e-14 (orders 25/35/45 agree); r_far = 1000 vs 1500 and step 0.04 vs 0.02
  agree to 2e-12 in R and 8e-15 in R'/R at r = 400. Cost 95 ms unoptimised (15000 steps; a few hundred
  would do for a 300 M wavelength).

  | outer data at r = 400 | |ΔR|/|R| | |Δ(R'/R)| |
  | --- | --- | --- |
  | series J_B = 7 (production) | 1.9e-2 | 3.8e-4 |
  | series J_B = 11 | 2.7e-3 | 7.8e-5 |
  | series J_B = 15 (optimal truncation) | 1.2e-3 | 4.6e-5 |
  | bare Coulomb H⁺₃(η, kr)/r, normalised at r_far = 1000 / 1500 / 3000 | 1.9e-2 / 2.4e-2 / 2.9e-2 | 7.1e-5 |

  So the bare Coulomb function is **not** a fix: the Kerr O(1/r²) terms (2amω, a²ω², the Δ expansion)
  accumulate a ≈ 2 % error between r = 400 and infinity, as large as the production series' error.
  The planned correction factor Σ c_j r^{−j} on top of H⁺ has no polynomial recursion either
  (substituting R = H⁺ T/r gives T'' + 2(H⁺'/H⁺) T' + δV T = 0, and H⁺'/H⁺ is not rational in 1/r).
  Exact outer data therefore means **numerical inward integration of the linear outer system** from the
  series-converged radius (2kr ≳ |η|² + 40) to the mesh endpoint, per propagating component; the series
  stays as the far data and the Coulomb prototype serves only as a unit test of the integrator.
  Decision pending (Catherine).

### Running / next steps
1. Decide the outer-boundary treatment for propagating components (Coulomb waves vs. energy
   normalization); until then production domains are limited to r ≈ 450 and the secondary's rate
   carries the endpoint systematic quantified above (≈ 40 % between r*_max = 200 and 400).
2. Spindown run of the same length with a physical f_a/M_Pl (needs the cheaper ∂_b, ∂_j derivatives first).
3. Settle the outer endpoint for the a = 0.99 secondary with stretched meshes (`stretched.log`);
   then convergence tables for (J_H, J_B), step and endpoints as in the outline's validation section.
   Decide (Catherine) how to treat the near-threshold q = (−1, 2) wave when k → 0 along an evolution.
4. Decide (Catherine) the metric-work convention; if wanted, derive w in the horizon-comoving chart.
5. Further speed if needed (measured on the 751-node two-mode problem, `validation/logs/assembly_profile.jl`):
   Jacobian 0.8 s = per-node source Jacobians ≈ 0.45 s (1500 evaluations of the phase-grid contraction;
   could be batched into one BLAS product per node) + interval products ≈ 0.1 s + sparse assembly 0.04 s;
   LU 0.4 s; residual 0.07 s. A preallocated sparse pattern is therefore *not* worth it.

### Scripts
The session scratchpad (and its logs) was wiped when the connection dropped, so the diagnostic
scripts now live in `validation/` (logs in `validation/logs/`); run them from the package root with
`julia --project=. validation/<name>.jl`. Earlier references to `scratchpad/*.log` below are to lost files;
the numbers quoted are preserved here.

#### Original scratchpad list (for the record)
- `continuation.jl`: linear → joint → single-mode solves at E = 650 for r*_max = 200 and 300.
- `evolution_linear.jl`: linear-regime evolution check.
- `check_jacobian.jl`, `profile_jacobian.jl`, `baseline.jl`: Jacobian verification and timing.
- Julia must be invoked by full path from non-interactive shells:
  `/global/common/software/nersc9/julia/1.11.7/bin/julia`.

## Day 1 (2026-11-07, end of day)

### Done and verified
- Full package skeleton in `src/` (see README for the file map), all with docstrings,
  descriptive ASCII identifiers (no Greek, no terse names), structs instead of tuple returns.
- Verified numerically (scratch scripts, now encoded in `test/runtests.jl`):
  - tortoise inversion exact to the horizon (`s ≈ 1.7e-11` at `r* = −80`);
  - polar harmonics normalised/orthogonal, `∂_θ` identity; mixing matrix `C^1_{11} = 1/5`;
  - nonlinear source matches the analytic cubic projection to all digits;
  - horizon and outer series satisfy the radial ODE; AD flows through them;
  - structured sparse Jacobian equals naive ForwardDiff Jacobian to `1e-14`;
  - single |211⟩ linear solve on the production mesh (a = 0.9, α = 0.3):
    ω = 0.2962955 + 1.478e-8 i, 6 Newton steps, warm restart converges in 1 step,
    `J/E = m/ω`, `2γE = −L_E` to 3e-6, `∂E/∂A = 2E/A` to 1e-13, mesh halving shifts γ by 1e-6.
  - Energy normalisation consistent with Catherine's plot: E ≈ 650 f_a²M ↔ peak |R| ≈ 0.12.

### Still not exercised
- `evolve_with_spindown`, `balance_terms`, `metric_work` (BL-fixed ∂_τ g^{αβ} via ForwardDiff;
  near-horizon cancellation argued benign but must be checked numerically).
- `scripts/profile_comparison.jl` (CairoMakie; env in `scripts/`).

### Amplitude scales (a = 0.99, α = 0.3)
- E ≈ 650 needs A₂₁₁ ≈ 4e-2 but A₃₂₂ ≈ 1.3e-4: the |322⟩ profile peaks at ≈ 600·A because of the
  r^{n−1} outer normalisation (κ₃₂₂ ≈ 0.03 so e^{−κr} ≈ e^{−6} at r* = 200 — outer-endpoint sensitivity).
- Ramping both amplitudes equally (yesterday) drove |322⟩ to E ≈ 11000, far outside the regime of interest.
