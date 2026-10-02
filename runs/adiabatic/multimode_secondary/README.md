# Secondary cloud on the nonlinear primary (Omiya et al. 2211.01949, Sec. III)

The l = m = 2 secondary of Eq. (61), solved as an eigenproblem for ω₂ on each
saved point of the cosine primary branch
(`runs/adiabatic/fig2_20260922/results/cosine_branch.jls`, a/M = 0.99, Mμ = 0.42).
The secondary keeps (l, m₂) = (2,2), (4,2), (6,2), uses the primary's radial domain
(r* from −100 to 150, matching at 5) at twice its step (0.2), and fixes the (2,2)
outer amplitude to 1. Library code: `src/adiabatic/multimode/`; driver:
`src/adiabatic/drivers/continue_secondary.jl`.

## Running

```bash
sbatch runs/adiabatic/multimode_secondary/scripts/submit_secondary.sh
```

The driver saves after every point and resumes an existing output prefix; one
writer per prefix. `AXION_MAX_POINTS=n` stops after n primary points (local tests).

Plots (written to `plots/`, ignored):

```bash
julia --project=. runs/adiabatic/multimode_secondary/scripts/export_profiles.jl
PYTHONPATH=src python runs/adiabatic/multimode_secondary/scripts/plot_secondary.py
```

`multimode_secondary_frequency.png` shows ω₂ along the branch against the
small-A₁ limits; `multimode_secondary_environment.png` shows the equatorial
⟨1 − cos ϕ₁⟩ and the secondary's |R₂₂| at six primary amplitudes;
`multimode_secondary_growth_vs_energy.png` shows Γ₁ and Γ₂ against E₁; and
`multimode_secondary_fig7_analogue.png` compares nonlinear and linear |R₁₁|,
|R₂₂| at E₁ = E₂ ≈ 967 F_a²M (A₁ = 0.672), the analogue of the paper's Fig. 7
(which is at Mμ = 0.3).

## Results

| File | Background term | Status |
|---|---|---|
| `results/cosine_secondary_fixed_time.{jls,csv}` | fixed-time projection of (1 − cos ϕ₁)ϕ₂ | **not physical beyond small A₁**; kept as a record (job 26624675, 2026-10-01) |
| `results/cosine_secondary.{jls,csv}` | period average ⟨1 − cos ϕ₁⟩ | complete, 232 points (local points 1–4, then job 26638955, 2026-10-02); preliminary until the convergence checks below |

Along the averaged branch every solve converged (≤ 4 Newton iterations,
residual ≤ 9.4e-9). Re ω₂ decreases monotonically from 0.415628 to 0.395101;
Im ω₂ increases monotonically from 6.37e-11 to 1.96e-6 at A₁ = 13 (×3.1e4) and
is positive throughout. Im ω₂ exceeds the primary's Im ω₁ above A₁ ≈ 12.2.

### Why the fixed-time run is not physical

Evaluating the integral in Eq. (61) at one instant also keeps a term at
(4ω₁ − ω₂, m₂) proportional to conj(R). That couples R to R*, so ω₂ depends on
the secondary's normalization phase. Fixing the (2,2) outer amplitude to e^{iθ}
instead of 1:

| A₁ | Im ω₂ at θ = 0, π/4, π/2, 3π/4 (fixed-time) |
|---|---|
| 2.05 | +5.8e-8, −1.4e-4, −2.7e-8, +1.4e-4 |
| 9.89 | −1.2e-6, −2.6e-3, −2.0e-5, +2.6e-3 |
| 13.0 | +1.4e-4, −3.3e-3, +1.5e-3, +3.6e-3 |

With the period-averaged term the equation is complex-linear, ω₂ is independent
of θ, and Im ω₂ is positive at all three points (+2.29e-8, +1.32e-6, +1.96e-6).
The average is what the paper's (ω, m) classification of components (Sec. II B,
after Eq. (22)) implies for the (ω₂, m₂) equation. At order A₁² the two readings
agree; the conjugate term starts at A₁⁴.

`results/cosine_secondary_energy.csv` (2026-10-02, from
`src/adiabatic/drivers/secondary_energies.jl`): the averaged secondary's energy
Ê₂(A₁) per unit (2,2) outer amplitude, with horizon and outer fluxes. Energy
includes the interaction term (μ²/2)⟨cos ϕ₁⟩ϕ₂². The balance
−(F_H + F_out)/(2 Γ₂ E₂) lies between 0.9999965 and 1.0000001 at all 232 points.
Ê₂ falls from 7.9e8 to 2.9e3; this reflects the normalization at r* = 150 as
the secondary contracts, not a physical energy change. About 0.15% of the
linear secondary's energy lies beyond r* = 150 (outer flux at small A₁).

## Checks so far

- Linear limit (A₁ = 0): ω₂ = 0.4156284036 + 6.3703e-11 i, matching spheroidal
  Leaver (`--N 2000`) in Re to 10 digits; Im is high by 7e-6 relative at step 0.2,
  r*_out = 150, and agrees to 5e-7 at step 0.1, r*_out = 250.
- Small A₁: Re shift = −0.0101 A₁², Γ₂/Γ₂(linear) − 1 = 42 A₁², constant to 5e-4
  over the three smallest saved amplitudes.
- Radial step (`tests/adiabatic/check_secondary_step_convergence.jl`): at A₁ = 2.05,
  re-solving primary and secondary at half and quarter step changes ω₂ by
  7.0e-11 then 4.4e-12 (ratio 16, RK4-dominated); Γ₂ at the production step is
  accurate to about 6e-6 relative.
- At A₁ = 13 the same check fails in the **primary**: re-solving it at step 0.05
  stagnates. The saved primary solutions (step 0.1), re-evaluated at steps 0.05
  and 0.025, have scaled matching residual ≈ 8e-8 at A₁ = 2.05, 3e-7 at 6.26,
  4e-5 at 9.14, 6e-4 at 10.2, 2e-2 at 12.4, and 4e-2 at 13. The primary branch is
  therefore not step-converged above A₁ ≈ 9, and neither is the secondary built on it.
- Not yet checked: outer-radius, angular, and l-truncation convergence.
