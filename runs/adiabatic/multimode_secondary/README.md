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

### Mμ = 0.3: the paper's Fig. 7

Same settings as Mμ = 0.42 except outer r* = 300 (both linear frequencies are
unchanged for outer r* from 200 to 400 and match spheroidal Leaver: primary to
all printed digits, secondary Γ₂ to 8e-6).

```bash
julia --project=. runs/adiabatic/multimode_secondary/scripts/prepare_alpha030_seed.jl
sbatch runs/adiabatic/multimode_secondary/scripts/submit_alpha030_primary.sh   # job 26640571
AXION_SECONDARY_ULP_TOLERANCE_FACTOR=2 julia --project=. src/adiabatic/drivers/continue_secondary.jl \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_branch.jls \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_secondary
julia --project=. src/adiabatic/drivers/secondary_energies.jl \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_branch.jls \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_secondary.jls \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_secondary_energy.csv
julia --project=. runs/adiabatic/multimode_secondary/scripts/export_profiles.jl alpha030 660
PYTHONPATH=src python runs/adiabatic/multimode_secondary/scripts/plot_fig7_analogue.py alpha030 660
```

The primary has 83 points to A₁ = 0.06 (E₁ = 651 at A₁ = 0.0364, the nearest
to 660). The secondary needs `AXION_SECONDARY_ULP_TOLERANCE_FACTOR=2`: its
residual floor is the residual's change from one ulp of Re ω₂ (5.9e-8 at
A₁ = 0.0014, 1.2e-6 at 0.016), which grows along the branch, so a fixed
tolerance fails partway. Energy balance is 1.0000046 at all 83 points.

Fig. 7 comparison (`plots/multimode_alpha030_fig7_analogue.png`), peaks of |R|:

| | l=m=1 linear | l=m=1 nonlinear | l=m=2 linear | l=m=2 nonlinear |
|---|---|---|---|---|
| paper, read from figure (E ≈ 660) | 0.116, r* ≈ 22 | 0.126, r* ≈ 22 | 0.030, r* ≈ 72 | 0.034, r* ≈ 62 |
| ours, E = 651 | 0.1152, r* = 24.4 | 0.1244, r* = 24.0 | 0.0301, r* = 71.8 | 0.0338, r* = 64.0 |
| ours × √(660/651) | 0.116 | 0.125 | 0.030 | 0.034 |

The near-horizon plateau (|R₁₁| ≈ 0.006) and the tails also agree by eye.

### Piece 4: the driven $m = 0$ field

The $m = 0$ field at $\omega_0 = 2\,\mathrm{Re}\,\omega_1 - \mathrm{Re}\,\omega_2$
(Eqs. (51)–(52)), driven by the $e^{-2i\psi_1}$ harmonic of $1 - \cos\phi_1$
times the conjugate secondary, with modes $(l, m) = (0,0), (2,0), (4,0)$. It
carries energy into the horizon: this is the dissipation of the perturbative
$F_0$ (Eq. (38), (A11)). Method, notation, and every check are explained in
[docs/multimode-piece4-guide.md](../../../docs/multimode-piece4-guide.md).
Library code: `src/adiabatic/multimode/dissipative_mode.jl`
(and `mass_deficit_harmonic` in `background.jl`).

```bash
julia --project=. src/adiabatic/drivers/dissipative_fluxes.jl \
    runs/adiabatic/fig2_20260922/results/cosine_branch.jls \
    runs/adiabatic/multimode_secondary/results/cosine_secondary.jls \
    runs/adiabatic/multimode_secondary/results/cosine_secondary_energy.csv \
    runs/adiabatic/multimode_secondary/results/cosine_dissipative_m0.csv
julia --project=. src/adiabatic/drivers/dissipative_fluxes.jl \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_branch.jls \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_secondary.jls \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_secondary_energy.csv \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_dissipative_m0.csv
```

About 3 s per point locally (12 min for $M\mu = 0.42$, 7 min for $0.3$). The
key column is `horizon_flux_per_energy` $= F_H/\hat E_2$: energy into the
horizon per unit time, per unit secondary energy, independent of the
secondary's normalization.

Results (2026-10-02), `results/cosine_dissipative_m0.csv` (232 points) and
`results/alpha030_cosine_dissipative_m0.csv` (83 points):

| | $M\mu = 0.42$ | $M\mu = 0.3$ |
|---|---|---|
| $\omega_0$ along the branch | 0.4020 → 0.3406 | 0.2942 → 0.2933 |
| $F_H/\hat E_2$ (monotonic) | $6.0\times10^{-14}$ → $0.040$ | $1.0\times10^{-11}$ → $1.8\times10^{-5}$ |
| power ratio $(F_H+F_\text{out})/\text{source work}$ | within $10^{-6}$ of 1 | within $8\times10^{-8}$ of 1 |
| residual after solve / before | $\le 8\times10^{-15}$ | $\le 3.5\times10^{-15}$ |
| $\lvert F_\text{out}/F_H\rvert$ | $\le 10^{-26}$ | $\le 5\times10^{-26}$ |
| share of $F_H$ in $l = 2$, $l = 4$ | $\le 3.4\times10^{-5}$, $\le 4\times10^{-10}$ | $\le 4\times10^{-8}$, $\le 10^{-16}$ |
| scaled condition number of the matching matrix | 4–227 | 12 → 4100 |

Checks (scripts in `tests/adiabatic/`, logs in `logs/`):

- `check_mass_deficit_harmonic.jl`, `check_dissipative_current.jl`,
  `check_dissipative_boundary.jl`, `check_dissipative_solve.jl`: steps 4a–4d,
  results in the guide, Section 8.
- **Small $A_1$:** $d\ln(F_H/\hat E_2)/d\ln A_1 = 3.9997, 3.9990, 3.9960$ as
  $A_1$ doubles from 0.0011 ($M\mu = 0.42$), and 3.997 at $M\mu = 0.3$. This is
  the $E_1^2$ law of the perturbative $F_0 E_1^2 E_2$. It bends below 4 as the
  nonlinearity sets in.
- **Radial step** (`check_dissipative_convergence.jl`, $A_1 = 2.05$): primary at
  steps 0.1, 0.05, 0.025 changes $F_H/\hat E_2$ by $-1.6\times10^{-10}$, then
  $+1.7\times10^{-10}$ (relative $4\times10^{-8}$). The sign alternates: this
  is the Newton-tolerance floor, not RK4 error. Not checked above
  $A_1 \approx 9$, where the primary is not step-converged.
- **$l$-truncation** (same script): adding $(6,0)$ changes $F_H/\hat E_2$ by
  $10^{-15}$, $1.4\times10^{-4}$, $4.5\times10^{-4}$ at $A_1 = 0.001, 2.05, 13$.
- **Resonance at $M\mu = 0.3$** (`check_dissipative_resonance.jl`): scanning the
  drive frequency at fixed background shows a resonance in the $m = 0$ response
  near $\omega \approx 0.295$–$0.296$, with a sharp dip just above. This is
  consistent with the hydrogenic $m = 0$ level near $0.2966$. $\omega_0$ sits on
  its lower flank, at 60–75% of the peak flux, at both $A_1 = 0.001$ and $0.06$.
  So the enhancement is already present in the linear limit (it is part of the
  perturbative $F_0$ as well).
- **Condition number at $M\mu = 0.3$** (same script): the growth from 12 to
  4100 is *not* that resonance. $\sigma_\text{min}/\sigma_\text{max}$ is flat
  across the frequency band, and its singular vector at $A_1 = 0.06$ lies
  entirely on the outer $l = 2, 4$ amplitudes. Those two outer solutions line
  up as they are integrated inward through the cloud, where the background
  mixes even $l$. This is a shooting artifact, harmless at $4\times10^3$ (the
  residual and power identity confirm the solution). The driver warns above
  $10^8$; the fix would be re-orthogonalizing the outer solutions.
- Not checked: outer-radius convergence of the $m = 0$ field, and an
  independent evaluation of the perturbative $F_0$ (Eq. (38)) to compare
  coefficients.

## Checks so far (secondary)

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
