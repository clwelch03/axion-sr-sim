# Secondary cloud on the nonlinear primary (Omiya et al. 2211.01949, Sec. III)

## Start here

This campaign computes a linear secondary cloud on a nonlinear primary,
the driven $m = 0$ and $m = 3$ fluxes, and two-cloud evolution on fixed Kerr
at $a/M = 0.99$, $M\mu = 0.42$ and $0.3$. It supplies the nonperturbative
Fig. 6 trajectories and the Fig. 7 profile comparison. All commands run from
the repository root. [HANDOFF.md](../../../HANDOFF.md) owns project status
and priorities; the records below give the detailed numerical evidence.

| Input/result | $M\mu = 0.42$ | $M\mu = 0.3$ |
|---|---|---|
| Primary checkpoint/table | `../fig2_20260922/results/cosine_branch.{jls,csv}` | `results/alpha030_cosine_branch.{jls,csv}` |
| Averaged secondary checkpoint/table | `results/cosine_secondary.{jls,csv}` | `results/alpha030_cosine_secondary.{jls,csv}` |
| Secondary energy table | `results/cosine_secondary_energy.csv` | `results/alpha030_cosine_secondary_energy.csv` |
| Driven flux tables | `results/cosine_dissipative_m0.csv`, `results/cosine_radiative_m3.csv` | `results/alpha030_cosine_dissipative_m0.csv`, `results/alpha030_cosine_radiative_m3.csv` |
| Evolution result | `results/two_cloud_evolution_alpha042.csv` | `results/two_cloud_evolution_alpha030.csv` |

Paths in the table are relative to this campaign. The coefficient tables
contain 232 aligned points at $M\mu = 0.42$ and 83 at $0.3$. The evolution
files each contain 1001 saved times. Use the averaged secondary results;
`cosine_secondary_fixed_time.*` is a superseded phase-dependent experiment,
retained for the comparison below.

Redraw Fig. 6 directly from saved evolution CSVs:

```bash
PYTHONPATH=src python runs/adiabatic/multimode_secondary/scripts/plot_fig6.py
```

To regenerate those CSVs from the saved coefficient tables, then redraw:

```bash
julia --project=. src/adiabatic/primary_secondary/drivers/two_cloud_evolution.jl
PYTHONPATH=src python runs/adiabatic/multimode_secondary/scripts/plot_fig6.py
```

Output: `plots/multimode_fig6_analogue.png`. The evolution command replaces
both trajectory CSVs. For Fig. 7, reconstruct profiles from saved checkpoints
and plot at the comparison energy:

```bash
julia --project=. runs/adiabatic/multimode_secondary/scripts/export_profiles.jl alpha030 660
PYTHONPATH=src python runs/adiabatic/multimode_secondary/scripts/plot_fig7_analogue.py alpha030 660
```

This writes `work/alpha030_profiles.csv` and
`plots/multimode_alpha030_fig7_analogue.png`. Full primary/secondary and flux
recipes appear below. Continuation resumes an existing prefix; energy and
flux exporters recompute and replace their output CSVs. Run costly solves
through Slurm, with one writer per output.

**Validation limits:** Fig. 6 agreement is at reading precision at $M\mu = 0.42$;
the $0.3$ primary peak is about 10% lower and plateaus 5–10% higher than the
paper's plotted values. Independent perturbative flux prefactors and the
physical impact of the Sec. III approximations remain unchecked. The full
high-amplitude tables inherit the primary's radial convergence failure above
$A_1 \simeq 9$; the $M\mu = 0.42$ Fig. 6 trajectory stays below that range.
See the channel and secondary checks below for the tested domains and
remaining boundary/angular checks.

## Dilute self-gravity estimates — October 7, 2026

Apply published hydrogenic self- and cross-gravity coefficients to the saved
scalar trajectories:

```bash
julia --project=. runs/adiabatic/multimode_secondary/scripts/estimate_self_gravity.jl
```

The script uses [Baryakhtar et al., Appendix G](https://arxiv.org/html/2011.11646#A7)
with its rounded coefficients, writing $q_i=M_{c,i}/M_0=\epsilon_f E_i$ and
$\epsilon_f=(f_a/M_{\rm Pl})^2$ with the unreduced Planck mass:

$$
\frac{\delta\omega_1}{\mu}=-\alpha^2(0.19q_1+0.11q_2),\qquad
\frac{\delta\omega_2}{\mu}=-\alpha^2(0.11q_1+0.09q_2).
$$

The occupation-to-mass conversion uses $M_{c,i}\simeq\mu N_i$, at leading
nonrelativistic order. These are extrapolated dilute estimates using the saved
energies; no actual-profile Poisson calculation, nonlinear scalar response,
metric correction, or corrected evolution is performed. Frequencies used for
the comparison margins are interpolated from the existing nonlinear tables.
The $0.19$ coefficient is the published rounded value of $237/1280$.

All entries below are maxima over the 1001 saved samples; maxima in different
columns need not occur together. Define the splitting $s=\omega_2-\omega_1$
and outward $m=3$ scalar threshold margin $d_3=2\omega_2-\omega_1-\mu$.
Their corrections are $\delta s=\delta\omega_2-\delta\omega_1$ and
$\delta d_3=2\delta\omega_2-\delta\omega_1$.

| $\alpha$ | $f_a$ [GeV] | $\max|\delta\omega_1|/\mu$ | $\max|\delta\omega_2|/\mu$ | $\max|\delta s/s|$ | $\max|\delta d_3/d_3|$ |
|---|---:|---:|---:|---:|---:|
| 0.42 | $10^{16}$ | $7.70\times10^{-6}$ | $4.49\times10^{-6}$ | 0.0177% | 0.0219% |
| 0.30 | $10^{16}$ | $5.80\times10^{-6}$ | $3.38\times10^{-6}$ | 0.0317% | 0.0483% |
| 0.42 | $10^{17}$ | $7.70\times10^{-4}$ | $4.49\times10^{-4}$ | 1.77% | 2.19% |
| 0.30 | $10^{17}$ | $5.80\times10^{-4}$ | $3.38\times10^{-4}$ | 3.17% | 4.83% |

Both frequencies decrease; the splitting increases and $d_3$ decreases.
No threshold changes sign through $f_a=10^{17}$ GeV. At that decay constant,
the primary superradiant margin increases by at most 1.24% (0.42) and 0.126%
(0.30); the secondary margin changes by at most 0.0417% and 0.0178%.
The distance of the driven $m=0$ frequency below the positive scalar mass
threshold increases by at most 2.31% and 3.98%, respectively.
These are changes in detunings, not calculated fractional flux changes.
Bound-resonance detunings and nonlinear profile responses have not been checked.

All shifts scale as $f_a^2$ for these unchanged trajectories: results at
$10^{15}$ GeV are 100 times smaller than at $10^{16}$ GeV.
The peak cloud mass fractions at $10^{17}$ GeV are 0.0232 and 0.0343, while
scalar cumulative angular-momentum extraction is 0.149 and 0.353 of the
initial BH angular momentum. Frozen BH parameters therefore need a separate
assessment even where these gravity shifts look modest. At $10^{18}$ GeV
the inferred peak cloud mass exceeds the initial BH mass and scalar
extraction exceeds its reservoirs; those rows are flagged and are not
physical predictions.

Outputs (replaced on rerun):

- `results/self_gravity_estimate_summary.csv`: four decay constants per alpha,
  including superradiant/threshold comparisons and scalar reservoir flags.
- `results/self_gravity_estimate_alpha042_fa1e16.csv` and
  `results/self_gravity_estimate_alpha030_fa1e16.csv`: every saved time,
  mass fractions, separate self/cross shifts, and unperturbed margins at
  $f_a=10^{16}$ GeV. Frequencies/shifts are in code units $GM_0\omega$;
  divide by $\alpha$ to obtain $\omega/\mu$. Physical time assumes axion
  rest energy $3\times10^{-13}$ eV and covers 0.709 and 2.841 days.

Checks: independent CSV recomputation of the mass-to-shift conversion,
combination-frequency signs, maximum threshold ratios, $f_a^2$ scaling, and
absence of threshold crossings through $10^{17}$ GeV passed. The scalar
normalization script retains its direct-run results and can now be included
without executing its report loop.

## Secondary eigenproblem and normalization

The l = m = 2 secondary of Eq. (61), solved as an eigenproblem for ω₂ on each
saved point of the cosine primary branch
(`runs/adiabatic/fig2_20260922/results/cosine_branch.jls`, a/M = 0.99, Mμ = 0.42).
The secondary keeps (l, m₂) = (2,2), (4,2), (6,2), uses the primary's radial domain
(r* from −100 to 150, matching at 5) at twice its step (0.2), and fixes the (2,2)
outer amplitude to 1. Library code: `src/adiabatic/primary_secondary/`; driver:
`src/adiabatic/primary_secondary/drivers/continue_secondary.jl`.

## Running

```bash
mkdir -p runs/adiabatic/multimode_secondary/logs
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
`src/adiabatic/primary_secondary/drivers/secondary_energies.jl`): the averaged secondary's energy
Ê₂(A₁) per unit (2,2) outer amplitude, with horizon and outer fluxes. Energy
includes the interaction term (μ²/2)⟨cos ϕ₁⟩ϕ₂². The balance
−(F_H + F_out)/(2 Γ₂ E₂) lies between 0.9999965 and 1.0000001 at all 232 points.
Ê₂ falls from 7.9e8 to 2.9e3; this reflects the normalization at r* = 150 as
the secondary contracts, not a physical energy change. About 0.15% of the
linear secondary's energy lies beyond r* = 150 (outer flux at small A₁).

Regenerate the $M\mu = 0.42$ secondary energy table after continuation:

```bash
julia --project=. src/adiabatic/primary_secondary/drivers/secondary_energies.jl \
    runs/adiabatic/fig2_20260922/results/cosine_branch.jls \
    runs/adiabatic/multimode_secondary/results/cosine_secondary.jls \
    runs/adiabatic/multimode_secondary/results/cosine_secondary_energy.csv
```

### Mμ = 0.3: the paper's Fig. 7

Same settings as Mμ = 0.42 except outer r* = 300 (both linear frequencies are
unchanged for outer r* from 200 to 400 and match spheroidal Leaver: primary to
all printed digits, secondary Γ₂ to 8e-6).

```bash
julia --project=. runs/adiabatic/multimode_secondary/scripts/prepare_alpha030_seed.jl
sbatch runs/adiabatic/multimode_secondary/scripts/submit_alpha030_primary.sh   # job 26640571
AXION_SECONDARY_ULP_TOLERANCE_FACTOR=2 julia --project=. src/adiabatic/primary_secondary/drivers/continue_secondary.jl \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_branch.jls \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_secondary
julia --project=. src/adiabatic/primary_secondary/drivers/secondary_energies.jl \
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
Library code: `src/adiabatic/primary_secondary/dissipative_mode.jl`
(and `mass_deficit_harmonic` in `background.jl`).

```bash
julia --project=. src/adiabatic/primary_secondary/drivers/dissipative_fluxes.jl \
    runs/adiabatic/fig2_20260922/results/cosine_branch.jls \
    runs/adiabatic/multimode_secondary/results/cosine_secondary.jls \
    runs/adiabatic/multimode_secondary/results/cosine_secondary_energy.csv \
    runs/adiabatic/multimode_secondary/results/cosine_dissipative_m0.csv
julia --project=. src/adiabatic/primary_secondary/drivers/dissipative_fluxes.jl \
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

### Piece 5: the driven $m = 3$ field

The $m = 3$ field at $\omega_3 = 2\,\mathrm{Re}\,\omega_2 - \mathrm{Re}\,\omega_1$
(Eqs. (51)–(52)), driven by $\tfrac12\mu^2\,\overline{[\sin\phi_1]_1}\,P_2^2$:
the $e^{+i\psi_1}$ harmonic of $\sin\phi_1$ times the secondary's profile
squared, with modes $(l, m) = (3,3), (5,3), (7,3)$. It radiates to infinity:
this is the perturbative $F_3$ (Eq. (39), (A14)). Library code:
`src/adiabatic/primary_secondary/radiative_mode.jl` (and `sine_harmonic` in
`background.jl`).

Two differences from piece 4:

- **Far outer leg.** $\omega_3$ is above $\mu$ by only $2$–$9\times10^{-3}$
  ($M\mu = 0.42$) and $6$–$9\times10^{-4}$ ($0.3$), so the outgoing series does
  not converge at the primary's outer boundary. `radiative_outer_boundary_state`
  starts it at $r_* = 5000$ (order 16) and integrates the source-free Kerr
  equation inward to the boundary at step 0.25. The primary's averaged
  background (its radiating tail, $\lesssim 10^{-7}$) is left out of that leg;
  see the comment in the function.
- **Two fluxes.** $F_\infty$ is the outward current at the primary's outer
  boundary (equal to the flux to infinity, since the leg has no source and a
  real operator). $F_H$ is negative: $\omega_3 < 3\Omega_H$, so the black hole
  feeds the field. Eq. (65) counts both.

```bash
julia --project=. src/adiabatic/primary_secondary/drivers/radiative_fluxes.jl \
    runs/adiabatic/fig2_20260922/results/cosine_branch.jls \
    runs/adiabatic/multimode_secondary/results/cosine_secondary.jls \
    runs/adiabatic/multimode_secondary/results/cosine_secondary_energy.csv \
    runs/adiabatic/multimode_secondary/results/cosine_radiative_m3.csv
julia --project=. src/adiabatic/primary_secondary/drivers/radiative_fluxes.jl \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_branch.jls \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_secondary.jls \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_secondary_energy.csv \
    runs/adiabatic/multimode_secondary/results/alpha030_cosine_radiative_m3.csv
```

About 4–9 s per point locally. The key column is
`infinity_flux_per_energy_squared` $= F_\infty/\hat E_2^2$: flux per unit
secondary energy squared, independent of the secondary's normalization.
`AXION_RADIATIVE_MAX_L` (odd, default 7) changes the $l$ truncation.

Results (2026-10-05), `results/cosine_radiative_m3.csv` (232 points) and
`results/alpha030_cosine_radiative_m3.csv` (83 points):

| | $M\mu = 0.42$ | $M\mu = 0.3$ |
|---|---|---|
| $\omega_3$ along the branch | 0.42233–0.43118 | 0.30060–0.30086 |
| $F_\infty/\hat E_2^2$ | $7.0\times10^{-13}$, peak $1.96\times10^{-6}$ at $A_1 = 5.0$, $2.4\times10^{-7}$ at 13 | $3.1\times10^{-12}$ → $5.2\times10^{-9}$ (monotonic) |
| power ratio $(F_H+F_\infty)/\text{source work}$ | within $6.2\times10^{-8}$ of 1 | within $2.4\times10^{-10}$ of 1 |
| residual after solve / before | $\le 10^{-13}$ | $\le 4\times10^{-14}$ |
| $-F_H/F_\infty$ | $\le 2.8\times10^{-4}$ | $\le 1.5\times10^{-11}$ |
| share of $F_\infty$ in $l = 5$, $l = 7$ | $\le 5.7\times10^{-4}$, $\le 10^{-9}$ | $\le 5\times10^{-4}$, $\le 2.3\times10^{-10}$ |
| scaled condition number of the matching matrix | 1.7–159 | up to 1009 |

Checks (scripts in `tests/adiabatic/`):

- `check_sine_harmonic.jl` (5a): $[\sin\phi_1]_0$ and $[\sin\phi_1]_2$ vanish to
  round-off; $[\sin\phi_1]_1 \to P_1$ with relative difference $0.592\,A_1^2$;
  azimuthal grid doubling changes nothing.
- `check_radiative_rhs.jl` (5b): with the source off the radial current is
  conserved up to RK4 error (ratio per halving reaches 15.5 inward by step
  0.0125); even-$l$ projections of the source vanish.
- `check_radiative_outer.jl`, `check_radiative_far_background.jl` (5c): the
  outgoing series satisfies the equation (16 per halving); the boundary state
  is independent of series order and far start; the series used directly at
  the boundary would be off by up to 14% ($0.42$) and 160% ($0.3$). Leaving
  out the primary's background changes the solution by $\le 10^{-7}$.
- `check_radiative_solve.jl` (5d): residual, power identity, $l$ shares; the
  small-$A_1$ log-slope of $F_\infty/\hat E_2^2$ is $1.9999$, the $E_1$ law of
  the perturbative $F_3E_1E_2^2$ (1.9989 at $M\mu = 0.3$). It bends below 2
  sooner than piece 4's slope bends below 4.
- `check_radiative_convergence.jl` (5e, log in `logs/`), relative changes of
  $F_\infty/\hat E_2^2$:
  - radial step at $A_1 = 2.05$: $7.9\times10^{-8}$ from step 0.1 to 0.05,
    ratio of successive changes 15.8 (RK4 error);
  - adding $l = 9$: $0$, $6.6\times10^{-7}$, $5.8\times10^{-6}$ at
    $A_1 = 0.001, 2.05, 13$;
  - far start 2500 or 10000, far step 0.125: $\le 5\times10^{-11}$. The flux
    is far less sensitive than the boundary state itself ($4\times10^{-7}$),
    because the far-leg error is mostly a phase that the outer amplitudes absorb.
- Not checked: the coefficient against an independent perturbative $F_3$
  (Eq. (A14)).

### Piece 6: two-cloud evolution, the paper's Fig. 6

Energy and angular-momentum balance (Eqs. (62)–(63)) written in the cloud
energies $E_1$ and $E_2 = \hat E_2|A_2|^2$ (the paper's unit-energy
normalization of the secondary); equations and choices in
[multimode-adiabatic-design.md](../../../docs/multimode-adiabatic-design.md),
"Piece 6 design". The coefficients come from the piece 1–5 tables, splined in
$E_1$. Library code: `src/adiabatic/primary_secondary/evolution.jl`; driver
and table reader: `src/adiabatic/primary_secondary/drivers/two_cloud_evolution.jl`.

```bash
julia --project=. src/adiabatic/primary_secondary/drivers/two_cloud_evolution.jl
PYTHONPATH=src python runs/adiabatic/multimode_secondary/scripts/plot_fig6.py
```

A few seconds in all. Writes `results/two_cloud_evolution_alpha042.csv` and
`_alpha030.csv` (time in units of $\omega_{1,I}$, energies in $F_a^2M$) and
`plots/multimode_fig6_analogue.png` (nonperturbative curves only).

Results (2026-10-05), against Fig. 6 read by eye on its log axes:

| | $M\mu = 0.42$ | paper | $M\mu = 0.3$ | paper |
|---|---|---|---|---|
| $E_1$ peak (value at $\omega_{1,I}t$) | 341 at 3.48 | ~350 at ~3.5 | 503 at 1.82 | ~550 at ~1.85 |
| $E_2$ peak | 125 at 3.61 | ~130 at ~3.6 | 172 at 1.95 | ~170 at ~1.95 |
| $E_1$, $E_2$ at $\omega_{1,I}t = 10$ | 38, 63 | ~38, ~65 | 60, 84 | ~55, ~78 |

Agreement is at reading precision at $M\mu = 0.42$; at $0.3$ our peak looks
about 10% lower and our plateaus 5–10% higher than the figure, at the limit of
what can be read off it.

Checks (scripts in `tests/adiabatic/`):

- `check_two_cloud_splines.jl` (6a): the not-a-knot spline reproduces a cubic
  to $10^{-14}$ and converges at fourth order; splines from every other saved
  point miss the omitted points by $\le 5.5\times10^{-3}$ (relative to the
  coefficient's maximum) at $M\mu = 0.42$, $E_1 \le 1000$, and
  $\le 4\times10^{-6}$ at $0.3$. Splines in $\ln E_1$ were 5–100 times worse.
- `check_two_cloud_rates.jl` (6b): one channel at a time, the rates obey the
  exact limits to $\le 6\times10^{-12}$: the primary alone grows at $-F_1$; the
  secondary alone leaves $E_1$ fixed; the $m = 0$ channel moves energy as
  $F/\omega_0$ quanta, each taking two primary quanta and giving one
  secondary quantum; the $m = 3$ channel the reverse at $F/\omega_3$. At the
  first point $F_1/E_1 = -2\,\mathrm{Im}\,\omega_1$ to $10^{-6}$.
- `check_two_cloud_evolution.jl` (6c, log in `logs/`): time step 2e-3 → 5e-4
  changes the final energies by $\le 4\times10^{-12}$ and the peaks by
  $\le 5\times10^{-5}$ (peak times move by one sample); splines from every
  other point change everything by $\le 2\times10^{-6}$; the initial slope of
  $\ln E_1$ is 2.001 ($0.42$) and 2.017 ($0.3$, where the growth rate
  already rises with amplitude).
- Not checked: how much each approximation of Sec. III matters (for
  example the secondary's own nonlinearity at $E_2 \sim 100$), and the
  perturbative curves of the figure.

## Checks so far (secondary)

- Linear limit (A₁ = 0): ω₂ = 0.4156284036 + 6.3703e-11 i, matching spheroidal
  Leaver (`--N 2000`) in Re to 10 digits; Im is high by 7e-6 relative at step 0.2,
  r*_out = 150, and agrees to 5e-7 at step 0.1, r*_out = 250.
- Small A₁: Re shift = −0.0101 A₁², Γ₂/Γ₂(linear) − 1 = 42 A₁², constant to 5e-4
  over the three smallest saved amplitudes.
- Radial step (`tests/adiabatic/primary_secondary/check_secondary_step_convergence.jl`): at A₁ = 2.05,
  re-solving primary and secondary at half and quarter step changes ω₂ by
  7.0e-11 then 4.4e-12 (ratio 16, RK4-dominated); Γ₂ at the production step is
  accurate to about 6e-6 relative.
- At A₁ = 13 the same check fails in the **primary**: re-solving it at step 0.05
  stagnates. The saved primary solutions (step 0.1), re-evaluated at steps 0.05
  and 0.025, have scaled matching residual ≈ 8e-8 at A₁ = 2.05, 3e-7 at 6.26,
  4e-5 at 9.14, 6e-4 at 10.2, 2e-2 at 12.4, and 4e-2 at 13. The primary branch is
  therefore not step-converged above A₁ ≈ 9, and neither is the secondary built on it.
- Not yet checked: outer-radius, angular, and l-truncation convergence.
