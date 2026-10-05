# Project handoff

Updated 2026-10-05. Active work: reproducing Omiya, Takahashi, Tanaka & Yoshino,
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

All six pieces are done at $M\mu = 0.42$ and $0.3$, $a/M = 0.99$, and the
nonperturbative curves of the paper's **Fig. 6 are reproduced** (piece 6):
peak and plateau energies agree with the figure to reading precision
(`plots/multimode_fig6_analogue.png`; table in the campaign README, "Piece 6").

1. Period average ⟨1 − cos ϕ₁⟩ of the primary (`src/adiabatic/multimode/background.jl`).
2. Secondary eigenproblem, Eq. (61), continued along all 232 primary points
   (`multimode/secondary_mode.jl`, driver `drivers/continue_secondary.jl`).
3. Secondary energy Ê₂(A₁) and boundary fluxes (`multimode/secondary_energy.jl`,
   driver `drivers/secondary_energies.jl`).
4. The driven $m = 0$ field at $\omega_0 = 2\,\mathrm{Re}\,\omega_1 - \mathrm{Re}\,\omega_2$
   and its horizon flux (`multimode/dissipative_mode.jl`, `mass_deficit_harmonic`
   in `background.jl`, driver `drivers/dissipative_fluxes.jl`). Explained step by
   step in [docs/multimode-piece4-guide.md](docs/multimode-piece4-guide.md).
   Results: `results/cosine_dissipative_m0.csv`, `results/alpha030_cosine_dissipative_m0.csv`.
5. The driven $m = 3$ field at $\omega_3 = 2\,\mathrm{Re}\,\omega_2 - \mathrm{Re}\,\omega_1$,
   radiating to infinity (`multimode/radiative_mode.jl`, `sine_harmonic` in
   `background.jl`, driver `drivers/radiative_fluxes.jl`). $\omega_3$ is just
   above $\mu$, so its outgoing series is started at $r_* = 5000$ and integrated
   in to the primary's boundary. Results: `results/cosine_radiative_m3.csv`,
   `results/alpha030_cosine_radiative_m3.csv`.
6. Evolution of $E_1, E_2$ by energy and angular-momentum balance, written in
   the cloud energies (`multimode/two_cloud_evolution.jl`, driver
   `drivers/two_cloud_evolution.jl`, plot `scripts/plot_fig6.py`).

Validated: linear limits against spheroidal Leaver (both clouds, both α);
small-A₁ coefficients (Re ω₂ shift −0.0101 A₁², Γ₂/Γ₂^lin − 1 = 42 A₁²);
exact energy balance 2Γ₂E₂ = −F to 3e-6 at every point; radial-step
convergence of the secondary at A₁ = 2.05 (RK4 ratio 16).

Pieces 1–3 also run at Mμ = 0.3 and **reproduce the paper's Fig. 7
quantitatively**: at E₁ = E₂ ≈ 660, nonlinear/linear peak |R| are
0.125/0.116 (l = m = 1) and 0.034/0.030 (l = m = 2), matching the figure to its
reading precision, with the same peak positions. Commands and the table are in
the campaign README.

Piece 4 validation (2026-10-02):
- the power identity $F_H + F_\text{out} = $ source work holds to $10^{-6}$ at
  every point;
- $F_H/\hat E_2 \propto A_1^4$ (that is, $E_1^2$) at small $A_1$, with log-slope
  3.9997;
- radial step converged to $4\times10^{-8}$ at $A_1 = 2.05$;
- adding $l = 6$ changes the flux by at most $4.5\times10^{-4}$;
- $l = 0$ carries the flux to parts in $10^5$.

At $M\mu = 0.3$, $\omega_0$ lies just below an $m = 0$ resonance near 0.296.
The matching matrix's condition number grows to $4\times10^3$ there, but
that is a shooting artifact: the outer $l = 2, 4$ solutions line up as they
are integrated inward. It is not the resonance.

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
- **Piece 6 in cloud energies** $E_1$, $E_2 = \hat E_2|A_2|^2$ (the paper's
  unit-energy normalization), so $d\hat E_2/dA_1$ never appears; the
  $d\omega_2/dE_1$ term from $J_2 = 2E_2/\omega_2$ is kept. Coefficients are
  splined in $E_1$, not $\ln E_1$.
- **$m = 3$ horizon flux included** in $F_E$, $F_J$ (Eq. (65)); it is negative
  (superradiant) and at most $3\times10^{-4}$ of the flux to infinity.

## Known caveats

- The Mμ = 0.42 primary branch is **not step-converged above A₁ ≈ 9**: its
  step-0.1 solutions re-evaluated at step 0.05 have scaled residual 4e-5 at
  A₁ = 9.1 rising to 4e-2 at 13. This bears on the cosine flux zero near
  A₁ ≈ 12.96. It does not affect Fig. 6 at Mμ = 0.42, which stays below
  E₁ ≈ 1000 (A₁ ≲ 0.8, step-converged to ~1e-7).
- Unchecked: secondary l-truncation (l ≤ 6) and outer boundary at large A₁.
  The paper reports neither.
- Pieces 4–5: the outer-radius convergence of the $m = 0$ field is unchecked, and
  so is each coefficient against an independent evaluation of the perturbative
  $F_0$ and $F_3$ (Eqs. (38)–(39), (A11), (A14), e.g. by a Green's function).
  The power laws match; the prefactors have not been compared.

## Next steps

- **Commit.** Pieces 5 and 6 (library, drivers, tests, results CSVs, docs)
  are uncommitted as of 2026-10-05.
- **Perturbative coefficients.** The small-$E_1$ limits of the spline
  coefficients `m0_flux_coefficient` and `m3_flux_coefficient` are the paper's
  $F_0$ and $F_3$ (Eqs. (38)–(39)). An independent evaluation of those (for
  example by the Green's-function formula) would check the prefactors of
  pieces 4–5, which so far are checked only through power laws and the power
  identity. The same coefficients would also give the perturbative curves of
  Fig. 6.
- **$M\mu = 0.3$ details.** Our $E_1$ peak looks ~10% below the figure and the
  plateaus 5–10% above; this is at the limit of reading a log axis. Not
  investigated.
- Beyond Sec. III: see the design doc's last section.

## Working with the user

Build incrementally with a checkpoint per piece; the user types library code
from pasted snippets (paste only code that has been run), Claude may write
drivers, scripts, and runs. Exception: on 2026-10-05 the user asked Claude to
write the piece 6 library code directly to finish that day. Define new notation in place, and prefer naming
physical quantities over bookkeeping symbols. Size checks to what the next
result needs; keep checks that can change a conclusion. Heavy runs go through
Slurm; check `git status` first — this is a shared workspace.

The user asked for explanations built from the ground up: the piece 4 guide is
the model (define every object before using it, give the procedure of each
check, say what it does and does not establish). Math in markdown docs is
LaTeX, not Unicode.

The previous handoff (2026-09-24, amplitude-evolution interpolation checks) is
in git history; its state is superseded by
[adiabatic-branch-work-2026-09-28.md](docs/adiabatic-branch-work-2026-09-28.md).
