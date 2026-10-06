# Project handoff

Updated 2026-10-06. This is the authoritative project status: what works,
open caveats, and next steps. Campaign READMEs own exact commands, saved-result
paths, and detailed validation records; design documents own equations and
implementation decisions.

Active work: reproducing Omiya, Takahashi, Tanaka & Yoshino,
*Impact of multiple modes…* ([2211.01949v1](docs/papers/2211.01949v1.pdf)),
Sec. III — a nonlinear l = m = 1 primary with a linear l = m = 2 secondary on it,
evolved by energy and angular-momentum balance (Eqs. 62–63). The single-cloud
adiabatic work (OTT 2201.04382) is the primary's input; see
[adiabatic-repro.md](docs/adiabatic-repro.md).

Read first: [multimode-adiabatic-design.md](docs/multimode-adiabatic-design.md)
(equations, decisions, component map) and the
[secondary campaign README](runs/adiabatic/multimode_secondary/README.md)
(commands, results, checks). Cite equation numbers from the PDF in
`docs/papers/`, not the arXiv HTML.

## Where things stand

- **Multimode adiabatic evolution:** all six Sec. III components are implemented
  at $M\mu = 0.42$ and $0.3$, $a/M = 0.99$: averaged primary background,
  secondary eigenproblem and energy, driven $m = 0$ and $m = 3$ fluxes, and
  evolution from energy and angular-momentum balance. The nonperturbative
  Fig. 6 curves agree with the paper to reading precision at $M\mu = 0.42$;
  at $0.3$ the primary peak is about 10% lower and plateaus 5–10% higher.
  Fig. 7 profile peaks at $M\mu = 0.3$ agree to reading precision. See the
  [secondary campaign](runs/adiabatic/multimode_secondary/README.md) for the
  saved tables, comparisons, and checks. Pieces 5–6 are committed in `5b1f322`.
- **Single-cloud adiabatic work:** cosine, quartic, and sextic stationary
  branches are preserved in the [Figure 2 campaign](runs/adiabatic/fig2_20260922/README.md).
  The [amplitude campaign](runs/adiabatic/amplitude_evolution/README.md) has a
  merged 426-point derivative table and a saved cosine trajectory through
  $t/M = 4\times10^7$. Its approach to $A_0 \simeq 12.96268$ is a result on
  the saved spatial grid, subject to the high-amplitude caveat below.
- **Dolan time-domain work:** the coupled-tower solver, eigenvalue/Leaver
  comparison, and frequency-filter analysis have a
  [linear validation record](docs/dolan-repro.md) and
  [Fig. 13 campaign](runs/dolan/fig13/README.md). The nonlinear source and
  synchronized RK4 have [separate checks](docs/nonlinear-progress.md);
  those do not establish accurate collapse evolution.
- **Extensions:** coupled black-hole mass/spin evolution, GW losses, self-gravity,
  mutual nonlinear cloud response, additional independent clouds, and a PDE
  handoff remain proposals. The [October 6 research note](docs/further_work_1006.md)
  discusses their ordering and scientific motivation.

## Implementation decisions

The [multimode design](docs/multimode-adiabatic-design.md#decisions-and-their-reasons)
records the averaged reading of Eq. (61), restored $\mu^2$, secondary domain
and normalization, and the energy variables used for the balance equations.
The campaign records the $M\mu = 0.3$ tolerance setting, phase test, resonance
scan, and convergence evidence. Use its averaged secondary results;
`cosine_secondary_fixed_time.*` is retained only as a superseded experiment.

## Known caveats

- The Mμ = 0.42 primary branch is **not step-converged above A₁ ≈ 9**: its
  step-0.1 solutions re-evaluated at step 0.05 have scaled residual 4e-5 at
  A₁ = 9.1 rising to 4e-2 at 13. This bears on the cosine flux zero near
  A₁ ≈ 12.96. It does not affect Fig. 6 at Mμ = 0.42, which stays below
  E₁ ≈ 1000 (A₁ ≲ 0.8, step-converged to ~1e-7).
- Unchecked: secondary l-truncation (l ≤ 6) and outer boundary at large A₁.
  The paper reports neither.
- Full high-amplitude angular-quadrature and mode-truncation convergence of
  the primary remains unestablished. Interpolation and time-step checks on
  saved tables do not supply those spatial checks.
- Pieces 4–5: the outer-radius convergence of the $m = 0$ field is unchecked, and
  so is each coefficient against an independent evaluation of the perturbative
  $F_0$ and $F_3$ (Eqs. (38)–(39), (A11), (A14), e.g. by a Green's function).
  The power laws match; the prefactors have not been compared.
- The physical effect of the Sec. III approximations, including neglect of
  secondary self-interaction and its deformation of the primary, is unchecked.
- The saved two-cloud trajectories include scalar fluxes only. GW power,
  torque, and radiation reaction have not been implemented; physical energy
  normalization alone does not add them.

## Next steps

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
- **Single-cloud saturation:** refine the high-amplitude spatial solve before
  interpreting the saved cosine flux zero or saturation trajectory physically.
- **Proposed next research direction:** physical normalization of the saved
  two-cloud trajectories, with GW power/torque estimates and a frozen-profile
  self-gravity diagnostic, followed by a limited coupled-BH calculation
  including GW losses and a matched perturbative
  comparison. See the [October 6 proposal](docs/further_work_1006.md) for the
  rationale and conditional extensions; these are not completed capabilities.

## Working with the user

Follow [CLAUDE.md](CLAUDE.md) for incremental development, notation,
validation, and Slurm conventions. The
[piece 4 guide](docs/multimode-piece4-guide.md) is the model for explanations:
define each object, describe the check procedure, and state what it establishes.
Use LaTeX for math in Markdown documentation.

The September 24 handoff is in Git history. The
[September 28 branch-work record](docs/adiabatic-branch-work-2026-09-28.md)
preserves the following stage's experiments; its next-step list is historical.
