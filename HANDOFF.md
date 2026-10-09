# Project handoff

Updated 2026-10-08. This is the authoritative project status: what works,
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
- **Stationary background continuation:** ten primary-only solves at fixed
  $A_0=0.0011390625$ around $(\alpha,\chi)=(0.42,0.99)$ pass twenty observable
  comparisons. Offsets of $\pm0.001$ in alpha/spin converge in two or three
  Newton corrections; two routes to the joint downward offset agree in growth
  rate to $4.01\times10^{-10}$ relative. Halving the radial step changes growth
  by at most $1.91\times10^{-8}$ relative. See the
  [background-continuation campaign](runs/adiabatic/background_continuation/README.md).
  The same ten solves and twenty comparisons also pass at
  $A_0=2.051758324382913$, where the primary peak field is about 1.446 and
  $M\omega_I$ is 6.137 times its dilute value. Half-step changes in imaginary
  frequency and horizon flux are below $1.1\times10^{-8}$ relative; the small
  outward flux changes by $6.8\times10^{-5}$ relative. These are local primary
  stationary checkpoints, not a coupled BH trajectory; angular/endpoint checks
  remain separate work.
- **Secondary and driven background response:** the averaged secondary and
  $m=0,3$ fields now solve on all ten nonlinear-primary configurations above.
  The initial 133/135 passing checks include route agreement, secondary energy
  balance, and driven power identities. Two coarse-grid $m=3$ horizon-flux
  comparisons miss the $10^{-4}$ relative tolerance. Re-solving baseline and
  the joint offset at primary step 0.025 passes all 28 follow-up checks, with
  changes in that flux reduced from about $1.6\times10^{-4}$ to at most
  $5.10\times10^{-6}$ relative. Original failures are preserved. See the
  [response results and refinement](runs/adiabatic/background_continuation/README.md#secondary-and-driven-fields-on-the-finite-amplitude-backgrounds).
  This establishes local scalar response at fixed $A_0$, not BH–cloud evolution.
- **Coupled-balance prototype:** the leading small-cloud-mass scalar-only
  balance is derived in [the balance note](docs/coupled-bh-cloud-balance.md),
  with BH mass/angular momentum and both cloud energies in fixed initial-BH
  units. The campaign's `check_coupled_balance.jl` passes 44 checks against
  the existing nonlinear RHS, dilute channel counting, unit conversion, and
  a changing-background occupation-conservation example. The latter exposes
  the metric-work term required beyond the leading approximation. Actual
  nonlinear metric work and cloud angular-momentum background derivatives
  remain unevaluated; the four-variable RHS is a driver prototype, not a
  time-integrated coupled trajectory. See the
  [checkpoint recipe and illustrative RHS](runs/adiabatic/background_continuation/README.md#coupled-balance-checkpoint--october-8-2026).
- **Dolan time-domain work:** the coupled-tower solver, eigenvalue/Leaver
  comparison, and frequency-filter analysis have a
  [linear validation record](docs/dolan-repro.md) and
  [Fig. 13 campaign](runs/dolan/fig13/README.md). The nonlinear source and
  synchronized RK4 have [separate checks](docs/nonlinear-progress.md);
  those do not establish accurate collapse evolution.
- **Self-gravity benchmark:** a standalone Newtonian Poisson check on an
  analytic hydrogenic $211$ density recovers the monopole/quadrupole shift
  $\delta\omega/\mu=-(237/1280)\alpha^2 M_c/M_{\rm BH}$. All 16 checks pass,
  including radial refinement, angular/outer-radius checks, and an independent
  exact radial integral. See the [derivation and command](docs/adiabatic-self-gravity-review.md#211-benchmark-checkpoint--october-7-2026).
  Saved nonlinear profiles and scalar evolution still have no gravity correction.
- **Dilute self-gravity screening:** published hydrogenic self/cross shifts
  have been evaluated on both saved two-cloud energy trajectories. At
  $f_a=10^{16}$ GeV, peak $|\delta\omega_1|/\mu$ is $7.70\times10^{-6}$
  (alpha 0.42) and $5.80\times10^{-6}$ (0.30); the outward scalar threshold
  margin changes by at most 0.0219% and 0.0483%. At $10^{17}$ GeV these
  changes become 2.19% and 4.83%. See the
  [campaign recipe and limits](runs/adiabatic/multimode_secondary/README.md#dilute-self-gravity-estimates--october-7-2026).
  This is an analytic-coefficient estimate, not the nonlinear frequency
  response or a bound on scalar flux errors. BH evolution remains separate.
- **Extensions:** coupled black-hole mass/spin evolution, GW losses, self-gravity feedback,
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

- **Coupled-BH preparation:** local primary continuation in alpha and spin
  now works at both dilute amplitude and $A_0\simeq2.052$; the secondary and
  driven channels are checked on the nonlinear backgrounds, with the targeted
  refinement described above. Use primary step 0.05 or finer for the tested
  small $m=3$ horizon-flux accuracy requirement.
  The leading balance is now derived and its algebra checked. Its functions
  now live in `src/adiabatic/multimode/coupled_balance.jl`, with units and
  equations documented in the source. The user stopped the time-integrator
  draft to require small, explained implementation steps. Stationary-table
  job 26757461 was cancelled; `runs/adiabatic/coupled_local/` contains paused
  drafts, not a validated trajectory. Review the balance piece together before
  adding state conversion and one time step. Subsequently build
  a limited trajectory with refreshed stationary data, fixed reference units,
  and specified physical seeds/$f_a$, monitoring small cloud mass and regular
  response. The fixed outer amplitude used in the stationary check is not
  an evolution prescription. Beyond the leading approximation, evaluate
  metric work, background angular response, and consistent gravitational
  energy/boundary accounting as described in the balance note.
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
