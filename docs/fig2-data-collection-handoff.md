# Figure 2 data collection — Astra handoff

Status checked on 2026-09-23. Repository:
`/global/scratch/projects/pc_heptheory/clwelch/axion-sr-sim`.

## Objective and working agreement

Reproduce **both panels of Omiya, Takahashi & Tanaka, arXiv:2201.04382v3,
Figure 2**, at `spin=0.99`, `alpha=0.42`:

- Energy versus fundamental outer amplitude `A0`.
- Total flux versus energy, retaining amplitude order.
- Quartic (black), sextic (red), and cosine (blue) potentials.

The user wants useful plot data soon, with targeted numerical checks rather
than an exhaustive validation campaign before every step. They are learning
the method and have explicitly returned to **code snippets they can place
themselves**, most recently while adding the polynomial potentials. Give exact
placement instructions, match existing docstrings, and work incrementally.
Do not treat earlier one-off permission for autonomous implementation as
blanket permission for further code edits. Check in before adding code.
Writing this handoff was explicitly requested. Preserve the user's edits.

The cosine diagnosis was authorized for autonomous diagnostic runs. Those
runs are now finished. Do not silently promote their experimental methods
into production or launch another long continuation without explaining the
concrete next step. Read `CLAUDE.md` for project conventions.

## Current repository state

Latest commit at inspection: `d695fb5` — `feat: improve adiabatic continuation
and record validation results`. Earlier energy/continuation commit: `f66bb0b`.
At inspection, the only modified tracked file was `src/adiabatic/nonlinear.jl`:
the user has added the two potential helpers from the supplied snippets.
Recheck `git status`; this is a shared workspace.

`potential_value(field, potential=:cosine)` and
`potential_remainder(field, potential=:cosine)` now exist immediately after
`sine_remainder`. **They are not wired into the solver or energy yet.**
`nonlinear_source` still calls `sine_remainder.(field)`, and `energy.jl` still
uses the cosine potential directly. No quartic/sextic branch has been collected.

Older `docs/adiabatic-design.md` and `docs/adiabatic-repro.md` have useful
conventions but stale progress sections. For recent numerical results use
this handoff, the run README, and the saved artifacts.

## Data already available

All paths below are relative to `runs/adiabatic/fig2_20260922/`
unless stated otherwise.

| Artifact | Contents |
| --- | --- |
| `results/cosine_branch.csv` | Main preliminary cosine branch, 106 accepted points |
| `results/cosine_branch.jls` | Full solutions, boundary amplitudes, settings, and rows |
| `seeds/seed_before_residual_check.jls` | Preserved copy of the 106-point branch used for diagnosis |
| `diagnostics/high_amplitude/curved_prediction_diagnosis.jls` | Separate successful point at exactly `A0=10`; **not appended to the main branch** |
| `scripts/plot_branch.py` | Plots main cosine CSV only; does not include the separate A0=10 point |
| `README.md` | Existing recipe, energy/flux conventions, earlier speed measurements |

Main branch endpoint:

```text
A0                 9.949017521664452
E/Fa²              1678.6572279955956
Ftot/Fa²           -0.012573139175095345
omega              0.37391947581441826 + 3.996742112434326e-6im
scaled residual    5.548628126938361e-10
sampled peak field 2.9437412218936405
```

The sampled flux minimum is near `A0=9.14303`, `E/Fa²=1653.54807`,
`Ftot/Fa²=-0.01300184`. The curve has begun turning upward, but has **not**
reached zero flux or the desired cosine extent near `A0=13`.

Current images:
`plots/preliminary_cosine_branch_20260922.png` and `.pdf`.
These are preliminary six-mode, finite-boundary results, not a completed
or fully converged Figure 2 reproduction.

CSV columns:
`amplitude,energy,normalized_energy,omega_real,omega_imag,horizon_flux,outer_flux,total_flux,balance_ratio,peak_field,residual,iterations`.
`normalized_energy` means `energy/A0^2`; plot `energy` for Figure 2.

## Numerical settings and conventions to preserve

- Modes: `[(1,1),(3,1),(5,1),(3,3),(5,3),(5,5)]`; fundamental `m0=1`.
- Angular quadrature: 16 polar and 31 azimuthal points.
- Production radial RK4 step: `0.1` in tortoise coordinate.
- Endpoints: `inner=-100`, `outer=150`, `match=5`.
- Boundary series defaults: horizon order 3, outer order 7.
- Matching tolerance: `norm(residual/A0) <= 1e-9`.
- State: `[R; Q]`, with `Q=dR/dr_*`; physical field is
  `2real(sum(R_lm Y_lm exp(im*m*phi)))`. Polar harmonics integrate to one
  over `dcos(theta)`; azimuthal weights integrate `dphi/(2pi)`.
- Energy uses composite Simpson quadrature on each radial branch; descending
  outer profiles are reversed. Uniform grids and even interval counts required.
- The source already separates out the linear mass term. Potential changes
  must change **both** the nonlinear remainder and the energy potential.
- Horizon flux is minus the finite-radius outward current; total flux is
  horizon plus outer flux. For `m0=1`, the outward current is
  `4pi*(r^2+a^2)*sum(real(im*m*omega*R*conj(Q)))`.
- `balance_ratio=-Ftot/(2imag(omega)*energy)` is a linear-limit check, not an
  identity to enforce along the nonlinear branch.
- The paper's reported outer endpoint is 100; this preliminary run uses 150.
  The finite inner cutoff is part of the current energy convention.

Earlier checks established fourth-order energy convergence at small amplitude
(refinement ratio 15.99) and negligible outer-endpoint change from 150 to 200
there. They do not establish high-amplitude spatial/angular convergence.

## Why the long run slowed, and what the diagnosis actually found

Production job `26286366` completed normally after 3h51m27s, stopping at its
internal runtime budget. It reached A0≈5 in 14 minutes, then progressively
shrunk amplitude steps above roughly 8. Near 10 it repeatedly failed descent
around the residual tolerance. Do not extrapolate the early speed into a
promise for the remaining branch.

The production predictor is still a **two-point secant** in physical boundary
amplitudes and frequency. The driver still defaults to a uniform Jacobian
difference scale of `1e-8` (library default `1e-6`). Its Newton solve has no
fixed iteration cap, but retains descent and stagnation guards.

Two separate limitations were demonstrated:

### 1. Jacobian differences are too coarse locally

At the saved A0≈9.949 point, residual evaluations repeat exactly. Perturbing
each real-packed input upward by one floating-point increment changes the
residual by at most `4.66e-13`, well below `1e-9`. This argues against a basic
Float64 resolution floor at that state; it is not a general precision proof.

Starting from the *same* saved state (residual `5.55e-10`), one full Newton
correction gives:

| Difference scale | Residual after correction |
| --- | --- |
| `1e-8` | `2.49e-9` — worse |
| `1e-9` | `3.04e-11` |
| `1e-10` | `4.26e-12` |
| `1e-11` | `4.06e-12` |

The largest derivative discrepancy is in the real/imaginary `(5,1)` horizon
coefficient columns: about 0.8–1% at `1e-8`, about `1e-4` at `1e-9`, and
about `1e-6` at `1e-10`, relative to the `1e-11` calculation. Other columns
show roundoff sensitivity at the smallest differences. Smaller is not
universally better; the useful tested compromise here is `1e-10`.

### 2. The secant starting guess also becomes inadequate

From the same last two points, the secant prediction at A0=10 starts at
residual `0.0400167`. Even `1e-10` differences do not rescue this guess within
the diagnostic's 12-correction budget (residual still about `0.01069`).
This is a bounded diagnostic failure, not proof that further Newton
iterations could never converge. A0=10.5 was **not attempted** by that script
because A0=10 did not pass.

A cubic extrapolation of physical amplitudes and frequency, using four
more widely separated historical points, **did** work:

```text
history A0: [9.791649280567519, 9.830827653295797,
             9.899017025202976, 9.949017521664452]
target A0:  10.0
difference scale: 1e-10
initial residual: 9.012459395094593e-4
corrections: 3
final residual: 4.1597797160037053e-13
energy: 1680.238109177804
total flux: -0.012511883387398089
frequency: 0.37380575782568665 + 3.974191089535886e-6im
measured solve plus diagnostics: 38.36 seconds
```

The cubic script picks the latest history point at or below
`last_A0 - offset`, for offsets `(0.15, 0.10, 0.05, 0.0)`, and uses Lagrange
extrapolation weights. It fixes the first outer amplitude exactly to target.
This is a **local diagnostic**, not a general-purpose predictor: hard-coded
offsets need guards/fallbacks for short or differently spaced histories.

Moving the matching point to 10 or 20 with the old `1e-8` differences did
not solve the A0=10 problem; both failed descent. These tests do not isolate
all combinations of matching point and derivative scale. There is currently
no reason to change the production matching location on this evidence.

Artifacts and completed jobs:

| Script / result stem | Job | Outcome |
| --- | --- | --- |
| `diagnostics/high_amplitude/check_residual_floor.jl` / `diagnostics/high_amplitude/residual_floor_diagnosis.jls` | 26299471 | Completed; recorded unsuccessful solve trials and sensitivity probes |
| `diagnostics/high_amplitude/check_jacobian_increment.jl` / `diagnostics/high_amplitude/jacobian_increment_diagnosis.jls` | 26299535 | Completed; local derivative comparison |
| `diagnostics/high_amplitude/check_finer_correction.jl` / `diagnostics/high_amplitude/finer_correction_diagnosis.jls` | 26299562 | Completed; secant A0=10 did not converge within 12 corrections |
| `diagnostics/high_amplitude/check_curved_prediction.jl` / `diagnostics/high_amplitude/curved_prediction_diagnosis.jls` | 26299798 | Completed; cubic A0=10 converged |

Scheduler exit zero does not mean every numerical trial succeeded: the
diagnostic scripts intentionally catch and record numerical failures.
All these jobs, including the production job, were finished at handoff.

## Next steps: collect the remaining cosine curve

1. Explain the two findings above. Supply a small, reviewable predictor
   change as a snippet, respecting the user's current workflow. Keep the
   existing secant fallback, enforce distinct history nodes and normalization,
   and account for unequal spacing. Do not blindly use the last four tightly
   clustered points to estimate curvature.
2. Preserve the main 106-point checkpoint. The successful A0=10 solution is
   available through `deserialize(...).trial.solution`; its amplitude is
   `.trial.target`, and `.trial.integration_options` records the radial settings.
   The result file is **not** a production continuation checkpoint. To use
   it as a seed, explicitly construct an accepted entry and compute the
   complete diagnostic row before saving a separate continuation checkpoint.
3. Try a modest next target such as 10.05 or 10.1 with the curved predictor,
   `difference_step=1e-10`, and unchanged `1e-9` acceptance tolerance.
   Compare one modest jump against subdivided steps at a common endpoint,
   checking frequency, energy, flux, and boundary amplitudes for branch agreement.
   The A0=10 success alone has not validated subsequent large jumps.
4. **Watch the restart rule:** `scripts/continue_cosine.jl` currently forces its first
   increment into 10–25% of A0. At A0≈10 this jumps by at least 1, far beyond
   the successful diagnostic step. Do not simply relaunch `scripts/submit_long.sh`
   and assume it uses the tested cubic predictor/small restart step.
5. Once the next-step check passes, collect toward A0≈13 and through the
   upward flux zero crossing, checkpointing every accepted point. Keep failed
   attempts out of the plotted data. Record predictor, difference scale,
   tolerance, grid, boundaries, and potential; distinguish historical settings
   when extending a checkpoint with a new method.
6. Refine sampling around the flux minimum/zero crossing as needed for a
   readable figure. Do a focused high-amplitude resolution check before
   claiming quantitative reproduction. A small matching residual alone is
   not a spatial or angular convergence result.

## Next steps: wire in quartic and sextic, then collect independent branches

The user has added these helpers but has not yet connected them:

| Potential | `potential_value` = U | `potential_remainder` = phi − U′ |
| --- | --- | --- |
| cosine | `2sin(phi/2)^2` | stable `sine_remainder(phi)` |
| quartic | `phi^2/2 - phi^4/24` | `phi^3/6` |
| sextic | `phi^2/2 - phi^4/24 + phi^6/720` | `phi^3/6 - phi^5/120` |

Give the following changes as incremental snippets, not automatic edits:

1. Add a `potential=:cosine` keyword to `nonlinear_source`; replace its
   sine-only remainder by `potential_remainder.(field, Ref(potential))`.
   Keep its existing `alpha^2*rho_squared` factor and radial source sign.
2. Forward that keyword through `nonlinear_tortoise_rhs!`,
   `nonlinear_matching_residual`, and `solve_nonlinear_mode` as appropriate.
   The last function already forwards `integration_options...`; inspect
   signatures before adding redundant plumbing.
3. Add the same keyword to `energy_per_tortoise_radius` and
   `integrate_energy_branch`; compute `alpha^2 .* potential_value.(field,
   Ref(potential))`. Rename the existing local array named `potential` so it
   does not overwrite the selector. Forward the selection when rebuilding
   profiles and computing diagnostics, as well as during the solve.
4. Keep cosine as the default and verify existing checks still pass. Add a
   focused consistency check between U and its nonlinear remainder, and an
   analytic polynomial angular-source/energy check, following existing tests.
5. Generalize the run recipe so source, energy, metadata, and plot label use
   the same selection. Make separate quartic/sextic output prefixes and check
   potential compatibility on restart. Never append a new potential to the
   cosine CSV or merely relabel a cosine solution.
6. Start each new potential near the linear regime, using a small-amplitude
   solution as a **guess**, then re-solve with the selected potential before
   accepting it. The old `runs/adiabatic/small_amplitude/results/continuation.jls` and early cosine
   checkpoints contain small-amplitude seeds; inspect their settings first.
7. Collect the three branches independently. The polynomial source is cheaper
   to evaluate, but faster Newton convergence has not been established.
   The quartic energy curve turns over: preserve A0 order when plotting
   Ftot against E, rather than sorting by energy.

## Running and inspecting

Heavy solves go through SLURM: account `pc_heptheory`, partition `lr6`,
QOS `lr_normal`, one CPU for the existing recipes. Set
`OPENBLAS_NUM_THREADS=1` and `JULIA_NUM_THREADS=1`.
Scheduler access may require tool escalation because of sandbox networking.
Use shared run directories, not node-local `/tmp`, for checkpoints and scripts.
Never run two writers against the same checkpoint prefix.

Julia executable: `/global/home/users/clwelch/.juliaup/bin/julia`.
Run with `--project=.` from the repository root. Full existing suite:

```bash
/global/home/users/clwelch/.juliaup/bin/julia --project=. tests/adiabatic/runtests.jl
```

The full suite passed before the latest potential-helper additions; those
user additions were not tested during this handoff. Earlier speed comparison
also passed 16 assertions (see run README).

To inspect the successful separate point without re-running a solve:

```julia
using Serialization
result = deserialize("runs/adiabatic/fig2_20260922/diagnostics/high_amplitude/curved_prediction_diagnosis.jls")
result.trial.target
result.trial.solution
result.trial.energy
result.trial.total_flux
```

Some diagnostic scripts execute their experiment when included; inspect them
or load their `.jls` results instead of casually including all scripts.
The production `scripts/continue_cosine.jl` guards its `main()`, but loads its seed
and global settings at include time. Set environment overrides before inclusion.

To redraw the current main CSV:

```bash
MPLCONFIGDIR=/tmp/axion_mpl_cache \
  /global/home/users/clwelch/.conda/envs/BHsuper/bin/python \
  runs/adiabatic/fig2_20260922/scripts/plot_branch.py
```

The plotting script currently reads only `results/cosine_branch.csv`; a separate
extended branch or polynomial branch requires an explicit plotting update.
Save standalone PNG/PDF artifacts and retain the numerical tables.
