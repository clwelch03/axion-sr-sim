# Branch driver merge, dense derivatives, and cosine flux zero (2026-09-25 to 09-28)

This note covers one working session. It started from Codex's question "are we
ready for a physical simulation?" and ended with two completed Slurm runs.
Numbers are from the saved outputs; job IDs are listed where they apply.

## 1. Why the sampling had to be denser

The flux-derivative check (job 26400508, `amplitude_evolution/results/hermite_flux_checks.csv`)
tested cubic Hermite interpolation of the 178-point cosine table against 12
independently solved reference points:

| Quantity | Within 0.1% | Worst |
|---|---:|---:|
| Hermite total flux | 12/12 | 0.021% |
| Hermite energy derivative | 8/12 | 0.157% |
| Amplitude rate | 8/12 | 0.157% |

The flux was accurate enough, but the energy derivative was not. The failures were all in sparse
intervals. The fix was to resample so that `ΔA0 <= min(0.25, 0.05 A0)`. Codex
had added that rule to the cosine continuation driver; the existing branch did
not satisfy it in 29 intervals.

Rerunning the continuation from scratch was rejected: it took several 4-hour
jobs to build, it would overwrite tracked Fig. 2 results, and a derivative pass
would still be needed afterwards. Filling the gaps reuses the 178 solved points.

## 2. Continuation drivers moved and merged

### Move (commit `aa17cb5`)

`continue_cosine.jl` and `continue_polynomial.jl` were general drivers living
in `runs/adiabatic/fig2_20260922/scripts/`. They were moved with `git mv` to
`src/adiabatic/drivers/`. Their output defaults had been relative to their own
folder, so they were rewritten as explicit `fig2_20260922/...` paths. Eleven
references (submit scripts, campaign scripts, diagnostics) and the docs were
updated.

Check: both drivers resolved the same output and seed files as before (178
cosine points, 126 quartic points), and all rewritten `include` paths exist.

`aa17cb5` is the last commit where the seven `fig2_20260922/diagnostics/`
scripts and `solve_quartic_target.jl` / `fill_quartic_gap.jl` run unchanged.
Those finished experiments were not rewritten for the new driver; reproduce
them by checking out that commit.

### Merge (commit `ff1ce83`)

`src/adiabatic/drivers/continue_branch.jl` replaces both drivers:

```bash
AXION_TARGET_AMPLITUDE=13 julia --project=. src/adiabatic/drivers/continue_branch.jl \
    cosine|quartic|sextic  seed.jls  output_prefix
```

It was written in five pieces, each checked before the next:

| Piece | Function | Check | Result |
|---|---|---|---|
| 1 | `load_branch_setup` | settings, grid, couplings equal the old driver's globals | all equal |
| 2 | `outward_flux`, `branch_diagnostics` | recompute a saved cosine row and a saved quartic row | identical |
| 3 | `save_branch` | rewrite the cosine branch to scratch | CSV byte-identical; `.jls` contents identical |
| 4a | `start_branch` | seed history, resume, and quartic re-solve from the cosine seed | quartic first row identical to production |
| 4b | `continue_branch` | quartic run with old settings vs saved `quartic_branch` rows 1-5; cosine run vs old `continue_cosine.jl` | both identical |
| 5 | `main` | command line vs function call; missing target refused | identical; refused |

Behaviour that differs from the old drivers:

- The step limit `min(maximum_amplitude_step, maximum_fractional_step * A0)` is
  applied before every attempt, for all potentials (defaults 0.25 and 0.05).
  The old polynomial driver only limited growth to `0.25 A0`.
- An existing output is always resumed; the seed is used only for a fresh
  output. To start over from a different seed, use a new output prefix.
- `AXION_TARGET_AMPLITUDE` is required. The old defaults (13 for cosine, 0.1
  for polynomial) differed.
- Couplings are built up to the largest ℓ in the modes instead of a fixed 5.
  Same result for every saved branch.
- Checkpoints record the step limits.

## 3. Dense cosine seed

`runs/adiabatic/amplitude_evolution/scripts/prepare_dense_seed.jl` wrote
`seeds/cosine_dense_20260925.jls`: 372 points, the 178 originals plus 194
inserted, evenly spaced within each interval that was too wide. Inserted guesses are Hermite
interpolants of the normalized state using the saved state derivatives. The
seed also stores, for each inserted point, the Hermite prediction of energy,
energy derivative, and fluxes from the 178-point table.

Check: originals retained exactly, amplitudes strictly increasing, and the
largest step equals the allowed step (ratio 1.0; several old intervals already sat exactly on
the limit).

## 4. Runs

### Dense derivatives: job array 26434701

`branch_energy_derivatives.jl` on the dense seed, 4 shards, 1 CPU each.
All tasks completed (18-27 min each, 5024 CPU-seconds total). Output:
`amplitude_evolution/results/cosine_dense_derivatives_part{0..3}.{csv,jls}`.

- 372/372 points, residual <= 9.2e-11, at most 3 Newton corrections.
- Every energy derivative is positive and every output finite.
- The 178 original points reproduce the earlier table exactly (energy,
  derivatives, and fluxes).

How accurate the old 178-point table was, measured at the 194 new points:

| Quantity | Within 0.1% | Worst |
|---|---:|---:|
| Energy derivative | 174/194 | 0.156% at A0 = 2.46, interval [2.05, 2.56] |
| Total flux (Hermite) | 194/194 | 0.020% |
| Amplitude rate | 174/194 | 0.156% |

The 20 failures lie between A0 = 0.089 and 4.84, in the old sparse intervals.

How accurate the dense table is: each interior point was predicted by Hermite
interpolation from its two neighbours, which is twice the dense spacing.

| Quantity | Within 0.1% | Worst |
|---|---:|---:|
| Energy derivative | 370/370 | 0.0036% at A0 = 3.21 |
| Amplitude rate | 370/370 | 0.0036% |

The energy-derivative error of cubic Hermite scales as h³, so at the actual
spacing the error is roughly 8 times smaller than these doubled-spacing values.
The dense table meets the 0.1% target with a wide margin. This is interpolation
accuracy on a fixed spatial grid, not radial or angular convergence.

### Cosine continuation: job 26434700

`submit_long.sh`, now calling `continue_branch.jl`. Resumed from A0 = 11.985
with the cubic predictor (history spacing 0.025), ForwardDiff Jacobians, and
tolerance 1e-9. It completed in 1 h 36 min at the target A0 = 13 and redrew
`plots/preliminary_fig2_20260922.png`.

- 54 new points (232 total), residual <= 9.6e-10. There were 26 failed attempts,
  each followed by a halved step. Several failures stopped at residuals of
  2.6e-9 to 5.9e-9, close to the 1e-9 tolerance, which matches the earlier
  high-amplitude residual-floor observations.
- Energy stays monotonic in A0.
- **The net flux changes sign once**, between A0 = 12.959108
  (total_flux = -2.97e-5) and 12.962707 (+2.09e-7). A linear estimate puts the
  zero at A0 ≈ 12.9627. With `total_flux = horizon_flux + outer_flux` as the net
  loss, this is the point where superradiant gain and emission balance.

## 5. Filesystem note

Late in the session, node n0002 could not open some recently rewritten files
(`cosine_dense_derivatives_part2.jls` and several files under `.git/`):
"Cannot send after transport endpoint shutdown". The same files read normally
from compute nodes, and the cluster status page reported no problems. So
n0002's filesystem client has stale handles for those files; the data is
intact. The assessments above ran on compute nodes through `srun`. If `git`
fails on n0002, run it from a compute node or wait for the client to recover.

## 6. Not done yet

- **Commit** the run outputs and recipe changes (listed by `git status` from a
  compute node): the updated `cosine_branch.{csv,jls}`, `submit_long.sh`,
  `prepare_dense_seed.jl`, the dense seed, and the four dense derivative shards.
- **Energy derivatives above A0 = 11.985.** The dense derivative table comes
  from the frozen 2026-09-24 branch copy, so it stops at 11.985. The 54 new
  cosine points, including the flux zero, have no derivatives yet. Time
  evolution through saturation needs them.
- **Merge the shards** into one table, and save the assessment scripts, which
  are currently in the ignored `amplitude_evolution/work/` folder.
- **Switch the remaining callers** (`submit.sh`, `submit_polynomial.sh`,
  `submit_refinement.sh` with `refine_polynomial.jl`) to `continue_branch.jl`,
  then remove `continue_cosine.jl` and `continue_polynomial.jl`.
- **Move the Hermite state guess** into `src/adiabatic/amplitude_evolution/interpolation.jl`.
  It is currently duplicated in `prepare_dense_seed.jl` and
  `tests/adiabatic/check_branch_interpolation.jl`.
- **Update the READMEs.** The `amplitude_evolution` and `fig2_20260922` READMEs
  and the design doc do not yet describe `continue_branch.jl`, the dense run,
  or the `aa17cb5` freeze point.
- **Spatial convergence** at large amplitude (radial step, outer radius, mode
  truncation, angular grid) is still unchecked. Interpret the flux zero near
  A0 ≈ 12.96 only after that check.
