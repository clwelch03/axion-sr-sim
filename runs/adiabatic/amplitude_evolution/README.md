# Amplitude evolution checks

## Physical interpolation checks — 2026-09-24

`tests/adiabatic/check_branch_interpolation.jl` prepares quarter-, mid-, and
three-quarter points in four saved intervals: [0.0820125, 0.12301875],
[2.56469791, 3.20587238], [6.26146950, 7.82683687], and
[11.80053105, 11.87454326]. The last interval supplies a finely sampled
comparison. Off-center checks avoid relying on possible midpoint cancellation.

The interpolated predictions are saved before solving the new points. Initial
guesses interpolate the normalized endpoint states; the derivative driver then
corrects each guess to a stationary solution and computes fresh implicit
derivatives, energies, and fluxes. It retains the endpoint table's spatial
settings and matching tolerance. The new points are separate validation data;
they do not modify the original 178-point table.

```bash
julia --project=. tests/adiabatic/check_branch_interpolation.jl prepare
sbatch runs/adiabatic/amplitude_evolution/scripts/submit_interpolation_checks.sh
# After all four tasks finish:
julia --project=. tests/adiabatic/check_branch_interpolation.jl assess
```

Preparation refuses to overwrite its frozen seed. Job array **26351040** runs
four one-CPU tasks, each capped at 30 minutes. Its input is
`seeds/interpolation_checks_20260924.jls`, and its outputs are
`results/interpolation_checks_part{0,1,2,3}.{csv,jls}` with logs in
`logs/interpolation-26351040_{0,1,2,3}.out`. Assessment verifies the source hash
and complete coverage, then writes `results/interpolation_checks.csv`.
It reports the energy derivative against the requested 0.1% target and also
compares the growth rate to that threshold. These are interpolation checks
on a fixed spatial grid, not radial/angular convergence tests.

### Completed assessment

The original array accepted nine reference points. Newton corrections from
linear state guesses failed for the three points in the highest-amplitude
interval. `prepare-retry` generated guesses for just those missing points using
cubic Hermite interpolation of the saved states and implicit state derivatives:

```bash
julia --project=. tests/adiabatic/check_branch_interpolation.jl prepare-retry
```

This saved `seeds/interpolation_retry_20260924.jls`. Retry job **26353282** ran
the derivative driver on that seed with output prefix
`results/interpolation_retry`; all three points converged in two Newton
corrections, and the job completed successfully in 2m08s. The retry changed
only the solver guesses: the observable predictions under test were preserved.
Assessment checks both input hashes and combines all 12 reference points.

The completed [comparison table](results/interpolation_checks.csv) gives:

| Quantity | Points within 0.1% | Maximum relative discrepancy |
|---|---:|---:|
| Hermite energy derivative | 8/12 | 0.1570% |
| Amplitude rate with linear flux interpolation | 4/12 | 1.6354% |

All three checks in the finest interval pass both thresholds. Sparse intervals
show errors in both the energy derivative and the interpolated flux. A passing
rate can also result from cancellation between those errors, so it does not
independently validate either ingredient. The interpolation is not yet validated
for production at the requested accuracy.

Next experiment: compute flux amplitude derivatives from the already saved state
derivatives and test Hermite flux interpolation against these same references.
This will isolate the benefit of improving flux interpolation. Energy errors in
sparse intervals will still require attention. No time integration or spatial
convergence claim follows from the present checks.

## Cosine branch derivatives — 2026-09-24

The reusable driver `src/adiabatic/drivers/branch_energy_derivatives.jl`
computes the implicit energy derivative at each saved stationary solution.
It checks/corrects matching to tolerance 1e-10, then recomputes energy and
signed finite-radius fluxes with the same spatial settings. Outputs include
`energy_derivative`, `amplitude_rate = -total_flux/energy_derivative`, and
`fractional_growth = amplitude_rate/amplitude`. An exactly zero derivative
is retained with an undefined (`NaN`) amplitude rate; no turning point is
regularized or crossed. No time integration is performed.

Input is frozen in `seeds/cosine_branch_20260924.jls`, copied from the
Figure 2 cosine checkpoint. Job array **26337413** was submitted with:

```bash
sbatch runs/adiabatic/amplitude_evolution/scripts/submit_derivatives.sh
```

The four array tasks each use one CPU for at most one hour. Task `k` reads
entries `k+1:4:end`, so every input point is covered once. Per-point results
and corrected states are checkpointed to `results/cosine_derivatives_part{k}`
as CSV/JLS. Logs are under `logs/derivatives-26337413_{k}.out`. The JLS metadata
records the input SHA256, physical and numerical settings, and shard selection.
Restarting the same task resumes its accepted points after checking metadata;
never start two writers for the same shard. Failures stop that task and remain
visible in the scheduler/log, preserving its previously accepted results.

The source revision is `3a3cc8e` plus the working-tree sensitivity functions,
residual extraction, driver, and submission recipe. These calculations retain
the preliminary spatial discretization of the saved branch.

All four tasks completed successfully (exit 0); elapsed times ranged from
8m12s to 13m13s. The merged table contains all **178** input amplitudes exactly
once, from 0.0011390625 to 11.985191522619052. The maximum matching residual is
9.18428e-11; 135 points required additional corrections to meet 1e-10.
All computed energy derivatives and amplitude rates are finite and positive.
The largest sampled amplitude rate is 0.0004160882172 at A0=9.235896808;
at the final point it is 0.0002187388637. The saved branch has not reached
zero net flux.

Against the implicit derivative on the same grid, the maximum absolute
relative differences of the three-point quadratic, five-point cubic, and
seven-point cubic fits are respectively 3.8383%, 1.1638%, and 6.7032%.
All three maxima occur at A0=2.564697905. The independent centered-difference
validation below applies at that point; it has not been repeated at every
amplitude, and spatial convergence remains to be checked.

After all tasks finish, merge their CSVs and plot the growth rates and local-fit
derivative differences with:

```bash
python runs/adiabatic/amplitude_evolution/scripts/summarize_derivatives.py
```

Use the BHsuper Python environment. The script requires every amplitude in
`seeds/cosine_branch_20260924.csv` exactly once and verifies matching residuals
and finite outputs. It writes `results/cosine_derivatives.csv` and
`plots/amplitude_growth_implicit_comparison.{png,pdf}`. Local polynomial fits
use the newly recomputed energies, so their comparison isolates derivative
estimation rather than differences in matching tolerance.

## Implicit energy derivative — 2026-09-24

The check driver is `tests/adiabatic/check_branch_sensitivity.jl`. It reads
the saved cosine branch without modifying it and tests the implicit state
and energy derivatives against centered differences of separately corrected
stationary solutions. The corrected solutions use the baseline state as
their initial guess, without using the implicit derivative as a predictor.

From the repository root:

```bash
mkdir -p runs/adiabatic/amplitude_evolution/work
sbatch --job-name=axion_sensitivity --account=pc_heptheory --partition=lr6 --qos=lr_normal --nodes=1 --ntasks=1 --cpus-per-task=1 --time=00:20:00 --output=runs/adiabatic/amplitude_evolution/work/sensitivity-%j.out --wrap='export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1; /global/home/users/clwelch/.juliaup/bin/julia --project=. tests/adiabatic/check_branch_sensitivity.jl'
```

Input: `runs/adiabatic/fig2_20260922/results/cosine_branch.jls`, selecting the
accepted amplitude nearest 2.564697905478641. Uses the checkpoint's settings:
spin 0.99, alpha 0.42, six modes through l=m=5, 16 polar and 31 azimuthal
samples, radial step 0.1, inner -100, match 5, outer 150. Correction tolerance
is 1e-10 with ForwardDiff Jacobians. This is a single-point differentiation
check on a fixed spatial discretization, not a continuum convergence result.

Job 26335774 passed all 13 checks in 148.7 seconds of test execution. Code:
base revision `3a3cc8e` plus the working-tree residual extraction, branch
sensitivity function, and check driver. Output is the ignored log
`work/sensitivity-26335774.out`. An earlier job, 26335564, failed because of
a variable-scope collision in the check driver; that collision was corrected
before the successful run.

The implicit derivative is **d(E/Fa²)/dA0 = 80.73934435268751**.

| Centered-difference spacing | Energy derivative | Relative discrepancy |
|---|---:|---:|
| 0.025646979054786412 | 80.74557757340474 | 7.72018e-5 |
| 0.012823489527393206 | 80.74090259949655 | 1.92997e-5 |
| 0.006411744763696603 | 80.73973389692831 | 4.82471e-6 |

Discrepancy decreases by approximately four per halving of the spacing.
The finest finite-difference state derivative agrees with the implicit
state derivative to 3.42123e-6 in relative vector norm. Testing other branch
points and spatial convergence remains necessary before time evolution.
