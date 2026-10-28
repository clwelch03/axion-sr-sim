# Adaptive radial boundary-value experiment

An isolated comparison with the shooting investigation in
[`docs/multiphase-target-solve-investigation-2026-10-28.md`](../../../docs/multiphase-target-solve-investigation-2026-10-28.md).
The motivation and literature review are in
[`docs/radial-solver-methods-review-2026-10-28.md`](../../../docs/radial-solver-methods-review-2026-10-28.md).

**Result: both direct target-amplitude guesses converged in three Newton
corrections**, where the recorded shooting attempts did not converge. No
amplitude continuation was needed for these two runs. Their final frequencies
agree within `3.31e-13`, and the largest per-mode relative profile difference
is `2.74e-8` (including slopes). A tighter spatial refinement also converged.

## Formulation

The experiment uses fourth-order adaptive MIRK collocation from
`BoundaryValueDiffEqMIRK` 1.19.2. Its source equations, angular quadrature,
boundary series, and reference data come from the existing frozen
`multiphase_two_cloud/work/amplitude_continuation_20261028_v1` snapshot.
It retains all 14 complex modes, both complex eigenfrequencies, the cosine
potential, the 16×31×31 angular/phase grid, and the finite radial interval
`r* ∈ [-100,300]` with boundary-series orders 3 and 7.

All radial profile values are simultaneous unknowns. Physical states are
packed as `[Re R; Re Q; Im R; Im Q]`, where `Q=dR/dr*`. The current 56 real
frequency/boundary parameters are represented by 56 additional constant
states. There are 112 real boundary constraints, 56 at each end. The two
fundamental outer amplitudes, including their phase convention, remain fixed.
This first adapter therefore has 112 unknowns per mesh point. Eliminating
the linear boundary coefficients or treating parameters as a global border
could reduce this overhead in a later implementation.

State scales are fixed from the full reference profiles. Each mode has its
own field and slope scale, with an absolute floor for tiny modes. Target
guesses with zero driven fields use the same scales as the driven guess.
Boundary-parameter scales use their endpoint influence; each frequency
component has displacement scale `0.001`. Every run saves these scales.
The collocation tolerance is a scaled differential-defect tolerance; it is
not the old shooting matching norm.

Available cases:

| Case | Secondary amplitude | Initial profile |
| --- | ---: | --- |
| `weak` | `1.0896058341192524e-7` | Previously matched joint weak solution |
| `scaled_weak` | `1.0896058341192524e-4` | Weak solution scaled by harmonic power |
| `target_zero` | `1.0896058341192524e-4` | Full primary/secondary references, zero driven fields |
| `target_driven` | `1.0896058341192524e-4` | Full primary/secondary and forced driven references |

The primary amplitude is `0.036390480353373555` throughout. Reference
profiles are Hermite interpolated using both saved `R` and `Q`. A mismatch
at a shooting join is retained in the guess; interpolation does not certify
that a reference is a solution.

## Library correction and derivative checks

The installed MIRK package's default two-point Jacobian band is too narrow
for a general coupled system. With `M` states and `La` left boundary
constraints, the two-node collocation blocks require lower bandwidth
`La+M-1` and upper bandwidth `2M-La-1` in the square case. The installed
implementation assumes `M+1` on each side.

`scripts/mirk_coupling_patch.jl` corrects this method at runtime, inside this
experiment only. It requires the pinned MIRK version. No package files or
production solver files are modified. A six-equation coupled linear example
exposes a maximum Jacobian-entry error of `0.010304921996152405` without the
correction; the corrected Jacobian agrees exactly with uncompressed AD.
Its solved profile agrees with the matrix-exponential solution to
`2.582936087236476e-10` on independent sample points.

The complete physical collocation Jacobian on a three-node test mesh also
agrees exactly with uncompressed AD. This is a derivative test, not a
physical discretization. Separate adapter tests check saved boundary states,
frequency/boundary derivatives, full field/frequency directional derivatives
against finite differences at three radial locations, and interpolation.
All 25 adapter checks passed.

Large differentiation batches caused excessive compilation. The working
configuration uses chunks of 28 and explicit greedy sparse coloring; startup
and warmed solve costs must be distinguished in any performance claim.

The installed package also overwrites old nodal values while interpolating
new mesh guesses. `scripts/test_remesh.jl` demonstrates the problem even for
the exact interpolant `y=t`: the error becomes `0.4`. The experiment-local
`mirk_remesh_patch.jl` preserves the old values, reducing this test error to
`1.11e-16`. This correction is used in subsequent refinement runs. The two
successful direct target solves needed no remeshing and therefore were not
affected by that issue. These tested corrections must be carried forward or
replaced by verified upstream fixes before production adoption.

## Running

The separate `Project.toml` and `Manifest.toml` pin Julia dependencies. The
wrapper uses a writable experiment-local depot followed by the existing
user depot, with Julia and BLAS each restricted to one thread.

```bash
bash runs/adiabatic/radial_bvp/scripts/run_julia.sh \
  -e 'using Pkg; Pkg.instantiate()'

bash runs/adiabatic/radial_bvp/scripts/run_julia.sh \
  runs/adiabatic/radial_bvp/scripts/test_adapter.jl

bash runs/adiabatic/radial_bvp/scripts/run_julia.sh --compiled-modules=existing \
  runs/adiabatic/radial_bvp/scripts/test_mirk_jacobian.jl

bash runs/adiabatic/radial_bvp/scripts/run_julia.sh --compiled-modules=existing \
  runs/adiabatic/radial_bvp/scripts/check_radial_jacobian.jl

sbatch runs/adiabatic/radial_bvp/scripts/submit_case.sh \
  weak 1e-4 2.0 runs/adiabatic/radial_bvp/results/unique_prefix
```

Arguments are case, tolerance, initial mesh spacing, output prefix, and an
optional earlier `.jls` result to restart the same case. Existing output
prefixes are refused. A batch job freezes the experiment scripts and project
files under `work/job-JOBID`; the result records source and input hashes.
Each solve is bounded by 1,200 intervals, 20 nonlinear iterations per mesh,
and a 2,400-second callback wall-time limit. Slurm supplies a final 55-minute
limit, including compilation and validation.

Run independent comparisons after solver results are available:

```bash
sbatch runs/adiabatic/radial_bvp/scripts/submit_analysis.sh \
  runs/adiabatic/radial_bvp/results/coarse.jls \
  runs/adiabatic/radial_bvp/results/fine.jls
```

The final input is the comparison reference. Analysis reports per-mode
profile changes, frequency changes, finite-domain E/J and endpoint fluxes
at radial quadrature spacings 1 and 0.5, and short RK4 reintegrations at
steps 0.05 and 0.025. No exterior charge tails are included.

## Acceptance and current results

A library success return alone is insufficient. Saved assessments include
boundary residuals, parameter constancy, off-mesh differential defects at
two independent points per interval, and uniformly sampled physical
profiles. Adoption additionally requires refinement of the frequencies,
profiles, and observables, preservation of the intended branch, and useful
convergence at the target amplitudes.

| Initial guess | Newton residual history (scaled 2-norm) | Boundary residual (infinity norm) | Independent off-mesh defect | Points | Solve seconds |
| --- | --- | ---: | ---: | ---: | ---: |
| Zero driven fields | `7.296 → 5.853 → 0.004655 → 2.02e-7` | `4.46e-8` | `3.18e-5` | 401 | 604.7 |
| Full driven profiles | `16.868 → 5.009 → 0.003637 → 1.67e-7` | `3.37e-8` | `3.18e-5` | 401 | 592.4 |

Both requested differential-defect tolerance `1e-4`. The solve timings include
approximately 121 seconds of Julia compilation, and exclude package loading
and input preparation. Complete Slurm durations were 12m56s and 12m04s.
Peak memory was approximately 19 GiB for this first implementation. These are
convergence results, not a controlled warmed speed benchmark.

The common frequencies are approximately
`ω1 = 0.2960068237643 - 2.9475427e-6im` and
`ω2 = 0.2983684033854 + 8.4141072e-7im`.
Machine-readable results are in [`results/target-convergence.toml`](results/target-convergence.toml).
Core checkpoints are `results/target_{zero,driven}_mirk4_v1.jls`; logs are
`logs/bvp-27390854.out` and `logs/bvp-27390855.out`.

Independent short RK4 reintegrations of the zero-driven result disagree by
at most `2.88e-5` in the fixed component scales. Halving RK4 spacing changes
that check by less than `8.90e-9`. This supports the off-mesh defect check.
The driven result passes the same check. Full numerical outputs are in
[`results/target-validation.json`](results/target-validation.json).

Refinement to requested tolerance `1e-6` converged on **881 mesh points**:
the boundary residual is `1.62e-13`, the independently sampled differential
defect is `3.04e-7`, and parameter drift is zero. The largest frequency
change from the 401-point solution is `1.2461e-10`. Maximum per-mode profile
changes are `0.0257%` for fields and `0.103%` for slopes; these are global
profile changes, distinct from the local defect tolerance. See
[`results/target-refinement.toml`](results/target-refinement.toml).
The refinement used two corrections on its initial 801-point mesh and one
after adaptive redistribution, with the mesh-transfer correction active.
It took 635.0 solve seconds, including 112.2 seconds of compilation, and
11m33s including startup and assessment. Peak memory was approximately
7.4 GiB with line search disabled for this nearby restart.

The refined profile also passes short independent RK4 reintegrations, with
maximum scaled discrepancy `4.01e-7` and integration-step sensitivity
`8.89e-9`. Refinement changes finite-domain energy by `1.01e-7` relative,
angular momentum by `2.60e-8`, and the largest boundary-flux change is
`4.75e-5` relative (`0.00475%`). The coarse charge quadrature was already
stable to about `1.02e-11` between radial spacings 1 and 0.5. See
[`results/target-refined-validation.json`](results/target-refined-validation.json)
and [`results/target-observable-refinement.json`](results/target-observable-refinement.json).

## Speed investigation

`scripts/speed_benchmark.jl` compares the original `BackTracking` configuration
with `BackTracking(autodiff=AutoForwardDiff(chunksize=1))`. The latter computes
the merit-function directional derivative directly. In the pinned library,
the default instead allocates a dense matrix and calls the supplied full
Jacobian before multiplying by the direction. With 401 nodes and 112 states,
that matrix has 44,912 rows and columns and occupies 15.03 GiB by itself.
The Newton correction still uses the corrected banded Jacobian in both runs.

The benchmark repeats the direct driven-target solve twice in each process,
always resetting to the same initial guess. The first pass includes compilation;
the second measures warmed work and collects a sampling profile. Both use
401 initial nodes, differential tolerance `1e-4`, nonlinear tolerance `1e-6`,
the same backtracking constants, and the same physical source snapshot.
Every pass checks the physical boundaries and independent off-mesh defects.
Callback counters and timers separate ordinary source evaluations from each
derivative batch width. Timed solve work excludes loading, initial-guess
preparation, and the separate assessment.

```bash
sbatch runs/adiabatic/radial_bvp/scripts/submit_speed.sh baseline \
  runs/adiabatic/radial_bvp/results/speed_baseline_unique
sbatch runs/adiabatic/radial_bvp/scripts/submit_speed.sh directional \
  runs/adiabatic/radial_bvp/results/speed_directional_unique
sbatch --export=ALL,RADIAL_BVP_AD_CHUNKSIZE=56 \
  runs/adiabatic/radial_bvp/scripts/submit_speed.sh directional \
  runs/adiabatic/radial_bvp/results/speed_directional56_unique
```

The runs should use the same compute-node model and one Julia/BLAS thread.
Their `-cold.toml` and `-warm.toml` files report timings and residual checks;
`-warm-profile-{flat,tree}.txt` records the sampling profile. The batch wrapper
freezes scripts and dependencies under its job ID.

`scripts/test_directional_linesearch.jl` passed all seven checks. On the physical
three-node collocation system, direct Jacobian-vector products agree with the
full corrected Jacobian to `8.30e-16` relative. A scalar Newton step that
overshoots the root receives the same backtracked step, `0.1`, through both
derivative paths. This checks the changed derivative calculation and step rule;
the full target benchmarks check convergence and physical accuracy.

**Recommendation:** proceed toward an optional collocation backend. The
target-convergence advantage is established for this case. Production work
should package the two tested library corrections and reduce avoidable
parameter/line-search overhead; a warmed speed advantage has not yet been
measured. The production shooting implementation is unchanged.

The comparison holds equations, physical boundaries, retained modes, and
amplitudes fixed. Collocation has different unknowns and residual scaling,
so this does not isolate which numerical improvement is responsible. It
establishes a root of the existing finite-domain problem; separate boundary,
angular-grid, and mode-truncation systematics remain to be assessed.

Cancelled setup runs: 27390715 used the uncorrected Jacobian band; 27390794
used an excessively large differentiation batch; 27390815 omitted explicit
sparse coloring. Weak-case run 27390835 reached its coarse-mesh root in one
correction, but its subsequent mesh transfer was affected by the overwrite
issue and was stopped when the investigation focused on the target. None
of those cancelled runs is counted as a physical convergence failure.
