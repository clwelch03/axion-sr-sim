# Multiphase derivative verification — October 27, 2026

The user authorized completing the derivative verification autonomously after
building the input, weak-secondary guess, weighted residual, and ordinary-AD
column helper together. The target is a documented comparison at the original
weak-secondary guess and a freshly matched seed, plus finite differences and
independent source tests. This work does not increase the secondary amplitude.

## Fixed calculation and reproducibility

- Input: `runs/adiabatic/multiphase_two_cloud/results/direct_fig7_20261026_v5-input.jls`.
- SHA-256: `f0152fda5c7ecbaffb06be7ad334cacfe2a0c1fbef94c8c590c4372e237ee6fb`.
- Primary outer amplitude: `0.036390480353373555`; secondary:
  `1.0896058341192524e-7` (0.1% of the target).
- Same 14 modes, 16 × 31 × 31 quadrature, spin 0.99, alpha 0.3, radial step
  0.1, endpoints -100/300, match 5, boundary orders 3/7 as yesterday.
- The original coefficient scale stays fixed; every secondary-dependent
  residual row is divided by another 0.001.
- Frozen source: `runs/adiabatic/multiphase_two_cloud/work/ad_verification_20261027_v1`.
  All 199 manifest files verified. The input copy is `inputs/initial.jls`.
- Numerical library behavior is unchanged. The user's rename to
  `two_phase_real_matching_residual` is incorporated.

The standalone experiment is
[`autodiff_vs_optimized_test.jl` at its historical revision](https://github.com/clwelch03/axion-sr-sim/blob/e29b40c22c6e0ba518e0dfe2a21bae086e33566d/src/adiabatic/multiphase/autodiff_vs_optimized_test.jl).
That standalone diagnostic was subsequently retired in `5698692`; its frozen
source is also preserved with the campaign inputs. The stale test and runner
reference were removed during the November 3 checkout validation.
It is not loaded by the multiphase library loader. It has `initial`, `seed`,
and `matched` execution stages. The existing Newton solver and matching/source
functions supply the numerical work. The seed stage allows 20 corrections at
the unchanged `1e-5` tolerance, saves accepted iterates, and reintegrates a
converged state before saving both radial branches.

## Checks and acceptance criteria

The optimized 56-column matrix uses the production full-width ForwardDiff
route. The ordinary reference computes one column at a time with the diagnostic
tag that invokes generic source arithmetic, bypassing `local_source_jacobian`.
For each column, report the norm of the difference divided by the larger
column norm, requiring at most `1e-8`. Two exactly zero columns have zero error;
nonfinite comparisons fail. Every column is checked even after a discrepancy.
Matrices, per-column errors, elapsed times, and completion counts are saved
after each column. Timing this reference route does not constitute a controlled
benchmark of the original pre-optimization implementation.

Two directions, based on normalized cosine and sine patterns with the same
coordinate scaling as yesterday, are checked with centered differences at
steps `1e-9` through `1e-15`. The original `1e-5` relative criterion is retained.
The diagnostic also reports comparison against the actual representable
displacement and the error introduced by rounding the requested displacement;
these additional measurements do not silently replace the acceptance criterion.
Both directions and all step sizes are saved even when checks fail.

Independent tests compare source values and derivatives against direct
quadrature, including zero, weak, and nonlinear fields and three potentials.
The compact regression additionally exercises all 16 Jacobian columns through
a short four-mode nonlinear matching calculation. Its short domain tests the
implemented derivative pipeline; it is not a physical boundary/convergence test.

## Completed verification

Recovered after the n0003 SSH outage from shared artifacts and Slurm accounting.
All three jobs completed with exit `0:0` on **n0060.lr7**. The queue was empty
at the recovery check. The initial and matched serialized matrices and
finite-difference records were loaded and checked for all 56 columns, both
directions, finite errors, unchanged thresholds, and identical comparison
states. All 199 frozen-source manifest files still verify. The saved seed loads,
its hash matches the matched-stage setup, and its packed state matches the
verified matrix input; both saved radial branches are finite.

On resuming from “Recover chats after SSH timeout,” a read-only audit again
confirmed the three successful Slurm exits and all 199 snapshot hashes. All 62
live Julia source files, dependency manifests, and numerical tests still match
the verified snapshot; the changed chat-recovery test is unrelated. Recomputing
the per-column errors from both saved matrices reproduces every AD pass, and
both saved finite-difference sweeps retain their original passing thresholds.
The seed state equals the matched comparison state exactly, its provenance hash
matches, and the norm recomputed from its saved radial mismatch is
`5.306800040310494e-6`. Both radial branches are finite. Repacking the physical
profile differs from the saved solver state only by coefficient-scaling
roundoff (maximum absolute difference `2.7755575615628914e-17`). This audit ran
no new integrations and submitted no jobs; no unfinished derivative check was
found in the agreed scope.

| Check | Status |
|---|---|
| Comparison, failure detection, finite-difference helper, atomic save, accepted-state save at iteration cap | 14/14 passed |
| Full matrices through short nonlinear matching | 3/3 passed; all 16 columns below `1e-10` |
| Independent direct-quadrature source tests | 38/38 passed |
| Full domain, original weak guess | All 56 AD columns pass; worst relative error `8.87909e-9`. Both finite-difference directions pass. |
| Weak-seed solve | Job **27334006** COMPLETED, 22m45s; correction 13 and fresh reintegration both give `5.306800040310494e-6` |
| Full domain, matched seed | Job **27334124** COMPLETED, 33m19s; all 56 AD columns and both finite-difference directions pass |

The 55 local checks ran with one BLAS/Julia thread:

```bash
julia --startup-file=no --project=. -e 'using Test; include("tests/adiabatic/multiphase/test_autodiff_comparison.jl"); include("tests/adiabatic/multiphase/test_source_basis.jl")'
```

The saved output is
[`ad-verification-local-20261027.log`](../runs/adiabatic/multiphase_two_cloud/logs/ad-verification-local-20261027.log).
The new regression is also included in the multiphase test runner. The unrelated
previous evolution-test field mismatch has not been addressed or reclassified.

The physical run prefixes are `results/ad_initial_20261027_v1` and
`results/ad_seed_20261027_v1` under the multiphase campaign. Logs:
[`27333738`](../runs/adiabatic/multiphase_two_cloud/logs/ad-verify-27333738.out),
[`27334006`](../runs/adiabatic/multiphase_two_cloud/logs/ad-verify-27334006.out).
The matched-state prefix is `results/ad_matched_20261027_v1`; job 27334124
starts only if the seed job exits successfully after reintegrating and saving
the seed. An invalid dependency cancels it. All allocations use one CPU and
8 GB; comparisons have a two-hour cap and the seed a one-hour cap. The agreed
derivative verification is complete for this configuration at the original
guess and matched weak seed. This does not establish target-amplitude
convergence or spatial accuracy.

The matched seed is saved as
[`ad_seed_20261027_v1-seed.jls`](../runs/adiabatic/multiphase_two_cloud/results/ad_seed_20261027_v1-seed.jls),
including both radial branches. The last accepted-state artifact and setup
record are also retained. The seed job completed at 15:26:26 PDT; the matched
comparison started at 15:26:29 PDT and independently reproduced the same
weighted matching norm before building its Jacobian. This is a weak-secondary
seed at 0.1% of the target A2, not the target-amplitude solution for the plot.

## Initial-state full-domain comparison

All 56 columns satisfy the predeclared `1e-8` relative threshold. The
matrix-wide Frobenius relative difference is `6.88584e-15`, but the individual
columns expose a larger discrepancy in two small columns:

| Coordinate | Column norm (optimized) | Difference norm | Relative difference |
|---|---:|---:|---:|
| Real free outer coefficient, l=5, q=(5,0) | 0.452281 | `4.01584e-9` | `8.87909e-9` |
| Imaginary free outer coefficient, l=5, q=(5,0) | 0.452073 | `3.11650e-9` | `6.89380e-9` |

These are the primary's fifth-harmonic coefficients, in the scaled solver
coordinates. The first is close to the threshold; it has not been omitted or
rounded into a claim of machine-precision agreement for every column. All
other columns have relative difference at most `2.52820e-11`. The complete
named-coordinate table is
[`ad_initial_20261027_v1-columns.csv`](../runs/adiabatic/multiphase_two_cloud/results/ad_initial_20261027_v1-columns.csv).
Optimized construction took 94.97 s; the 56 ordinary one-direction evaluations
took 1628.71 s in total, including compilation effects.

The first finite-difference direction has relative errors `4.56835e-4`,
`4.62714e-6`, `4.63085e-8`, `7.44209e-10`, and `1.94898e-8` at steps
`1e-11` through `1e-15`. Direction 2 likewise reaches `5.29663e-6` at `1e-12`,
`5.30071e-8` at `1e-13`, and `2.24105e-9` at `1e-14`. This establishes a useful
range of agreement with approximately second-order improvement before
floating-point effects limit the smallest probes. The old steps `1e-6` through
`1e-8` had not reached that range. No acceptance threshold was relaxed.

## Matched-state full-domain comparison

Job **27334124** completed at **15:59:48 PDT** on October 27. All 56 columns
satisfy `1e-8`; the worst relative error is `9.587805056868246e-9` in column 21,
the real free outer coefficient for l=5, q=(5,0). Column 49, its imaginary
counterpart, has relative error `8.690135834245296e-9`. These remain close to
the threshold and are not described as machine-precision column agreement.
The matrix-wide Frobenius relative difference is `2.01415e-16`. Optimized
construction took 92.08 s; ordinary column evaluations totaled 1610.56 s,
including compilation effects and with different batching.

Both centered-difference directions pass the unchanged `1e-5` criterion.
At `h=1e-12`, their relative errors are `6.42847e-6` and `8.79351e-6`;
at `1e-13`, `6.42379e-8` and `8.80294e-8`; their minima occur at `1e-14`,
with `2.72366e-9` and `1.36802e-9`. Both increase again at `1e-15`.

The saved evidence is `results/ad_matched_20261027_v1-{setup,ad,fd}.jls` and
[`ad_matched_20261027_v1-columns.csv`](../runs/adiabatic/multiphase_two_cloud/results/ad_matched_20261027_v1-columns.csv).
The named-column CSV was exported from the completed saved matrix during
recovery, using the same coordinate ordering as the initial-state table.

**Decision:** local-Jacobian AD is no longer an open derivative blocker for
this configuration and these tested states. Resume from the saved weak seed;
review the next amplitude step with the user before launching continuation.
The seed remains at 0.1% of target A2. Newton convergence at higher amplitudes,
spatial refinement, and the architecture review remain separate work.
