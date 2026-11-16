# Benchmark and regression harness

Measures the cost of one spindown right-hand-side evaluation and checks that the numbers it
produces have not changed. It uses the library as it is (nothing under `src/` is modified).

## Files

| File | Purpose |
| --- | --- |
| `common.jl` | State definitions, loading, problem setup, thin wrappers around library internals (module `BenchmarkHarness`) |
| `quantities.jl` | The compared observables and the comparison measure (included by `common.jl`) |
| `regression.jl` | `--write-reference` or compare against `reference.toml`; PASS/FAIL and exit code |
| `timing.jl` | Min-of-N wall time and allocations of each piece; table on stdout, CSV in `logs/` |
| `noise_floor.jl` | Measures the noise floors behind `tolerances.toml` |
| `build_states.jl` | Rebuilds the stored states from scratch (see "Fixed inputs") |
| `submit.sh` | Slurm wrapper: regression, then timing |
| `states/*.bin` | The four reference states: raw `Float64` unknown vectors, tracked in git |
| `reference.toml` | Reference values, 17 significant digits |
| `tolerances.toml` | Relative tolerance and measured floor per quantity class |

## How to run

All commands from the repository (or worktree) root, `julia --project=.`. The first run
precompiles for 1 to 2 minutes.

```
julia --project=. benchmarks/regression.jl                  # compare: about 4 min on one core
julia --project=. benchmarks/regression.jl --skip-evolution # about 2 min
julia --project=. benchmarks/timing.jl --repeats 3          # about 10 min
benchmarks/submit.sh <cpus> <julia-threads> [blas-threads=1] [walltime]   # both, in one lr7 job
```

`submit.sh` writes `logs/bench-<jobid>.out` and the timing CSV `logs/timing-<commit>-<host>-j<N>b<M>.csv`
(columns: commit, `source_modified` flag for `src/`, hostname, CPU model, Julia and BLAS thread counts,
state, piece, minimum seconds, repeats, allocated bytes, note). Thread counts come from the
environment (`julia -t`, `OPENBLAS_NUM_THREADS`, set by `submit.sh` from its arguments). Options of
`regression.jl` and `timing.jl`: `--states linear,plateau,...`, `--skip-evolution`, `--reference PATH`
(regression), `--repeats N` or `BENCH_REPEATS` and `--csv PATH` (timing).

Do not use `--write-reference` on an optimised branch. The references belong to a baseline commit: they were
first written from b4c4100 and rewritten on 2026-11-14 from 6bec64a (the original solver with the
conjugate-frequency rule of 41fc2e7, which changes the equations) by checking out that commit with this
directory from `speedup/base` and running `regression.jl --write-reference`. The stored states are unchanged;
they were converged under the earlier rule, so they are fixed evaluation points, not solutions of the present
equations. When a later commit changes the equations again, rewrite the reference the same way, from the
unoptimised code if it still exists and otherwise from the last commit before the change plus a separate
argument for the change.
To prove a speed-up, run `regression.jl` (must PASS) and `timing.jl` on the same node type as the
baseline table below and compare the CSV rows.

## The states

Production problem: `KerrBackground(0.9)`, `ModeSet([BoundMode(2,1,1), BoundMode(3,2,2)])`,
`RadialMesh(background, α, mode_set)`, `CosinePotential()`, `decay_constant_ratio = 0.0466873`.
"Reference problem" below means the problem at the initial background (it fixes units and mesh).

| State | α | b, j | Amplitudes | Mesh nodes | Unknowns |
| --- | --- | --- | --- | --- | --- |
| `linear` | 0.3 | 1, 0.9 | (9.856e-7, 3.057e-9), the seed that `scripts/evolution_spindown.jl` prints | 690 | 22112 |
| `plateau` | 0.3 | 0.98147538879, 0.83749443414 | (2.8077e-3, 3.396e-9), t = 1.5e9 of log 61110567 | 690 | 22112 |
| `nonlinear` | 0.3 | 0.99546020550, 0.86714826011 | (1.4719e-4, 4.100e-6), t = 1e11 of the same log | 690 | 22112 |
| `large_mesh` | 0.187 | 1, 0.9 | (1e-6, 1e-6), linear | 1355 | 43392 |

### Fixed inputs

Each state is its converged unknown vector in `states/<name>.bin` (tracked). Amplitudes, `b`, `j`, α
and the mesh are in `StateSpec` in `common.jl`; the unknown scale is `max(amplitudes)`. `load_state`
rebuilds a `RadialSolution` from the vector with `solution_from_unknowns` (one function in
`common.jl`: fresh Jacobian, `NewtonScaling`, `ScaledFactorization`), so no serialized library object
is stored. If an optimisation changes `RadialSolution` or the factorisation type, that function is
the one to adapt.

`build_states.jl` rebuilds the vectors from scratch (a cold solve at the target amplitudes for
`linear` and `large_mesh`; a cold solve at 1e-3 times the amplitudes on the target background and 30
geometric amplitude steps for `plateau` and `nonlinear`, as the restart branch of
`scripts/evolution_spindown.jl` does; then a polish at tolerance 1e-12, or 1e-9 on the large mesh).
A rebuilt state agrees with the stored one only to the Newton tolerance, and a state built with a
different probe tolerance is a different state, so rebuilding requires rewriting `reference.toml`.
The rebuild takes 30 s to 2 min per state.

## What each check establishes

Quantities are grouped in two tiers.

**Evaluation tier** (pure evaluations at the stored unknown vector or at a fixed smooth displacement of it
by relative 1e-3, `fixed_vector` in `common.jl`):

* residual: 2-norm and 12 fixed-index entries at the displaced vector, and the max-norm at the
  converged vector (roundoff-sized, compared on the scale of the displaced residual);
* Jacobian at the displaced vector: stored-entry count, Frobenius norm, norm and 12 entries of `J v` for
  a fixed vector `v`;
* `cloud_charges`, `phase_actions`, `boundary_losses`, `action_losses`;
* `amplitude_response` (column norms and 12 entries per column), `charge_jacobian`, `action_jacobian`;
* every field of `background_charge_derivatives` and of `balance_terms`, `amplitude_rates` from the
  terms (spindown balance) and from the solution (fixed background).

**Re-solve tier**: `solve_radial` from the state after a fixed perturbation (second amplitude times 1.01,
`b` times `1 − 1e-6`, `j` times `1 − 1.5e-6`): mode frequencies, 12 unknown entries, cloud charges;
and `evolve_with_spindown` from `plateau` over t = 1.5e9 to 1.9e9: the saved times, both amplitudes,
`b` and `j`. Iteration counts, accepted-step count and right-hand-side call count (4 steps, 31 calls at
the baseline) are printed next to their reference values but never asserted.

What this does not establish: accuracy or convergence of the discretisation (the references are the
baseline's own numbers, not truth), results on other states, other potentials, `balance=:energy`,
or step rejections (none occur in the short segment). Only 12 sampled entries of long vectors are
kept, so a change confined to other entries of the unknowns, Jacobian or response vectors is seen only
through the norms and the derived physical quantities. The imaginary part of the `|211>` frequency on
the plateau is noise (about 1e-13) and is inside a 2-norm dominated by the real part.

Comparison measure for each quantity (a group of numbers):
`‖new − reference‖₂ / max(‖reference‖₂, anchor)`. `anchor` is 0 except for differences of nearly equal
terms. For those, the anchor, stored with the reference, is the size of the terms:
losses (outer minus inner flux; anchor = the larger flux, `edge_flux_scales`), metric work (sum of
absolute values of its products), amplitude rates (`|[∂N/∂A]⁻¹|` applied to the sizes of the loss
and drift terms; on the plateau the rates are about 1e-15 against inputs of 1e-8 to 1e-13), the
converged residual (anchor = displaced residual size), and sampled entries of long vectors
(anchor = 1e-3 of the vector's max-norm, so a near-zero entry is not judged relatively).

## Tolerances and measured noise floors

`noise_floor.jl --evolution` compares each quantity with itself under variations that must not change
the answer. The floor of a class is the largest relative difference over the four states.
Variations: `blas_threads` (BLAS 4 threads instead of 1), `roundoff_input` (every unknown multiplied
by `1 + 4e-16 × pattern`, about 2 ulp; this stands in for reordered floating-point sums in optimised
code and measures how much the computation amplifies input roundoff), and for the re-solve tier
`solver_path` (a ten times tighter Newton tolerance, except on the large mesh whose residual floor of
5e-10 makes that unreachable, and a different warm-start solution). The `blas_threads` variation is
exactly zero or below 5e-16 everywhere: the dense linear algebra here is too small for threading to
reorder anything. A different CPU model (Xeon Gold 6330 on lr7 for the reference, Gold 6230 on the
login node for the comparison) also gave differences of at most 5e-16 on every quantity. Tolerances are
set a factor of about 3 to 10 above the floor, with a lower limit of 1e-12.

| Class | Tier | Floor | Tolerance | Quantities, and where the floor comes from |
| --- | --- | --- | --- | --- |
| `exact` | evaluation | 0 | 0 | Jacobian stored-entry count |
| `charge` | evaluation | 4.2e-16 | 1e-12 | `cloud_charges`, `phase_actions` |
| `residual` | evaluation | 2.3e-11 | 3e-10 | residual norm and entries; large mesh: 2 ulp of input becomes 2e-11 (cause not investigated) |
| `residual_converged` | evaluation | 4.2e-13 | 1e-11 | max-norm at the converged vector, anchored |
| `jacobian` | evaluation | 2.4e-11 | 3e-10 | Frobenius norm, `J v` (large mesh) |
| `response` | evaluation | 2.3e-11 | 3e-10 | `amplitude_response`, both observable Jacobians (large mesh) |
| `loss` | evaluation | 6.3e-11 | 1e-9 | boundary and action losses, metric work (anchored; large mesh) |
| `rate` | evaluation | 2.2e-12 | 3e-11 | amplitude rates (anchored; plateau) |
| `rate_background` | evaluation | 1.4e-15 | 1e-12 | `b'`, `j'` |
| `derivative` | evaluation | 3.1e-7 | 3e-6 | `BackgroundDerivatives` fields: the explicit part is a central difference with relative step 1e-6, which turns 1e-16 input roundoff into 1e-10 and more (3e-7 for `∂N/∂b` on the linear state) |
| `solve_frequency` | re-solve | 1.1e-14 | 1e-12 | perturbed-solve mode frequencies |
| `solve_unknown` | re-solve | 2.9e-14 | 1e-12 | 12 unknowns of the perturbed solution |
| `solve_charge` | re-solve | 1.2e-12 | 1e-11 | `cloud_charges` of the perturbed solution |
| `trajectory_amplitude` | re-solve | 1.1e-8 | 1e-7 | evolution amplitudes (`evolution.amplitude_2` on the plateau); the integrator's own `reltol` is 1e-7 |
| `trajectory_background` | re-solve | 2.5e-13 | 1e-11 | `b(t)`, `j(t)` |
| `trajectory_time` | re-solve | 2.7e-8 | 3e-7 | accepted step times |
| `reported` | none | n/a | not asserted | iteration, step and right-hand-side counts (differ by up to 20% between solver paths) |

Points to know:

* The `derivative` tolerance (3e-6) is loose because of the finite-difference step in
  `background_charge_derivatives`, not because of the optimisation under test. A change that
  alters that step or its arithmetic will be caught only if it changes the result by more than this.
* The trajectory tolerances are of the order of the integrator's tolerance. A change that flips an
  accept or reject decision of a step changes `times` at once by O(1) relative; the step count is
  printed so such a flip can be seen.
* On the large mesh the Newton residual cannot go below about 5e-10 (scaled max-norm), so the default
  `tol = 1e-10` never converges there (30 iterations, about 90 s per solve). The large-mesh state is built at 1e-9
  and its perturbed solves use `tol = 1e-8` (`run_tolerance` in `common.jl`); its
  timing rows are for that tolerance.
* The reference was written on lr7 (Xeon Gold 6330). Compare runs on the login node (Xeon Gold 6230)
  and on lr7 both agreed to 5e-16 or better.

## Baseline timings

Slurm job 27815524 (the `large_mesh` right-hand-side row is from job 27815832, same node type, after a harness fix), lr7 node n0004 (Intel(R) Xeon(R) Gold 6330 CPU @ 2.00GHz), 1 Julia thread, 1 BLAS thread, commit b4c4100870, minimum of 3 calls after a warm-up (evolution segment: one run). Seconds, with allocated MB in parentheses. `large_mesh` solves use tol 1e-8.

| Piece | linear | plateau | nonlinear | large_mesh |
| --- | --- | --- | --- | --- |
| residual | 0.060 (220) | 0.061 (220) | 0.065 (220) | 0.114 (408) |
| jacobian assembly | 0.528 (1416) | 0.523 (1435) | 0.541 (1419) | 0.621 (1683) |
| factorisation | 0.124 (190) | 0.122 (190) | 0.137 (190) | 0.259 (374) |
| back-substitution | 0.003 (1) | 0.003 (1) | 0.003 (1) | 0.006 (2) |
| solve_radial, perturbed | 1.937 (4380); 5 iterations | 1.806 (4163); 4 iterations | 2.141 (4380); 5 iterations | 4.500 (8830); 4 iterations |
| solve_radial, same point | 1.852 (4800); 1 iterations | 1.808 (4801); 1 iterations | 2.013 (4800); 1 iterations | 3.676 (8443); 1 iterations |
| cloud_charges | 0.031 (90) | 0.032 (90) | 0.036 (90) | 0.065 (176) |
| phase_actions | 0.002 (8) | 0.002 (8) | 0.003 (8) | 0.005 (15) |
| boundary_losses | 0.000 (0) | 0.000 (0) | 0.000 (0) | 0.000 (0) |
| action_losses | 0.000 (0) | 0.000 (0) | 0.000 (0) | 0.000 (0) |
| amplitude_response | 0.032 (41) | 0.031 (41) | 0.033 (41) | 0.042 (44) |
| charge_jacobian | 0.211 (343) | 0.202 (343) | 0.207 (343) | 0.421 (673) |
| action_jacobian | 0.005 (25) | 0.005 (25) | 0.005 (25) | 0.010 (50) |
| background_charge_derivatives | 0.714 (1324) | 0.692 (1324) | 0.742 (1324) | 1.484 (2546) |
| balance_terms (total) | 0.987 (1733) | 0.984 (1733) | 1.105 (1733) | 2.066 (3313) |
| spindown rhs call | 3.022 (6631) | 2.992 (6632) | 3.419 (6631) | 6.103 (11947) |
| evolution 1.5e9 -> 1.9e9 |  | 94.288 (182065); 4 accepted steps, 31 rhs calls |  |  |

## Library internals the harness depends on

Functions reached through `AxionCloud.` that are not exported (`INTERNAL_INTERFACE` in `common.jl`):
`residual`, `jacobian`, `UnknownLayout`, `NewtonScaling`, `ScaledFactorization`, `RadialParameters`,
`parameter_range`, `frequency_range`, `n_unknowns`, `_complexify`, `_spindown_rhs!`, and the
`EvolutionCache` fields `accepted_solution` and `diagnostics`. Wrappers for them are in `common.jl`.
`ode_solution.stats.nf` is used as the right-hand-side call count.
