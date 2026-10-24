# Multiphase work — October 23, 2026

Today we prepared and launched a direct joint nonlinear profile calculation
for a Fig. 7-style comparison, diagnosed its unexpectedly expensive Jacobian,
and began implementing a substantially faster source calculation together.
The target joint profile has not converged yet, and the final overlay has not
been produced. Job status below was checked at **17:25 PDT on October 23**.

See the [campaign record](../runs/adiabatic/multiphase_two_cloud/README.md#direct-fig-7-profile-attempt-october-23)
for run recipes and the [learning map](multiphase-learning-map.md) for the
function walkthrough. This is a dated record; its job states and implementation
status will become historical as work continues.

## The deliverable and the controlled comparison

The requested deliverable is an overlay of the primary and secondary radial
profiles, following Omiya 2023 Fig. 7, with an additional pair of curves from
our joint nonlinear calculation. The joint source is evaluated from the total
retained field, so it includes interactions between the populated clouds.

We chose a **direct solve at the target amplitudes**. An incremental increase
of the secondary amplitude would be a continuation strategy to help Newton
converge; it is not required to define this comparison. There is no amplitude
scan or time evolution in today's direct-profile driver.

The control variables are the two prescribed outer coefficients, A1 and A2.
We use saved primary/secondary solutions to choose those normalizations and
construct a consistent initial guess. Their interior profiles are then allowed
to change when solving the joint equations. Using the reference profiles as a
seed does not constrain the final joint profiles to retain their shapes.

The saved reference nearest the paper's energy 660 is point 72 of 83:

| Reference quantity | Value |
|---|---:|
| Alpha, mu M | 0.3 |
| Spin, a/M | 0.99 |
| Reference E1 = E2, in Fa² M | 650.5419114886128 |
| A1 | 0.036390480353373555 |
| Unit-amplitude secondary energy | 5.479443353694717e10 |
| A2 | 0.00010896058341192524 |

The secondary is normalized using `A2 = sqrt(reference_energy / unit_energy)`.
The secondary table's `amplitude` field identifies its primary background; it
is not the physical secondary normalization.

**Equal energies describe the reference.** The joint solve holds A1 and A2
fixed and does not constrain the resulting joint energies to equal each other
or the reference energy. The eventual plot needs to state that normalization.

## Driver, settings, and saved outputs

We built [direct_profile.jl](../runs/adiabatic/multiphase_two_cloud/scripts/direct_profile.jl)
as sequential commented sections: load the reference data, select the target,
normalize the secondary, choose components, assemble the boundary guess,
construct the grids/problem, check packing and ordering, and run/save a solve.
The user worked through the construction before authorizing the remaining
driver sections to be added together.

The retained phase labels and angular harmonics are:

| Phase label q | Retained l values |
|---|---|
| (1, 0) | 1, 3, 5 |
| (3, 0) | 3, 5 |
| (5, 0) | 5 |
| (0, 1) | 2, 4, 6 |
| (2, -1) | 0, 2, 4 |
| (-1, 2) | 3, 5 |

This gives 14 complex radial functions and their 14 complex slopes. Newton
has **56 real unknowns**: two complex frequencies, 14 complex horizon
coefficients, and 12 complex free outer coefficients. The prescribed outer
coefficients occupy flattened positions 1 and 7. The mixed-label blocks
q=(2,-1) and q=(-1,2) start with zero coefficients and are solved for jointly.

| Numerical setting | Value |
|---|---|
| Polar × phase-1 × phase-2 nodes | 16 × 31 × 31 |
| Tortoise endpoints | -100, 300 |
| Matching radius in tortoise coordinates | 5 |
| Maximum radial RK4 step | 0.1 |
| Horizon/outer series orders | 3 / 7 |
| Coefficient/residual scale | min(A1, A2) = A2 |
| Scaled matching tolerance | 1e-5 |
| Maximum Newton corrections | 12 |
| Slurm allocation | 1 CPU, 8 GB, 4-hour cap |

The four-hour allocation is a limit, not a completion estimate. The matching
tolerance is also distinct from radial, angular, domain, or mode-truncation
accuracy. Those physical convergence checks remain outstanding.

The driver supports setup-only operation with no arguments, or
`--check INPUT_DIRECTORY`. Supplying `INPUT_DIRECTORY OUTPUT_PREFIX` runs the
solve. It checks the selected reference, mode ordering, fixed coefficients,
packing round trip, and initial admissibility before integration.

The [preparation script](../runs/adiabatic/multiphase_two_cloud/scripts/prepare_direct_profile.py)
freezes the source and three reference inputs, records provenance and hashes,
and does not submit a job. The [submission wrapper](../runs/adiabatic/multiphase_two_cloud/scripts/submit_direct_profile.sh)
validates its arguments and verifies the snapshot before launching Julia.
Snapshots v2 and v3 each verified 193 files. Historical snapshots are preserved.

Each fresh result prefix reserves a lock and owns:

- `-input.jls`: problem, initial guess, amplitudes, settings, and provenance.
- `-solution.jls`: successful solution and both radial integration histories.
- `-profiles.csv`: complex and absolute primary R11 and secondary R22 on both branches.
- `-failure.txt`: an exception report if the driver catches a failure.

Existing prefixes are not overwritten. The direct jobs currently have their
input checkpoints but no solution or profile CSV.

## Making the slow solve observable

The first running version only printed accepted Newton iterations. We added
`logged_profile_newton` to the driver, which times and flushes messages for:

- The initial residual evaluation.
- Each ForwardDiff pass used to build a Jacobian.
- Each trial step during backtracking.
- Accepted corrections, input saving, and final profile integration.

The wrapper preserves the residual and derivative-carrying values passed to
Newton; only the printed norm extracts ordinary values. Seventeen focused
checks covered unchanged results/evaluation sequences, damping, 56-variable
systems, and error propagation. The extra logging was frozen in v3 and
submitted as a separate job; v2 remained running.

At these settings, one residual requires about 4,000 RK4 steps, or 16,000
source evaluations. ForwardDiff's default width is 12 derivative directions,
so a 56-column Jacobian requires **five complete residual passes**. Repeated
residual norms during those passes describe the same Newton point; they do
not indicate accepted corrections.

The earlier dilute pilot used only four radial functions, an 8 × 9 × 9 grid,
and about 276 radial steps. The projection-loop work count for the current
setup is roughly 1,204 times larger. A warmed source benchmark estimated
about 532 seconds of source work per ordinary residual, close to the observed
537 seconds. Earlier live accounting also showed CPU time close to wall time
and less than 1 GB peak resident memory. The measurements point to expensive
computation, rather than a job waiting for I/O.

## The new optimization: evidence, installed code, and next steps

The first experiments found a type-specialization problem in field
reconstruction. Moving its inner loop into a specialized helper gave roughly
a threefold source speedup. Caching the angular/phase basis made ordinary
source evaluations much faster, but generic differentiation through those
matrix operations still limited the derivative speedup to about threefold.
The user explicitly requested an order-of-magnitude improvement.

The stronger prototype uses the source's **local derivative matrix** and
propagates all 56 incoming solver directions together. It keeps the same
retained modes, quadrature, and cosine potential.

For this explanation, define a real vector c containing the real parts of all
14 radial coefficients followed by their imaginary parts. Define B as the
cached real basis matrix, with one row per angular/phase sample and cosine
columns followed by sine columns. The reconstructed field samples are

\[
\Phi = 2 B c.
\]

Let h be the implemented `potential_remainder`, and let G be a diagonal
matrix containing the quadrature weights multiplied by
`radius^2 + spin^2*cos(theta)^2`. The packed real source s and its coefficient
Jacobian are

\[
s = \alpha^2 B^T G h(2Bc),\qquad
\frac{\partial s}{\partial c}
= 2\alpha^2 B^T G\,\operatorname{diag}(h'(2Bc))B.
\]

These are implementation notation, not new physical variables. The prototype
uses scalar automatic differentiation of the existing potential helper, so
its small-field polynomial branch is preserved. It also checks derivatives
with respect to radius, spin, and alpha.

The **28 × 28 source Jacobian at one radius** is different from the
**56 × 56 boundary-matching Jacobian** used by Newton. The local matrix
accelerates propagation of derivatives through the source. Boundary data,
linear operators, and RK4 still participate in the full differentiation.

Warmed medians from the standalone experiments, with one BLAS thread on the
interactive host:

| Calculation | Original | Optimized prototype | Speedup |
|---|---:|---:|---:|
| Ordinary source | 32.008 ms | 0.872 ms | 36.7× |
| Source work covering 56 directions: five width-12 passes versus one width-56 pass | 224.851 ms | 3.378 ms | 66.6× |
| Complete Jacobian through short radial endpoint integrations | 5.849 s | 0.105 s | 55.8× |

The short radial test uses the actual boundary helpers and target angular
grids, but integrates only 0.2 in tortoise radius from each endpoint. Its
endpoint difference is a test observable, not physical matching across the
full domain. The boundary fields in that test are dilute; separate source
checks include nonlinear fields. **A 56× full-target solve speedup has not
been established.** Full-domain matching, Newton convergence, and runtime
still need measurement.

The standalone scripts live under
[`runs/adiabatic/multiphase_two_cloud/work/`](../runs/adiabatic/multiphase_two_cloud/work/),
which is ignored by Git:

| Script | Evidence |
|---|---|
| `measure_direct_profile_source_20261023.jl` | Old/current source dimensions and source-only timing estimates |
| `cached_source_experiment_20261023.jl` | 91 checks: cached/specialized reconstruction, source values, analytic quartic coefficients, and derivatives |
| `batched_source_experiment_20261023.jl` | 96 checks: batched products/local source Jacobians and incoming derivative widths 12, 28, 56 |
| `short_radial_jacobian_experiment_20261023.jl` | 58 checks, including all 56 Jacobian columns through short radial integrations |
| `source_cache_adapter_check_20261023.jl` | 21 adapter checks with the experimental helper, including missing modes and fallback output requests |
| `source_matrix_step_check_20261023.jl` | 36 checks of the ordinary matrix source and its generic derivatives on two grids and three potentials |

These are scoped checks, not a rerun of the full project suite. The 17
existing nonlinear-source checks also passed after installing the adapter.

### Where the collaborative implementation actually stands

The working [projection.jl](../src/adiabatic/multiphase/projection.jl) now contains:

1. `nonlinear_two_phase_source(...; source_cache=nothing)`: delegates to the
   cached helper when a cache is supplied; its ordinary path retains the
   original reconstruction and projection.
2. `cached_two_phase_source`: orders dictionary coefficients by the cache,
   fills absent inputs with zero, packs the parameters, calls the numerical
   helper, and restores the requested output dictionary. Requests for modes
   outside the cache fall back to the original calculation.
3. The ordinary `local_chain_source` method: reconstructs with a matrix
   product, applies potential/geometry/weights, and projects with the transpose.

The user explicitly authorized direct insertion of the first adapter, then
added the next two functions from tested snippets. Continue that incremental
workflow and preserve the user's learning-map checkboxes.

**Still missing from the library:** the cache builder, the optimized local
source-Jacobian helper and derivative-specific method, cache forwarding through
the radial/matching calls, and configurable derivative batching in Newton.
The cache stores basis data; it does not store a previously computed source
value. It must correspond to the same component ordering and quadrature.

For the direct solve, the planned existing-function edits are
`nonlinear_two_phase_source`, `nonlinear_two_phase_rhs!`,
`two_phase_matching_residual`, `two_phase_real_residual`,
`correct_two_phase_profile`, `newton_solve`, and `residual_jacobian`, plus
the driver's cache setup and logging wrapper. The general reconstruction and
projection helpers can remain as readable references. Reusing the same work
in amplitude response later also touches `two_phase_trial_observables`,
`two_phase_profile_response`, `implicit_branch_response`, and
`implicit_observable_derivatives`.

The next implementation steps are to build the matching cache, add the fast
derivative method, connect the direct-solve path, and validate a full-domain
residual/Jacobian before preparing a new optimized snapshot. The current
running jobs use frozen v2/v3 source and receive none of these changes.

## Running jobs and their results at this checkpoint

Slurm and logs were checked at **17:25 PDT, October 23, 2026**. All three
allocations have a four-hour limit.

| Job | Snapshot | State / elapsed time | Latest result |
|---|---|---|---|
| 27241322 | `direct_fig7_20261023_v1` | FAILED after 4 s | A whitespace snapshot-directory argument caused `cd` to fail before Julia started. The wrapper now validates arguments. |
| 27241421 | `direct_fig7_20261023_v2` | RUNNING, 2 h 13 m 7 s | Accepted Newton iteration 1 at scaled residual 6.832061250638738e8. This log reports accepted iterates only. |
| 27241725 | `direct_fig7_20261023_v3` | RUNNING, 1 h 56 m 38 s | Same accepted iteration-1 residual; building correction 2, Jacobian pass 2. |

The detailed v3 log records:

| Stage | Duration | Scaled residual / outcome |
|---|---:|---|
| Initial residual | 537.04 s | 1.9064826994960344e9 |
| First Jacobian, five passes | About 4,238 s (70.6 min) | Pass durations 832.78, 839.90, 835.22, 871.10, 859.09 s |
| First full-step trial | 554.21 s | 2.253960839827775e19; rejected |
| Next trial, half step | 547.12 s | 6.832061250638738e8; accepted |
| Second Jacobian, pass 1 | 865.25 s | Completed; pass 2 started |

The first accepted correction occurred at 5,877.72 seconds on the matching
timer, approximately 98 minutes after that timer started. The accepted
residual is **64.2% lower** than the initial value, but it is still far above
the requested 1e-5 tolerance. One improving step does not establish convergence
or a useful completion estimate. The rejected full step explains why the
first Jacobian was followed by two expensive residual evaluations.

Logs are
[`direct-profile-27241421.out`](../runs/adiabatic/multiphase_two_cloud/logs/direct-profile-27241421.out)
and
[`direct-profile-27241725.out`](../runs/adiabatic/multiphase_two_cloud/logs/direct-profile-27241725.out).
The two result prefixes are `results/direct_fig7_20261023_v2` and
`results/direct_fig7_20261023_v3` within the campaign. Only their `-input.jls`
files exist at this check; neither has produced a converged solution or the
new plot curves. Both jobs were left running. On resuming, inspect their
latest logs/results alongside completing the optimization; do not treat the
short-benchmark speedup as an improvement already applied to these jobs.
