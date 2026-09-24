# Radial solver experiments

See [RESULTS.md](RESULTS.md) for the completed initial comparison and conclusions.

Isolated comparisons for the Figure 2 cosine branch. Production solver files
and checkpoints are not changed. Code lives here; frozen inputs, saved iterates,
results, and SLURM logs live in `runs/adiabatic/radialsolvertest/`.

## Common problem

`seed_cosine.jls` freezes the production checkpoint at A0=10.235655300252361.
Every method receives the same cubic extrapolation (history spacing 0.025) to
A0=10.33565530025236. Preserve spin 0.99, alpha 0.42, six harmonics, angular
quadrature 16 by 31, endpoints -100 and 150, matching location 5, and boundary
series orders 3 and 7. Shooting uses the saved RK4 spacing 0.1. No finite
difference Jacobian increment is used by the AD methods.

## Cases

1. Two-sided shooting with production centered differences, increment 1e-8.
2. Same with increment 1e-10.
3. Same with ForwardDiff through the real-packed discrete RK4 map and boundary series.
4. AD shooting with Deuflhard-inspired simplified-correction damping.
5. Four-segment multiple shooting, local AD and sparse Jacobian.
6. Eight-segment multiple shooting, local AD and sparse Jacobian.
7. Fourth-order Hermite–Simpson/Lobatto radial collocation, spacing 0.5.
8. Same at spacing 0.25.
9. Same at spacing 0.125.

All use a 1e-9 native residual tolerance, at most 20 corrections, and the same
physical bound-frequency admissibility constraint. Multiple shooting integrates
inward on the outer half and outward on the inner half. Collocation divides
each defect by interval length. This normalization is not a shooting mismatch;
native convergence and independent replay agreement are reported separately.

The Deuflhard experiment factors the Jacobian, computes a trial iterate, and
uses the same factorization to compute a simplified correction at that trial.
It accepts contraction of this correction norm and estimates reduced damping
from departure from the linear model. This is an intentionally limited
implementation of the idea, not the complete NLEQ-ERR algorithm, and its
performance must not be attributed to all Deuflhard solvers.

## Verification and timing

`check_methods.jl` compares the generic boundary and angular/radial kernel with
the production functions, checks local AD derivatives against directional
differences, verifies correction on a known root, and checks the assembled
sparse spatial Jacobians. The experiments reuse production boundary recurrences
and energy quadrature. A separate real-packed angular kernel avoids hard-coded
Float64 storage without changing production functions.

Setup time includes initial Jacobian evaluation and compilation. Solve time
starts afterward; independent replay and diagnostics are outside it. Initial
guess construction for spatial methods is included in setup. Every iteration
saves the current vector to `<method>.jls.iterate`, so a failed correction can
be diagnosed directly. Residual-call counts in saved results count the Newton
driver's explicit evaluations; they exclude evaluations internal to Jacobian
construction and are not total RHS-work counts.

Replay uses the original production integrator, not the optimized benchmark
kernel. Collocation also reports energy integrated on its native profile.
Compare both frequency components, energy, and signed flux, and inspect mesh
convergence before ranking methods with different radial discretizations.
One difficult target tests local behavior; a winner still needs subdivision
and additional-amplitude checks before promotion to production.

## Run

Use the isolated Project/Manifest, not the root environment:

```bash
julia --project=src/adiabatic/radialsolvertest -e 'using Pkg; Pkg.instantiate()'
sbatch src/adiabatic/radialsolvertest/submit_checks.sh
sbatch --dependency=afterok:<checks-job> src/adiabatic/radialsolvertest/submit.sh
python src/adiabatic/radialsolvertest/summarize.py
```

The array allows three simultaneous single-core jobs, one hour and 8 GB each.
Generated result binaries are stored under `runs/`, not among source files.

Initial submission (2026-09-24): checks **26321963** passed 22 assertions;
comparison array **26321971** follows those checks, and report job **26321990**
follows the array, including failed cases. A scheduler completion code alone
does not imply Newton convergence; inspect the recorded `status` field.

Follow-ups: Jacobian audit **26322002**; multiple-shooting tolerance 1e-12
**26322029**; AD shooting at RK4 spacing 0.05 **26322030** and 0.025
**26322046**; collocation spacing 0.0625 **26322031**; AD subdivision validation
**26322044** (17 checks passed). Eight-segment shooting **26322061** aligns its
segment endpoints with the original RK4 grid and uses tolerance 1e-12.
Without alignment, uniformly spaced shooting segments slightly change the RK4
step lengths because each subinterval must contain an integer number of steps.
That changes the discrete problem and can leave a nonzero long-shooting replay
mismatch even after local residuals converge. The aligned experiment separates
that discretization effect from accumulated local matching error.

Follow-up settings are controlled by `AXION_NATIVE_TOLERANCE`,
`AXION_RK_STEP`, `AXION_RADIAL_RESOLUTION`, `AXION_ALIGN_MESH=true`, and
`AXION_LABEL_SUFFIX` (use a distinct suffix to preserve earlier results).
The original nine cases remain reproducible with these unset. To repeat the
Jacobian audit, run `audit_jacobian.jl` in the same environment; to repeat the
subdivision check, run `validate_ad.jl` after the AD cases have completed.

## Method references

- [ForwardDiff Jacobians](https://juliadiff.org/ForwardDiff.jl/stable/user/api/).
- [Shooting and variational equations](https://bifurcationkit.github.io/BifurcationKitDocs.jl/stable/bvpShooting/).
- [Ontañón and Alcubierre: global Newton methods, especially Appendix A](https://arxiv.org/pdf/2103.13993).
- [Fourth-order collocation and residual control](https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_bvp.html).
- [Original axion calculation, Appendix A](https://arxiv.org/html/2201.04382v3).
