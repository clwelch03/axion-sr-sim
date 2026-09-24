# Initial solver comparison, 2026-09-24

All four requested approaches have been implemented and run on the frozen
cosine problem at A0=10.33565530025236, a 0.1 jump beyond the checkpoint where
production continuation stalled. Fifteen case summaries cover the original
nine cases and six follow-ups. This is a local benchmark, not a general solver
ranking. Production continuation has not been changed by these experiments.

## Main results

Times below measure corrections after initial Jacobian warmup; setup adds
roughly 10–26 seconds and is recorded separately. All timings use the same
cached experimental RHS. They should not be interpreted as an isolated speedup
relative to the older production RHS, which recomputes more intermediate data.

| Method | Corrections | Solve seconds | Outcome |
|---|---:|---:|---|
| Shooting, finite differences 1e-8 | 20 | 152.6 | Did not converge; residual 4.70e-6 |
| Shooting, finite differences 1e-10 | 4 | 30.5 | Converged; original-integrator replay 3.54e-13 |
| Shooting, automatic differentiation | 4 | 20.6 | Converged; replay 2.00e-12 |
| AD plus Deuflhard-inspired controller | 5 | 15.8 | Converged; replay 2.72e-12 |
| Four-segment shooting, tight local tolerance | 3 | 31.1 | Converged; replay 1.15e-12 |
| Eight-segment shooting, aligned grid and tight tolerance | 3 | 18.5 | Converged; replay 1.59e-13 |
| Radial collocation, spacing 0.25 | 2 | 5.4 | Converged discrete problem; mesh error remains |
| Radial collocation, spacing 0.0625 | 2 | 17.7 | Converged and agrees with refined shooting observables |

The Deuflhard controller is a simplified-correction contraction implementation,
not the complete NLEQ-ERR algorithm. Single-run timings on different compute
nodes do not establish a statistically significant difference between methods
that are only a few seconds apart.

## What caused the finite-difference difficulty?

The AD matching Jacobian has condition number about 6.63e9 in the current
coordinates. Relative to AD, the 1e-8 finite-difference Jacobian differs by
about 0.9% in matrix norm, and its resulting Newton correction differs by
211%. At 1e-10 these differences fall to 1.09e-6 and 1.02e-4 respectively.
The actual solve succeeds with 1e-10 and with AD at the same fixed target.
This identifies an actionable derivative-accuracy problem in the tested
configuration, not merely an inadequate predictor. Naively scaling variables
by their magnitudes worsened the measured condition number; this experiment
does not justify that particular scaling.

The audit saves the complete AD Jacobian and the prediction. See
`runs/adiabatic/radialsolvertest/JACOBIAN_AUDIT.md` and `jacobian_audit.jls`.
The audit concerns the discrete RK4 map, not the continuum differential equation.

## Verification and radial accuracy

- Expanded boundary/RHS/derivative/root checks passed all 24 assertions,
  including checks for the aligned multiple-shooting mesh (job 26322075).
- AD subdivision and agreement with the Deuflhard-style correction passed
  17 assertions. Both routes reach the same energy, flux, frequency, and
  boundary coefficients. Results are in `ad_validation.jls`.
- Multiple shooting needs a tighter local matching tolerance than the original
  two-sided residual tolerance to guarantee a comparably small long-shooting
  replay residual. Segment boundaries also must align with the original RK4
  grid for an identical discretization. Both effects were isolated in follow-ups.
- Collocation spacings 0.5, 0.25, 0.125, and 0.0625 show approximately fourth-order
  convergence: the final three meshes give successive-difference ratios 15.92
  for energy, 15.93 for flux, and 15.91 for the real frequency.
- At spacing 0.0625, native collocation energy differs from AD shooting at
  RK4 spacing 0.025 by about 4.74e-9 fractionally; flux differs by 3.01e-7.
  The 0.025 shooting reference itself has finite discretization error.
- A collocation solution need not have a tiny shooting replay residual on a
  different radial discretization. Use its native profile for energy and test
  mesh convergence; do not rank it by long-shooting mismatch alone.

Angular resolution, asymptotic boundary truncation, and harmonic truncation
were held fixed. These experiments do not establish convergence in those choices.

## Recommendation

AD shooting is the smallest well-supported change to the current production
method: it preserves the matching formulation and passed step subdivision.
The 1e-10 finite-difference setting is an even smaller interim change, although
AD avoids having to retune that increment. Keep radial collocation as a strong
alternative: it converged in two corrections throughout the mesh study and
its refined observables agree with refined shooting. Before selecting a
long-term solver, repeat at other amplitudes and potentials and record repeated
timings; do not infer the complete branch behavior from this one target.

Raw results, full iterates, and regenerated comparison tables are in
`runs/adiabatic/radialsolvertest/`. Run `summarize.py` to regenerate `SUMMARY.md`.
