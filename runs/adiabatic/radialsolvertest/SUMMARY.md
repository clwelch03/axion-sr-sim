# Radial solver comparison

All cases use the same frozen cosine seed, cubic prediction, and A0 increment 0.1.
Setup includes profile construction and initial residual/Jacobian compilation.
Solve time excludes that warmup and the independent production replay.
These are single-target timings, not broad performance claims.

| Method | Status | Setup s | Solve s | Corrections | Native residual | Shooting replay |
|---|---|---:|---:|---:|---:|---:|
| collocation_0125 | converged | 22.8 | 15.0 | 2 | 1.37e-11 | 0.000419 |
| collocation_025 | converged | 12.6 | 5.4 | 2 | 5.58e-12 | 0.00449 |
| collocation_050 | converged | 11.4 | 3.3 | 2 | 1.91e-12 | 0.0858 |
| collocation_050_h00625 | converged | 17.4 | 17.7 | 2 | 2.19e-11 | 0.000711 |
| multiple_4 | converged | 16.0 | 16.2 | 2 | 7.95e-10 | 6.94e-06 |
| multiple_4_tight | converged | 25.9 | 31.1 | 3 | 3.47e-15 | 1.15e-12 |
| multiple_8 | converged | 16.0 | 13.3 | 2 | 3.77e-10 | 6.93e-06 |
| multiple_8_aligned_tight | converged | 15.7 | 18.5 | 3 | 1.24e-15 | 1.59e-13 |
| multiple_8_tight | converged | 25.9 | 31.6 | 3 | 3.95e-15 | 3.99e-06 |
| shooting_ad | converged | 17.0 | 20.6 | 4 | 2.1e-12 | 2e-12 |
| shooting_ad_deuflhard | converged | 10.6 | 15.8 | 5 | 2.73e-12 | 2.72e-12 |
| shooting_ad_rk0025 | converged | 18.5 | 35.6 | 3 | 2.52e-11 | 2.52e-11 |
| shooting_ad_rk005 | converged | 22.1 | 30.9 | 3 | 2.35e-11 | 2.35e-11 |
| shooting_fd_1e10 | converged | 10.7 | 30.5 | 4 | 3.69e-13 | 3.54e-13 |
| shooting_fd_1e8 | iteration_limit | 10.2 | 152.6 | 20 | 4.7e-06 | 4.7e-06 |

15 case summaries available (nine original cases plus any follow-ups).

Collocation residuals measure radial defects; shooting residuals measure matching.
A small native residual does not establish discretization accuracy.

## Observable agreement with AD shooting

| Method | Relative energy difference | Relative flux difference | Relative Im(omega) difference |
|---|---:|---:|---:|
| collocation_0125 | 8.58e-08 | 2.24e-05 | 6.84e-07 |
| collocation_025 | 1.08e-06 | 6.71e-05 | 1.09e-05 |
| collocation_050 | 2.05e-05 | 0.00149 | 0.000205 |
| collocation_050_h00625 | 1.59e-07 | 2.8e-05 | 1.34e-06 |
| multiple_4 | 5.38e-09 | 1.42e-09 | 4.08e-07 |
| multiple_4_tight | 1.47e-12 | 1.5e-11 | 1.11e-08 |
| multiple_8 | 1.66e-09 | 1.71e-07 | 4e-07 |
| multiple_8_aligned_tight | 1.47e-12 | 1.5e-11 | 1.11e-08 |
| multiple_8_tight | 7.01e-09 | 1.73e-07 | 2.1e-08 |
| shooting_ad | 0 | 0 | 0 |
| shooting_ad_deuflhard | 1.46e-12 | 1.49e-11 | 1.11e-08 |
| shooting_ad_rk0025 | 1.64e-07 | 2.83e-05 | 1.16e-06 |
| shooting_ad_rk005 | 1.53e-07 | 2.7e-05 | 1.09e-06 |
| shooting_fd_1e10 | 1.12e-14 | 3.6e-13 | 3.07e-10 |

## Collocation versus refined shooting (RK4 spacing 0.025)

| Method | Relative native energy difference | Relative flux difference | Relative Im(omega) difference |
|---|---:|---:|---:|
| collocation_0125 | 7.79e-08 | 5.92e-06 | 4.76e-07 |
| collocation_025 | 1.24e-06 | 9.54e-05 | 1.2e-05 |
| collocation_050 | 2.07e-05 | 0.00151 | 0.000207 |
| collocation_050_h00625 | 4.74e-09 | 3.01e-07 | 1.75e-07 |

Deuflhard denotes a simplified-correction contraction controller inspired by NLEQ-ERR;
it is not a full reproduction of that algorithm. Default native residual tolerance is
1e-9 with a 20-correction budget. Suffixed cases override settings as listed below.

| Method | Native tolerance | RK4 spacing | Spatial resolution |
|---|---:|---:|---:|
| collocation_0125 | 1e-9 | 0.1 | see method name |
| collocation_025 | 1e-9 | 0.1 | see method name |
| collocation_050 | 1e-9 | 0.1 | see method name |
| collocation_050_h00625 | 1.0e-9 | 0.1 | 0.0625 |
| multiple_4 | 1e-9 | 0.1 | see method name |
| multiple_4_tight | 1.0e-12 | 0.1 | 2.0 |
| multiple_8 | 1e-9 | 0.1 | see method name |
| multiple_8_aligned_tight | 1.0e-12 | 0.1 | 4.0 |
| multiple_8_tight | 1.0e-12 | 0.1 | 4.0 |
| shooting_ad | 1e-9 | 0.1 | see method name |
| shooting_ad_deuflhard | 1e-9 | 0.1 | see method name |
| shooting_ad_rk0025 | 1.0e-9 | 0.025 | 0.0 |
| shooting_ad_rk005 | 1.0e-9 | 0.05 | 0.0 |
| shooting_fd_1e10 | 1e-9 | 0.1 | see method name |
| shooting_fd_1e8 | 1e-9 | 0.1 | see method name |

## Collocation refinement

Ratios compare successive differences on meshes whose spacing halves.
Fourth-order asymptotic convergence predicts a ratio near 16.

- collocation_050 → collocation_025 → collocation_0125: native_energy ratio 16.679.
- collocation_050 → collocation_025 → collocation_0125: total_flux ratio 15.854.
- collocation_050 → collocation_025 → collocation_0125: omega_real ratio 16.773.
- collocation_025 → collocation_0125 → collocation_050_h00625: native_energy ratio 15.916.
- collocation_025 → collocation_0125 → collocation_050_h00625: total_flux ratio 15.93.
- collocation_025 → collocation_0125 → collocation_050_h00625: omega_real ratio 15.906.
