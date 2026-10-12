# Stationary continuation in BH parameters

## Start here

This campaign tests whether the existing nonlinear primary solver can follow
small changes in $\alpha=M_{\rm BH}\mu$ and $\chi=a/M_{\rm BH}$ at fixed outer
amplitude $A_0$. It is a local stationary-solve checkpoint toward coupled
BH–cloud evolution. The subsequent checks below solve the secondary and
driven scalar fields on the nonlinear primary backgrounds. This campaign does
not evolve the BH or any cloud population.

Current result: the primary, averaged secondary, and driven $m=0,3$ fields
can all be solved on the nearby nonlinear backgrounds. The ten-point response
run passed 133/135 checks; two coarse-grid $m=3$ horizon-flux comparisons
required another radial refinement. That targeted refinement passed all
28 checks at the original tolerances. Details and preserved failures are below.

The inputs are selected cosine solutions in
`../fig2_20260922/results/cosine_branch.jls`, at $\alpha=0.42$, $\chi=0.99$:
the dilute point (branch index 1, $A_0=0.0011390625$), and a nonlinear point
(branch index 50, $A_0=2.051758324382913$). The retained $(l,m)$ sectors
are $(1,1),(3,1),(5,1),(3,3),(5,3),(5,5)$, with 16 polar and 31 azimuthal
nodes. The reference tortoise-coordinate domain is $[-100,150]$, matching
at 5, with radial step 0.1 in instantaneous $M_{\rm BH}=1$ units.

`scripts/check_nearby_backgrounds.jl` calls the existing stationary solver
and observable routines. Newton corrects the frequency and all free boundary
amplitudes on each new Kerr background. Its initial frequency guess preserves
$\omega_R/\alpha$ from the preceding solution; this is only a predictor.
The scaled matching tolerance is $10^{-11}$, with a ForwardDiff Jacobian and
at most 12 corrections per point. No library changes are required.

## Cases and checks

The ten solves are:

1. Reproduce the input background.
2. Change $\alpha$ to 0.419 and 0.421 separately, keeping $\chi=0.99$.
3. Change $\chi$ to 0.989 and 0.991 separately, keeping $\alpha=0.42$.
4. Reach $(0.419,0.989)$ by changing $\alpha$ first or $\chi$ first.
5. Return from the joint offset to the input background.
6. Re-solve the input background and joint offset at radial step 0.05.

The comparison file checks real and imaginary frequency, energy, and horizon
energy flux: relative tolerances are $10^{-7}$ for reproducing the saved
point, $10^{-5}$ for route/return agreement, and $10^{-4}$ for radial-step
refinement. These are acceptance thresholds, not claimed measured errors.
Every accepted solution must satisfy the matching tolerance and remain
superradiant in this small neighborhood.

Output columns also include outer energy flux, peak field, Newton iterations,
and the superradiant margin $\Omega_H-\omega_R$ for the primary $m=1$ pattern.
The primary horizon torque is $F_J^H=F_E^H/\omega_R$: all retained primary
harmonics have the same ratio $m_n/\omega_n=1/\omega_R$. Fluxes are positive
into the BH or outward at the finite outer endpoint. No $f_a$ normalization
has been applied.

Holding the outer coefficient $A_0$ fixed across backgrounds does **not** hold
physical cloud energy, occupation, or the radial profile fixed. These cases
are neighboring configurations, not points selected by a physical trajectory.
The first two stages solve only the primary; secondary and driven responses
are recorded in the final stage below. Radial-step comparison alone
does not establish angular, endpoint, or boundary-series convergence.

## Run and restart

From the repository root:

```bash
mkdir -p runs/adiabatic/background_continuation/logs
sbatch --parsable runs/adiabatic/background_continuation/scripts/submit.sh \
    runs/adiabatic/background_continuation/results/amplitude2052 50
```

The submission reserves one CPU on `lr6` for up to 45 minutes. Only one job
may write a given output directory at a time. Supply the output directory and
one-based saved branch index as script arguments. The script resumes
accepted cases from that directory's `nearby_backgrounds.jls` after requiring identical
input, driver, and solver metadata. The checkpoint records the seed and driver
SHA-256 hashes and Git revision; `.jls` and `.csv` are replaced atomically
after each accepted point. `checks.csv` is written after all cases.
For a distinct run, supply a different output directory.

The original dilute driver is preserved byte-for-byte as
`scripts/check_dilute_backgrounds.jl`, matching the original checkpoint's
driver hash. To re-read that completed checkpoint and rerun its inexpensive
observable comparisons:

```bash
julia --project=. runs/adiabatic/background_continuation/scripts/check_dilute_backgrounds.jl
```

The general driver's syntax is `check_nearby_backgrounds.jl OUTPUT_DIRECTORY
BRANCH_INDEX`; its extra branch-index metadata intentionally belongs to a
separate checkpoint from the original dilute run.

## October 8, 2026 dilute run

- Source revision: `2259c7bc2bd82d1f5904b94c165f1cdfaf943d01`, with the new
  campaign driver uncommitted; its exact hash is recorded in the checkpoint.
- Historical submission command: `sbatch --parsable runs/adiabatic/background_continuation/scripts/submit.sh`
  (before adding the branch-index and output-directory arguments).
- Slurm job: `26747446`.
- Log: `logs/nearby-26747446.out`.
- Results: `results/nearby_backgrounds.{jls,csv}` and `results/checks.csv`.
- Status: completed successfully, exit `0:0`, elapsed 2 minutes 30 seconds.
  All ten solves and twenty observable comparisons passed. A subsequent
  restart skipped all accepted solves and reproduced the comparison checks.

### Results

Frequencies below are $M_{\rm BH}\omega$; energy and flux use the existing
dimensionless scalar normalization with instantaneous $M_{\rm BH}=1$.
These six distinct backgrounds all use radial step 0.1 and the same $A_0$.

| Case | $\alpha$ | $\chi$ | $M\omega_R$ | $M\omega_I$ | Code energy | Horizon energy flux |
|---|---:|---:|---:|---:|---:|---:|
| Baseline | 0.420 | 0.990 | 0.4088054935 | 1.50388869e-7 | 0.0274811590 | -8.26572006e-9 |
| Lower alpha | 0.419 | 0.990 | 0.4078970748 | 1.50206395e-7 | 0.0282872606 | -8.49785408e-9 |
| Higher alpha | 0.421 | 0.990 | 0.4097132895 | 1.50435026e-7 | 0.0266989175 | -8.03290393e-9 |
| Lower spin | 0.420 | 0.989 | 0.4088040958 | 1.45568480e-7 | 0.0274599296 | -7.99459966e-9 |
| Higher spin | 0.420 | 0.991 | 0.4088068863 | 1.53538650e-7 | 0.0275023499 | -8.44534656e-9 |
| Both lower | 0.419 | 0.989 | 0.4078957086 | 1.45904078e-7 | 0.0282656549 | -8.24814784e-9 |

All scaled matching residuals are below $1.45\times10^{-12}$. The perturbed
backgrounds require two or three Newton corrections; each half-step solve
requires one. Re-evaluating the saved baseline reproduces its frequency,
energy, and flux exactly, without a Newton correction.

| Comparison | Relative frequency error, real part | Relative growth-rate error | Relative energy error | Relative horizon-flux error |
|---|---:|---:|---:|---:|
| Two routes to joint offset | 1.36e-16 | 4.01e-10 | 1.33e-14 | 2.65e-14 |
| Return to baseline | 2.72e-16 | 1.14e-11 | 3.71e-14 | 1.60e-13 |
| Baseline, step 0.1 vs 0.05 | 1.88e-12 | 1.85e-8 | 1.10e-9 | 4.69e-10 |
| Joint offset, step 0.1 vs 0.05 | 1.86e-12 | 1.90e-8 | 1.08e-9 | 4.93e-10 |

At fixed alpha, decreasing spin from 0.990 to 0.989 reduces the dimensionless
field-amplitude growth rate by 3.205% and the magnitude of the horizon energy
flux by 3.280%. Energy changes across the alpha offsets partly reflect holding
the outer normalization coefficient fixed, so they are not adiabatic energy
changes of a conserved cloud population.

This establishes local continuation and radial-step stability for a dilute
primary. It does not establish strongly nonlinear continuation, a varying-background
cloud balance law, or the secondary/driven-channel response.

## October 8, 2026 finite-amplitude run

The same ten cases and acceptance thresholds are applied to branch point 50,
$A_0=2.051758324382913$. The saved baseline has peak dimensionless field
$\max|\Phi|\simeq1.446$, compared with 0.00306 in the dilute run. Its growth
rate is about 6.14 times the dilute rate, so this tests a profile substantially
affected by self-interaction. Earlier radial residual checks at this amplitude
are recorded in [the reproduction notes](../../../docs/adiabatic-repro.md#high-amplitude-radial-check--2026-10-02).
The new run includes actual re-solves at half the radial step.

- Source revision: `2259c7bc2bd82d1f5904b94c165f1cdfaf943d01`, with the
  generalized campaign driver uncommitted; its exact hash is in the checkpoint.
- Submission: the two-argument `sbatch` command in **Run and restart** above.
- Slurm job: `26747597`.
- Log: `logs/nearby-26747597.out`.
- Results: `results/amplitude2052/nearby_backgrounds.{jls,csv}` and `checks.csv`.
- Status: completed successfully, exit `0:0`, elapsed 3 minutes 15 seconds.
  All ten solves and twenty comparisons passed.

### Finite-amplitude results

The baseline frequency moves from $M\omega_R=0.40880549$ in the dilute run to
$0.39742366$ here. This shift is about 1.02 times the dilute binding gap
$\alpha-M\omega_R$; the nonlinear effect is therefore substantial even though
the fractional change in the total real frequency is only 2.78%.
The imaginary frequency is 6.137 times its dilute value. This comparison
is between stationary configurations at the same BH parameters.

| Case | $\alpha$ | $\chi$ | $M\omega_R$ | $M\omega_I$ | Code energy | Horizon energy flux |
|---|---:|---:|---:|---:|---:|---:|
| Baseline | 0.420 | 0.990 | 0.3974236550 | 9.22983036e-7 | 1343.7387151 | -2.43302183e-3 |
| Lower alpha | 0.419 | 0.990 | 0.3965506083 | 9.20968385e-7 | 1353.3680296 | -2.44508064e-3 |
| Higher alpha | 0.421 | 0.990 | 0.3982961863 | 9.24731255e-7 | 1334.1665541 | -2.42029413e-3 |
| Lower spin | 0.420 | 0.989 | 0.3974192736 | 9.18254691e-7 | 1343.1951052 | -2.41955727e-3 |
| Higher spin | 0.420 | 0.991 | 0.3974280127 | 9.23755208e-7 | 1344.2814003 | -2.43606461e-3 |
| Both lower | 0.419 | 0.989 | 0.3965463141 | 9.17145491e-7 | 1352.8267438 | -2.43393573e-3 |

All scaled matching residuals are below $10^{-11}$. The nearby backgrounds
need three or four Newton corrections; the two half-step solves need one
each. The saved baseline already meets the tighter requested tolerance and
requires no correction.

| Comparison | Relative frequency error, real part | Relative imaginary-frequency error | Relative energy error | Relative horizon-flux error |
|---|---:|---:|---:|---:|
| Two routes to joint offset | 5.60e-16 | 3.74e-14 | 3.66e-14 | 3.63e-14 |
| Return to baseline | 4.64e-14 | 1.74e-8 | 1.46e-12 | 9.72e-12 |
| Baseline, step 0.1 vs 0.05 | 9.71e-11 | 9.81e-9 | 8.89e-10 | 1.08e-8 |
| Joint offset, step 0.1 vs 0.05 | 9.51e-11 | 7.86e-9 | 8.65e-10 | 1.07e-8 |

The baseline return comparison is with the original accepted solution, whose
scaled residual is $9.53\times10^{-12}$; the route comparison uses two newly
corrected solutions. These agreement numbers are not bounds on discretization
or physical-model errors.

The additional outward-flux comparison, evaluated from the saved rows, is
less precise: half-step changes are $6.79\times10^{-5}$ at baseline and
$6.71\times10^{-5}$ at the joint offset, relative to the outward flux itself.
That flux is small here ($6.68\times10^{-7}$ at baseline, compared with
$|F_H|=2.43\times10^{-3}$); net scalar flux changes by at most
$2.95\times10^{-8}$ relative under refinement. These outward/net comparisons
are diagnostics in addition to the twenty acceptance checks.

Lowering spin from 0.990 to 0.989 at fixed $\alpha$ and $A_0$ reduces
$M\omega_I$ by 0.5123% and the magnitude of horizon extraction by 0.5534%.
The corresponding dilute changes were 3.205% and 3.280%. This is a response
of the stationary family at fixed outer amplitude, not a predicted response
along a conserved-population or coupled-BH trajectory. In particular, the
nonlinear cloud evolution must use the branch energy derivative and boundary
fluxes; $\dot E=2\omega_I E$ is only a dilute-limit check (the saved nonlinear
`balance_ratio` here is about 0.98059).

Both finite-amplitude and frozen dilute restart checks skip all accepted
solves and reproduce their observable comparisons. These results establish
local primary continuation at this nonlinear amplitude, with the stated
radial checks. Angular/endpoint convergence and the secondary/driven-field
response on the shifted backgrounds remain separate checks at this stage.

## Secondary and driven fields on the finite-amplitude backgrounds

`scripts/check_secondary_backgrounds.jl` takes the completed primary checkpoint
and solves the averaged secondary, the driven $m=0$ field, and the driven
$m=3$ field for each of its ten entries. It uses the saved secondary at
primary branch index 50 as the baseline seed and follows the same source-point
routes as the primary checks. The secondary's real frequency divided by alpha
is preserved only in the initial guess; Newton then corrects the frequency
and free boundary amplitudes on the new primary profile and Kerr background.

The secondary modes are $(2,2),(4,2),(6,2)$; driven modes are
$(0,0),(2,0),(4,0)$ and $(3,3),(5,3),(7,3)$. The 16-by-31 angular grid is
unchanged. Secondary and driven near-zone integrations use twice the primary
step, as in the existing solver. Thus the half-step entries refine both the
primary and these near-zone solves. The radiative far leg retains its checked
start at $r_*=5000$ and spacing 0.25; it is not refined in this experiment.
Secondary matching tolerance is $10^{-8}$ with a ForwardDiff Jacobian.

For the secondary reference normalization $A_2=1$, let $\hat E_2$ be its
computed energy. The saved coefficients include

$$
\frac{F_{2,H}}{\hat E_2},\qquad
\frac{F_{0,H}}{\hat E_2},\qquad
\frac{F_{3,\infty}}{\hat E_2^2},\qquad
\frac{F_{3,H}}{\hat E_2^2}.
$$

These remove the arbitrary secondary normalization. At a chosen physical
population in code units, multiply the first two by $E_2$ and the latter two
by $E_2^2$. They still include the primary-amplitude dependence: these are not
the additionally $E_1$-factored spline coefficients of the two-cloud ODE.
The driven frequencies are $\omega_0=2\omega_{1,R}-\omega_{2,R}$ and
$\omega_3=2\omega_{2,R}-\omega_{1,R}$. The recorded radiative threshold margin
is $\omega_3-\alpha$ in instantaneous $M=1$ units.

The 135 acceptance checks comprise seven observable comparisons for each of
five reference/route/refinement pairs, plus ten checks on each configuration:
primary reconstruction, secondary matching and energy balance, each driven
field's power identity and relative matching residual, each matching matrix's
condition number, and the expected channel signs/thresholds. Comparison
tolerances are $10^{-6}$ for the saved baseline, $10^{-5}$ for route and return,
and $10^{-4}$ for radial refinement. Secondary energy balance must agree within
$10^{-4}$; driven source-work/flux balances within $10^{-5}$. Both relative
driven matching residuals must be below $10^{-9}$ and column-scaled matrix
condition numbers below $10^8$.

Run from the repository root:

```bash
sbatch --parsable runs/adiabatic/background_continuation/scripts/submit_secondary.sh
```

One CPU, 45-minute time limit. The driver checkpoints each completed background
and refuses a resume with changed input/driver metadata. It saves the secondary
and driven boundary amplitudes, but not the large reconstructed profile arrays;
these can be regenerated from the checkpoint and primary inputs. To rerun the
completed checkpoint's comparisons:

```bash
julia --project=. runs/adiabatic/background_continuation/scripts/check_secondary_backgrounds.jl \
    runs/adiabatic/background_continuation/results/amplitude2052/nearby_backgrounds.jls \
    runs/adiabatic/background_continuation/results/amplitude2052/secondary
```

### October 8, 2026 submission

- Source revision: `0ef376427ce99b8d9bbb271b1b7f6899c9a3c6bd`; new driver
  uncommitted, with its hash and all input hashes recorded in the checkpoint.
- Slurm job: `26749758`; log: `logs/secondary-26749758.out`.
- Primary input: `results/amplitude2052/nearby_backgrounds.jls`.
- Secondary seed: `../multimode_secondary/results/cosine_secondary.jls`, index 50.
- Outputs: `results/amplitude2052/secondary/secondary_backgrounds.{jls,csv}`
  and `secondary_checks.csv` in the same directory.
- Status: all ten backgrounds solved; 133 of 135 checks passed. The job
  exited `1:0` after 2 minutes 52 seconds because the two $m=3$ horizon-flux
  refinement comparisons exceeded their $10^{-4}$ relative tolerance.
  The original failure records are preserved in `secondary_checks.csv`.

### Secondary and driven-field results

At all ten configurations the secondary grows, $m=0$ carries energy into the
BH, and $m=3$ radiates outward above its mass threshold. The positive margin
$M\omega_3-\alpha$ ranges from 0.009492 to 0.009684. The matching matrices
have column-scaled condition numbers 188–200 for $m=0$ and 154–164 for $m=3$.
No nearly singular matching matrix appears in these local samples.

| Case | $\alpha$ | $\chi$ | $M\omega_{2,R}$ | $M\omega_{2,I}$ | $F_{0,H}/\hat E_2$ | $F_{3,\infty}/\hat E_2^2$ |
|---|---:|---:|---:|---:|---:|---:|
| Baseline | 0.420 | 0.990 | 0.4135056505 | 2.29141014e-8 | 4.42216142e-3 | 7.67755296e-7 |
| Lower alpha | 0.419 | 0.990 | 0.4125211958 | 2.28490507e-8 | 4.48536708e-3 | 7.75976895e-7 |
| Higher alpha | 0.421 | 0.990 | 0.4144898835 | 2.29751033e-8 | 4.35905591e-3 | 7.59372768e-7 |
| Lower spin | 0.420 | 0.989 | 0.4135097803 | 2.28599708e-8 | 4.41346620e-3 | 7.65763311e-7 |
| Higher spin | 0.420 | 0.991 | 0.4135015393 | 2.29682182e-8 | 4.43081064e-3 | 7.69735840e-7 |
| Both lower | 0.419 | 0.989 | 0.4125252945 | 2.27955105e-8 | 4.47676082e-3 | 7.74023721e-7 |

All table entries use primary step 0.1 and secondary/driven near-zone step 0.2.
The saved baseline observables reproduce within $7.0\times10^{-13}$ relative.
Two routes to the joint offset agree in all seven observables within
$2.35\times10^{-12}$; returning to baseline agrees within $2.78\times10^{-11}$.
These are solver consistency checks, not independent accuracy estimates.

Across all ten cases, the secondary energy-balance discrepancy is at most
$2.56\times10^{-6}$, while source-work/flux discrepancies are below
$7.46\times10^{-8}$ for $m=0$ and $2.99\times10^{-8}$ for $m=3$.
Driven matching residuals relative to the source residual are below
$3.74\times10^{-15}$ and $2.37\times10^{-14}$, respectively.

Primary-step refinement 0.1 to 0.05 (also refining the secondary and driven
near-zone steps from 0.2 to 0.1) gives:

| Observable | Baseline relative change | Joint-offset relative change |
|---|---:|---:|
| $\omega_{2,R}$ | 1.70e-10 | 1.71e-10 |
| $\omega_{2,I}$ | 5.89e-6 | 5.90e-6 |
| $\hat E_2$ | 9.10e-8 | 9.17e-8 |
| $F_{2,H}/\hat E_2$ | 3.78e-6 | 3.61e-6 |
| $F_{0,H}/\hat E_2$ | 3.72e-8 | 4.06e-8 |
| $F_{3,\infty}/\hat E_2^2$ | 7.92e-8 | 7.79e-8 |
| $F_{3,H}/\hat E_2^2$ | **1.590e-4** | **1.515e-4** |

Only the final row exceeds the predeclared $10^{-4}$ tolerance. At baseline,
$F_{3,H}/F_{3,\infty}=-2.00\times10^{-7}$, so this concerns a very small
component of the radiative channel; it remains part of the BH flux budget.
Its failing tolerance is not relaxed or replaced by normalization to the
larger infinity flux.

### Targeted further radial refinement

`scripts/refine_secondary_backgrounds.jl` re-solves the primary, secondary,
and both driven fields at primary step 0.025 at baseline and the joint offset.
It starts from their accepted step-0.05 solutions and applies the same
$10^{-4}$ relative observable tolerance to this next refinement. The
original ten-point results and failing checks are left intact.

```bash
sbatch --parsable runs/adiabatic/background_continuation/scripts/submit_secondary_refinement.sh
```

- Job: `26749926`; log: `logs/secondary-refine-26749926.out`.
- One CPU, 30-minute limit; source revision as above, with the exact new driver
  hash and input hashes in its checkpoint.
- Output directory: `results/amplitude2052/secondary/quarter_step/`, containing
  `refined_responses.{jls,csv}` and `refinement_checks.csv`.
- Status: completed successfully, exit `0:0`, elapsed 3 minutes 15 seconds.
  Both backgrounds pass all 28 checks. Resuming the completed refinement
  reproduces the checks without rerunning any stationary solves.

The 0.05-to-0.025 primary-step comparison (secondary/driven near-zone steps
0.1 to 0.05) gives:

| Observable | Baseline relative change | Joint-offset relative change |
|---|---:|---:|
| $\omega_{2,R}$ | 1.07e-11 | 1.08e-11 |
| $\omega_{2,I}$ | 3.57e-7 | 3.58e-7 |
| $\hat E_2$ | 5.74e-9 | 5.79e-9 |
| $F_{2,H}/\hat E_2$ | 1.43e-7 | 1.38e-7 |
| $F_{0,H}/\hat E_2$ | 3.76e-8 | 3.69e-8 |
| $F_{3,\infty}/\hat E_2^2$ | 5.00e-9 | 4.92e-9 |
| $F_{3,H}/\hat E_2^2$ | **5.10e-6** | **4.87e-6** |

The previously failing $m=3$ horizon-flux changes shrink by factors 31.21 and
31.13. This supports convergence at the finer spacing; it does not establish
pure fourth-order error for that observable. The secondary imaginary-frequency
and $m=3$ infinity-flux changes shrink by approximately 16. The $m=0$ coefficient
changes stay near $4\times10^{-8}$, consistent with the previously observed
solver-tolerance floor rather than a measured RK4 convergence order.
On the finest grid, the secondary energy balance agrees within
$1.74\times10^{-8}$, and driven power balances within $4.10\times10^{-9}$
($m=0$) and $1.15\times10^{-10}$ ($m=3$).

The refinement rows retain the source labels `baseline_half_step` and
`combined_half_step`; their `primary_step` column is **0.025**. For the next
calculation, use primary spacing 0.05 or finer when requiring $10^{-4}$
relative accuracy in the small $m=3$ horizon contribution at these tested
points. Only baseline and the joint offset have this extra refinement;
the other backgrounds retain their original spatial settings.

The original 135-check script intentionally still exits with its two recorded
coarse-grid failures on replay. Replay the successful targeted follow-up with:

```bash
julia --project=. runs/adiabatic/background_continuation/scripts/refine_secondary_backgrounds.jl
```

This completes the local stationary-response checkpoint for all retained scalar
channels at $A_0\simeq2.052$. It does not establish angular/endpoint convergence
on these backgrounds, independently validate the transfer prefactors, or
include secondary feedback on the primary. The next coupled-evolution step
is the cloud energy/angular-momentum balance on the changing BH background,
using fixed reference units and a specified $f_a$ and seed prescription.

## Coupled-balance checkpoint — October 8, 2026

The balance equations now live in
[`src/adiabatic/primary_secondary/coupled_balance.jl`](../../../src/adiabatic/primary_secondary/coupled_balance.jl).
The check script below loads them through the library and supplies only the
saved-data/analytic checkpoint. Coupled time integration is paused for
piece-by-piece review; see the [draft status](../coupled_local/README.md).

The derivation is in [Coupled BH–cloud balance](../../../docs/coupled-bh-cloud-balance.md).
The four state variables are BH mass ratio $b=M/M_0$, BH angular momentum
$j=J/M_0^2$, and cloud energies $q_i=E_i^{\rm physical}/(\epsilon_f M_0)$.
The existing code energies are $e_i=q_i/b$. The prototype combines the
instantaneous two-cloud balance matrix with separate horizon mass and torque
updates. It is the leading small-cloud-mass model: bulk metric work and cloud
angular-momentum background-response terms are written down but not evaluated
on the nonlinear profiles. No coupled trajectory has been integrated.

Run the inexpensive saved-data/algebra checkpoint on the login node:

```bash
julia --compiled-modules=existing --project=. \
    runs/adiabatic/background_continuation/scripts/check_coupled_balance.jl
```

The compiled-module flag permits existing Julia caches without attempting to
write missing optional test-extension caches in the read-only home directory.
The checks perform no stationary solves or Slurm submissions. They verify:

- the existing nonlinear two-cloud rates at the saved baseline;
- dilute primary/secondary growth and $m=0,3$ transfer multiplicities;
- horizon-only BH updates and the mass contribution to the spin derivative;
- physical/code-energy and torque conversion when $M/M_0\ne1$;
- the energy and angular-momentum audit identities of the stated approximation;
- an independent two-mode model with conserved occupations under externally
  prescribed slow frequency changes, including a counterexample to omitting
  metric work while retaining the background angular response.

All **44 checks pass**. This validates the algebra, units, and stated limits;
it does not calculate the actual metric work or establish its smallness on
our nonlinear branch. Library code and the existing fixed-BH integrator are
unchanged. The driver is a tested prototype for the next implementation step.

`results/balance/leading_rhs_example.csv` and `balance_checkpoint.jls` contain
a local RHS example with the saved coarse baseline $A_0=2.051758324382913$,
$e_1=1343.738715$, illustrative $e_2=0.1$, $b=1$, $\chi=0.99$, and
$f_a=10^{16}$ GeV. The secondary frequency's energy derivative comes from the
existing branch spline, not the fixed-amplitude background differences.
The JLS metadata records the code revision and hashes of the driver, primary,
secondary responses, and five branch-table inputs. These physical choices
are not initial conditions selected for a production trajectory and do not
impose a QCD-axion mass–decay-constant relation.

| Quantity | Local example |
|---|---:|
| $\epsilon_f=(f_a/M_{\rm Pl})^2$ | 6.70883075e-7 |
| Cloud mass fraction $\epsilon_f(e_1+e_2)$ | 9.01558649e-4 |
| Balance-matrix condition number | 13.5982 |
| $db/d\tau$ | -1.33560091e-9 |
| $dj/d\tau$ | -4.10715128e-9 |
| $d\chi/d\tau$ | -1.46266147e-9 |
| $d\alpha/d\tau$ | -5.60952384e-10 |
| $dq_1/d\tau$ | 1.51063808e-3 |
| $dq_2/d\tau$ | 4.79496844e-4 |

Here $\tau=t/M_0$; no physical duration is inferred without specifying $M_0$.
The apparent difference between BH and cloud rate magnitudes is the
$\epsilon_f$ in the definition of $q_i$. The next work is a trajectory driver
that refreshes stationary data as $(e_1,\alpha,\chi)$ changes and monitors the
light-cloud and regular-response assumptions. More accurate balance also
requires the work/background-response calculation described in the derivation.
