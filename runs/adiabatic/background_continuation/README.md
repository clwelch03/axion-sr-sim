# Stationary continuation in BH parameters

## Start here

This campaign tests whether the existing nonlinear primary solver can follow
small changes in $\alpha=M_{\rm BH}\mu$ and $\chi=a/M_{\rm BH}$ at fixed outer
amplitude $A_0$. It is a local stationary-solve checkpoint toward coupled
BH–cloud evolution. It does not evolve the BH or any cloud population.

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
Only the primary is solved here; the secondary and driven responses on these
new backgrounds remain a subsequent checkpoint. Radial-step comparison alone
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
response on the shifted backgrounds remain separate checks.
