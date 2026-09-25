# Figure 2 experiment history

Preserved campaign notes; see [README.md](README.md) for current entry points.
Paths were updated during the layout migration; numerical claims and historical
settings below are retained from the original notes.

## Directory layout

- `scripts/`: cosine and polynomial drivers, submit scripts, and plotting.
  (2026-09-25: the two drivers moved to `src/adiabatic/drivers/`; references
  to `scripts/continue_*.jl` below describe the layout at the time.)
- `results/`: accepted branch CSVs and checkpoints. Figures are still written to the repository's `plots/`.
- `seeds/`: preserved checkpoints used as inputs to diagnostic experiments.
- `diagnostics/extrapolation/`: early rescaled-versus-secant comparison.
- `diagnostics/speed/`: Jacobian increment and step-size performance checks.
- `diagnostics/high_amplitude/`: residual-floor, Jacobian, and cubic-prediction
  experiments near A0=10, together with their result files and submit scripts.
- `logs/`: historical and future SLURM output.

Run submit commands from the repository root. Serialized results and old logs
retain their original provenance paths; use the locations above to find moved
files. Numerical artifacts were moved without changing their contents.

The cosine production driver still uses the secant predictor and its original restart
step rule. The successful cubic A0=10 result is a separate diagnostic; the next
bounded step comparison belongs in `diagnostics/high_amplitude/` before any
production restart.

This run extends the existing small-amplitude adiabatic continuation toward
finite amplitude. It is a preliminary cosine-potential branch, not a complete
or converged reproduction of OTT Fig. 2. Quartic and sextic short-run support
is described below; their full nonlinear branches remain unfinished.

## Quartic and sextic short runs

`scripts/continue_polynomial.jl` selects `quartic` or `sextic` with
`AXION_POTENTIAL` (default `quartic`). It passes the same selection to the
matching solve, profile reconstruction, and energy diagnostics. Each branch
has its own `results/quartic_branch.*` or `results/sextic_branch.*` CSV and JLS.
Existing checkpoints resume only if their stored potential matches.

A fresh run uses only the first entry (A0=0.0011390625) of
`seeds/preliminary_cosine_20260922.jls` as an initial guess, re-solves it for
the selected potential, and recomputes its diagnostics. No cosine rows are
copied into the polynomial branch. The radial/angular settings and secant
predictor are the same as in the cosine driver.

From the repository root, submit the two independent short runs:

```bash
sbatch --job-name=axion_quartic --export=ALL,AXION_POTENTIAL=quartic runs/adiabatic/fig2_20260922/scripts/submit_polynomial.sh
sbatch --job-name=axion_sextic --export=ALL,AXION_POTENTIAL=sextic runs/adiabatic/fig2_20260922/scripts/submit_polynomial.sh
```

Each targets A0=0.1 with a 1000-second driver budget and a 20-minute SLURM
limit. `AXION_TARGET_AMPLITUDE`, `AXION_RUN_SECONDS`, and
`AXION_DIFFERENCE_STEP` can override the defaults. Logs go to
`logs/polynomial-<job-name>-<job-id>.out`. Submit at most one job per potential
at a time because each potential uses a fixed checkpoint path.

Polynomial continuation now limits each attempt to 12 Newton corrections
(`AXION_NEWTON_MAX_ITERATIONS`). Failure halves the amplitude increment and
retries from the last accepted point; the 1e-9 residual tolerance is unchanged.
Restarts use 1.5 times the last accepted spacing, capped at 10% of A0,
without a forced minimum jump. `AXION_INITIAL_STEP` can set the first increment.
`AXION_PREDICTOR=cubic` selects the cubic predictor with secant fallback;
the default remains `secant`. `AXION_HISTORY_SPACING` defaults to 0.05.
Checkpoints record the predictor, history spacing, derivative increment,
and correction limit. Already-running jobs retain their loaded settings.

The bounded quartic comparison uses
`diagnostics/high_amplitude/submit_quartic_check.sh` (job 26305457) and
`seeds/quartic_before_bounded_correction.jls`, preserved through A0=5.58376018.
It retries A0=5.60131013 using the cubic predictor and difference scale 1e-10,
and compares the result against two half-sized steps. Results are saved to
`diagnostics/high_amplitude/quartic_correction_validation.jls`.

The original quartic job eventually accepted A0=5.60131013 after 56
corrections (709 seconds), then advanced to A0=5.61447259. The cubic/1e-10
trial did not show fast convergence (residual still about 1e-3 after 11
corrections; it failed the 12-correction benchmark), so it was not adopted for production. The unbounded job
26301502 was stopped and replaced with **26305515**, retaining secant/1e-8
while limiting attempts to 12 corrections and starting with a 0.0065 step:

```bash
sbatch --job-name=axion_quartic_bounded --time=04:00:00 --export=ALL,AXION_POTENTIAL=quartic,AXION_TARGET_AMPLITUDE=6,AXION_RUN_SECONDS=13800,AXION_DIFFERENCE_STEP=1e-8,AXION_PREDICTOR=secant,AXION_NEWTON_MAX_ITERATIONS=12,AXION_INITIAL_STEP=0.0065 runs/adiabatic/fig2_20260922/scripts/submit_polynomial.sh
```

The old job saved one further point before stopping, so replacement job
26305515 actually resumed at A0=5.6276350465.

### Isolated quartic method comparison

Array job **26305583** runs seven cases from the same preserved A0=5.58376018
seed toward A0=5.7: secant with 1e-8 or 1e-10 Jacobian differences; cubic
with 0.01 history spacing and either increment; cubic with 0.025 or 0.05
spacing and 1e-10 differences; and secant/1e-8 with a quarter-sized initial
step. All use the same grids and 1e-9 tolerance, a 12-correction cap, and a
600-second budget checked between attempts. Each has a 15-minute SLURM
limit; at most three one-core cases run concurrently.

Each writes its own CSV summary and full JLS attempt history under
`diagnostics/high_amplitude/quartic_methods/`. Production checkpoints are
never written by these cases. The automatic summary job **26305596** depends
on completion of the whole array (including failures). It writes
`diagnostics/high_amplitude/quartic_methods/SUMMARY.md`, ranking amplitude
progress per minute and comparing energy, frequency, and flux at the common
endpoint when at least two cases reach it. Rankings without that comparison
are provisional; the report does not change production settings.

Winner validation uses `diagnostics/high_amplitude/submit_quartic_validation.sh`.
The shared A0=5.60131013 point agrees with the original production solution
to 9.05e-12 in energy and 1.85e-9 in flux fractionally; radial profile arrays
agree to 4.25e-11 or better in relative norm. The first validation job
(26306865) exposed a history-selection issue during step subdivision: two
target offsets selected the same saved point, triggering secant fallback
despite enough distinct earlier points being available. That final solve
exhausted its 12-correction budget.

The cubic predictor now retains the original selection when valid, otherwise
walks backward by the requested spacing from each selected node to obtain
four distinct inputs. It still uses secant when history is insufficient.
Regression checks cover the repeated-node case (23 predictor checks pass).
Validation job **26307003** passed all 31 existing-solution and subdivided-step
checks in 1m35.9s with this fix. The inserted midpoint A0=5.6871589523 converged
in three corrections; the endpoint A0=5.7 took four, reaching residual 1.04e-12.
Compared with the original direct step, endpoint energy agrees to 2.37e-10,
net flux to 2.04e-8, and profiles to 6.30e-10 in relative norm.
This supports the winning quartic settings (cubic predictor, history spacing
0.025, Jacobian difference scale 1e-10, at most 12 corrections) at the current
grid and mode truncation. It does not establish high-amplitude discretization
convergence or validate these settings for the sextic potential.
Results are saved in `diagnostics/high_amplitude/quartic_winner_validation.jls`.

From the repository root, the reproducible commands are:

```bash
sbatch runs/adiabatic/fig2_20260922/diagnostics/high_amplitude/submit_quartic_methods.sh
sbatch --dependency=afterany:<array-job-id> runs/adiabatic/fig2_20260922/diagnostics/high_amplitude/submit_quartic_summary.sh
```

Initial short jobs: quartic **26300924**, sextic **26300925**. Both driver
configurations loaded successfully before submission, and all eight analytic
polynomial source/energy checks passed. Job logs and checkpoints record the
continuation outcomes.

The quartic short run completed at A0=0.1 with 23 saved points and final
scaled residual 4.15e-13. Its extension toward A0=6 was submitted as job
**26301502**, resuming the same checkpoint with a four-hour allocation:

```bash
sbatch --job-name=axion_quartic_to6 --time=04:00:00 --export=ALL,AXION_POTENTIAL=quartic,AXION_TARGET_AMPLITUDE=6,AXION_RUN_SECONDS=13800,AXION_DIFFERENCE_STEP=1e-8 runs/adiabatic/fig2_20260922/scripts/submit_polynomial.sh
```

The sextic short run ended normally at its runtime budget, with 22 saved
points through A0=0.0958576968 and scaled residual 6.41e-10. Its extension
toward A0=10 was submitted as job **26305114**, resuming the same checkpoint:

```bash
sbatch --job-name=axion_sextic_to10 --time=04:00:00 --export=ALL,AXION_POTENTIAL=sextic,AXION_TARGET_AMPLITUDE=10,AXION_RUN_SECONDS=13800,AXION_DIFFERENCE_STEP=1e-8 runs/adiabatic/fig2_20260922/scripts/submit_polynomial.sh
```

## Cosine run settings and history

Settings: spin=0.99, alpha=0.42; modes [(1,1),(3,1),(5,1),(3,3),(5,3),(5,5)];
16 polar and 31 azimuthal samples; radial RK4 step=0.1; rstar endpoints -100
and 150; matching coordinate 5; scaled Newton residual tolerance 1e-9.
Energy uses Simpson integration over the same finite domain. The sampled
peak field uses every tenth radial sample and is not a precise extremum.

The current driver predicts frequency and physical boundary amplitudes by
secant extrapolation in A0 from the last two accepted points, accounting for
unequal amplitude steps. It fixes the fundamental outer amplitude exactly.
Newton correction continues until the scaled residual is at most 1e-9, with
no fixed iteration cap. Nonfinite residuals, failed descent, and eight
successive corrections improving the residual by less than 0.1% are failure
conditions. Failed corrections halve the amplitude increment and retry from
the last accepted point. Successful correction counts govern step growth.
Every Newton residual is printed in the SLURM log.

The driver uses centered Jacobian differences with scale 1e-8 (overridable
with `AXION_DIFFERENCE_STEP`), while the library default remains 1e-6.
On restart, the initial amplitude increment is 1.5 times the last accepted
spacing, bounded between 10% and 25% of A0; later increments are capped at 25%.
Failed steps can shrink below 10%. The original driver used rescaled
amplitudes and a ten-correction cap and stopped at A0=0.270585259 with 21
accepted points. Those points remain in the checkpoint.

The initial five points were produced before a login-node change, in
`seeds/preliminary_cosine_20260922.jls`. The shared driver resumes that checkpoint
if its own `results/cosine_branch.jls` does not exist; otherwise it resumes its own
latest solution. Checkpoints are replaced atomically after every accepted point.
All successful points, frequencies, boundary amplitudes, and settings are retained.

Run from the repository root:

```bash
sbatch runs/adiabatic/fig2_20260922/scripts/submit.sh
```

Do not start a second job using the same output while the first is running.
The driver attempts further continuation for about 1000 seconds (checked
between solves); SLURM imposes a 20-minute wall-time limit. An unfinished
solve may be discarded at that limit, but completed points remain available.
The original resumed job is 26265137.

A four-hour extension ran as job 26265520 after the short job, and completed
when the original continuation stopping rule was reached. The tested secant
driver resumed as job 26283970 on 2026-09-23 UTC. It was stopped after
55 minutes to resume its saved checkpoint with the validated finer Jacobian
differences and larger restart spacing as job **26286366**.
The long recipe uses `scripts/submit_long.sh`, sets `AXION_RUN_SECONDS=13800`, targets A0=13,
and refreshes the plot on normal completion. Failed Newton solves reduce
the amplitude increment; repeated failures can end continuation before the
target. A time-limit kill retains the last complete checkpoint.

To repeat the long run after all existing jobs have finished:

```bash
sbatch runs/adiabatic/fig2_20260922/scripts/submit_long.sh
```

The initial delivered plot contains nine points through A0=0.12301875,
E/Fa²=216.4900641, sampled peak field=0.2890454. At that point E/A0² is
about 32.5% below the small-amplitude value. Later checkpoints may extend
this range; the CSV and SLURM logs are authoritative.

To compare the two predictors at the original rejected amplitude
A0=0.2782139804271685, run `diagnostics/extrapolation/submit_extrapolation_check.sh` through `sbatch`.
This standalone check leaves the branch unchanged and writes
`diagnostics/extrapolation/extrapolation_validation.csv`. It uses the preserved 21-point checkpoint
`seeds/seed_before_secant.jls`; `AXION_SEED_PATH` can select another saved seed.
`AXION_OUTPUT_PREFIX` optionally selects a separate continuation output.
`AXION_TARGET_AMPLITUDE` optionally changes the default target of 13.

The saved-branch comparison passed: the old scaled guess required 11 Newton
corrections (136 seconds), while the secant guess required 8 (95 seconds).
Initial residuals were 1.0566 and 0.05235, respectively; both reached the
unchanged 1e-9 tolerance. Energies agreed to 1.6e-11 fractionally. The first
run's rejection at ten corrections was reproduced before the eleventh
correction converged. This is a local performance comparison, not a promise
of the same speedup throughout the branch.

## Continuation speed check

Run `sbatch runs/adiabatic/fig2_20260922/diagnostics/speed/submit_speed_check.sh` to
repeat the comparison from `seeds/seed_before_speed_check.jls` without modifying
the main branch. Job 26286005 passed all 16 assertions. Results are in
`diagnostics/speed/continuation_speed_validation.csv`.

At A0=0.4018594396, replaying the same predictor with Jacobian difference
scale 1e-8 instead of 1e-6 reduced correction count from 10 to 3. Energy
agreed with the saved solution to 1.0e-11 fractionally; the replay residual
was 7.3e-13. The finer difference was motivated by normalized horizon
coefficients as small as 2e-11, for which the former 1e-6 perturbation was
large compared with the solution coefficient.

A direct 25% amplitude jump to A0=0.5023242995 converged in five corrections
(58 seconds). Two half-sized steps reached the same point with energy
agreement of 1.8e-13 fractionally. Frequency, flux, and boundary-amplitude
checks also passed. The residual tolerance, radial grid, angular grid,
physical equations, and six retained modes were unchanged. This validates
the larger step locally; it does not guarantee equally fast convergence
through the high-amplitude branch. Checkpoint `difference_step` records
the current solver setting; historical rows retain their original results.

After loading the BHsuper environment, redraw the latest saved points:

```bash
python runs/adiabatic/fig2_20260922/scripts/plot_branch.py
```

This writes both Figure 2 panels for quartic (black), sextic (red), and cosine
(blue) to `plots/preliminary_fig2_20260922.png` and `.pdf`. Lines connect saved
solutions in amplitude order, without smoothing or extrapolation. The plot
uses the paper's approximate axis ranges, leaving unfinished portions empty.
The plot is a snapshot; rerun this command after additional points complete.
Numerical data are in `results/{quartic,sextic,cosine}_branch.csv`, and full
checkpoints have the corresponding `.jls` suffix. High-amplitude grid and mode
convergence remain to be established.

### Fill sparse polynomial sampling

Run `sbatch runs/adiabatic/fig2_20260922/scripts/submit_refinement.sh` to
fill quartic and sextic amplitude gaps in two independent jobs. Array
**26307340** is the initial run; dependent job **26307355** redraws the plot
after both tasks finish. Each task has four hours allocated (13800 seconds
of driver budget, checked between solves).

`refine_polynomial.jl` bisects the widest gaps until spacing is at most 0.1,
the runtime budget expires, or remaining midpoint solves have failed. It
interpolates bracketed boundary data only to obtain initial guesses; each
new point must solve the nonlinear equations to residual <=1e-9 with at most
12 corrections. Failed midpoints are logged and skipped for that invocation.
Quartic uses difference scale 1e-10; sextic retains 1e-8. Solutions and
diagnostics are saved after each success to `results/{potential}_refined.jls`
and `.csv`, leaving production checkpoints untouched. Rerunning resumes the
combined saved history and retries unresolved gaps. The plotter merges new
amplitudes from these files with current production data, keeping production
values at shared amplitudes. This refines sampling in amplitude, not the
radial/angular discretization.

## Flux convention and checks

Cosine predictor comparison job **26307889** uses
`diagnostics/high_amplitude/submit_cosine_predictors.sh` and the frozen
`seeds/cosine_before_predictor_comparison.jls`. It screens secant and cubic
history spacings 0.01, 0.025, 0.05, and 0.1 at amplitude increments 0.05 and
0.1. The secant and two best cubic guesses at increment 0.1 receive bounded
12-correction solves with fixed Jacobian difference scale 1e-8. If any solve
succeeds, the fastest is checked against two half steps using energy,
frequency, flux, and boundary-amplitude agreement. Results are saved to
`diagnostics/high_amplitude/cosine_predictor_comparison.jls`. Failed solves
are recorded; submission or successful screening alone does not validate a
method. This experiment does not write production checkpoints.

Branch extension status on 2026-09-23: quartic job **26305515** continues
toward A0=6 and sextic job **26305114** toward A0=10. Cosine was restarted
from A0=9.9490175217 in job **26307661**, targeting A0=13 with four hours
allocated. Its driver now preserves the last accepted spacing on restart
and limits each correction attempt to 12 iterations before halving the
amplitude step. It retains the secant predictor and difference scale 1e-8.
The cosine submission script redraws all three branches on normal completion.
Targets are requested endpoints, not guarantees of completion within one job.

At the user's request, cosine job 26307661 was cancelled and replaced by
**26307933**, dependent on the old job ending to prevent concurrent checkpoint
writes. This trial uses `AXION_PREDICTOR=cubic`, history spacing 0.025, initial
step 0.1, and the unchanged difference scale 1e-8 and residual tolerance 1e-9.
Each attempt permits 12 corrections, with failed attempts halving the step.
It targets A0=13 with four hours allocated. The independent predictor
comparison 26307889 continues; its validation was not yet complete when this
production trial was launched. Predictor choice and spacing are checkpointed.

For m0=1, the outward flux at a finite radius is evaluated as

    F_out/Fa² = 4π (r²+a²) sum_lm Re[i m omega R_lm conj(Q_lm)].

The CSV records horizon_flux=-F_out(inner), outer_flux=F_out(outer), and
net loss total_flux=outer_flux+horizon_flux. These are finite-boundary
currents, not extrapolated asymptotic fluxes. The plotting script preserves
amplitude ordering. The balance_ratio is -F_tot/(2 Im(omega) E); it is a
linear-limit check, not an identity imposed on nonlinear configurations.
Before initial continuation, the diagnostic was checked against prescribed
outgoing/ingoing waves. The first point has balance_ratio=0.999999895.

Prior energy checks established fourth-order radial refinement at A0=1e-4
on [-100,100], and a 1.24e-8 fractional energy change extending the outer
endpoint from 150 to 200. These are small-amplitude checks; mode truncation,
angular sampling, and radial/boundary errors must be checked at larger
amplitudes before interpreting saturation. The finite inner cutoff also
changes the signed near-horizon energy and remains part of the convention.
