# Preliminary cosine continuation — 2026-09-22

## Directory layout

- `branch/`: production driver, submit scripts, accepted cosine CSV/checkpoint,
  and plotting script. Figures are still written to the repository's `plots/`.
- `seeds/`: preserved checkpoints used as inputs to diagnostic experiments.
- `diagnostics/extrapolation/`: early rescaled-versus-secant comparison.
- `diagnostics/speed/`: Jacobian increment and step-size performance checks.
- `diagnostics/high_amplitude/`: residual-floor, Jacobian, and cubic-prediction
  experiments near A0=10, together with their result files and submit scripts.
- `logs/`: historical and future SLURM output.

Run submit commands from the repository root. Serialized results and old logs
retain their original provenance paths; use the locations above to find moved
files. Numerical artifacts were moved without changing their contents.

The production driver still uses the secant predictor and its original restart
step rule. The successful cubic A0=10 result is a separate diagnostic; the next
bounded step comparison belongs in `diagnostics/high_amplitude/` before any
production restart.

This run extends the existing small-amplitude adiabatic continuation toward
finite amplitude. It is a preliminary cosine-potential branch, not a complete
or converged reproduction of OTT Fig. 2. Quartic and sextic branches are absent.

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
`../preliminary_cosine_20260922.jls`. The shared driver resumes that checkpoint
if its own `branch/cosine_branch.jls` does not exist; otherwise it resumes its own
latest solution. Checkpoints are replaced atomically after every accepted point.
All successful points, frequencies, boundary amplitudes, and settings are retained.

Run from the repository root:

```bash
sbatch runs/adiabatic/preliminary_20260922/branch/submit.sh
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
The long recipe uses `branch/submit_long.sh`, sets `AXION_RUN_SECONDS=13800`, targets A0=13,
and refreshes the plot on normal completion. Failed Newton solves reduce
the amplitude increment; repeated failures can end continuation before the
target. A time-limit kill retains the last complete checkpoint.

To repeat the long run after all existing jobs have finished:

```bash
sbatch runs/adiabatic/preliminary_20260922/branch/submit_long.sh
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

Run `sbatch runs/adiabatic/preliminary_20260922/diagnostics/speed/submit_speed_check.sh` to
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
python runs/adiabatic/preliminary_20260922/branch/plot_branch.py
```

This writes PNG and PDF files to `plots/preliminary_cosine_branch_20260922.*`.
The plot is a snapshot; rerun this command after additional points complete.
Numerical data are in `branch/cosine_branch.csv`, and full checkpoints in
`branch/cosine_branch.jls`.

## Flux convention and checks

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
