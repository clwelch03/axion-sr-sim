# Preliminary cosine continuation — 2026-09-22

This run extends the existing small-amplitude adiabatic continuation toward
finite amplitude. It is a preliminary cosine-potential branch, not a complete
or converged reproduction of OTT Fig. 2. Quartic and sextic branches are absent.

Settings: spin=0.99, alpha=0.42; modes [(1,1),(3,1),(5,1),(3,3),(5,3),(5,5)];
16 polar and 31 azimuthal samples; radial RK4 step=0.1; rstar endpoints -100
and 150; matching coordinate 5; scaled Newton residual tolerance 1e-9.
Energy uses Simpson integration over the same finite domain. The sampled
peak field uses every tenth radial sample and is not a precise extremum.

The initial five points were produced before a login-node change, in
`../preliminary_cosine_20260922.jls`. The shared driver resumes that checkpoint
if its own `cosine_branch.jls` does not exist; otherwise it resumes its own
latest solution. Checkpoints are replaced atomically after every accepted point.
All successful points, frequencies, boundary amplitudes, and settings are retained.

Run from the repository root:

```bash
sbatch runs/adiabatic/preliminary_20260922/submit.sh
```

Do not start a second job using the same output while the first is running.
The driver attempts further continuation for about 1000 seconds (checked
between solves); SLURM imposes a 20-minute wall-time limit. An unfinished
solve may be discarded at that limit, but completed points remain available.
The original resumed job is 26265137.

A four-hour extension is queued as job 26265520 with
`--dependency=afterany:26265137`, so it starts only after the short job ends.
It uses `submit_long.sh`, sets `AXION_RUN_SECONDS=13800`, targets A0=13,
and refreshes the plot on normal completion. Failed Newton solves reduce
the amplitude increment; repeated failures can end continuation before the
target. A time-limit kill retains the last complete checkpoint.

To repeat the long run after all existing jobs have finished:

```bash
sbatch runs/adiabatic/preliminary_20260922/submit_long.sh
```

The initial delivered plot contains nine points through A0=0.12301875,
E/Fa²=216.4900641, sampled peak field=0.2890454. At that point E/A0² is
about 32.5% below the small-amplitude value. Later checkpoints may extend
this range; the CSV and SLURM logs are authoritative.

After loading the BHsuper environment, redraw the latest saved points:

```bash
python runs/adiabatic/preliminary_20260922/plot_branch.py
```

This writes PNG and PDF files to `plots/preliminary_cosine_branch_20260922.*`.
The plot is a snapshot; rerun this command after additional points complete.
Numerical data are in `cosine_branch.csv`, and full checkpoints in
`cosine_branch.jls`.

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
