# Small-amplitude adiabatic checks

## Start here

This is the preserved small-amplitude validation campaign. Its seven-point
`results/continuation.{csv,jls}` branch and `results/profiles.csv` are reference
inputs for checking the stationary solver; `results/nonlinear_convergence.csv`
records five convergence samples. Finite-amplitude production inputs are in
the [Figure 2 campaign](../fig2_20260922/README.md).

Redraw the saved profiles from the repository root:

```bash
python src/adiabatic/plot_profiles.py
```

This writes `plots/adiabatic_profiles.png`. These checks establish behavior in
the small-amplitude regime; they do not validate high-amplitude branches or
two-cloud evolution. See the [validation record](../../../docs/adiabatic-repro.md)
for interpretation and [HANDOFF.md](../../../HANDOFF.md) for project status.

## Regenerate the reference calculations

`results/continuation.{csv,jls}` contains the initial saved continuation;
`results/profiles.csv` contains reconstructed radial profiles;
`results/nonlinear_convergence.csv` records the nonlinear convergence check.
These are preserved reference results, separate from the finite-amplitude
[Figure 2 campaign](../fig2_20260922/README.md).

From the repository root:

```bash
julia --project=. src/adiabatic/drivers/continuation.jl
julia --project=. src/adiabatic/drivers/inspect_profile_agreement.jl
python src/adiabatic/plot_profiles.py
julia --project=. tests/adiabatic/check_nonlinear_convergence.jl
julia --project=. tests/adiabatic/check_energy_convergence.jl
```

The first two commands replace the corresponding saved results; plotting only
reads profiles and writes `plots/adiabatic_profiles.png`. Run costly calculations
through Slurm. See [the reproduction notes](../../../docs/adiabatic-repro.md)
for settings and interpretation.
