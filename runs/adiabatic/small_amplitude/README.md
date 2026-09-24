# Small-amplitude adiabatic checks

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
