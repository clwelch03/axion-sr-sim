# Nonlinear primary with a secondary on its background

This is the existing Omiya et al. Sec. III approximation: the primary deforms
nonlinearly, the secondary responds linearly to its averaged background, and
driven fields carry energy and angular momentum. The secondary does not reshape
the primary. Two mutually deforming clouds will belong to `multiphase/`.

Load `src/adiabatic/primary_secondary/load.jl`; it loads `single_cloud/` and
`core/` automatically. It defines library functions without running solves or
writing campaign output.

| Files | Responsibility |
|---|---|
| `background.jl`, `secondary_mode.jl`, `secondary_energy.jl` | Averaged primary background, secondary eigenproblem, energy, and fluxes |
| `dissipative_mode.jl`, `radiative_mode.jl` | Driven m=0 and m=3 responses |
| `evolution.jl`, `coupled_balance.jl` | Fixed-BH cloud rates and leading coupled BH/cloud balances |
| `scalar_energy_accounting.jl` | Reconstruct boundary fluxes and torques from saved trajectories |
| `continuation.jl` | Read/reconstruct primary checkpoints and continue secondary solutions |
| `secondary_tables.jl`, `dissipative_tables.jl`, `radiative_tables.jl` | Build observable tables and evaluate source-work/flux diagnostics |
| `stationary_point.jl` | Compose primary, secondary, and driven solves for the QCD campaign |
| `reference_tables.jl` | Read the saved reference tables and reproduce Fig. 6 evolution |
| `qcd_spindown.py` | Python interpolation and leading BH/cloud evolution |
| `drivers/` | Argument handling and explicit campaign execution |

Use `from adiabatic.primary_secondary import qcd_spindown` in Python. Scripts
that configure its module-level physical constants must still run each case in
its own process; this refactor does not change that convention.

```bash
julia --project=. tests/adiabatic/primary_secondary/runtests.jl
julia --project=. tests/adiabatic/primary_secondary/check_two_cloud_rates.jl
PYTHONPATH=src python -m unittest discover -s tests/adiabatic/primary_secondary
```

See the [model design](../../../docs/multimode-adiabatic-design.md),
[secondary campaign](../../../runs/adiabatic/multimode_secondary/README.md), and
[QCD campaign](../../../runs/adiabatic/qcd_spindown/README.md).
