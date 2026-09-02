# Adiabatic Kerr backreaction

This directory contains the first Tier-1 backreaction baseline: one scalar
quasi-bound level grows while the black hole moves through a sequence of Kerr
backgrounds.  The model is deliberately lossless and single-level.  It includes
the horizon-flux relation but not level transitions, gravitational radiation,
axion radiation, nonlinear detuning, or dynamical metric evolution.

## Endpoint calculation

For a mode with azimuthal number `m`, each incremental transfer obeys

```text
dM_BH / dN = -omega_R,
dJ_BH / dN = -m,
```

or, equivalently,

```text
dJ_BH / dM_BH = m / omega_R.
```

The integration stops at the superradiant threshold

```text
omega_R = m * Omega_H.
```

The evolving-frequency calculation keeps the boson mass `mu` fixed, updates
`alpha = M_BH * mu` and `chi = J_BH / M_BH**2`, and solves the Leaver problem at
each point on the trajectory.  Cloud energy and angular momentum are the
cumulative losses from the black hole.

## First result: spectral-feedback control

Parameters: initial `M_BH=1`, `chi=0.99`, `alpha=0.42`, and `(n,l,m)=(2,1,1)`.
The initial cloud fraction is `1e-12`; it affects the bookkeeping offset but not
the lossless endpoint at the displayed precision.

| frequency model | final `M_BH` | final `chi` | `M_cloud/M_BH` |
|---|---:|---:|---:|
| fixed `omega_R=0.4087971` control | 0.96869435 | 0.97341273 | 0.032317 |
| evolving Leaver `omega_R(M_BH,J_BH)` | 0.96965450 | 0.97406098 | 0.031295 |

The fixed-frequency row is an exact-answer control for the ODE and event
detection.  Its analytic trajectory is

```text
J_BH(M_BH) = J_BH,0 + (m / omega_R,0) * (M_BH - M_BH,0).
```

The numerical control reproduces the analytic endpoint to machine precision.
Allowing the Leaver frequency to evolve reduces the extracted cloud energy by
about `3.1%` relative to that control.  This isolates spectral feedback from
the omitted transition and radiation channels.

## Reproduction

Load the Lawrencium environment, then run:

```bash
source env/lawrencium.sh
GWENV=/global/scratch/projects/pc_heptheory/clwelch/software/envs/gwenv/bin/python
$GWENV Examples/AxionSuperradiance/radial1d/backreaction/adiabatic_endpoint.py \
  --validate-fixed-frequency
$GWENV Examples/AxionSuperradiance/radial1d/backreaction/adiabatic_endpoint.py
```

The first command checks the exact constant-frequency control.  The second
performs the evolving-Leaver endpoint calculation and checks threshold,
energy/angular-momentum conservation, and nondecreasing Kerr horizon area.
