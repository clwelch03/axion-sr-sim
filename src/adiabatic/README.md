# Adiabatic cloud models

Choose the physical approximation first. Each model keeps its profile solves,
observables, continuation, and slow evolution together.

| Folder | What it implements | Julia entry point |
|---|---|---|
| [core](core/README.md) | Shared Kerr geometry, angular bases, potentials, explicit-frequency operators, boundary data, and numerical methods | `core/load.jl` |
| [single_cloud](single_cloud/README.md) | One nonlinear cloud, including harmonics of its single fundamental phase | `single_cloud/load.jl` |
| [primary_secondary](primary_secondary/README.md) | A nonlinear primary, a secondary on its averaged background, and driven radiation; no secondary deformation of the primary | `primary_secondary/load.jl` |
| [multiphase](multiphase/README.md) | Reserved for mutually deforming independent clouds; no multiphase solver yet | `multiphase/load.jl` currently loads only core |

Run Julia from the repository root with `julia --project=.`. The existing
`src/adiabatic/load.jl` entry point loads both implemented models. For a narrower
dependency, include the model's loader directly. Loaders are idempotent within
the calling Julia module; they define functions there, without introducing a
new package or changing serialized solution types.

Dependencies are `core -> single_cloud -> primary_secondary`. Multiphase starts
from `core`; future seed/comparison drivers can explicitly load the other models.
Core code does not select cloud populations or depend on a model. Command-line
drivers live inside the model they run; reusable functions live in the model's
library and are available through its loader.

## Commands

```bash
# Complete Julia unit suite, or one group independently.
julia --project=. tests/adiabatic/runtests.jl
julia --project=. tests/adiabatic/core/runtests.jl
julia --project=. tests/adiabatic/single_cloud/runtests.jl
julia --project=. tests/adiabatic/primary_secondary/runtests.jl

# Python balance and adapter tests.
PYTHONPATH=src python -m unittest discover -s tests/adiabatic

# Continue one nonlinear cloud (explicit seed and output required).
AXION_TARGET_AMPLITUDE=13 julia --project=. \
  src/adiabatic/single_cloud/drivers/continue_branch.jl cosine seed.jls output_prefix
```

Use [campaign READMEs](../../runs/README.md) for actual inputs, outputs, restart
rules, and Slurm recipes. Longer `check_*.jl` files live under the corresponding
model in `tests/adiabatic/`; they are separate from the unit suite.

## October 12 layout migration

The old global `stationary/`, `observables/`, `amplitude_evolution/`,
`multimode/`, and `drivers/` source folders have been replaced. Command and import
paths changed; the root Julia loader and numerical function names remain usable.

- Single-cloud drivers: `src/adiabatic/single_cloud/drivers/`.
- Primary/secondary drivers: `src/adiabatic/primary_secondary/drivers/`.
- Python balance import: `adiabatic.primary_secondary.qcd_spindown`.
- Fig. 13 adapter import: `adiabatic.single_cloud.fig13_points`.
- Superseded cosine/polynomial drivers: the Figure 2 campaign's
  `scripts/legacy_continue_cosine.jl` and `legacy_continue_polynomial.jl`.

Campaign data and saved checkpoints retain their locations. Frozen source
manifests intentionally reject the new layout until a separately validated
migration or new campaign records its sources; do not bypass their hashes.
See the [implementation plan](../../docs/multiphase-implementation-plan.md)
for the next physics milestones and [HANDOFF.md](../../HANDOFF.md) for results.
