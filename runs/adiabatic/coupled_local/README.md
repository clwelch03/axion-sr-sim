# Coupled evolution: paused draft

On October 8, 2026, development was stopped at the user's request to rebuild
the implementation in small, explained pieces under `src/adiabatic/`.
Stationary-table Slurm job **26757461 was cancelled**. Any partial stationary
output is exploratory; no coupled trajectory has been run or validated.
The scripts here preserve the unfinished draft, not a reproduction recipe.
Their source hashes have changed since job submission, so their existing
resume checks intentionally reject any older checkpoint.

The first piece is the instantaneous balance evaluator in
[`coupled_balance.jl`](../../../src/adiabatic/multimode/coupled_balance.jl).
It takes the current cloud response and boundary fluxes and returns BH/cloud
rates in fixed initial-BH units. The existing
[`check_coupled_balance.jl`](../background_continuation/scripts/check_coupled_balance.jl)
driver checks that piece against saved stationary data and analytic balances.
After the move, all **44 checks passed**; verification output was written to a
temporary directory, preserving the earlier saved balance checkpoint.

Review that piece before adding the state-to-stationary-data conversion and
then a single time step. Interpolation, trajectory output, checkpointing, and
new stationary runs remain deferred.
