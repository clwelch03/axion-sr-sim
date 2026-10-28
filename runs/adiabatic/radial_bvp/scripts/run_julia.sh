#!/bin/bash
# Keep dependency installation and compiler caches inside this experiment.
set -euo pipefail
experiment_directory=$(cd "$(dirname "$0")/.." && pwd)
export JULIA_DEPOT_PATH="$experiment_directory/work/depot:/global/home/users/clwelch/.julia"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1
exec /global/home/users/clwelch/.juliaup/bin/julia --startup-file=no \
    --project="$experiment_directory" "$@"
