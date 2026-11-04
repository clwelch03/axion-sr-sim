#!/bin/bash
# Run after activating the recreated BHsuper environment and selecting Julia.
set -euo pipefail
repository=$(cd "$(dirname "$0")/../.." && pwd)
cd "$repository"
python_executable=${AXION_PYTHON:-python}
julia_executable=${AXION_JULIA:-julia}
"$python_executable" -c 'import sys, tomllib; assert sys.version_info >= (3,11), "Use Python >=3.11"'
"$julia_executable" --startup-file=no -e \
    'VERSION == v"1.13.0" || error("Select Julia 1.13.0; found $VERSION")'
export JULIA_DEPOT_PATH="${AXION_JULIA_DEPOT:-$HOME/.julia}:"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 JULIA_NUM_PRECOMPILE_TASKS=1
"$python_executable" -m pip install -e . --no-deps
"$python_executable" -c 'import numpy, scipy, mpmath, matplotlib, pandas'
for project in . runs/adiabatic/radial_bvp; do
    "$julia_executable" --startup-file=no --project="$project" -e \
        'using Pkg; Pkg.instantiate(); Pkg.precompile()'
done
mkdir -p runs/adiabatic/multiphase_two_cloud/logs
echo 'Environment ready. Submit scripts/hpc/smoke_perlmutter.jl with perlmutter_job.sh.'
