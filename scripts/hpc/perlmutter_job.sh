#!/bin/bash -l
#SBATCH --job-name=axion_perlmutter
#SBATCH --constraint=cpu
#SBATCH --qos=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/perlmutter-%j.out
# Supply --account=YOUR_NERSC_PROJECT to sbatch.
set -euo pipefail
[[ $# -ge 2 ]] || { echo 'usage: SNAPSHOT_OR_CHECKOUT RELATIVE_DRIVER [DRIVER_ARGUMENTS...]' >&2; exit 2; }
snapshot_directory=$(realpath "$1")
driver=$2
shift 2
driver_arguments=()
for argument in "$@"; do
    if [[ $argument == '$SLURM_ARRAY_TASK_ID' ]]; then
        driver_arguments+=("${SLURM_ARRAY_TASK_ID:?Submit this driver with --array}")
    else
        driver_arguments+=("$argument")
    fi
done
[[ $driver != /* && -f "$snapshot_directory/$driver" ]] || {
    echo 'Driver must be a file relative to the snapshot or checkout.' >&2; exit 2;
}
export JULIA_DEPOT_PATH="${AXION_JULIA_DEPOT:-$HOME/.julia}:"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
export JULIA_NUM_PRECOMPILE_TASKS=1 MPLBACKEND=Agg
cd "$snapshot_directory"
if [[ -f source-manifest.sha256 ]]; then
    sha256sum --quiet -c source-manifest.sha256
fi
case "$driver" in
    *.jl)
        julia_executable=${AXION_JULIA:-julia}
        "$julia_executable" --startup-file=no -e \
            'VERSION == v"1.13.0" || error("Use Julia 1.13.0 for these checkpoints; found $VERSION")'
        exec "$julia_executable" --startup-file=no --compiled-modules=existing \
            --project=. "$driver" "${driver_arguments[@]}"
        ;;
    *.py)
        exec "${AXION_PYTHON:-python}" "$driver" "${driver_arguments[@]}"
        ;;
    *) echo 'Supported drivers end in .jl or .py.' >&2; exit 2 ;;
esac
