#!/bin/bash
# Run the regression check and then the timing script in one Slurm job.
#
#   benchmarks/submit.sh <cpus-per-task> <julia-threads> [blas-threads=1] [walltime=02:00:00] [extra timing.jl arguments...]
#
# Logs go to benchmarks/logs/bench-<jobid>.out. Run from any directory; the script
# locates the repository from its own path. The job's exit status is that of the
# regression check (the timing script runs afterwards in any case). To benchmark
# a branch, check it out (or use its worktree) and run this script from there.
set -u
benchmarks_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(dirname "$benchmarks_directory")"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    cpus="${1:?usage: submit.sh <cpus-per-task> <julia-threads> [blas-threads] [walltime] [timing.jl args...]}"
    julia_threads="${2:?usage: submit.sh <cpus-per-task> <julia-threads> [blas-threads] [walltime] [timing.jl args...]}"
    blas_threads="${3:-1}"
    walltime="${4:-02:00:00}"
    shift $(( $# < 4 ? $# : 4 ))
    mkdir -p "$benchmarks_directory/logs"
    exec sbatch --account=pc_heptheory --partition=lr7 --qos=lr_normal --ntasks=1 --cpus-per-task="$cpus" \
        --time="$walltime" --job-name=axion-bench \
        --output="$benchmarks_directory/logs/bench-%j.out" \
        --export=ALL,BENCH_REPOSITORY="$repository",BENCH_JULIA_THREADS="$julia_threads",BENCH_BLAS_THREADS="$blas_threads" \
        "${BASH_SOURCE[0]}" "$@"
fi

# --- inside the job ---
# Slurm runs a copy of this script, so the repository path comes from the submit step.
cd "$BENCH_REPOSITORY"
export OPENBLAS_NUM_THREADS="${BENCH_BLAS_THREADS:-1}"
julia="${JULIA:-$HOME/.juliaup/bin/julia}"
echo "job $SLURM_JOB_ID on $(hostname), commit $(git rev-parse --short HEAD), julia threads $BENCH_JULIA_THREADS, BLAS threads $OPENBLAS_NUM_THREADS"

"$julia" --project=. -t "$BENCH_JULIA_THREADS" benchmarks/regression.jl
regression_status=$?
echo "regression exit status $regression_status"

"$julia" --project=. -t "$BENCH_JULIA_THREADS" benchmarks/timing.jl "$@"
echo "timing exit status $?"
exit $regression_status
