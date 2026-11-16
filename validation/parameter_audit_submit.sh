#!/bin/bash
# Slurm wrapper for validation/parameter_audit.jl: solves and times the listed studies at 8 Julia threads, then
# retimes them at 1 thread, in one lr7 job.
#
#   validation/parameter_audit_submit.sh <output-dir> <walltime> <study> [study ...]   (study "tolerance" allowed)
set -u
validation_directory="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repository="$(dirname "$validation_directory")"

if [[ -z "${SLURM_JOB_ID:-}" ]]; then
    output_directory="$(mkdir -p "$1" && cd "$1" && pwd)"
    walltime="$2"
    shift 2
    mkdir -p "$output_directory/logs"
    exec sbatch --account=pc_heptheory --partition=lr7 --qos=lr_normal --ntasks=1 --cpus-per-task=8 --time="$walltime" \
        --job-name=axion-audit --output="$output_directory/logs/audit-%j.out" \
        --export=ALL,AUDIT_REPOSITORY="$repository",AUDIT_OUTPUT="$output_directory" "${BASH_SOURCE[0]}" "$@"
fi

cd "$AUDIT_REPOSITORY"
export OPENBLAS_NUM_THREADS=1
julia="${JULIA:-$HOME/.juliaup/bin/julia}"
echo "job $SLURM_JOB_ID on $(hostname) ($(lscpu | grep 'Model name' | cut -d: -f2 | xargs)), commit $(git rev-parse --short HEAD), studies: $*"
for study in "$@"; do
    if [[ "$study" == tolerance ]]; then
        "$julia" --project=. -t 8 validation/parameter_audit.jl tolerance "$AUDIT_OUTPUT"
    else
        "$julia" --project=. -t 8 validation/parameter_audit.jl solve "$study" "$AUDIT_OUTPUT"
    fi
done
for study in "$@"; do
    if [[ "$study" == tolerance ]]; then
        "$julia" --project=. -t 1 validation/parameter_audit.jl tolerance "$AUDIT_OUTPUT"
    else
        "$julia" --project=. -t 1 validation/parameter_audit.jl time "$study" "$AUDIT_OUTPUT"
    fi
done
echo "done"
