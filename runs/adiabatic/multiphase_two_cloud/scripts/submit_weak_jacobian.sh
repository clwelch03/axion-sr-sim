#!/bin/bash
#SBATCH --job-name=axion_weak_jacobian
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/weak-jacobian-%j.out
set -euo pipefail
[[ $# -eq 2 ]] || { echo "usage: submit_weak_jacobian.sh SNAPSHOT OUTPUT_JLS" >&2; exit 2; }
snapshot_directory=$1
output_path=$2
[[ ! -e "$output_path" ]] || { echo "Output already exists" >&2; exit 2; }
export OPENBLAS_NUM_THREADS=1 JULIA_NUM_THREADS=1
cd "$snapshot_directory"
python3 runs/adiabatic/multiphase_two_cloud/scripts/verify_snapshot.py .
# Reuse exactly the tested setup and derivative diagnostic, without Newton.
/global/home/users/clwelch/.juliaup/bin/julia --startup-file=no --project=. -e '
include("runs/adiabatic/multiphase_two_cloud/scripts/weak_secondary_profile.jl")
setup = weak_secondary_setup(deserialize("inputs/initial.jls"))
println("Checking the original weak-secondary guess; no Newton solve is run here.")
admissible = state -> two_phase_profile_admissible(state, setup.problem)
check_matching_jacobian(setup.residual, setup.state, admissible, ARGS[1])
' "$output_path"
