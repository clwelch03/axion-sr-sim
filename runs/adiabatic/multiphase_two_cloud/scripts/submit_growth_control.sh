#!/bin/bash
#SBATCH --job-name=axion_growth_control
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=08:00:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/growth-control-%j.out
set -euo pipefail
[[ $# -eq 3 ]] || { echo 'usage: SNAPSHOT chiral|cosine|omiya OUTPUT_PREFIX' >&2; exit 2; }
snapshot_directory=$(realpath "$1")
control_model=$2
output_prefix=$(realpath -m "$3")
[[ "$control_model" == chiral || "$control_model" == cosine || "$control_model" == omiya ]] || exit 2
export JULIA_DEPOT_PATH="$SLURM_SUBMIT_DIR/runs/adiabatic/radial_bvp/work/depot:/global/home/users/clwelch/.julia"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/axion-growth-${SLURM_JOB_ID}"
cd "$snapshot_directory"
python3 runs/adiabatic/multiphase_two_cloud/scripts/verify_snapshot.py .
evolution_status=0
/global/home/users/clwelch/.julia/juliaup/julia-1.13.0+0.x64.linux.gnu/bin/julia \
    --startup-file=no --compiled-modules=existing --project=. \
    "runs/adiabatic/multiphase_two_cloud/scripts/${control_model}_growth.jl" \
    inputs/initial.jls inputs/reference.toml "$output_prefix" || evolution_status=$?
if [[ "$control_model" != omiya && -f "${output_prefix}-steps1000.csv" && -f "${output_prefix}-linear-reference.csv" ]]; then
    PYTHONPATH=src /global/home/users/clwelch/.conda/envs/BHsuper/bin/python \
        runs/adiabatic/multiphase_two_cloud/scripts/plot_chiral_growth.py \
        "$output_prefix" --potential "$control_model"
fi
exit "$evolution_status"
