#!/bin/bash
#SBATCH --job-name=axion_qcd_threeway
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=02:00:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/qcd-threeway-%A_%a.out
set -euo pipefail
[[ $# -eq 2 ]] || { echo 'usage: SNAPSHOT OUTPUT_DIRECTORY' >&2; exit 2; }
snapshot_directory=$(realpath "$1")
output_directory=$(realpath -m "$2")
export JULIA_DEPOT_PATH="$SLURM_SUBMIT_DIR/runs/adiabatic/radial_bvp/work/depot:/global/home/users/clwelch/.julia"
export JULIA_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/axion-qcd-mpl-${SLURM_JOB_ID}"
cd "$snapshot_directory"
sha256sum --quiet -c source-manifest.sha256
/global/home/users/clwelch/.julia/juliaup/julia-1.13.0+0.x64.linux.gnu/bin/julia \
    --startup-file=no --compiled-modules=existing --project=. \
    runs/adiabatic/multiphase_two_cloud/scripts/qcd_211_comparison.jl \
    "$snapshot_directory/inputs/cases.toml" "${SLURM_ARRAY_TASK_ID:?Submit an array}" "$output_directory"
/global/home/users/clwelch/.conda/envs/BHsuper/bin/python - \
    "$snapshot_directory/inputs/cases.toml" "$SLURM_ARRAY_TASK_ID" "$output_directory" <<'PY'
import subprocess,sys,tomllib
from pathlib import Path
config=tomllib.loads(Path(sys.argv[1]).read_text())
case=config['cases'][int(sys.argv[2])-1]
prefix=Path(sys.argv[3])/case['name']/'comparison'
subprocess.run([sys.executable,'runs/adiabatic/multiphase_two_cloud/scripts/plot_qcd_both_clouds.py',str(prefix)],check=True)
PY
