#!/bin/bash
#SBATCH --job-name=axion_qcd_plots
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr7
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=runs/adiabatic/multiphase_two_cloud/logs/qcd-collect-%j.out
set -euo pipefail
[[ $# -eq 3 ]] || { echo 'usage: COLLECTOR_SCRIPT CAMPAIGN OUTPUT_DIRECTORY' >&2; exit 2; }
/global/home/users/clwelch/.conda/envs/BHsuper/bin/python \
    "$1" "$2" "$3"
