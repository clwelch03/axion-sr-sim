#!/bin/bash
#SBATCH --job-name=axion_quartic_validation
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --output=runs/adiabatic/fig2_20260922/logs/quartic-validation-%j.out
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
export JULIA_NUM_THREADS=1
/global/home/users/clwelch/.juliaup/bin/julia --project=. runs/adiabatic/fig2_20260922/diagnostics/high_amplitude/validate_quartic_winner.jl
