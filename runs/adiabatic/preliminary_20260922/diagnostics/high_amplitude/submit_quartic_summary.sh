#!/bin/bash
# Submit with --dependency=afterany:<method-array-job-id>.
#SBATCH --job-name=axion_quartic_summary
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:02:00
#SBATCH --output=runs/adiabatic/preliminary_20260922/logs/quartic-summary-%j.out
set -euo pipefail
/global/home/users/clwelch/.conda/envs/BHsuper/bin/python runs/adiabatic/preliminary_20260922/diagnostics/high_amplitude/summarize_quartic_methods.py
