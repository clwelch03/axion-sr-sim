#!/bin/bash
#SBATCH -A pc_heptheory
#SBATCH -p lr7
#SBATCH -q lr_normal
#SBATCH -n 1
#SBATCH -t 06:00:00
# Usage: sbatch submit_witte_two_level.sh <f_a/M_Pl> <end time in M0> <tag for output name>
WITTE=/global/scratch/projects/pc_heptheory/clwelch/witte-axion-sr
HERE=/global/scratch/projects/pc_heptheory/clwelch/axion-sr-sim/runs/witte_two_level
cd $HERE || exit 1
julia --project=$WITTE/env $HERE/scripts/run_witte_two_level.jl $WITTE/repo/src $1 $2 $HERE/scripts/two_channel_list.txt $HERE/results/witte_two_level_fa$1_t$3_raw.csv
