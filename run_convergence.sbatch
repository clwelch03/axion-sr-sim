#!/bin/bash
#SBATCH --job-name=r1d_conv
#SBATCH --account=pc_heptheory
#SBATCH --partition=lr6
#SBATCH --qos=lr_normal
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --output=convergence_%j.out
#SBATCH --error=convergence_%j.err
source /global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/env/lawrencium.sh
cd /global/scratch/projects/pc_heptheory/clwelch/software/GRChombo-Axion-SR/Examples/AxionSuperradiance/radial1d
export OMP_NUM_THREADS=1
python3 convergence_test.py
