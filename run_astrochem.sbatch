#!/bin/bash

# Comments can go here with hashtags
#SBATCH --time=30:00:00 # estimated wall time
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=25
#SBATCH --mem-per-cpu=1000
#SBATCH --job-name=r030_iso
#SBATCH --output=%x-%j.out
#SBATCH --error=%x-%j.err
#SBATCH --partition=broadwl
#SBATCH --mail-user=ericvc@uchicago.edu
#SBATCH --mail-type=ALL

# load modules needed here
module load python

# Run code here
source activate candy
mpiexec python run_growth_parallel.py
