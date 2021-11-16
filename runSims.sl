#!/bin/bash -l
#SBATCH --job-name=smallSim
#SBATCH --account=def-edgrant
#SBATCH --time=00:05:00
#SBATCH --mem=100000

module load matlab
srun matlab -nodisplay -singleCompThread -r "start_sim" 