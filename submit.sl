#!/bin/bash -l
#SBATCH --job-name=results/job
#SBATCH --output=%x-%j.out
#SBATCH --account=def-edgrant
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=512G


module load matlab
srun matlab -nodisplay -r "start_sim_soc" 