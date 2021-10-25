#!/bin/bash
#SBATCH --job-name=horizontalDisaster
#SBATCH -p physical
#SBATCH --constraint=physg5
#SBATCH --time=48:00:00
#SBATCH --mem=6G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1040


# Load the required modules
module load gcc/8.3.0
module load cmake/3.18.4 
module load cgal/4.14.1-python-3.7.4
module load gmp/6.1.2
module load mpfr/4.0.2
# Move into folder and run
cd ../cpp_code

# Call job. Note first bool specifies print info, second whether graph should be plotted

# This one has 1040 jobs
./augmentation run preliminaryRun ${SLURM_ARRAY_TASK_ID} 1 0

