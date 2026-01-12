#!/bin/bash
#SBATCH --job-name=SCION-test
#SBATCH --time=48:00:00
#SBATCH --mem=80G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=XX@leeds.ac.uk

module load matlab/R2023a
matlab -nodisplay -r SCION_sens