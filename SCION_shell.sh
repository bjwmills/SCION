# Run with current environment
#$ -V

# Set a time limit
#$ -l h_rt=48:00:00

#Request some memory per core
#$ -l h_vmem=1G

# Ask for lots of smp cores (insert $ after # to comment this in)
#$ -pe smp 40

#Get email at start and end of the job
#$ -m be

# Load matlab module
module add matlab

# run matlab using command file
# -nodisplay flag should be given to suppress graphics
matlab -nodisplay < SCION_sens.m