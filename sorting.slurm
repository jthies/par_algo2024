#!/bin/bash
# account for participants of the MasterMath course "parallel algorithms"
#SBATCH --account=education-eemcs-courses-mastermath
# number of processes (will be spread arbitrarily over nodes)
#SBATCH --ntasks=4
# memory per process (note that the underlying MPI needs at least 500MB or so!)
#SBATCH --mem-per-cpu=1GB
# request up to two minutes
#SBATCH --time 00:02:00
# output will bw written to this file:
#SBATCH --output="sorting.out"
# load environment
module load 2024r1 openmpi opencoarrays

# run the code with "srun" to make sure it uses the allocated resources.
# Sort 40 million elements in total.(10M per process)
# We "pipe" the output to the grep command to get rid of
# pointless warnings about UCX tag-matching
srun ./main_sorting.x 10000000 | grep -v "unexpected tag"

# to investigate scaling, let's also run with 2 and 1 process:
srun -n 2 ./main_sorting.x 20000000 | grep -v "unexpected tag"
srun -n 1 ./main_sorting.x 40000000 | grep -v "unexpected tag"
