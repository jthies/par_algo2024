#!/bin/bash
#SBATCH --output="benchmarks.out"
# run on one (phase-1) node
#SBATCH -N 1
# using all 48 CPU cores
#SBATCH -n 48
# maximum runtime (estimate)
#SBATCH -t 00:30:00
# make sure no other jobs are running on the node
# to avoid resource contention
#SBATCH --exclusive
#SBATCH --partition="compute"
# request all RAM on the node (as we request it exclusively)
#SBATCH --mem=0
# make all cores use the same clock speed so it looks more like a BSP computer
#SBATCH --cpu-freq=medium
#SBATCH --account=education-eemcs-courses-paf24

# the driver writes matlab/octave executable output,
# we pipe all of it into an m-file for post-processing
echo "" > bench_data.m

# Force process s to run on a specific core (again, to make it
# look more like a BSP computer)
NUMA_DOMAINS="0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,3,3"

# load necessary modules for running Fortran coarray code
module load 2024rc1 openmpi opencoarrays
for np in 12 24 48; do
    srun  -n $np --cpu-bind=map_ldom:${NUMA_DOMAINS} ./main_benchmarks.x |& grep -v UCX >> bench_data.m
done
