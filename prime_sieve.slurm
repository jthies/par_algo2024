#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=3800MB
#SBATCH --time=00:01:00
#SBATCH --partition=compute-p1

module load 2024r1 openmpi opencoararys

git show HEAD |grep commit
git status -s

for N in 1000000 10000000 100000000; do
    srun ./prime_sieve.x $N
done
