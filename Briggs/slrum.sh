#!/bin/bash
#SBATCH --job-name=briggs
#SBATCH --output=textout_%j.txt
#SBATCH --error=texterr_%j.txt
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem-per-cpu=4G

module load julia/1.12.1
cd /work/home/ng66sume/Eigenvalue/Briggs || exit 1

echo "Job started on $(date)"
echo "Running in $(pwd)"

stdbuf -oL -eL julia briggs.jl

echo "Job ended on $(date)"