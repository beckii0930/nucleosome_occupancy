#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:2
#SBATCH --time=1:00:00
#SBATCH --mem=64GB
#SBATCH --output=test.out
#SBATCH --partition=rohs

module load cuda/10.1.243

# python3 test.py

for f in *.hdf5
do
    python3 -i $f test.py
done

