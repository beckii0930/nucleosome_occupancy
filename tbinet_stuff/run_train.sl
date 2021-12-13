#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:2
#SBATCH --time=10:00:00
#SBATCH --mem=64GB
#SBATCH --output=train_model2.out
#SBATCH --partition=rohs
# --partition=gpu
#SBATCH --mail-user=yibeijia@usc.edu

. ~/.bashrc
conda activate gpuenv2

module load usc
module load cuda/10.1.243

python3 train.py
