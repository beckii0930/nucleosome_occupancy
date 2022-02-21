#!/bin/bash

#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --gres=gpu:2
#SBATCH --time=24:00:00
#SBATCH --mem=64GB
#SBATCH --output=train_dense10_kernel10_run2.out
#--partition=rohs
#SBATCH --partition=gpu
#SBATCH --mail-user=yibeijia@usc.edu

. ~/.bashrc
conda activate gpuenv2

module load usc
module load cuda/10.1.243

folder=/project/rohs_108/yibeijia/nucleosome_occupancy/tbinet_stuff/train_dense10_kernel10_run2/
mkdir $folder
python3 train.py $folder
#python3 train_model3_nonneg.py
#python3 train_model3_removeLSTM.py
