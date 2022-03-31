#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=64GB
#SBATCH --gres=gpu:2
#SBATCH --partition=gpu
#SBATCH --mail-user=yibeijia@usc.edu

#SBATCH --output=train_dense10_kernel9_run4.out

. ~/.bashrc
conda activate gpuenv2
module load usc
module load cuda/10.1.243
folder=/project/rohs_108/yibeijia/nucleosome_occupancy/tbinet_stuff/train_dense10_kernel9_run4/

mkdir $folder
python3 train.py $folder
