#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=64GB
#SBATCH --gres=gpu:2
#SBATCH --partition=gpu
#SBATCH --mail-user=yibeijia@usc.edu

#SBATCH --output=train_ShapeDNN_run6.out

. ~/.bashrc
conda activate gpuenv2
module load usc
module load cuda/10.1.243
#python3 trainShapeDNN.py
# dense 50 with val_accuracy as monitor 
folder=/project/rohs_108/yibeijia/nucleosome_occupancy/ShapeDNN/train_debug/train_run6/

mkdir $folder
python3 trainShapeDNN.py $folder
