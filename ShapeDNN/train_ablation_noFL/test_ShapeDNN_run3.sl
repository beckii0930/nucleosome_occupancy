#!/bin/bash
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --time=24:00:00
#SBATCH --mem=64GB
#SBATCH --gres=gpu:2
#SBATCH --partition=gpu
#SBATCH --mail-user=yibeijia@usc.edu

#SBATCH --output=test_ShapeDNN_run3.out

. ~/.bashrc
conda activate gpuenv2
module load usc
module load cuda/10.1.243
folder=/project/rohs_108/yibeijia/nucleosome_occupancy/ShapeDNN/train_ablation_noFL/train_run3/
mkdir $folder
python3 testShapeDNN.py -i $folder
