#!/bin/bash
#SBATCH --ntasks-per-node=8
#SBATCH --time=30:00:00
#SBATCH --mem=150G
#SBATCH -N 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yibeijia@usc.edu
#SBATCH --output=scaleInput.out

. ~/.bashrc
conda activate gpuenv

module load usc
python3 ScaleInputData.py Test

#declare -a TrainOrTest=('Train' 'Test')

#for i in ${TrainOrTest[@]}; do
#   echo $i
#   python3 ScaleInputData.py $i &
#done
