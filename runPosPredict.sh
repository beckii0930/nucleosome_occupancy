#!/bin/bash
#SBATCH --ntasks-per-node=8
#SBATCH --time=30:00:00
#SBATCH --mem=150G
#SBATCH -N 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yibeijia@usc.edu
#SBATCH --array=1-50 # job array index

# python3  nucPosPredict.py 50 ${SLURM_ARRAY_TASK_ID}
python3  nucPosPredict2.py 50 ${SLURM_ARRAY_TASK_ID}
#python3  nucPosPredict.py 10 $(for i in $(seq 3 10); do echo $i; done)
