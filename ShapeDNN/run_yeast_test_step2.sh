#!/bin/bash
#SBATCH --ntasks-per-node=8
#SBATCH --time=30:00:00
#SBATCH --mem=150G
#SBATCH -N 1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yibeijia@usc.edu
#SBATCH --array=1
#SBATCH --output=yease_test.out
. ~/.bashrc
conda activate gpuenv
#python3  PreprocessShape_deplete.py yeast

#python3 -m cProfile cluster2shape.py

python3 -m cProfile cluster2shape_makeTest.py 5 ${SLURM_ARRAY_TASK_ID} yeast
#python3  nucPosPredict.py 10 $(for i in $(seq 3 10); do echo $i; done)
