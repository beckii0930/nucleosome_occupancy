#!/bin/bash
#SBATCH --ntasks-per-node=8
#SBATCH --time=30:00:00
#SBATCH --mem=150G
#SBATCH -N 2
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yibeijia@usc.edu
#SBATCH --partition=rohs
#python3 getInputFromMat.py
#cd-hit-est -i Eseqs.txt -o Eseqs80 -c 0.8
#cd-hit-est -i Dseqs.txt -o Dseqs80 -c 0.8

python3 cluster2seq.py
