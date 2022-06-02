#!/bin/bash
#SBATCH --ntasks-per-node=8
#SBATCH --time=30:00:00
#SBATCH --mem=150G
#SBATCH -N 8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=yibeijia@usc.edu
#SBATCH --array=1
#SBATCH --output=yeast_deplete_therest.out
. ~/.bashrc
conda activate gpuenv

declare -a All_Shapes=('Buckle-FL' 'Buckle' 'EP' 'HelT-FL' 'HelT' 'MGW-FL' 'MGW'
			 'Opening-FL' 'Opening' 'ProT-FL' 'ProT' 'Rise-FL' 'Rise' 'Roll-FL'
			 'Roll' 'Shear-FL' 'Shear' 'Shift-FL' 'Shift' 'Slide-FL' 'Slide'
			 'Stagger-FL' 'Stagger' 'Stretch-FL' 'Stretch' 'Tilt-FL' 'Tilt')

declare -a All_Shapes=('Buckle' 'EP' 'HelT-FL' 'HelT' 'MGW-FL' 'MGW')

# declare -a All_Shapes=('Buckle-FL' 'Rise' 'Roll-FL')
# declare -a All_Shapes=('Buckle-FL')
# declare -a All_Shapes=('Buckle-FL' 'Buckle' 'EP' 'HelT-FL' 'HelT' 'MGW-FL' 'MGW'
# 			 'Opening-FL' 'Opening' )

# for shape in ${All_Shapes[@]}; do
#    echo $shape
#    python3  PreprocessShape.py yeast $shape depleted_
# done
for shape in ${All_Shapes[@]}; do
   echo $shape
   python3  PreprocessShape.py yeast $shape depleted_ &
done


#python3 -m cProfile cluster2shape.py

#python3 -m cProfile cluster2shape2.py 5 ${SLURM_ARRAY_TASK_ID} yeast
#python3  nucPosPredict.py 10 $(for i in $(seq 3 10); do echo $i; done)
