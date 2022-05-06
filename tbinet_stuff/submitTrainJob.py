import os
nodes="#!/bin/bash\n#SBATCH --nodes=2\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=8\n"
time="#SBATCH --time=24:00:00\n"
mem="#SBATCH --mem=64GB\n"
gpu="#SBATCH --gres=gpu:2\n#SBATCH --partition=gpu\n"
mail="#SBATCH --mail-user=yibeijia@usc.edu\n"
header = nodes+time+mem+gpu+mail+"\n"

for i in range(1, 11):
	fname="train_dense10_kernel9_run"+str(i)
	out="#SBATCH --output="+fname+".out\n"
	env=". ~/.bashrc\nconda activate gpuenv2\nmodule load usc\nmodule load cuda/10.1.243\nfolder="

	fpath="/project/rohs_108/yibeijia/nucleosome_occupancy/tbinet_stuff/train_dense10_kernel9_run"+str(i)+"/\n"
	run="mkdir $folder\npython3 train.py $folder\n"
	
	script = header+out+"\n"+env+fpath+"\n"+run
	f=open(fname+".sl","w")

	f.write(script)
	os.system("chmod 777 "+fname+".sl")
	f.close()
	os.system("sbatch "+fname+".sl")

