import os
import time

nodes="#!/bin/bash\n#SBATCH --nodes=2\n#SBATCH --ntasks-per-node=1\n#SBATCH --cpus-per-task=8\n"
time="#SBATCH --time=24:00:00\n"
mem="#SBATCH --mem=64GB\n"
gpu="#SBATCH --gres=gpu:2\n#SBATCH --partition=gpu\n"
mail="#SBATCH --mail-user=yibeijia@usc.edu\n"
header = nodes+time+mem+gpu+mail+"\n"
proj_path = "/project/rohs_108/yibeijia/nucleosome_occupancy/tbinet_stuff/"
train_path="train_dense10_kernel9_run"
test_path="test_dense10_kernel9_run"

for i in range(1, 11):
	fname=test_path+str(i)
	out="#SBATCH --output="+fname+".out\n"
	env=". ~/.bashrc\nconda activate gpuenv2\nmodule load usc\nmodule load cuda/10.1.243\nfolder="

	fpath=proj_path+train_path+str(i)+"/\n"
	run="mkdir $folder\npython3 test.py $folder\n"
	
	script = header+out+"\n"+env+fpath+"\n"+run
	f=open(fname+".sl","w")

	f.write(script)
	os.system("chmod 777 "+fname+".sl")
	f.close()
	#os.system("sbatch "+fname+".sl")


def getSeconds(time):
    if(len(time.split(":"))<2): return 0
    units = list(map(int,time.split(":")))
    total = 60 * units[-2] + units[-1]
    if(len(units) == 3):
        total += 360 * units[-3]
    return total

import time
while True:
	wait = 900
	for i in range(1,11): 
		fpath=proj_path+train_path+str(i)+"/"
		if (os.path.isdir(fpath)):
			for filename in os.listdir(fpath):
				if (filename.endswith(".h5")):
					os.system("sbatch "+test_path+str(i)+".sl")
	time.sleep(wait)
	 
