import pandas as pd
import numpy as np
import scipy.io as sio
from os.path import dirname, join as pjoin
import os.path
import time
import random

tic = time.perf_counter()
Edata_list = np.array([])
Ddata_list = np.array([])

for i in range(1, 51):
    print(i)
    N = 1000 #take every 100000 points
    fname = '/scratch2/yibeijia/data/train_test_data/nucleosome_occupancy_' + str(i) + '.mat'
    #fname = '/project/rohs_108/yibeijia/nucleosome_occupancy/nucleosome_occupancy_' + str(i) + '.mat'
    if os.path.isfile(fname):
        mat_fname = pjoin(fname)
        mat_contents = sio.loadmat(mat_fname)
        
        edata=mat_contents['EnrichedData']
        ddata=mat_contents['DepletedData']

        if edata.size >0:
            if Edata_list.size == 0:
                print("Init Elist")
                #print(len(edata))
                
                #print("subsection of Edata:")
                Edata_list=edata[::N]
                #print(len(Edata_list))
                #Edata_list = edata[1:N:len(edata),1:N:len(edata[0])];
            else:
                Edata_list = np.concatenate((Edata_list, edata[::N]), axis=0)
          
        if ddata.size >0:
            if Ddata_list.size == 0:
                print("Init Dlist")
                Ddata_list = ddata[::N];
            else:

                Ddata_list = np.concatenate((Ddata_list, ddata[::N]), axis=0)
    else:
#         print("not exist, below are subfiles")
        for j in range(0, 20):
#             print(j)
            fname = '/scratch2/yibeijia/data/train_test_data/nucleosome_occupancy_' + str(i) + '_' + str(j) + '.mat'

            #fname = '/project/rohs_108/yibeijia/nucleosome_occupancy/nucleosome_occupancy_' + str(i) + '_' + str(j) + '.mat'
            if os.path.isfile(fname):
                mat_fname = pjoin(fname)
                mat_contents = sio.loadmat(mat_fname)
                edata=mat_contents['EnrichedData']
                ddata=mat_contents['DepletedData']

                if edata.size >0:
                    if Edata_list.size == 0:
                        print("Init Elist")
                        Edata_list = edata[::N];
                    else:
                        Edata_list = np.concatenate((Edata_list, edata[::N]), axis=0)
                  
                if ddata.size >0:
                    if Ddata_list.size == 0:
                        print("Init Dlist")
                        Ddata_list = ddata[::N];
                    else:
                        Ddata_list = np.concatenate((Ddata_list, ddata[::N]), axis=0)

print("Enriched list size")
print(Edata_list.shape)    
print("Depleted list size")
print(Ddata_list.shape)    
toc=time.perf_counter()
print(f"Create master list took {toc - tic:0.4f} seconds")
Edf = pd.DataFrame(Edata_list)
Ddf = pd.DataFrame(Ddata_list)
print("Edf")
print(Edf)
print("Ddf")
print(Ddf)
# Edf = pd.DataFrame(Edata_list, columns = ['A','C','G','T'])
# Ddf = pd.DataFrame(Ddata_list, columns = ['A','C','G','T'])
Etraining_data = Edf.sample(frac=0.6, random_state=25)
print("Etraining_data")
print(Etraining_data)
Etesting_data = Edf.drop(Etraining_data.index)


Dtraining_data = Ddf.sample(frac=0.6, random_state=25)
Dtesting_data = Ddf.drop(Dtraining_data.index)

print("Dtraining_data")
print(Dtraining_data)

print("Etesting_data")
print(Etesting_data)

print("Dtesting_data")
print(Dtesting_data)


print(f"No. of training sequeces: {Etraining_data.shape[0]}")
print(f"No. of testing sequences: {Dtesting_data.shape[0]}")


# Data = {"EnrichedData": np.array(allEnrichSeqArr), "DepletedData": np.array(allDepleteSeqArr)};
Data = {"Etraining_data" : np.array(Etraining_data)}
sio.savemat('/scratch2/yibeijia/data/Etraining_data.mat', Data, do_compression=True)

Data = {"Dtraining_data" : np.array(Dtraining_data)}
sio.savemat('/scratch2/yibeijia/data/Dtraining_data.mat', Data,  do_compression=True)

Data = {"Etesting_data" : np.array(Etesting_data)}
sio.savemat('/scratch2/yibeijia/data/Etesting_data.mat', Data,  do_compression=True)

Data = {"Dtesting_data" : np.array(Dtesting_data)}
sio.savemat('/scratch2/yibeijia/data/Dtesting_data.mat', Data,  do_compression=True)
