import math
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

N=5
print(f"N is: {N}")
for i in range(1, 51):
# for i in range(1, 5):
    print(i)
    # fname = '/scratch2/yibeijia/data/nucleosome_occupancy_' + str(i) + '.mat'
    fname = '/Users/yibeijia/Downloads/nucleosome_occupancy/data/nucleosome_occupancy_' + str(i) + '.mat'
    if os.path.isfile(fname):
        mat_fname = pjoin(fname)
        mat_contents = sio.loadmat(mat_fname)
        
        edata=mat_contents['EnrichedData']
        ddata=mat_contents['DepletedData']

        if edata.size >0:
            if Edata_list.size == 0:
                print("Init Elist")
                Edata_list=edata[::N]
            else:
                Edata_list = np.concatenate((Edata_list, edata[::N]), axis=0)
          
        if ddata.size >0:
            if Ddata_list.size == 0:
                print("Init Dlist")
                Ddata_list = ddata[::N];
            else:

                Ddata_list = np.concatenate((Ddata_list, ddata[::N]), axis=0)
    else:
        for j in range(0, 20):
            fname = '/Users/yibeijia/Downloads/nucleosome_occupancy/data/nucleosome_occupancy_' + str(i) + '_' + str(j) + '.mat'
            # fname = '/scratch2/yibeijia/data/nucleosome_occupancy_' + str(i) + '_' + str(j) + '.mat'
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

from sklearn.model_selection import train_test_split

data_list =np.empty((1, 1, 1))
if (Ddata_list.shape[0] > 0 and Edata_list.shape[0] > 0):
    print("both non empty")
    data_list = np.concatenate((Ddata_list, Edata_list), axis=0)
elif (Ddata_list.shape[0] > 0):
    print("D non empty, E is empty")
    data_list = Ddata_list
elif (Edata_list.shape[0] > 0):
    print("E non empty, D is empty")
    data_list = Edata_list

print(f" All data shape is: {data_list.shape}")

Ddata_labels = np.zeros(Ddata_list.shape[0])
Edata_labels = np.ones(Edata_list.shape[0])
data_labels =  np.concatenate((Ddata_labels, Edata_labels), axis=0)
print(f" All data label shape is: {data_labels.shape}")
Train_data, Test_data, Train_labels, Test_labels = train_test_split(data_list, data_labels, test_size=0.40, random_state=42)

print(f"No. of training sequeces: {Train_data.shape[0]}")
print(f"No. of testing sequences: {Test_data.shape[0]}")

Data = {"Train_data" : np.array(Train_data), "Train_labels" : Train_labels}
sio.savemat('/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/Train_data.mat', Data, do_compression=True)
# sio.savemat('/scratch2/yibeijia/data/train_test_data/Train_data.mat',Data, do_compression=True)

Data = {"Test_data" : np.array(Test_data), "Test_labels" : Test_labels}
sio.savemat('/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/Test_data.mat', Data,  do_compression=True)
# sio.savemat('/scratch2/yibeijia/data/train_test_data/Test_data.mat',Data, do_compression=True)
