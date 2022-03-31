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
Ndata_list = np.array([])
N=1
print(f"N is: {N}")
# for i in range(1, 51):
for i in range(1, 6):
    print(i)
    # fname = '/scratch2/yibeijia/data/nucleosome_occupancy_' + str(i) + '.mat'
    # fname = '/Users/yibeijia/Downloads/nucleosome_occupancy/data/nucleosome_occupancy_' + str(i) + '.mat'
    fname = '/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/wormSeqs_'+str(i)+'.mat'
    if os.path.isfile(fname):
        mat_fname = pjoin(fname)
        mat_contents = sio.loadmat(mat_fname)
        
        edata=mat_contents['EnrichedData']
        ddata=np.array([])# when the seqs are all enriched
        ndata=np.array([])
        # ddata=mat_contents['DepletedData']
        # ndata=mat_contents['NeutralData']

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

        if ndata.size >0:
            if Ndata_list.size == 0:
                print("Init Nlist")
                Ndata_list = ndata[::N];
            else:

                Ndata_list = np.concatenate((Ndata_list, ndata[::N]), axis=0)
    # else:
    #     for j in range(0, 20):
    #         fname = '/Users/yibeijia/Downloads/nucleosome_occupancy/data/nucleosome_occupancy_' + str(i) + '_' + str(j) + '.mat'
    #         # fname = '/scratch2/yibeijia/data/nucleosome_occupancy_' + str(i) + '_' + str(j) + '.mat'
    #         if os.path.isfile(fname):
    #             mat_fname = pjoin(fname)
    #             mat_contents = sio.loadmat(mat_fname)
    #             edata=mat_contents['EnrichedData']
    #             ddata=mat_contents['DepletedData']
    #             ndata=mat_contents['NeutralData']
    #             if edata.size >0:
    #                 if Edata_list.size == 0:
    #                     print("Init Elist")
    #                     Edata_list = edata[::N];
    #                 else:
    #                     Edata_list = np.concatenate((Edata_list, edata[::N]), axis=0)
                  
    #             if ddata.size >0:
    #                 if Ddata_list.size == 0:
    #                     print("Init Dlist")
    #                     Ddata_list = ddata[::N];
    #                 else:
    #                     Ddata_list = np.concatenate((Ddata_list, ddata[::N]), axis=0)

    #             if ndata.size >0:
    #                 if Ndata_list.size == 0:
    #                     print("Init Dlist")
    #                     Ndata_list = ndata[::N];
    #                 else:
    #                     Ndata_list = np.concatenate((Ndata_list, ndata[::N]), axis=0)

print("Enriched list size")
print(Edata_list.shape)    
print("Depleted list size")
print(Ddata_list.shape)  
print("Neutral list size")
print(Ndata_list.shape)    
toc=time.perf_counter()
print(f"Create master list took {toc - tic:0.4f} seconds")

D_empty = True
E_empty = True
Ddata_labels=np.array([])
Edata_labels=np.array([])
data_labels=np.array([])
data_list =np.empty((1, 1, 1))
start = True
if Ddata_list.size >0:
    if (Ddata_list.shape[0] > 0):
        print("D non empty")
        D_empty =False
        Ddata_labels = np.zeros(Ddata_list.shape[0])

if Edata_list.size >0:
    if (Edata_list.shape[0] > 0):
        if start:
            # print(Edata_list[0])
            start=False
        print("E non empty")
        E_empty =False
        Edata_labels = np.ones(Edata_list.shape[0])

if (E_empty==False and D_empty==False):
    print("both non empty")
    Test_labels =  np.concatenate((Edata_labels, Ddata_labels), axis=0)
    Test_data = np.concatenate((Edata_list, Ddata_list), axis=0)
elif (E_empty and D_empty==False):
    print("D non empty, E is empty")
    Test_labels =  Ddata_labels
    Test_data=Ddata_list
elif (E_empty==False and D_empty):
    print("E non empty, D is empty")
    Test_labels = Edata_labels
    Test_data=Edata_list

# print(data_list[2][0:4])
## Decide to not mix neutral data in the dataset because training is not with the neutral data.
print(f" All data shape is: {Test_data.shape}")
print(f" All data label shape is: {Test_labels.shape}")
# from sklearn.model_selection import train_test_split
# Train_data, Test_data, Train_labels, Test_labels = train_test_split(data_list, data_labels, test_size=0.40, random_state=42)

# print(f"No. of training sequeces: {Train_data.shape[0]}")
# print(f"No. of testing sequences: {Test_data.shape[0]}")

# Data = {"Train_data" : np.array(Train_data), "Train_labels" : Train_labels}
# sio.savemat('/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/Train_data.mat', Data, do_compression=True)
# # sio.savemat('/scratch2/yibeijia/data/train_test_data/Train_data.mat',Data, do_compression=True)

# Data = {"Test_data" : np.array(Test_data), "Test_labels" : Test_labels}
# sio.savemat('/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/Test_data.mat', Data,  do_compression=True)
# Test_data = np.concatenate((Edata_list, Ddata_list), axis=0)
# Test_labels = data_labels 
Data = {"Test_data" : np.array(Test_data), "Test_labels" : Test_labels}   
print(f"No. of test sequeces: {Test_data.shape[0]}")
# sio.savemat('/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/wormTest.mat', Data
sio.savemat('/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/wormTest.mat', Data, do_compression=True)
# print(Data)
# sio.savemat('/scratch2/yibeijia/data/train_test_data/Test_data.mat',Data, do_compression=True)
