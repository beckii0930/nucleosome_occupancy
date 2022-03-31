### Turn one hot encoded .mat to fasta seq
import math
import pandas as pd
import numpy as np
import scipy.io as sio
from os.path import dirname, join as pjoin
import os.path
import time
import random

def readInputAsArray(fileName):
    with open(fileName, 'r') as myfile:
        data = myfile.readlines()

    # Strip newline
    for i in range(0, len(data)):
        data[i] = data[i].rstrip()
    # print(data)
    return data


def oneHotEncode(seq):
    import numpy as np
    seq2=list()
    mapping = {"A":[1., 0., 0., 0.], "C": [0., 1., 0., 0.], "G": [0., 0., 1., 0.], "T":[0., 0., 0., 1.]};
    for i in seq:
    	seq2.append(mapping[i]  if i in mapping.keys() else [0., 0., 0., 0.]);
    return seq2;


def cluster2seq(Ecluster, Eseqs,species):
    debug=100
    all_curr_seq_index = []
    for i in range(len(Ecluster)):
        # header
    #     if (Ecluster[i][0] == '>'and i < debug):
        if (Ecluster[i][0] == '>'):
            curr_seq_index = []
            all_curr_seq_index += [curr_seq_index]
    #     elif (Ecluster[i][0] != '>' and i < debug):
        elif (Ecluster[i][0] != '>'):
            idx = int(Ecluster[i].split()[2][1:-3])
            curr_seq_index.append(idx)

####################### Separate clusters for train and test split #######################
    print(f"num clusters {len(all_curr_seq_index)}") # 0 encoded
    train_seq_num = int(math.ceil(0.6*len(all_curr_seq_index)))
    print(f"target train seqs {train_seq_num}")
    test_seq_num = len(all_curr_seq_index) - train_seq_num
    print(f"target test seqs {test_seq_num}")
    
    ## train test data split
    E_train_data = []
    E_test_data = []
    seq_num = 0
    max_seq_in_clstr = 1
    for clstr in all_curr_seq_index:
        
        selection_encode = []
        if (len(clstr) <= max_seq_in_clstr):
            selection = clstr
        else:
            # print(clstr)
            # print("large cluster")
            random.shuffle(clstr)
            selection = clstr[-max_seq_in_clstr:] # Only want x seqs from large clusters
            # print(selection)
        for seq_idx in selection:
            seq = oneHotEncode(Eseqs[seq_idx*2+1])
            selection_encode += [seq]
            seq_num +=1

        if seq_num > train_seq_num:
            E_test_data += selection_encode
        else:
            E_train_data += selection_encode
#     print(f"train_data.shape {np.array(E_train_data).shape}")
#     print(f"test_data.shape {np.array(E_test_data).shape}")
    return np.array(E_train_data), np.array(E_test_data)

import sys
# Eseqs = readInputAsArray('/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/cluster_seq/Eseqs.txt')[1:]
# Eseqs = readInputAsArray('/scratch2/yibeijia/data/Eseqs.txt')[1:]
# Ecluster = readInputAsArray('/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/cluster_seq/Eseqs80.clstr')
# Ecluster = readInputAsArray('/scratch2/yibeijia/data/Eseqs80.clstr')
species = sys.argv[1]
Eseqs = readInputAsArray('/project/rohs_102/share/nucleosome_occupancy_data/'+species+'Eseqs.txt')[1:] # avoid the first \n
# Eseqs = readInputAsArray('/Users/yibeijia/Downloads/data/'+species+'Eseqs.txt')[1:] # avoid the first \n
# Eseqs = readInputAsArray('/scratch2/yibeijia/data/Eseqs.txt')[1:]
# Ecluster = readInputAsArray('/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/cluster_seq/Eseqs80.clstr')
if species =='human':
    Ecluster = readInputAsArray('/project/rohs_102/share/nucleosome_occupancy_data/'+species+'Eseqs80_every20.clstr')    
else:
    Ecluster = readInputAsArray('/project/rohs_102/share/nucleosome_occupancy_data/'+species+'Eseqs80_every5.clstr')
# Ecluster = readInputAsArray('/Users/yibeijia/Downloads/data/'+species+'Eseqs80_every5.clstr')
# Dseqs = readInputAsArray('/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/cluster_seq/Dseqs.txt')[1:]
# Dseqs = readInputAsArray('/scratch2/yibeijia/data/Dseqs.txt')[1:]

# Dcluster = readInputAsArray('/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/cluster_seq/Dseqs80.clstr')
# Dcluster = readInputAsArray('/scratch2/yibeijia/data/Dseqs80.clstr')
# Nseqs = readInputAsArray('/Users/yibeijia/Downloads/data/Nseqs_labeled.txt')
# Nseqs = readInputAsArray('/Users/yibeijia/Downloads/data/Nseqs_labeled.txt')
print(f">>>>>>>>Preprocessing Enriched data for {species}")
E_train_data, E_test_data = cluster2seq(Ecluster, Eseqs,species)
# print(">>>>>>>>Preprocessing Depleted data")
# D_train_data, D_test_data = cluster2seq(Dcluster, Dseqs)

## Experiemnt 1: combine E & D to create training, testing dataset & labels
## train test from different cluster
# print(">>>>>>>>Preparing Train Test dataset")
# Train_data = np.concatenate((E_train_data, D_train_data), axis=0)
# Edata_labels = np.ones(E_train_data.shape[0])
# Ddata_labels = np.zeros(D_train_data.shape[0])
# Train_labels = np.concatenate((Edata_labels, Ddata_labels), axis=0)
# print(f"train data  shape is: {Train_data.shape}")
# print(f"train label shape is: {Train_labels.shape}")

# Test_data = np.concatenate((E_test_data, D_test_data), axis=0)
# Edata_labels = np.ones(E_test_data.shape[0])
# Ddata_labels = np.zeros(D_test_data.shape[0])
# Test_labels = np.concatenate((Edata_labels, Ddata_labels), axis=0)
# print(f"test data shape is: {Test_data.shape}")
# print(f"test label shape is: {Test_labels.shape}")

## Experiemnt 2: combine E & D to create training, testing dataset & labels
## train test from same cluster
from sklearn.model_selection import train_test_split

print(">>>>>>>>Preparing Train Test dataset")
Edata_list = np.concatenate((E_train_data, E_test_data), axis=0)
# Ddata_list = np.concatenate((D_train_data, D_test_data), axis=0)
Edata_labels = np.ones(Edata_list.shape[0])
# Ddata_labels = np.zeros(Ddata_list.shape[0])

# data_list = np.concatenate((Ddata_list, Edata_list), axis=0)
# data_labels =  np.concatenate((Ddata_labels, Edata_labels), axis=0)
Test_data=Edata_list
Test_labels=Edata_labels

# print(f" All data label shape is: {data_labels.shape}")
# Train_data, Test_data, Train_labels, Test_labels = train_test_split(data_list, data_labels, test_size=0.40, random_state=42)

# print(f"No. of training sequeces: {Train_data.shape[0]}")
print(f"No. of testing sequences: {Test_data.shape[0]}")

## Write to .mat file
# Data = {"Train_data" : np.array(Train_data), "Train_labels" : Train_labels}
# #sio.savemat('/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/Train_data.mat', Data, do_compression=True)
# sio.savemat('/scratch2/yibeijia/data/train_test_data/Train_data.mat',Data, do_compression=True)

Data = {"Test_data" : np.array(Test_data), "Test_labels" : Test_labels}
sio.savemat('/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/'+species+'_Test_data.mat', Data,  do_compression=True)
# sio.savemat('/scratch2/yibeijia/data/train_test_data/wormTest_data.mat',Data, do_compression=True)
