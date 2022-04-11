#!/usr/bin/env python
# coding: utf-8

# In[160]:


### Turn one hot encoded .mat to fasta seq
import math
import pandas as pd
import numpy as np
import scipy.io as sio
from os.path import dirname, join as pjoin
import os.path
import time
import random
import sys
import gzip

def readGZInputAsArray(fileName):
    with gzip.open(fileName, 'rb') as myfile:
        data = myfile.readlines()

    # Strip newline
    for i in range(0, len(data)):
        data[i] = data[i].rstrip()
    # print(data)
    return data

def readInputAsArray(fileName):
    with open(fileName, 'r') as myfile:
        data = myfile.readlines()

    # Strip newline
    for i in range(0, len(data)):
        data[i] = data[i].rstrip()
    # print(data)
    return data

def list_str_to_float(arr):
    out = []
    for i in arr:
        out+=[float(i)]
    return out

def my_print(name, val):
    print(f'{name}: {val}')
    return

def append_list_to_np_arr(lst, np_arr):
    bs_lst=[]
    for i in range(len(lst)):
#         bs_lst+= [[lst[i], lst[i]]]
        bs_lst+= [[lst[i]]]
    
    np_lst=np.array([bs_lst])
#     my_print('np_lst.shape', np_lst.shape)

    if np_arr.shape[0] == 0:
        np_arr = np_lst 
    else:
        np_arr = np.concatenate((np_arr, np_lst), axis=0) 
#     my_print('np_lst.shape', np_arr.shape)
    return np_arr

def append_bplist_to_np_arr(lst, np_arr):
    bp_lst = []
    for i in range(1, len(lst)):
        bp_lst+= [[lst[i-1], lst[i]]]
    np_bp=np.array([bp_lst])
    
#     my_print('np_bp', np_bp.shape)
#     my_print('np_arr', np_arr.shape)
    if np_arr.shape[0] == 0:
        np_arr = np_bp
    else:
        np_arr = np.concatenate((np_arr, np_bp), axis=0) 
#     print(np_arr.shape)
    return np_arr


# In[162]:


def makeSingleShapeArr(seq_file, index):
#     print(">>>>>>>>> In making makeSingleShapeArr ")
    np_all_shape_arr = np.array([])
    
    Seqs = readGZInputAsArray(seq_file)
    seq_count=1
    bp_shape=False
    line = Seqs[index]
    
    l = line.split()

    curr_seq_shape_list = list_str_to_float(l[2:]) # get the shape values only
    if len(l[0]) > len(curr_seq_shape_list):
#         print("base step parameter")
        bp_shape = False
    else:
#         print("base pair parameter")
        bp_shape = True

    ## 1. for base pair shapes, there are 147 vals for each seq 
    if bp_shape: 
        if len(curr_seq_shape_list) < 147:
            while len(curr_seq_shape_list) < 147: # if the sequence is shorter than 147
                curr_seq_shape_list+=[0]

            np_all_shape_arr = append_bplist_to_np_arr(curr_seq_shape_list, np_all_shape_arr)
        else:
            start = 0
             # is a base pair shape, each 147bp seq has 147 shape vals, we use 147-2+1
            while (start < len(curr_seq_shape_list)-146):
                curr_start = start;
                curr_end = 146 + start;
                sub_seq_shape_list = curr_seq_shape_list[curr_start: curr_end+1];
                start+=1
                np_all_shape_arr = append_bplist_to_np_arr(sub_seq_shape_list, np_all_shape_arr)

    ## 2. for base step shapes, there are 146 vals for each seq 
    else:
        if len(curr_seq_shape_list) < 146:
            while len(curr_seq_shape_list) < 146: # if the sequence is shorter than 146
                curr_seq_shape_list+=[0]

            np_all_shape_arr = append_list_to_np_arr(curr_seq_shape_list, np_all_shape_arr)
        else:
#           print('>>>>>> Found longer sequences')
            start = 0
            # is a base step shape, each 147bp seq has 146 shape vals, we use 146-2+1
            while (start < len(curr_seq_shape_list)-145):
                curr_start = start;
                curr_end = 145 + start;
                sub_seq_shape_list = curr_seq_shape_list[curr_start: curr_end+1];
                start+=1
                np_all_shape_arr = append_list_to_np_arr(sub_seq_shape_list, np_all_shape_arr)
    seq_count+=1
#     my_print('np_all_shape_arr.shape', np_all_shape_arr.shape)
    return np_all_shape_arr;


###### Now merge the single shape arrs to a whole np array
def shapeEncode(seq_file, seq_index):
    All_Shapes=['Buckle-FL', 'Buckle', 'EP', 'HelT-FL', 'HelT', 'MGW-FL', 'MGW', 
                'Opening-FL', 'Opening', 'ProT-FL', 'ProT', 'Rise-FL', 'Rise', 'Roll-FL',
                'Roll', 'Shear-FL', 'Shear', 'Shift-FL', 'Shift', 'Slide-FL', 'Slide', 
                'Stagger-FL', 'Stagger', 'Stretch-FL', 'Stretch', 'Tilt-FL', 'Tilt']
#     All_Shapes=['Buckle-FL', 'Buckle']
    
    np_shape_E = np.array([]) # enriched seqs output

    for i in range(len(All_Shapes)):
        my_print('\ncurr shape is ', All_Shapes[i])
        np_all_shape_arr = makeSingleShapeArr(seq_file + All_Shapes[i] + '.gz', seq_index)
        if np_shape_E.shape[0] == 0:
            np_shape_E = np_all_shape_arr
        else:
            ## merge: consider the base step as separate features
            np_shape_E = np.concatenate([np_shape_E, np_all_shape_arr], axis=2) #(218, 146, 3)
    
    return np_shape_E


# In[163]:


def cluster2seq(Cluster, seq_file, species):
    debug=200
    all_curr_seq_index = []
    for i in range(len(Cluster)):
        if (Cluster[i][0] == '>'and i < debug):
#         if (Cluster[i][0] == '>'):
            curr_seq_index = []
            all_curr_seq_index += [curr_seq_index] # a 2d list of cluster indexes
        elif (Cluster[i][0] != '>' and i < debug):
#         elif (Cluster[i][0] != '>'):
            idx = int(Cluster[i].split()[2][1:-3]) #seq indexes
            curr_seq_index.append(idx)

######################## train test from different clusters ########################
    print(f"num clusters {len(all_curr_seq_index)}") # 0 encoded
    train_seq_num = int(math.ceil(0.6*len(all_curr_seq_index)))
    print(f"target train seqs {train_seq_num}")
    test_seq_num = len(all_curr_seq_index) - train_seq_num
    print(f"target test seqs {test_seq_num}")
    
    ## train test data split
    Train_data = np.array([])
    Test_data = np.array([])
    train_clsts = []
    test_clsts = []
    seq_num = 0
    max_seq_in_clstr = 1 # **** can change this numbe ****
    clstr_id = 0
    test_start = True
    for clstr in all_curr_seq_index:
        
        # 1. select the seqs from clstr
        np_selection_encode = np.array([])
        if (len(clstr) <= max_seq_in_clstr):
            selection = clstr
        else:
            # print(clstr)
            # print("large cluster")
            random.shuffle(clstr)
            selection = clstr[-max_seq_in_clstr:] # Only want x seqs from large clusters
            # print(selection)
        
        # 2. encode the seqs from clstr   
        for seq_idx in selection:
            np_seq = shapeEncode(seq_file, seq_idx)
            if np_selection_encode.shape[0] == 0:
                np_selection_encode = np_seq
            else:
                np_selection_encode = np.concatenate([np_selection_encode, np_seq], axis=2)
            seq_num +=1
        
        # 3. assign the encoded seqs to Train/Test
        if seq_num > train_seq_num: # Now enough train seqs, encode test seqs
#             print('>>>> test seq')
            if test_start: # make sure test is in a different cluster
#                 print('skipping to next clustr')
                test_start = False
                np_selection_encode = np.array([])
            test_clsts+=[clstr_id]
#             print(f"np_selection_encode.shape {np.array(np_selection_encode).shape}")
            if Test_data.shape[0] == 0:                
#                 print('****** Init Test data')
                Test_data = np_selection_encode
            else:
                Test_data = np.concatenate([Test_data, np_selection_encode], axis=0)
            seq_num +=1
#             print(f"test_data.shape {np.array(Test_data).shape}")
            
        else: # first encode train seqs
#             print('>>>> train seq')
            train_clsts+=[clstr_id]
#             print(f"np_selection_encode.shape {np.array(np_selection_encode).shape}")
            if Train_data.shape[0] == 0:
#                 print('**** Init Train data')
                Train_data = np_selection_encode
            else:
                Train_data = np.concatenate([Train_data, np_selection_encode], axis=0)
#             print(f"train_data.shape {np.array(Train_data).shape}")
        clstr_id+=1
    print("test_clsts")
    print(test_clsts)
    print("train_clsts")
    print(train_clsts)
    
    return Train_data, Test_data


# In[164]:


###### First, preprocess the shape files
# seq_file='/home/yibei/Downloads/yibei_predictions/'
# preprocessShapeFile(seq_file, All_Shapes)

# species = sys.argv[1]
species ='yeast'

# path='/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/cluster_seq/'
# path='/project/rohs_102/share/nucleosome_occupancy_data/'
# path='/home/yibei/Projects/nucleosome_occupancy/tbinet_stuff/cluster_seq/'
path='/project/rohs_108/yibeijia/nucleosome_occupancy/tbinet_stuff/cluster_seq/'
path='/project/rohs_102/share/nucleosome_occupancy_data/'
if species =='human': 
    Dseqs = readInputAsArray(path+'Dseqs.txt')[1:] # There is no Dseqs
    Dcluster = readInputAsArray(path+'Dseqs80.clstr')
    
    Eseqs = readInputAsArray(path+species+'Eseqs.txt')[1:] # avoid the first \n
    Ecluster = readInputAsArray(path+species+'Eseqs80_every20.clstr')    
    
elif species == 'yeast':
    Eseqs = readInputAsArray(path+'Eseqs.txt')[1:] # avoid the first \n
    Ecluster = readInputAsArray(path+'Eseqs80.clstr')
    
    Dcluster = readInputAsArray(path+'Dseqs80.clstr')
    Dseqs = readInputAsArray(path+'Dseqs.txt')[1:]

else:
    Dseqs = readInputAsArray(path+'Dseqs.txt')[1:] # There is no Dseqs
    Dcluster = readInputAsArray(path+'Dseqs80.clstr')
    
    Eseqs = readInputAsArray(path+species+'Eseqs.txt')[1:] # avoid the first \n
    Ecluster = readInputAsArray(path+species+'Eseqs80_every5.clstr')

print(f">>>>>>>>Preprocessing Enriched data for {species}")
seq_file='/project/rohs_108/yibeijia/data/yibei_predictions/processed_enriched_'
# seq_file='/home/yibei/Downloads/yibei_predictions/processed_enriched_'
E_train_data, E_test_data = cluster2seq(Ecluster, seq_file, species)

print(">>>>>>>>Preprocessing Depleted data")
seq_file='/project/rohs_108/yibeijia/data/yibei_predictions/processed_depleted_'
# seq_file='/home/yibei/Downloads/yibei_predictions/processed_depleted_'
D_train_data, D_test_data = cluster2seq(Dcluster, seq_file, species)

print(f"E train  shape is: {E_train_data.shape}")
print(f"E TEST  shape is: {E_test_data.shape}")
print(f"D train  shape is: {D_train_data.shape}")
print(f"D TEST  shape is: {D_test_data.shape}")
print(">>>>>>>>Preparing Train Test dataset")
Train_data = np.concatenate((E_train_data, D_train_data), axis=0)
Edata_labels = np.ones(E_train_data.shape[0])
Ddata_labels = np.zeros(D_train_data.shape[0])
Train_labels = np.concatenate((Edata_labels, Ddata_labels), axis=0)
print(f"train data  shape is: {Train_data.shape}")
print(f"train label shape is: {Train_labels.shape}")

Test_data = np.concatenate((E_test_data, D_test_data), axis=0)
Edata_labels = np.ones(E_test_data.shape[0])
Ddata_labels = np.zeros(D_test_data.shape[0])
Test_labels = np.concatenate((Edata_labels, Ddata_labels), axis=0)
print(f"test data shape is: {Test_data.shape}")
print(f"test label shape is: {Test_labels.shape}")

# data_path='/home/yibei/Projects/data/train_test_data/'
data_path='/project/rohs_108/yibeijia/nucleosome_occupancy/data/train_test_data/'
Data = {"Test_data" : np.array(Test_data), "Test_labels" : Test_labels}
sio.savemat(data_path+species+'ShapeDNNSeqs_Test.mat', Data,  do_compression=True)

Data = {"Train_data" : np.array(Train_data), "Train_data" : Train_labels}
sio.savemat(data_path+species+'ShapeDNNSeqs_Train.mat', Data,  do_compression=True)

