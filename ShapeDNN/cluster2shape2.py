
### Turn one hot encoded .mat to fasta seq
### This is a faster version of cluster2shape
### We encoded only max 5 seqs from each cluster using PreprocessShape.py and get the 
### train and test sequences in each shape files.
### Here, we encode/concat all the shape values for the sequences we selected into .mat
### input files for training.

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

def makeSingleShapeArr(Seqs, index):
#     print(">>>>>>>>> In making makeSingleShapeArr ")
    np_all_shape_arr = np.array([])
    
    
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

def SeqToMat(path, All_Shapes):
    np_E_train_data = np.array([])

    for shape in All_Shapes:

        my_print('\nLoading shape file for ', shape)

        # E_path= seq_file+'all_processed_'+'enriched_'+shape+".txt"
        E_train_path = path+shape+".txt"
        E_train_Seqs = readInputAsArray(E_train_path)

        E_train_data = []

        if len(E_train_Seqs[0].split()[0]) > len(E_train_Seqs[0].split()[2:]):
            print("base step parameter")
            bp_shape = False
        else:
            print("base pair parameter")
            bp_shape = True

        for seqs in E_train_Seqs:
            seq = seqs.split()[2:]

            curr_seq = []
            if bp_shape:
                for i in range(1, len(seq)):
                    curr_seq+=[[float(seq[i-1]), float(seq[i])]]
                E_train_data +=[curr_seq]
            else:
                for  item in seq:

                    curr_seq+=[[float(item)]]
                E_train_data +=[curr_seq]

        print(f"Finished loading {shape}")

        if np_E_train_data.shape[0] == 0:
            np_E_train_data = np.array(E_train_data)
        else:
            np_all_shape_arr = np.array(E_train_data)
            np_E_train_data = np.concatenate([np_E_train_data, np_all_shape_arr], axis=2)
        print(f"np_data.shape: {np_E_train_data.shape}")
    return np_E_train_data



seq_file='/project/rohs_108/yibeijia/data/yibei_predictions/'
seq_file='/home/yibei/Downloads/yibei_predictions/'
species='yeast'

All_Shapes=['Buckle-FL', 'Buckle', 'EP', 'HelT-FL', 'HelT', 'MGW-FL', 'MGW',
              'Opening-FL', 'Opening', 'ProT-FL', 'ProT', 'Rise-FL', 'Rise', 'Roll-FL',
              'Roll', 'Shear-FL', 'Shear', 'Shift-FL', 'Shift', 'Slide-FL', 'Slide',
			  'Stagger-FL', 'Stagger', 'Stretch-FL', 'Stretch', 'Tilt-FL', 'Tilt']
# All_Shapes=['Buckle-FL', 'Stretch','EP']
print(All_Shapes)

E_train_path= seq_file+'train_processed_'+'enriched_'
E_test_path= seq_file+'test_processed_'+'enriched_'

print("Gathering E train data")
np_E_train_data = SeqToMat(E_train_path, All_Shapes)

print("Gathering E test data")
np_E_test_data = SeqToMat(E_test_path, All_Shapes)
print(f"np_E_test_data.shape: {np_E_test_data.shape}")
print(f"np_E_train_data.shape: {np_E_train_data.shape}")

D_train_path= seq_file+'train_processed_'+'depleted_'
D_test_path= seq_file+'test_processed_'+'depleted_'

print("Gathering D train data")
np_D_train_data = SeqToMat(D_train_path, All_Shapes)

print("Gathering D test data")
np_D_test_data = SeqToMat(D_test_path, All_Shapes)
print(f"np_D_test_data.shape: {np_D_test_data.shape}")
print(f"np_D_train_data.shape: {np_D_train_data.shape}")


Train_data = np.concatenate((np_E_train_data, np_D_train_data), axis=0)
Edata_labels = np.ones(np_E_train_data.shape[0])
Ddata_labels = np.zeros(np_D_train_data.shape[0])
Train_labels = np.concatenate((Edata_labels, Ddata_labels), axis=0)
print(f"train data  shape is: {Train_data.shape}")
print(f"train label shape is: {Train_labels.shape}")

Test_data = np.concatenate((np_E_test_data, np_D_test_data), axis=0)
Edata_labels = np.ones(np_E_test_data.shape[0])
Ddata_labels = np.zeros(np_D_test_data.shape[0])
Test_labels = np.concatenate((Edata_labels, Ddata_labels), axis=0)
print(f"test data shape is: {Test_data.shape}")
print(f"test label shape is: {Test_labels.shape}")

data_path='/project/rohs_108/yibeijia/nucleosome_occupancy/data/train_test_data/'
data_path='/home/yibei/Projects/data/train_test_data/'

section=int(sys.argv[1])
total_sections=int(sys.argv[2])
print(f"Current regions is: {section}\n");
print(f"Total # of regions is: {total_sections}\n");

total_lines = Test_data.shape[0]
section_length = math.floor(total_lines / total_sections);
print(f"test curr region line is {section_length}")

start_line = (section-1) * section_length;
end_line = section * section_length-1;
if (end_line > total_lines):
	end_line = total_lines-1;
Test_Data = {"Test_data" : Test_data[start_line:end_line,:,:], "Test_labels" : Test_labels[start_line:end_line]}
sio.savemat(data_path+species+'All_Shapes_Test_5seqsPerClustr_'+str(section)+'_'+str(total_sections)+'.mat', Test_Data,  do_compression=True)

print(f"Current regions is: {section}\n");
total_lines = Train_data.shape[0]
section_length = math.floor(total_lines / total_sections);
print(f"train curr region line is {section_length}")
start_line = (section-1) * section_length;

end_line = section * section_length-1;
if (end_line > total_lines):
	end_line = total_lines-1;
Train_Data = {"Train_data" : Train_data[start_line:end_line,:,:], "Train_data" : Train_labels[start_line:end_line]}
sio.savemat(data_path+species+'All_Shapes_Train_5seqsPerClustr_'+str(section)+'_'+str(total_sections)+'.mat', Train_Data,  do_compression=True)
