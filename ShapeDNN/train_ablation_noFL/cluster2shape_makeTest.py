
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

def SeqToMat(path, All_Shapes):
    np_E_train_data = np.array([])

    for shape in All_Shapes:

        # E_path= seq_file+'all_processed_'+'enriched_'+shape+".txt"
        E_train_path = path+shape+".txt"
        E_train_Seqs = readInputAsArray(E_train_path)
        my_print('\nLoading shape file for ', E_train_path)

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
        # print(E_train_data[0])
        print(f"Finished loading {shape}")
        print(f"b4 merge, np_E_train_data.shape: {np_E_train_data.shape}")
        np_curr_shape_arr = np.array([])
        if np_E_train_data.shape[0] == 0:
            np_E_train_data = np.array(E_train_data)
            print("In if")
            # print(np_E_train_data)
        else:
            print("########## ############################## in else")
            np_curr_shape_arr = np.array(E_train_data)
            print(f"np_curr_shape_arr.shape: {np_curr_shape_arr.shape}")
            # print(np_curr_shape_arr)
            # np_E_train_data = np.concatenate([np_E_train_data, np_curr_shape_arr], axis=1)
            np_E_train_data = np.concatenate((np_E_train_data, np_curr_shape_arr), axis=2)
            # print(np_E_train_data)
            print(f"after merging, np_E_train_data.shape.shape: {np_E_train_data.shape}")
        # print(np_E_train_data[0])
    return np_E_train_data


def main(seq_file, species, All_Shapes, E_data_path, D_data_path, out_data_path, ablation):

    print(">>>>>>>>>>>>.. Processing " + D_data_path)
    np_D_data = SeqToMat(D_data_path, All_Shapes)
    print(f"np_D_data.shape: {np_D_data.shape}")

    print(">>>>>>>>>>>>.. Processing " + E_data_path)
    np_E_data = SeqToMat(E_data_path, All_Shapes)
    print(f"np_E_data.shape: {np_E_data.shape}")

    Train_Test_data = np.concatenate((np_E_data, np_D_data), axis=0)
    Edata_labels = np.ones(np_E_data.shape[0])
    Ddata_labels = np.zeros(np_D_data.shape[0])
    Train_Test_labels = np.concatenate((Edata_labels, Ddata_labels), axis=0)
    print(f"train/test data  shape is: {Train_Test_data.shape}")
    print(f"train/test label shape is: {Train_Test_labels.shape}")

    section=int(sys.argv[2])
    total_sections=int(sys.argv[1])
    print(f"Current regions is: {section}\n");
    print(f"Total # of regions is: {total_sections}\n");
    print(f"Current regions is: {section}\n");
    total_lines = Train_Test_data.shape[0]
    section_length = math.floor(total_lines / total_sections);
    print(f"train/test curr region line is {section_length}")
    start_line = (section-1) * section_length;

    end_line = section * section_length-1;
    if (end_line > total_lines):
    	end_line = total_lines-1;

    data_label="Train_data"
    labels_label="Train_labels"
    mat_out=out_data_path+species+ablation+'_Shapes_Train_5seqsPerClustr_'
    if 'test_processed_' in E_data_path:
        data_label="Test_data"
        labels_label="Test_labels"
        mat_out=out_data_path+species+ablation+'_Shapes_Test_5seqsPerClustr_'

    DataToSave = {data_label: Train_Test_data[start_line:end_line,:,:], 
                labels_label: Train_Test_labels[start_line:end_line]}
    sio.savemat(mat_out+str(section)+'_'+str(total_sections)+'.mat', 
        DataToSave,  do_compression=True)

######################## ######################## Main ######################## ########################
seq_file='/project/rohs_108/yibeijia/data/yibei_predictions2/'
#seq_file='/Users/yibeijia/Downloads/data/yibei_predictions2/'
#seq_file='/home/yibei/Downloads/yibei_predictions/'
species='yeast'

All_Shapes=['Buckle-FL', 'Buckle', 'EP', 'HelT-FL', 'HelT', 'MGW-FL', 'MGW',
              'Opening-FL', 'Opening', 'ProT-FL', 'ProT', 'Rise-FL', 'Rise', 'Roll-FL',
              'Roll', 'Shear-FL', 'Shear', 'Shift-FL', 'Shift', 'Slide-FL', 'Slide',
              'Stagger-FL', 'Stagger', 'Stretch-FL', 'Stretch', 'Tilt-FL', 'Tilt']

All_Shapes_FL=['Buckle-FL','HelT-FL','MGW-FL','Opening-FL','ProT-FL', 'Rise-FL','Roll-FL',
			'Shear-FL','Shift-FL','Slide-FL','Stagger-FL','Stretch-FL',  'Tilt-FL']

All_Shapes_noFL=['Buckle', 'EP','HelT','MGW','Opening', 'ProT','Rise','Roll','Shear', 'Shift',
			'Slide','Stagger','Stretch', 'Tilt']

print(All_Shapes)
print(len(All_Shapes_FL))
print(len(All_Shapes_noFL))

out_data_path='/project/rohs_108/yibeijia/nucleosome_occupancy/data/train_test_data/'
#out_data_path='/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/'
#out_data_path='/home/yibei/Projects/data/train_test_data/'

E_train_path= seq_file+'train_processed_'+'enriched_'
D_train_path= seq_file+'train_processed_'+'depleted_'
#main(seq_file, species, All_Shapes, E_train_path, D_train_path, out_data_path,'All')
#main(seq_file, species, All_Shapes, E_train_path, D_train_path, out_data_path,'FL')
#main(seq_file, species, All_Shapes, E_train_path, D_train_path, out_data_path,'noFL')
# out_data_path='/project/rohs_108/yibeijia/nucleosome_occupancy/data/train_test_data/'
#out_data_path='/home/yibei/Projects/data/train_test_data/'

#E_train_path= seq_file+species+'_train_processed_'+'enriched_'
#D_train_path= seq_file+species+'_train_processed_'+'depleted_'

#main(seq_file, species, All_Shapes, E_train_path, D_train_path, out_data_path)

E_test_path= seq_file+species+'_test_processed_'+'enriched_'
D_test_path= seq_file+species+'_test_processed_'+'depleted_'

#main(seq_file, species, All_Shapes_FL, E_test_path, D_test_path, out_data_path,'All')
main(seq_file, species, All_Shapes_FL, E_test_path, D_test_path, out_data_path,'FL')
main(seq_file, species, All_Shapes_noFL, E_test_path, D_test_path, out_data_path,'noFL')

