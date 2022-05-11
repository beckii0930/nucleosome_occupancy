#!/usr/bin/env python
# coding: utf-8

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



def preprocessShapeFile(seq_file, All_Shapes, test_clsts, train_clsts, file_list ):
	out = ''
	test_out = ''
	train_out = ''
	seq_count = -1
#     debug = 0
    
	for enriched in file_list:
		print(enriched)
		for shape in All_Shapes:
			print(shape)
			Seqs = readInputAsArray(seq_file + enriched+'regions_seqOnly_'+shape+'.txt')
			#Seqs = readInputAsArray(seq_file +'humanEseqs_seqOnly_200000_'+shape+'.txt')
			#Seqs = readInputAsArray(seq_file +'wormEseqs_seqOnly_200000_'+shape+'.txt')
			#Seqs = readInputAsArray(seq_file +'flyEseqs_seqOnly_200000_'+shape+'.txt')
			for line in Seqs:
				l = line.split()
				if l[1] != '3': # Only consider flanking region of size 3
					continue;
				curr_seq_shape_list = list_str_to_float(l[2:])
				curr_seq = l[0]

				if len(curr_seq) > len(curr_seq_shape_list):
					bp_shape = False
				else:
					bp_shape = True

				if bp_shape:
					if len(curr_seq_shape_list) < 147:
						while len(curr_seq_shape_list) < 147: # if the sequence is shorter than 147
							curr_seq_shape_list+=[0]
							curr_seq += 'N'
						seq_count+=1
						if seq_count in test_clsts:
							# print(f'seq {seq_count} is in test')
							# print(curr_seq)
							out+=curr_seq + ' ' + l[1] + ' '
							out+=' '.join(str(e) for e in curr_seq_shape_list)+'\n'
							test_out+=curr_seq + ' ' + l[1] + ' '
							test_out+=' '.join(str(e) for e in curr_seq_shape_list)+'\n'
						elif seq_count in train_clsts:
							# print(f'seq {seq_count} is in train')
							# print(curr_seq)
							out+=curr_seq + ' ' + l[1] + ' '
							out+=' '.join(str(e) for e in curr_seq_shape_list)+'\n'
							train_out+=curr_seq + ' ' + l[1] + ' '
							train_out+=' '.join(str(e) for e in curr_seq_shape_list)+'\n'

					else:
						start = 0
						 # is a base pair shape, each 147bp seq has 147 shape vals, we use 147-2+1
						while (start < len(curr_seq_shape_list)-146):
							curr_start = start;
							curr_end = 146 + start;
							out_seq = curr_seq[curr_start:curr_end+1] 
							seq_count+=1

							out_seq_shape_list = curr_seq_shape_list[curr_start: curr_end+1];
							start+=1
							if seq_count in test_clsts:
								# print(f'seq {seq_count} is in test')
								# print(out_seq)
								out+=out_seq + ' ' + l[1] + ' '
								out+=' '.join(str(e) for e in out_seq_shape_list)+'\n'
								test_out+=out_seq + ' ' + l[1] + ' '
								test_out+=' '.join(str(e) for e in out_seq_shape_list)+'\n'
							elif seq_count in train_clsts:
								# print(f'seq {seq_count} is in train')
								# print(out_seq)
								out+=out_seq + ' ' + l[1] + ' '
								out+=' '.join(str(e) for e in out_seq_shape_list)+'\n'
								train_out+=out_seq + ' ' + l[1] + ' '
								train_out+=' '.join(str(e) for e in out_seq_shape_list)+'\n'
								
				else:
					if len(curr_seq_shape_list) < 146:
						while len(curr_seq_shape_list) < 146: # if the sequence is shorter than 146
							curr_seq_shape_list+=[0]
							curr_seq += 'N'
						seq_count+=1
						if seq_count in test_clsts:
							# print(f'seq {seq_count} is in test')
							# print(curr_seq)
							out+=curr_seq + ' ' + l[1] + ' '
							out+=' '.join(str(e) for e in curr_seq_shape_list)+'\n'
							test_out+=curr_seq + ' ' + l[1] + ' '
							test_out+=' '.join(str(e) for e in curr_seq_shape_list)+'\n'
						elif seq_count in train_clsts:
							# print(f'seq {seq_count} is in train')
							# print(curr_seq)
							out+=curr_seq + ' ' + l[1] + ' '
							out+=' '.join(str(e) for e in curr_seq_shape_list)+'\n'
							train_out+=curr_seq + ' ' + l[1] + ' '
							train_out+=' '.join(str(e) for e in curr_seq_shape_list)+'\n'
					else:
#						  print('>>>>>> Found longer sequences')
#						  print(l[0])
#						  print((curr_seq_shape_list))
						start = 0
						# is a base step shape, each 147bp seq has 146 shape vals, we use 146-2+1
						while (start < len(curr_seq_shape_list)-146):
							curr_start = start;
							curr_end = 146 + start;
							out_seq = curr_seq[curr_start:curr_end+1] 
							seq_count+=1

							out_seq_shape_list = curr_seq_shape_list[curr_start: curr_end];
							start+=1

							if seq_count in test_clsts:
								# print(f'seq {seq_count} is in test')
								# print(out_seq)
								out+=out_seq + ' ' + l[1] + ' '
								out+=' '.join(str(e) for e in out_seq_shape_list)+'\n'
								test_out+=out_seq + ' ' + l[1] + ' '
								test_out+=' '.join(str(e) for e in out_seq_shape_list)+'\n'
							elif seq_count in train_clsts:
								# print(f'seq {seq_count} is in train')
								# print(out_seq)
								out+=out_seq + ' ' + l[1] + ' '
								out+=' '.join(str(e) for e in out_seq_shape_list)+'\n'
								train_out+=out_seq + ' ' + l[1] + ' '
								train_out+=' '.join(str(e) for e in out_seq_shape_list)+'\n'
			f = open(seq_file+'all_processed_'+enriched+shape+".txt", "w+")
			f.write(out)
			f.close()
			f = open(seq_file+'train_processed_'+enriched+shape+".txt", "w+")
			f.write(train_out)
			f.close()
			f = open(seq_file+'test_processed_'+enriched+shape+".txt", "w+")
			f.write(test_out)
			f.close()
            # with gzip.open(seq_file+'processed_'+enriched+shape+'.gz', 'wb') as f:
            #     f.write(out.encode())
            # with gzip.open(seq_file+'train_processed_'+enriched+shape+'.gz', 'wb') as f:
            #     f.write(train_out.encode())
            # with gzip.open(seq_file+'test_processed_'+enriched+shape+'.gz', 'wb') as f:
            #     f.write(test_out.encode())


def getClusterIndex(cluster_file, seq_file, species):
    all_curr_seq_index = []
    curr_seq_index = []
    train_clsts = []
    test_clsts = []

    seq_num =1
    test_start = True
    total_select_seq_num = 0
    max_seq_in_clstr = 5
    ## First, get all the selected index for each cluster_file:
    # for i in range(2):
    for i in range(len(cluster_file)):
       # if (cluster_file[i][0] == '>'and i < debug):
        if (cluster_file[i][0] == '>'):  # a new cluster
            if len(curr_seq_index) > max_seq_in_clstr:
                total_select_seq_num += max_seq_in_clstr
            else:
                total_select_seq_num += len(curr_seq_index)

            curr_seq_index = []
            all_curr_seq_index += [curr_seq_index] # a 2d list of cluster_file indexes
        #elif (cluster_file[i][0] != '>' and i < debug):
        elif (cluster_file[i][0] != '>'):
            idx = int(cluster_file[i].split()[2][1:-3]) #seq indexes
            curr_seq_index.append(idx)
            

    print(f"num clusters {len(all_curr_seq_index)}") # 0 encoded
    print(f"Total seq selected is {total_select_seq_num}")

    train_seq_num = int(math.ceil(0.6*total_select_seq_num))
    print(f"target train seqs {train_seq_num}")
    test_seq_num = total_select_seq_num - train_seq_num
    print(f"target test seqs {test_seq_num}")
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

        # 2. assign the seq indexes to train/test
        if len(train_clsts) > train_seq_num: # Now enough train seqs, encode test seqs
            # print('>>>> test seq')
            if test_start: # make sure test is in a different cluster
                test_start = False
            test_clsts+=selection
        else: 
            # print('>>>> train seq')
            train_clsts+=selection
                
    print(f"train_data.shape {len(train_clsts)}")
    print(f"test_data.shape {len(test_clsts)}")

    train_clsts.sort() 
    test_clsts.sort() 
    # print(test_clsts)
    # print(train_clsts)
    return test_clsts, train_clsts

All_Shapes=['Buckle-FL', 'Buckle', 'EP', 'HelT-FL', 'HelT', 'MGW-FL', 'MGW',
             'Opening-FL', 'Opening', 'ProT-FL', 'ProT', 'Rise-FL', 'Rise', 'Roll-FL',
             'Roll', 'Shear-FL', 'Shear', 'Shift-FL', 'Shift', 'Slide-FL', 'Slide',
             'Stagger-FL', 'Stagger', 'Stretch-FL', 'Stretch', 'Tilt-FL', 'Tilt']

All_Shapes=['Buckle', 'EP', 'HelT-FL', 'HelT', 'MGW-FL', 'MGW',
'Opening-FL', 'Opening', 'ProT-FL', 'ProT', 'Rise-FL', 'Rise', 'Roll-FL',
'Roll', 'Shear-FL', 'Shear', 'Shift-FL', 'Shift', 'Slide-FL', 'Slide',
'Stagger-FL', 'Stagger', 'Stretch-FL', 'Tilt-FL', 'Tilt']
All_Shapes=['Buckle-FL', 'Stretch']

seq_file='/project/rohs_108/yibeijia/data/yibei_predictions2/'
#seq_file='/home/yibei/Downloads/yibei_predictions/'
species='yeast'

#path='/project/rohs_108/yibeijia/nucleosome_occupancy/tbinet_stuff/cluster_seq/'
path='/project/rohs_102/share/nucleosome_occupancy_data/'
#path='/home/yibei/Projects/data/'

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

# file_list = ['enriched_']
# E_test_clsts, E_train_clsts = getClusterIndex(Ecluster, seq_file, species)
# preprocessShapeFile(seq_file, All_Shapes, E_test_clsts, E_train_clsts,file_list )

file_list = ['depleted_']
D_test_clsts, D_train_clsts = getClusterIndex(Dcluster, seq_file, species)
preprocessShapeFile(seq_file, All_Shapes, D_test_clsts, D_train_clsts,file_list )

