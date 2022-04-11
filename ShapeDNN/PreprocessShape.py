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

def preprocessShapeFile(seq_file, All_Shapes):
    out = ''
#     debug = 0
    file_list = ['enriched_', 'depleted_']
    for enriched in file_list:
        print(enriched)
        for shape in All_Shapes:

            print(shape)
            Seqs = readInputAsArray(seq_file + enriched+shape+'.txt')
            for line in Seqs:
    #             if debug > 200: break;
    #             debug+=1
                l = line.split()
                if l[1] != '3': # Only consider flanking region of size 3
                    continue;
                curr_seq_shape_list = list_str_to_float(l[2:])
                curr_seq = l[0]

                if len(curr_seq) > len(curr_seq_shape_list):
    #                 print(shape)
    #                 print("base step parameter")
                    bp_shape = False
                else:
    #                 print("base pair parameter")
                    bp_shape = True

                if bp_shape:
                    if len(curr_seq_shape_list) < 147:
                        while len(curr_seq_shape_list) < 147: # if the sequence is shorter than 147
                            curr_seq_shape_list+=[0]
                            curr_seq += 'N'
                        out+=curr_seq + ' ' + l[1] + ' '
                        out+=' '.join(str(e) for e in curr_seq_shape_list)+'\n'

                    else:
                        start = 0
                         # is a base pair shape, each 147bp seq has 147 shape vals, we use 147-2+1
                        while (start < len(curr_seq_shape_list)-146):
                            curr_start = start;
                            curr_end = 146 + start;
                            out_seq = curr_seq[curr_start:curr_end+1] 
                            out_seq_shape_list = curr_seq_shape_list[curr_start: curr_end+1];
                            start+=1
                            out+=out_seq + ' ' + l[1] + ' '
                            out+=' '.join(str(e) for e in out_seq_shape_list)+'\n'
                else:
                    if len(curr_seq_shape_list) < 146:
                        while len(curr_seq_shape_list) < 146: # if the sequence is shorter than 146
                            curr_seq_shape_list+=[0]
                            curr_seq += 'N'
                        out+=curr_seq + ' ' + l[1] + ' '
                        out+=' '.join(str(e) for e in curr_seq_shape_list)+'\n'
                    else:
    #                     print('>>>>>> Found longer sequences')
    #                     print(l[0])
    #                     print((curr_seq_shape_list))
                        start = 0
                        # is a base step shape, each 147bp seq has 146 shape vals, we use 146-2+1
                        while (start < len(curr_seq_shape_list)-146):
                            curr_start = start;
                            curr_end = 146 + start;
                            out_seq = curr_seq[curr_start:curr_end+1] 
                            start+=1
                            out+=out_seq + ' ' + l[1] + ' '

                            out_seq_shape_list = curr_seq_shape_list[curr_start: curr_end];
                            out+=' '.join(str(e) for e in out_seq_shape_list)+'\n'
#             f = open(".txt", "w+")
            with gzip.open(seq_file+'processed_'+enriched+shape+'.gz', 'wb') as f:
                f.write(out.encode())
#             f.write(out)
#             f.close()

#seq_file='/home/yibei/Downloads/yibei_predictions/'
All_Shapes=['Buckle-FL', 'Buckle', 'EP', 'HelT-FL', 'HelT', 'MGW-FL', 'MGW',
             'Opening-FL', 'Opening', 'ProT-FL', 'ProT', 'Rise-FL', 'Rise', 'Roll-FL',
             'Roll', 'Shear-FL', 'Shear', 'Shift-FL', 'Shift', 'Slide-FL', 'Slide',
             'Stagger-FL', 'Stagger', 'Stretch-FL', 'Stretch', 'Tilt-FL', 'Tilt']

seq_file='/project/rohs_108/yibeijia/data/yibei_predictions/'
preprocessShapeFile(seq_file, All_Shapes)


# In[ ]:




