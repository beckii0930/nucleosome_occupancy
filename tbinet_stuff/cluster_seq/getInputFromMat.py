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

N=1
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

def oneHotDecode(seq):
    import numpy as np
    seq2=""
    mapping = {"A":[1., 0., 0., 0.], "C": [0., 1., 0., 0.], "G": [0., 0., 1., 0.], "T":[0., 0., 0., 1.]};
    for i in seq:
        if i[0] == 1.:
            seq2+='A'
        elif i[1] == 1.:
            seq2+='C'
        elif i[2] == 1.:
            seq2+='G'
        elif i[3] == 1.:
            seq2+='T'
        else:
            seq2+='N'
    return seq2;


print("Enriched list size")
print(Edata_list.shape)    
print("Depleted list size")
print(Ddata_list.shape)    
toc=time.perf_counter()
print(f"Create master list took {toc - tic:0.4f} seconds")
seq=""

# Opening a file

file1 = open('Eseqs.txt', 'w')
i = 0;
for seq_arr in Edata_list:
    seq = oneHotDecode(seq_arr)
    file1.write("\n>" + str(i) + "\n")
    i+=1;
    file1.write(seq)
    
# Closing file
file1.close()

# Opening a file
file2 = open('Dseqs.txt', 'w')
i = 0;
for seq_arr in Ddata_list:
    seq = oneHotDecode(seq_arr)
    file2.write("\n>" + str(i) + "\n")
    i+=1;
    file2.write(seq)
    
# Closing file
file2.close()