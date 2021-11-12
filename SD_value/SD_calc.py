import pandas as pd
import numpy as np
import os.path
import time
import random
import math

def readInputAsArray(fileName):
    with open(fileName, 'r') as myfile:
        data = myfile.readlines()

    # Strip newline
    for i in range(0, len(data)):
        data[i] = data[i].rstrip()
    # print(data)
    return data

def readInputAsString(fileName):
    with open(fileName, 'r') as myfile:
        data=myfile.read().replace('\n', ' ')
        # print 
    return data

def reverse_complement(seq):    
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    reverse_complement = ''.join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement

def read_SD_table(fileName):
    tb_data = readInputAsArray(fileName)
    tb = {}
    DNAshape = tb_data[0].split("\t")
    # print(DNAshape)
    for i in range(1, len(tb_data)):
        line = tb_data[i].split("\t")
        key = line[0]
        value = line[1:]
        for j in range(0, len(value)):
            value[j] = float(value[j])
        tb[key] = value
    # print(tb)
    return(tb)

def SD_calculation(seq, tb, DNAshape):
    # print(tb.keys())
    shape_values = []

    #Initial the first pentamer
    pentamer = seq[0:5]
    if pentamer not in tb.keys():
        pentamer = reverse_complement(pentamer)
    for i in range(len(DNAshape)):
        shape_values.append(tb.get(pentamer)[i])

    for i in range(1, len(seq)-5):
        pentamer = seq[i:i+5]
        if pentamer not in tb.keys():
            pentamer = reverse_complement(pentamer)
        else:
            for j in range(0, len(DNAshape)):
                shape_values[j] = shape_values[j] + tb.get(pentamer)[j]
    return shape_values;

def main(fname, out):

    # Open and read file
    data = readInputAsArray(fname);

    fileName = "/Users/yibeijia/Downloads/ProjectStuff/mc_pentamer_06042020.txt"
    tb = read_SD_table(fileName);

    shapes = []
    for d in data:
        if 'N' not in d:
            shape_values = SD_calculation(d, tb, DNAshape)
            shapes.append(shape_values)

    string = ""
    for line in shapes:
        curr_str = " ".join(str(s) for s in line)
        string = string + curr_str + "\n"
    
    f = open(out, 'w+');
    f.write(string)
    f.write("\n")
    f.close()

DNAshape = ['Shear', 'Stretch', 'Stagger', 'Buckle', 'Propel', 'Opening', 
'Shift_-1', 'Shift_1', 'Slide_-1', 'Slide_1', 'Rise_-1', 'Rise_1', 'Tilt_-1', 'Tilt_1', 
'Roll_-1', 'Roll_1', 'Twist_-1', 'Twist_1', 'mingw', 'mingd', 'majgw', 'majgd', 
'Shear-SD', 'Stretch-SD', 'Stagger-SD', 'Buckle-SD', 'Propel-SD', 'Opening-SD', 
'Shift-SD', 'Slide-SD', 'Rise-SD', 'Tilt-SD', 'Roll-SD', 'Twist-SD', 'mingw-SD', 
'mingd-SD', 'majgw-SD', 'majgd-SD']

E_out = "Edata_SD_values.txt"
E_fname = '/Users/yibeijia/Downloads/ProjectStuff/E_sequences_to_save.txt'
D_out = "Ddata_SD_values.txt"
D_fname = '/Users/yibeijia/Downloads/ProjectStuff/D_sequences_to_save.txt'
# main(E_fname, E_out)
main(D_fname, D_out)