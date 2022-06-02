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

path='~/Downloads/data/'
Eseqs = readInputAsArray(path+'Eseqs.txt')[1:] # avoid the first \n
Ecluster = readInputAsArray(path+'Eseqs80.clstr')

Dcluster = readInputAsArray(path+'Dseqs80.clstr')
Dseqs = readInputAsArray(path+'Dseqs.txt')[1:]