#!/usr/bin/env python
# coding: utf-8

# In[35]:


import numpy as np
import scipy.io
from sklearn import metrics
import pandas as pd
import os
import scipy.io as sio
import os.path
import time
import random
import math
# os.environ['THEANO_FLAGS'] = "device=cuda0,force_device=True,floatX=float32"
# from aesara_theano_fallback import aesara as theano
#import theano

# print(theano.config.device)
# import sys, getopt
import tensorflow.keras as keras
import tensorflow as tf
import matplotlib.pyplot as plt

from keras.layers import Embedding
from keras.models import Model
from keras.layers import Dense, Dropout, Activation, Flatten, Layer, merge, Input, Concatenate, Reshape
from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.layers.pooling import GlobalMaxPooling1D
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import Bidirectional, TimeDistributed
from keras.models import load_model
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras import backend as K
import tensorflow.keras as keras
import tensorflow as tf
import matplotlib.pyplot as plt


# In[36]:


data_folder = "/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/model2/"

data_folder = "/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/model2_reduced_dense10/"
data_folder = "/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/model3/"
data_folder = "/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/train_dense10_kernel10_run10/"
model = load_model(data_folder+"tbinet.h5")
# model.load_weights("/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff//model2_reduced_dense10/tbinet.38-0.08.hdf5")
# model.load_weights("/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/model2/tbinet.22-0.21.hdf5")


# In[37]:


model.summary()


# In[38]:


attention_wt_model = tf.keras.models.Model(inputs=model.inputs, outputs=[model.output, model.get_layer('multiply').output])


# In[39]:


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

def oneHotEncode(seq):
    import numpy as np
    seq2=list()
    mapping = {"A":[1., 0., 0., 0.], "C": [0., 1., 0., 0.], "G": [0., 0., 1., 0.], "T":[0., 0., 0., 1.]};
    for i in seq:
    	seq2.append(mapping[i]  if i in mapping.keys() else [0., 0., 0., 0.]);
    return seq2;

def ReverseOneHotEncode(arr):
    seq = ''
    nuc_map = ['A', 'C', 'G', 'T']
    for subarr in arr:
        nuc = nuc_map[np.argmax(subarr)]
        seq+=nuc
    return seq;

def writeSeqToFile(d,ones,zeros):
    fname = '/Users/yibeijia/Downloads/nucleosome_occupancy/data/SeqsForFeature.txt'
    out = 'Enriched:\n'
    for i in ones:
        seq = ReverseOneHotEncode(d[i])
        out = out + str(1) + ' ' + seq + '\n'
    out += 'Depleted:\n'
    for i in zeros:
        seq = ReverseOneHotEncode(d[i])
        out = out + str(0) + ' ' + seq + '\n'
    with open(fname, 'w+') as f:
        f.write(out)


# In[40]:


## Set up a small set of data for predict call
data_folder = "/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/"
testmat = scipy.io.loadmat(data_folder+'Test_data.mat')
d = testmat['Test_data'][::1000]

labels = testmat['Test_labels']
new_labels = labels[0][::1000]
print(f"new_labels.shape: {new_labels.shape}")
ones = np.flatnonzero(new_labels == np.max(new_labels))
zeros = np.flatnonzero(new_labels == np.min(new_labels))
print(f'ones.shape {ones.shape}')
print(f'zeros.shape {zeros.shape}')

writeSeqToFile(d,ones,zeros)
# print(f"Total Test Seq is: {len(testmat['Test_data'])}, we take: {len(d)}")
attention_wts = attention_wt_model.predict(np.transpose(d,axes=(0,1,2)))


# In[41]:


model_outputs = attention_wts[0]
attention_outputs = attention_wts[1]


# In[42]:


r = len(attention_outputs)
c = len(attention_outputs[0])
l = len(attention_outputs[0][0])
print(f"The attention_outputs Array is {r} x {c} x {l}")


# In[43]:


attention_outputs[0]


# In[44]:


attention_outputs[0][-len(np.transpose(d,axes=(0,1,2))):]


# In[45]:


len(d[0])


# In[46]:


ReverseOneHotEncode(d[0])


# In[47]:


len(attention_outputs[0][-len(np.transpose(d,axes=(0,1,2))):])


# In[58]:


def tokenize(str):
    out = []
    for i in range(len(str)-4):
        out.append(str[i:i+5])
    return out


# In[76]:


class CharVal(object):
    def __init__(self, char, val):
        self.char = char
        self.val = val

    def __str__(self):
        return self.char

def rgb_to_hex(rgb):
    return '#%02x%02x%02x' % rgb

def color_charvals(s):
    r = 255-int(s.val*255)
    color = rgb_to_hex((255, r, r))
    return 'background-color: %s' % color


def plot_attention(token, currseq_attention_output):
    # match each 5mer with sum of the 320 dim array attention weight array
    char_vals = [CharVal(c, v) for c, v in zip(token, currseq_attention_output)]
    char_df = pd.DataFrame(char_vals).transpose()
    # apply coloring values
    char_df = char_df.style.applymap(color_charvals)
    return char_df


# In[78]:


# if you are using batches the outputs will be in batches
# get exact attentions of chars
i=0 #get attention map for ith sequence
currseq_attention_output = sum(attention_outputs[i])
token = tokenize(ReverseOneHotEncode(d[i]))
char_df = plot_attention(token, currseq_attention_output)
char_df


# In[82]:


mapped = dict(zip(token, currseq_attention_output))
print(mapped)


# In[115]:


# if you are using batches the outputs will be in batches
# get exact attentions of chars
mapped = {}
## Lots of sequences seq=1:151
for seq in range(len(attention_outputs)):
    tokens = tokenize(ReverseOneHotEncode(d[seq]))
    print(seq)
#     print(tokens)
    
    # for 147bp each sequence, there are 320 feature scores for each bp
    for i in range(len(attention_outputs[0])-4):
        currseq_attention_output = sum(attention_outputs[seq][i:i+5])
        print(currseq_attention_output)
        print("currseq_attention_output")
#         print(len(attention_outputs[0][i]))        
        for token in tokens:
            if token not in mapped.keys():
                mapped[token] = currseq_attention_output
            else:
                new_val = mapped[token] + currseq_attention_output
                mapped[token] = new_val
# print(mapped['ACGCG'])
            # char_df = plot_attention(token, currseq_attention_output)
# char_df


# In[ ]:




