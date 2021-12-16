import numpy as np
import scipy.io
from sklearn import metrics
import pandas as pd
import os
# os.environ['THEANO_FLAGS'] = "device=cuda0,force_device=True,floatX=float32"
# from aesara_theano_fallback import aesara as theano
#import theano

# print(theano.config.device)
# import sys, getopt

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

data_folder = "/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/model2/"
model = load_model(data_folder+"tbinet.h5")
model.load_weights("/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/model2/tbinet.12-0.12.hdf5")

print("model.summary()")
print(model.summary())

layer_outputs = [layer.output for layer in model.layers]

# Run the model to predict once
feature_map_model = tf.keras.models.Model(inputs=model.inputs, outputs=layer_outputs)

# data_folder = "./data/"
data_folder = "../data/train_test_data/"
testmat = scipy.io.loadmat(data_folder+'Test_data.mat')
feature_maps = feature_map_model.predict(np.transpose(testmat['Test_data'],axes=(0,1,2)))

labels = testmat['Test_labels']
ones = np.flatnonzero(labels == np.max(labels)).T
zeros = np.flatnonzero(labels == np.min(labels)).T
print(f'ones.shape {ones.shape}')
print(f'zeros.shape {zeros.shape}')

def plot_feature_map(y, plt_title):
    plt.rcParams["figure.figsize"] = 5,2
    x = range(1,10)

    fig, (ax,ax2) = plt.subplots(nrows=2, sharex=True)

    extent = [x[0]-(x[1]-x[0])/2., x[-1]+(x[1]-x[0])/2.,0,1]
    ax.imshow(y[np.newaxis,:], cmap="plasma", aspect="auto", extent=extent)
    ax.set_yticks([])
    ax.set_xlim(extent[0], extent[1])

    ax2.plot(x,y)
    fig.suptitle(plt_title)
    plt.tight_layout()
    plt.show()

def plot_both_feature_maps(y0, y1, plt_title):
    plt.rcParams["figure.figsize"] = 5,2
    x = range(1,10)
    data = np.concatenate(([y0], [y1]))
    print(data)
    fig, axis = plt.subplots()

    extent = [x[0]-(x[1]-x[0])/2., x[-1]+(x[1]-x[0])/2.,0,1]
    heatmap = axis.pcolor(data, cmap=plt.cm.Blues) 
    axis.set_yticks(np.arange(data.shape[0])+0.5, minor=False)
    axis.set_xticks(np.arange(data.shape[1])+0.5, minor=False)
    
    row_labels = ["Depleted - 0", "Enriched - 1"]
    axis.set_yticklabels(row_labels, minor=False)
    axis.set_xticklabels(x, minor=False)
    fig.suptitle(plt_title)
    plt.colorbar(heatmap)
    plt.tight_layout()
    plt.show()


for layer_name, feature_map in zip(layer_names, feature_maps):

    if layer_name == 'attention':
        print(feature_map)
        feature_map -= feature_map.mean()
        feature_map /= feature_map.std ()
        print(f'layer_name: {layer_name}')
        print(f'feature_map.shape: {feature_map.shape}')
        feature_map_ones = feature_map[ones]
        feature_map_zeros = feature_map[zeros]
        print(f'feature_map_ones.shape: {feature_map_ones.shape}')
        print(f'feature_map_zeros.shape: {feature_map_zeros.shape}')
        feature_map_avg_ones = feature_map_ones.sum(axis=0)
        feature_map_avg_zeros = feature_map_zeros.sum(axis=0)
        plot_feature_map(feature_map_avg_ones, "Feature Map Avg Ones")
        plot_feature_map(feature_map_avg_zeros, "Feature Map Avg Zeros")
        plot_both_feature_maps(feature_map_avg_ones, feature_map_avg_zeros,"Feature Maps") 
