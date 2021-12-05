## TODO: get the layers from the model.
## then print and visualize
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

def layer_to_visualize(model, layer):
    inputs = [K.learning_phase()] + model.inputs

    _convout1_f = K.function(inputs, [layer.output])
    def convout1_f(X):
        # The [0] is to disable the training phase flag
        return _convout1_f([0] + [X])

    convolutions = convout1_f(img_to_visualize)
    convolutions = np.squeeze(convolutions)

    print ('Shape of conv:', convolutions.shape)

    n = convolutions.shape[0]
    n = int(np.ceil(np.sqrt(n)))

    # Visualization of each filter of the layer
    fig = plt.figure(figsize=(12,8))
    for i in range(len(convolutions)):
        ax = fig.add_subplot(n,n,i+1)
        ax.imshow(convolutions[i], cmap='gray')

data_folder = "/Users/yibeijia/Downloads/nucleosome_occupancy/tbinet_stuff/model/"
#data_folder = "/scratch2/yibeijia/nucleosome_occupancy/tbinet_stuff/model/"
model = load_model(data_folder+"tbinet.h5")
keras.utils.plot_model(model, show_shapes=True, dpi=90)
# choose any image to want by specifying the index
img_to_visualize = X_train[65]

# To understand how the attention has helped improve the performance 
# and interpretability of TF-DNA binding prediction, 
# we visualized the attention scores generated from TBiNet

# Conv1D layer

# attention layer: get attentino scores from tbinet
# attention vector ,size 75, for each DNA seq
# layer_to_visualize(model, convout1)