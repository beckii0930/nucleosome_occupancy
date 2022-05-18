#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import h5py
import scipy.io
from sklearn import metrics
from sklearn import utils
import pandas as pd
import os
import sys
# from aesara_theano_fallback import aesara as theano
# os.environ['THEANO_FLAGS'] = "device=cuda0,force_device=True,floatX=float32,gpuarray.preallocate=\0.3"

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from keras.layers import Embedding
from keras.models import Sequential
from keras.models import Model
from keras.layers import Dense, Dropout, Activation, Flatten, Layer, merge, Input, Concatenate, Reshape, concatenate,Lambda,multiply,Permute,Reshape,RepeatVector
from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.layers.pooling import GlobalMaxPooling1D
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import Bidirectional, TimeDistributed
from keras.models import load_model
from keras.callbacks import ModelCheckpoint, EarlyStopping
from keras import optimizers
from keras import backend as K
from keras import regularizers
import tensorflow as tf
print("all packages loaded")
if len(sys.argv) != 2:
    print("python3 trainShapeDNN.py model_folder")
model_folder = sys.argv[1]

# In[2]:


import keras;
print(keras.__version__)


# ### Load data (training and validation)

# In[3]:


#data_folder = "/scratch2/yibeijia/data/train_test_data/"
# data_folder = "/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/"
data_folder = '/home/yibei/Projects/data/train_test_data/'
data_folder='/project/rohs_108/yibeijia/nucleosome_occupancy/data/train_test_data/'

# In[4]:


total_sections = 5
X_train_og = np.array([]) 
y_train_og = np.array([])
for i in range(total_sections):
    train_fn = 'yeastAll_Shapes_Train_5seqsPerClustr_'+str(i+1)+'_'+str(total_sections)+'.mat'
    trainmat = scipy.io.loadmat(data_folder+train_fn)
    if X_train_og.shape[0] == 0:
        X_train_og = np.array(trainmat['Train_data'])
        y_train_og = np.array(trainmat['Train_labels']).T
    else:
        curr_X_train = np.array(trainmat['Train_data'])
        curr_y_train = np.array(trainmat['Train_labels']).T
        X_train_og = np.concatenate([X_train_og, curr_X_train], axis=0)
        y_train_og = np.concatenate([y_train_og, curr_y_train], axis=0)
X_train, y_train = utils.shuffle(X_train_og, y_train_og)
#X_train, y_train = sklearn.utils.shuffle(X_train_og, y_train_og)
array_sum = np.sum(X_train)
array_has_nan = np.isnan(array_sum)
if array_has_nan == True:
	print("NaN in train")
	exit()

array_sum = np.sum(y_train)
array_has_nan = np.isnan(array_sum)
if array_has_nan == True:
	print("NaN in test")
	exit()
print("No NaN in train or test data")	

# first test mat is val and rest are test mats        
valid_fn = 'yeastAll_Shapes_Test_5seqsPerClustr_1_'+str(total_sections)+'.mat'
print(valid_fn)
validmat = scipy.io.loadmat(data_folder+valid_fn)
print("valid mat shape")
print(validmat['Test_data'].shape)
print(f"X_train.shape {X_train.shape}")
print(f"y_train.shape {y_train.shape}")

# ### Run TBiNet
sequence_input = Input(shape=(146,54))
sequence_input0 = tf.transpose(sequence_input, perm=[0, 2, 1])

# Convolutional Layer - shape
output0 = Conv1D(320,kernel_size=1,padding="same",activation="relu")(sequence_input0)
output0 = MaxPooling1D(pool_size=1, strides=1)(output0)
output0 = Dropout(0.2)(output0)

#Attention Layer 1 - shape
attention0 = Dense(1)(output0)
attention0 = Permute((2, 1))(attention0)
attention0 = Activation('softmax')(attention0)
attention0 = Permute((2, 1))(attention0)
attention0 = Lambda(lambda x: K.mean(x, axis=2), name='shape_attention',output_shape=(54,))(attention0)
attention0 = RepeatVector(146)(attention0)
attention0 = Permute((2,1))(attention0)
output = multiply([sequence_input0, attention0])

# Convolutional Layer 2 - seq
output= tf.transpose(output, perm=[0, 2, 1])
output = Conv1D(320,kernel_size=5,padding="same",activation="relu")(output)
output = MaxPooling1D(pool_size=5, strides=1)(output)
output = Dropout(0.2)(output)

#Attention Layer 2 - seq
attention = Dense(1)(output)
attention = Permute((2, 1))(attention)
attention = Activation('softmax')(attention)
attention = Permute((2, 1))(attention)
attention = Lambda(lambda x: K.mean(x, axis=2), name='seq_attention',output_shape=(75,))(attention)
attention = RepeatVector(320)(attention)
attention = Permute((2,1))(attention)
output = multiply([output, attention])

#BiLSTM Layer
output = Bidirectional(LSTM(320,return_sequences=True))(output)
output = Dropout(0.5)(output)
flat_output = Flatten()(output)

#FC Layer
FC_output = Dense(50)(flat_output)
FC_output = Activation('relu')(FC_output)

#Output Layer
output = Dense(1)(FC_output)
output = Activation('sigmoid')(output)

model = Model(inputs=sequence_input, outputs=output)

print('compiling model')
model.compile(loss='binary_crossentropy', optimizer='adam', metrics=[tf.keras.metrics.Precision(),
																	tf.keras.metrics.Recall(),
																	tf.keras.metrics.Accuracy()])

print('model summary')
model.summary()

checkpointer = ModelCheckpoint(filepath=model_folder+"tbinet.{epoch:02d}-{val_loss:.2f}.hdf5", verbose=1, save_best_only=False)
#earlystopper = EarlyStopping(monitor='val_loss', patience=10, verbose=1)
earlystopper = EarlyStopping(monitor='val_accuracy',patience=10, verbose=1)
reduce_lr = ReduceLROnPlateau(monitor='val_loss', factor=0.2,
                              patience=5, min_lr=0.001)

b=validmat['Test_labels'].T
X_val, y_val = utils.shuffle(validmat['Test_data'],validmat['Test_labels'])
array_sum = np.sum(X_val)
array_has_nan = np.isnan(array_sum)
if array_has_nan:
	print('Nan in X validation')
	exit()
array_sum = np.sum(y_val)
array_has_nan = np.isnan(array_sum)
if array_has_nan:
	print('Nan in y validation')
	exit()
print('No Nan in validation')
	
# model.fit(X_train, y_train, batch_size=100, epochs=60, shuffle=True, verbose=1, validation_data=(np.transpose(validmat['Train_data'],axes=(0,2,1)),validmat['Train_vals'][:,125:815]), callbacks=[checkpointer,earlystopper])

# added reduce learning rate callback to slow down learning rate
history = model.fit(X_train, y_train, batch_size=100, epochs=60, shuffle=True, verbose=1, validation_data=(X_val,y_val.T), callbacks=[checkpointer,earlystopper, reduce_lr])

model.save(model_folder+'tbinet.h5')

print("listing all data in history")
print(history.history.keys())
try:
	print("Loss")
	print(history.history['loss'])
	print("\n")
	print(history.history['val_loss'])
	
	print("accuracy")
	print(history.history['accuracy'])
	print("\n")
	print(history.history['val_accuracy'])

	print("precision")
	print(history.history['precision'])
	print("\n")
	print(history.history['val_precision'])

	print("recall")
	print(history.history['recall'])
	print("\n")
	print(history.history['val_recall'])

	fig, axs = plt.subplots(2, 2)
	axs[0, 0].plot(history.history['accuracy'], label='train')
	axs[0, 0].plot(history.history['val_accuracy'], label='test')
	axs[0, 0].set_title('Accuracy')
	axs[0, 0].legend(['train', 'test'], loc='upper right',prop={'size': 8})
	axs[0, 0].set_ylim([0, 1.1])
	axs[0, 0].xaxis.set_major_locator(MaxNLocator(integer=True))

	axs[0, 1].plot(history.history['loss'], label='train')
	axs[0, 1].plot(history.history['val_loss'], label='test')
	axs[0, 1].set_title('Loss')
	axs[0, 1].legend(['train', 'test'], loc='upper right',prop={'size': 8})
	axs[0, 1].set_ylim([0, 1.1])
	axs[0, 1].xaxis.set_major_locator(MaxNLocator(integer=True))

	axs[1, 0].plot(history.history['precision'], label='train')
	axs[1, 0].plot(history.history['val_precision'], label='test')
	axs[1, 0].set_title('Precision')
	axs[1, 0].legend(['train', 'test'], loc='lower right',prop={'size': 8})
	axs[1, 0].set_ylim([0, 1.1])
	axs[1, 0].xaxis.set_major_locator(MaxNLocator(integer=True))

	axs[1, 1].plot(history.history['recall'], label='train')
	axs[1, 1].plot(history.history['val_recall'], label='test')
	axs[1, 1].set_title('Recall')
	axs[1, 1].legend(['train', 'test'], loc='lower right',prop={'size': 8})
	axs[1, 1].set_ylim([0, 1.1])
	axs[1, 1].xaxis.set_major_locator(MaxNLocator(integer=True))
	
	fig.supxlabel('Epoch')
	fig.supylabel('Metric Value')
	plt.tight_layout()
	plt.savefig(model_folder+"model_metrics.png")

	# Hide x labels and tick labels for top plots and y ticks for right plots.
	#for ax in axs.flat:
	#    ax.label_outer()
except:
	print("An exception occurred") 

#plt.plot(history.history['accuracy'])
#plt.plot(history.history['val_accuracy'])
#plt.title('model accuracy')
#plt.ylabel('accuracy')
#
#plt.plot(history.history['loss'], label='train')
#plt.plot(history.history['val_loss'], label='test')
#plt.title('model loss')
#plt.ylabel('loss')
#plt.xlabel('epoch')
#plt.legend(['train', 'test'], loc='upper left')
#plt.savefig(model_folder+"model_loss.png")
