import numpy as np
import h5py
import scipy.io
from sklearn import metrics
import pandas as pd
import os
os.environ['THEANO_FLAGS'] = "device=cuda0,force_device=True,floatX=float32,gpuarray.preallocate=0.3"
import theano
import matplotlib.pyplot as plt

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

data_folder = "/scratch2/yibeijia/data/train_test_data/"
# data_folder = "/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/"

# trainmat = h5py.File(data_folder+'Etraining_data.mat')
trainmat = scipy.io.loadmat(data_folder+'Train_data.mat')
validmat = scipy.io.loadmat(data_folder+'Test_data.mat')

# X_train is the one hot encode #seq X 4
X_train = np.array(trainmat['Train_data'])
# y_train is
y_train = np.array(trainmat['Train_labels']).T

# trainmat.close()

# RUN TBINET
# sequence_input = Input
sequence_input = Input(shape=(147,4))

# Convolutional Layer
output = Conv1D(320,kernel_size=26,padding="valid",activation="relu")(sequence_input)
output = MaxPooling1D(pool_size=13, strides=13)(output)
output = Dropout(0.2)(output)

#Attention Layer
attention = Dense(1)(output)
attention = Permute((2, 1))(attention)
attention = Activation('softmax')(attention)
attention = Permute((2, 1))(attention)
attention = Lambda(lambda x: K.mean(x, axis=2), name='attention',output_shape=(75,))(attention)
attention = RepeatVector(320)(attention)
attention = Permute((2,1))(attention)
output = multiply([output, attention])

#BiLSTM Layer
output = Bidirectional(LSTM(320,return_sequences=True))(output)
output = Dropout(0.5)(output)

flat_output = Flatten()(output)

#FC Layer
FC_output = Dense(695)(flat_output)
FC_output = Activation('relu')(FC_output)

#Output Layer
output = Dense(1)(FC_output)
output = Activation('sigmoid')(output)

model = Model(inputs=sequence_input, outputs=output)

print('compiling model')
model.compile(loss='binary_crossentropy', optimizer='adam')

print('model summary')
model.summary()

# checkpointer = ModelCheckpoint(filepath="./model/tbinet_loss.temp_check.hdf5", verbose=1, save_best_only=False)
# checkpointer2 = ModelCheckpoint(filepath="./model/tbinet_acc.{epoch:02d}-{val_acc:.2f}.hdf5", monitor='val_acc',verbose=1, save_best_only=False)
checkpointer = ModelCheckpoint(filepath="./model/tbinet.{epoch:02d}-{val_loss:.2f}.hdf5", verbose=1, save_best_only=False)
earlystopper = EarlyStopping(monitor='val_loss', patience=10, verbose=1)

history = model.fit(X_train, y_train, batch_size=100, epochs=60, shuffle=True, verbose=1, validation_data=(validmat['Test_data'],validmat['Test_labels'].T), callbacks=[checkpointer,earlystopper])

# model.save('./model/tbinet_tmp.h5')
model.save('./model/tbinet.h5')

####### Learning curve
print(history.history.keys())
# summarize history for accuracy
# plt.plot(history.history['accuracy'])
# plt.plot(history.history['val_accuracy'])
# plt.title('model accuracy')
# plt.ylabel('accuracy')
# plt.xlabel('epoch')
# plt.legend(['train', 'test'], loc='upper left')
# plt.savefig("model_accuracy.png")
# summarize history for loss
plt.plot(history.history['loss'])
plt.plot(history.history['val_loss'])
plt.title('model loss')
plt.ylabel('loss')
plt.xlabel('epoch')
plt.legend(['train', 'test'], loc='upper left')
plt.savefig("model_loss.png")
