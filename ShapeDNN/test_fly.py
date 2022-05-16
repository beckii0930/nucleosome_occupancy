import numpy as np
import scipy.io
from sklearn import metrics
import pandas as pd
import os
os.environ['THEANO_FLAGS'] = "device=cuda0,force_device=True,floatX=float32"
from aesara_theano_fallback import aesara as theano
#import theano

print(theano.config.device)
import sys, getopt

from keras.layers import Embedding
from keras.models import Sequential
from keras.models import Model
from keras.layers import Dense, Dropout, Activation, Flatten, Layer, merge, Input, Concatenate, Reshape
from keras.layers.convolutional import Conv1D, MaxPooling1D
from keras.layers.pooling import GlobalMaxPooling1D
from keras.layers.recurrent import LSTM
from keras.layers.wrappers import Bidirectional, TimeDistributed
from keras.models import load_model
from keras.callbacks import ModelCheckpoint, EarlyStopping

def get_auroc(preds, obs):
    fpr, tpr, thresholds  = metrics.roc_curve(obs, preds, drop_intermediate=False)
    auroc = metrics.auc(fpr,tpr)
    return auroc

def get_aupr(preds, obs):
    precision, recall, thresholds  = metrics.precision_recall_curve(obs, preds)
    aupr = metrics.auc(recall,precision)
    return aupr

def get_aurocs_and_auprs(tpreds, tobs):
    tpreds_df = pd.DataFrame(tpreds)
    tobs_df = pd.DataFrame(tobs)
    
    task_list = []
    auroc_list = []
    aupr_list = []
    for task in tpreds_df:
        pred = tpreds_df[task]
        obs = tobs_df[task]
        auroc=round(get_auroc(pred,obs),5)
        aupr = round(get_aupr(pred,obs),5)
        task_list.append(task)
        auroc_list.append(auroc)
        aupr_list.append(aupr)
    return auroc_list, aupr_list

### Load data (test)
#data_folder = "/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/"
data_folder = "/project/rohs_108/yibeijia/data/"
#data_folder = "/scratch2/yibeijia/data/train_test_data/"
testmat = scipy.io.loadmat(data_folder+'fly_Test_data.mat')

### Load model
if len(sys.argv) != 2:
	print("python3 test.py model_folder")
model_folder = sys.argv[1]

#argv = sys.argv[1:]
#try:
#    opts, args = getopt.getopt(argv,"hi:",["ifile="])

#except getopt.GetoptError:
#    print ('test.py -i <inputfile>')
    # print ('test.py -i <inputfile> -o <outputfile>')
#    sys.exit(2)
#for opt, arg in opts:
#    if opt == '-h':
#        print ('test.py -i <inputfile>')
        # print ('test.py -i <inputfile> -o <outputfile>')
#        sys.exit()
#    elif opt in ("-i", "--ifile"):
#        inputmodel = arg
    # elif opt in ("-o", "--ofile"):
    #     output = arg
#print ('Input file is "', inputmodel)
# print ('Output file is "', output)

#model = load_model("./model/tbinet.h5")
#model = load_model("./model2_reduced/tbinet.h5")
#model = load_model("./model2_reduced_dense10/tbinet.h5")
model = load_model(model_folder+"tbinet.h5")
#model.load_weights(inputmodel)
print('model summary')
model.summary()

### Calculate averaged AUROC and AUPR
#tpreds = model.predict_classes(testmat['Test_data'])
tpreds = model.predict(testmat['Test_data'],verbose=1)
tpreds_temp = np.copy(tpreds)
reverse_start_id = int(testmat['Test_labels'].shape[0]/2)

for i in range(len(tpreds)):
    print("TrueLabel=%s, Predicted=%s" % (testmat['Test_labels'].T[i], tpreds_temp[i]))

for i in range(reverse_start_id):
    tpreds_avg_temp = (tpreds_temp[i] + tpreds_temp[reverse_start_id+i])/2.0
    tpreds_temp[i] = tpreds_avg_temp
    tpreds_temp[reverse_start_id+i] = tpreds_avg_temp


aurocs, auprs = get_aurocs_and_auprs(tpreds_temp,testmat['Test_labels'].T)
print("Averaged AUROC:",np.nanmean(aurocs))
print("Averaged AUPR:", np.nanmean(auprs))
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score 
print("accuracy_score", accuracy_score(tpreds_temp,testmat['Test_labels'].T))
print("precision_score", precision_score(tpreds_temp,testmat['Test_labels'].T))
print("recall_score", recall_score(tpreds_temp,testmat['Test_labels'].T))
