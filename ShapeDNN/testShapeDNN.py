import numpy as np
import scipy.io
from sklearn import metrics
from sklearn import utils
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
data_folder = "/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/"
data_folder = "/scratch2/yibeijia/data/train_test_data/"
data_folder = "/Users/yibeijia/Downloads/nucleosome_occupancy/data/train_test_data/"
data_folder = "/project/rohs_108/yibeijia/nucleosome_occupancy/data/train_test_data/"


#testmat = scipy.io.loadmat(data_folder+'Test_data.mat')
total_sections = 5
X_test_og = np.array([])
y_test_og = np.array([])
for i in range(1, total_sections):
	test_fn = 'yeastAll_Shapes_Test_5seqsPerClustr_'+str(i+1)+'_'+str(total_sections)+'.mat'
	print(test_fn)
	testmat = scipy.io.loadmat(data_folder+test_fn)
	if X_test_og.shape[0] == 0:
		X_test_og = np.array(testmat['Test_data'])
		y_test_og = np.array(testmat['Test_labels']).T
		print('Xtest')
		print(X_test_og)
		print('y test')
		print(y_test_og)
	else:
		curr_X_test = np.array(testmat['Test_data'])
		curr_y_test = np.array(testmat['Test_labels']).T
		X_test_og = np.concatenate([X_test_og, curr_X_test], axis=0)
		y_test_og = np.concatenate([y_test_og, curr_y_test], axis=0)
X_test, y_test = utils.shuffle(X_test_og,y_test_og)

print(f"X_test.shape {X_test.shape}")


# In[7]:


print(f"y_test.shape {y_test.shape}")

### Load model
argv = sys.argv[1:]
try:
    opts, args = getopt.getopt(argv,"hi:",["ifile="])

except getopt.GetoptError:
    print ('test.py -i <inputfile>')
    # print ('test.py -i <inputfile> -o <outputfile>')
    sys.exit(2)
for opt, arg in opts:
    if opt == '-h':
        print ('test.py -i <inputfile>')
        # print ('test.py -i <inputfile> -o <outputfile>')
        sys.exit()
    elif opt in ("-i", "--ifile"):
        inputmodel = arg
    # elif opt in ("-o", "--ofile"):
    #     output = arg
print ('Input file is ', inputmodel)
# print ('Output file is "', output)

model = load_model(inputmodel+"tbinet.h5")
#model = load_model("./model2_reduced/tbinet.h5")
#model = load_model("./model2_reduced_dense10/tbinet.h5")
#model = load_model("./model3_removeLSTM_addDense/tbinet.h5")
#model.load_weights(inputmodel)
print('model summary')
model.summary()

### Calculate averaged AUROC and AUPR
#tpreds = model.predict_classes(testmat['Test_data'])
tpreds = model.predict(X_test,verbose=1)

#tpreds = model.predict(np.transpose(X_test,axes=(0,2,1)),verbose=1)
tpreds_temp = np.copy(tpreds)
print(f"tpreds_temp.shape {tpreds_temp.shape}")
print(f"y_test.shape {y_test.shape}")


with open(inputmodel+'tpreds_temp_run1.npy', 'wb') as f:
	np.save(f, tpreds_temp)

with open(inputmodel+'y_test_run1.npy', 'wb') as f:
	np.save(f, y_test)

for i in range(len(tpreds)):
    if tpreds[i] <= 0.1:
#        print(f"original {tpreds[i]}")
        tpreds[i] = 0
#        print(f"changed to {tpreds[i]}")
    if tpreds[i] >= 0.9:
#        print(f"original {tpreds[i]}")
        tpreds[i] = 1
#        print(f"changed to {tpreds[i]}")

reverse_start_id = int(y_test.shape[0]/2)

#for i in range(len(tpreds)):
#    print("TrueLabel=%s, Predicted=%s" % (y_test.T[i], tpreds_temp[i]))
#	print("TrueLabel=%s, Predicted=%s" % (y_test[i], tpreds_temp[i]))

for i in range(reverse_start_id):
    tpreds_avg_temp = (tpreds_temp[i] + tpreds_temp[reverse_start_id+i])/2.0
    tpreds_temp[i] = tpreds_avg_temp
    tpreds_temp[reverse_start_id+i] = tpreds_avg_temp

aurocs, auprs = get_aurocs_and_auprs(tpreds_temp,y_test)
#aurocs, auprs = get_aurocs_and_auprs(tpreds_temp,y_test.T)
print("Averaged AUROC:",np.nanmean(aurocs))
print("Averaged AUPR:", np.nanmean(auprs))
from sklearn.metrics import accuracy_score
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score 
print("accuracy_score", accuracy_score(y_test,tpreds))
print("precision_score", precision_score(y_test,tpreds))
print("recall_score", recall_score(y_test,tpreds))
