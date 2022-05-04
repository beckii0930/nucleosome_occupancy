import numpy as np
import scipy.io
from sklearn import metrics
import pandas as pd
import os


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

data_folder = "/project/rohs_108/yibeijia/nucleosome_occupancy/data/train_test_data/"

#	#testmat = scipy.io.loadmat(data_folder+'Test_data.mat')
#	total_sections = 5
#	y_test = np.array([])
#	for i in range(1, total_sections):
#		test_fn = 'yeastAll_Shapes_Test_5seqsPerClustr_'+str(i+1)+'_'+str(total_sections)+'.mat'
#		print(test_fn)
#		testmat = scipy.io.loadmat(data_folder+test_fn)
#		if y_test.shape[0] == 0:
#			y_test = np.array(testmat['Test_labels']).T
#			print('y test')
#			print(y_test)
#		else:
#			curr_y_test = np.array(testmat['Test_labels']).T
#			y_test = np.concatenate([y_test, curr_y_test], axis=0)


tpreds = np.load('tpreds_temp_run1.npy')
tpreds_temp = np.copy(tpreds)
for i in range(len(tpreds)):
	if tpreds[i] <= 0.1:
		print(f"original {tpreds[i]}")
		tpreds[i] = 0
		print(f"changed to {tpreds[i]}")
	if tpreds[i] >= 0.9:
		print(f"original {tpreds[i]}")
		tpreds[i] = 1
		print(f"changed to {tpreds[i]}")

y_test=np.load('y_test_run1.npy')
print(f"y_test.shape {y_test.shape}")

#for i in range(len(tpreds_temp)):
#    print("TrueLabel=%s, Predicted=%s" % (y_test.T[i], tpreds_temp[i]))
# 	print("TrueLabel=%s, Predicted=%s" % (y_test[i], tpreds_temp[i]))
reverse_start_id = int(y_test.shape[0]/2)
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
