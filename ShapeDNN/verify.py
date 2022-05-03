import numpy as np
import h5py
import scipy.io
from sklearn import metrics
import pandas as pd
import os
import sys
data_folder='/project/rohs_108/yibeijia/nucleosome_occupancy/data/train_test_data/'
total_sections = 5
X_test = np.array([])
y_test = np.array([])
for i in range(total_sections):
	test_fn = 'yeastAll_Shapes_Test_5seqsPerClustr_'+str(i+1)+'_'+str(total_sections)+'.mat'
	print(test_fn)
	testmat = scipy.io.loadmat(data_folder+test_fn)
	if X_test.shape[0] == 0:
		X_test = np.array(testmat['Test_data'])
		y_test = np.array(testmat['Test_labels']).T
	else:
		curr_X_test = np.array(testmat['Test_data'])
		curr_y_test = np.array(testmat['Test_labels']).T
		X_test = np.concatenate([X_test, curr_X_test], axis=0)
		y_test = np.concatenate([y_test, curr_y_test], axis=0)

print('\n')
print(f"X_test.shape {X_test.shape}")
print(f"y_test.shape {y_test.shape}")

zero_cnt=0
one_cnt=0
for i in y_test:
	if i ==0:
		zero_cnt +=1
	else:
		one_cnt+=1
print("IN test")
print(f"zero cnt is {zero_cnt}")
print(f"one cnt is {one_cnt}")
print(f"total cnt is {zero_cnt+one_cnt}")


#		total_sections = 5
#		X_train = np.array([]) 
#		y_train = np.array([])
#		for i in range(total_sections):
#			train_fn = 'yeastAll_Shapes_Train_5seqsPerClustr_'+str(i+1)+'_'+str(total_sections)+'.mat'
#			trainmat = scipy.io.loadmat(data_folder+train_fn)
#			if X_train.shape[0] == 0:
#				X_train = np.array(trainmat['Train_data'])
#				y_train = np.array(trainmat['Train_labels']).T
#			else:
#				curr_X_train = np.array(trainmat['Train_data'])
#				curr_y_train = np.array(trainmat['Train_labels']).T
#				X_train = np.concatenate([X_train, curr_X_train], axis=0)
#				y_train = np.concatenate([y_train, curr_y_train], axis=0)
#		zero_cnt=0
#		one_cnt=0
#		for i in y_train:
#			if i ==0:
#				zero_cnt +=1
#			else:
#				one_cnt+=1
#		print("IN train")
#		print(f"zero cnt is {zero_cnt}")
#		print(f"one cnt is {one_cnt}")
#		print(f"total cnt is {zero_cnt+one_cnt}")

