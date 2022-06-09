import numpy as np
import scipy.io
import time

data_folder='/project/rohs_108/yibeijia/nucleosome_occupancy/data/train_test_data/'
total_sections = 5
X_train_og = np.array([])
y_train_og = np.array([])

print("loading train data")
tic = time.perf_counter()

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
toc = time.perf_counter()
print(f"loaded train data and took {toc - tic:0.4f} seconds")
print("loading test  data")
tic = time.perf_counter()

X_test_og = np.array([])
y_test_og = np.array([])
for i in range(total_sections):
    test_fn = 'yeastAll_Shapes_Test_5seqsPerClustr_'+str(i+1)+'_'+str(total_sections)+'.mat'
    testmat = scipy.io.loadmat(data_folder+test_fn)
    if X_test_og.shape[0] == 0:
        X_test_og = np.array(testmat['Test_data'])
        y_test_og = np.array(testmat['Test_labels']).T
    else:
        curr_X_test = np.array(testmat['Test_data'])
        curr_y_test = np.array(testmat['Test_labels']).T
        X_test_og = np.concatenate([X_test_og, curr_X_test], axis=0)
        y_test_og = np.concatenate([y_test_og, curr_y_test], axis=0)
toc = time.perf_counter()
print(f"loaded test data and took {toc - tic:0.4f} seconds")

print(f"X_train_og.shape {X_train_og.shape}")
print(f"y_train_og.shape {y_train_og.shape}")

for i in range(len(X_train_og[0])):
	print(X_train_og[0][i][0])

for i in range(len(X_test_og[0])):
	print(X_test_og[0][i][0])

#print(min(X_train_og))
#print(max(X_train_og))
#print(mean(X_train_og))
print(f"X_test_og.shape {X_test_og.shape}")
print(f"y_test_og.shape {y_test_og.shape}")
#print(min(X_test_og))
#print(max(X_test_og))
#print(mean(X_test_og))
