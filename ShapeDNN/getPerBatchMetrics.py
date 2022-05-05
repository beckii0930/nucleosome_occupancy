import numpy as np
import re
 
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

def readInputAsArray(fileName):
	with open(fileName, 'r') as myfile:
		data = myfile.readlines()
	for i in range(0, len(data)):
		data[i] = data[i].rstrip()
	return data
    
train_out = '/project/rohs_108/yibeijia/nucleosome_occupancy/ShapeDNN/train_ShapeDNN_run4.out'

train=readInputAsArray(train_out)   
start = False
Acc = []
Loss = []
Prec = []
Rec = []

for i in train:
	line = i.split()
	if start and len(line) > 0:
		if '1944' in line[0]:
			#print(line)
			for j in range(1,len(line)):
				if 'loss' in line[j-1]:
					Loss[epoch-1].append(float(line[j]))
				elif 'precision' in line[j-1]:
					Prec[epoch-1].append(float(line[j]))
				elif 'recall' in line[j-1]:
					Rec[epoch-1].append(float(line[j]))
				elif 'accuracy' in line[j-1]:
					acc_val = re.sub(r'[^\x20-\x7e]', '', line[j])
					Acc[epoch-1].append(float(acc_val))
					
	if len(line) == 2:
		if "Epoch" in line[0] and '/60' in line[1]:
			epoch = int(line[1].split('/')[0])
			print(f"epoch {epoch}")
			start=True
			Acc.append([])
			Loss.append([])
			Prec.append([])
			Rec.append([])

for epoch in range(len(Acc)):
	acc = Acc[epoch]
	loss = Loss[epoch]
	prec = Prec[epoch]
	rec = Rec[epoch]
	#print(len(acc))
	fig, axs = plt.subplots(2, 2)

	axs[0, 0].plot(acc)
	axs[0, 0].set_title('Accuracy')
	axs[0, 0].set_ylim([0, 1.1])
	axs[0, 0].xaxis.set_major_locator(MaxNLocator(integer=True))
	axs[0, 0].tick_params(axis='x', labelrotation = 45)

	axs[0, 1].plot(loss)
	axs[0, 1].set_title('Loss')
	axs[0, 1].set_ylim([0, 1.1])
	axs[0, 1].xaxis.set_major_locator(MaxNLocator(integer=True))
	axs[0, 1].tick_params(axis='x', labelrotation = 45)

	axs[1, 0].plot(prec)
	axs[1, 0].set_title('Precision')
	axs[1, 0].set_ylim([0, 1.1])
	axs[1, 0].xaxis.set_major_locator(MaxNLocator(integer=True))
	axs[1, 0].tick_params(axis='x', labelrotation = 45)

	axs[1, 1].plot(rec)
	axs[1, 1].set_title('Recall')
	axs[1, 1].set_ylim([0, 1.1])
	axs[1, 1].xaxis.set_major_locator(MaxNLocator(integer=True))
	axs[1, 1].tick_params(axis='x', labelrotation = 45)

	fig.supxlabel('Batch')
	fig.supylabel('Metric Value')
	fig.suptitle('Epoch '+str(epoch))
	plt.tight_layout()

	model_folder = '/home/yibei/Projects/nucleosome_occupancy/ShapeDNN/'
	fig_name="./batchwise_model_metrics_epoch_"+str(epoch)+".png"
	print("saving figure to:")
	print(fig_name)
	plt.savefig(fig_name)

#Acc_np = np.array(Acc)
#Loss_np = np.array(Loss)
#Prec_np = np.array(Prec)
#Rec_np = np.array(Rec)

#print(Acc_np)
#print(Loss_np)
#print(Prec_np)
#print(Rec_np)

#             print(score[0][1])
#            if (float(score[1][1:len(score[1])-1]) > 0.5):
#                pred += [1]
#            else:
#                pred += [0]  
#
#		from sklearn.metrics import accuracy_score
#		acc = accuracy_score(true, pred)
#		from sklearn.metrics import precision_score
#		prec = precision_score(true, pred)
#		from sklearn.metrics import recall_score
#		rec = recall_score(true, pred)
#
#		print(f"The accuracy_score: {acc}")
#		print(f"The precision_score: {prec}")
#		print(f"The recall_score: {rec}")
