import numpy as np

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
			for j in range(1,len(line)):
				if 'loss' in line[j-1]:
					print('loss')
					Loss[epoch-1].append(float(line[j].split('x08')[0]))
				elif 'precision' in line[j-1]:
					print('precision')
					Prec[epoch-1].append(float(line[j].split('\x')[0]))
				elif 'recall' in line[j-1]:
					print('recall')
					Rec[epoch-1].append(float(line[j].split('\x')[0]))
				elif 'accuracy' in line[j-1]:
					print('accuracy')
					Acc[epoch-1].append(float(line[j].split('\x')[0]))
					
	if len(line) == 2:
		if "Epoch" in line[0] and '/60' in line[1]:
			epoch = int(line[1].split('/')[0])
			print(f"epoch {epoch}")
			start=True
			Acc.append([])
			Loss.append([])
			Prec.append([])
			Rec.append([])

Acc_np = np.array(Acc)
Loss_np = np.array(Loss)
Prec_np = np.array(Prec)
Rec_np = np.array(Rec)

print(Acc_np)
print(Loss_np)
print(Prec_np)
print(Rec_np)

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
