import sys

def readInputAsArray(fileName):
	with open(fileName, 'r') as myfile:
		data = myfile.readlines()
	for i in range(0, len(data)):
		data[i] = data[i].rstrip()
	return data

test_out = sys.argv[1]
test=readInputAsArray(test_out) 

true = []
pred = []

for i in test:
	line = i.split()
	for k in line:
		score = k.split("=")
		if score[0]=="TrueLabel":
			true += [int(score[1][1])]
		if score[0]=="Predicted":
			if (float(score[1][1:len(score[1])-1]) > 0.5):
				pred += [1]
			else:
				pred += [0]
print(pred[:5])
print(true[:5])


from sklearn.metrics import accuracy_score
acc = accuracy_score(true, pred)
from sklearn.metrics import precision_score
prec = precision_score(true, pred)
from sklearn.metrics import recall_score
rec = recall_score(true, pred)
print(f"accuracy_score is: {acc}")
print(f"precision_score is: {prec}")
print(f"recall_score is: {rec}")
