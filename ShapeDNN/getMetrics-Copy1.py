def readInputAsArray(fileName):
    with open(fileName, 'r') as myfile:
        data = myfile.readlines()

    # Strip newline
    for i in range(0, len(data)):
        data[i] = data[i].rstrip()
    # print(data)
    return data
    
test_out = '/Users/yibeijia/Downloads/data/test_model2_REDUCED_dense10.out'
test_out = '/Users/yibeijia/Downloads/data/test_model3_removeLSTM.out'
test=readInputAsArray(test_out)   

true = []
pred = []

for i in test:
    line = i.split()
    for k in line:
        score = k.split("=")
        if score[0]=="TrueLabel":
#             print(score[1][1])
            true += [int(score[1][1])]
        if score[0]=="Predicted":
#             print(score[0][1])
            if (float(score[1][1:len(score[1])-1]) > 0.5):
                pred += [1]
            else:
                pred += [0]  

from sklearn.metrics import accuracy_score
acc = accuracy_score(true, pred)
from sklearn.metrics import precision_score
prec = precision_score(true, pred)
from sklearn.metrics import recall_score
rec = recall_score(true, pred)

print(f"The accuracy_score: {acc}")
print(f"The precision_score: {prec}")
print(f"The recall_score: {rec}")