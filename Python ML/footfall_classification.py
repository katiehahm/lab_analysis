# 6/6/22
# used to classify if footfall falls into person 1 or person 2 category
# decision tree works the best

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn import tree
from sklearn.tree import DecisionTreeClassifier
import matplotlib.image as pltimg
from sklearn import svm
from sklearn.ensemble import GradientBoostingClassifier

data = pd.read_csv('C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/April 3/ProcessedData/ExcelData/alltakes_localization_p1_withta.csv',header=None)
testdata = pd.read_csv('C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/April 3/ProcessedData/ExcelData/both_weight2_testing_060622.csv',header=None)
data2 = pd.read_csv('C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/April 3/ProcessedData/ExcelData/alltakes_localization_p2_withta.csv',header=None)
savefile = "both_weight2_testing_060622_results.csv"


data = data.apply(pd.to_numeric, errors='coerce')
dataM = pd.DataFrame(data).to_numpy()
sh = dataM.shape
print(sh)
nrows = sh[0]
ncols = sh[1]
xData = dataM[:,4:27]
print(xData.shape)
yData = np.full((1, nrows), 1)

testdata = testdata.apply(pd.to_numeric, errors='coerce')
testdataM = pd.DataFrame(testdata).to_numpy()
testsh = testdataM.shape
print(testsh)
nTestrows = testsh[0]
nTestcols = testsh[1]
X_test = testdataM[:,1:24]
print(X_test.shape)
y_test = np.array([1,2,1,1,2,1,2,1,2,1,2,1,2,1,1,1,2])
y_test = y_test.flatten()

data2 = data2.apply(pd.to_numeric, errors='coerce')
dataM2 = pd.DataFrame(data2).to_numpy()
sh2 = dataM2.shape
print(sh2)
nrows2 = sh2[0]
ncols2 = sh2[1]
xData2 = dataM2[:,4:27]
print(xData2.shape)
yData2 = np.full((1, nrows2), 2)

X_train = np.concatenate((xData,xData2),axis=0)
y_train = np.concatenate((yData,yData2),axis=1)
y_train = y_train.flatten()

# Feature Scaling
# from sklearn.preprocessing import StandardScaler
# sc = StandardScaler()
# X_train = sc.fit_transform(X_train)
# X_test = sc.transform(X_test)

# Training the Naive Bayes model on the Training set
from sklearn.naive_bayes import GaussianNB
classifier = GaussianNB()
classifier.fit(X_train, y_train)

# Predicting the Test set results
y_pred = classifier.predict(X_test)

print("naive bayes")
print(y_pred)
print(y_test)

# Making the Confusion Matrix
# from sklearn.metrics import confusion_matrix, accuracy_score
# ac = accuracy_score(y_test,y_pred)
# cm = confusion_matrix(y_test, y_pred)

print("decision tree")
dtree = DecisionTreeClassifier()
dtree = dtree.fit(X_train, y_train)
print(dtree.predict(X_test))
print(y_test)


clf = svm.SVC()
clf.fit(X_train,y_train)
print("svm")
print(clf.predict(X_test))
print(y_test)

lr_list = [0.05, 0.075, 0.1, 0.25, 0.5, 0.75, 1]
for learning_rate in lr_list:
    gb_clf = GradientBoostingClassifier(n_estimators=20, learning_rate=learning_rate, max_features=2, max_depth=2, random_state=0)
    print("GBC learning rate: ", learning_rate)
    gb_clf.fit(X_train, y_train)
    print(gb_clf.predict(X_test))
    print(y_test)