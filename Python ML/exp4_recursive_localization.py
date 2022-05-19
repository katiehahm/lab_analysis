import numpy as np  
import pandas as pd
from sklearn.model_selection import KFold
from sklearn import datasets, svm, metrics, preprocessing
from scipy import stats 
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.ensemble import RandomForestRegressor
from tryGBR import GBR
from tryClassifierCV import tryClassifierCV, tryClassifierSingleCV
from GBRtest import GBRtest
from RFtest import RFtest
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import KFold

# used to run localization prediction from experiment 3
# takes in csv with localization features
# outputs a csv with predicted locations for every impact

# read csv file
data = pd.read_csv('C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/Jenny 1/ProcessedData/ExcelData/both_regular1_localization_p1_withta.csv',header=None)
# data
# data_labels = data.iloc[0]
# take out the first row
# data = data.drop(data.index[0])
data = data.apply(pd.to_numeric, errors='coerce')
# convert from pandas dataframe to numpy array
dataM = pd.DataFrame(data).to_numpy()
sh = dataM.shape
print(sh)
nrows = sh[0]
ncols = sh[1]

takes = dataM[:,0]
trainVal = dataM[:,1]

# GBR, 20 iter
model_g = GradientBoostingRegressor(n_estimators=50)
param_distributions_g = {
    'n_estimators': stats.randint(low=10,high=1000),
    'max_depth': stats.randint(low=2, high=6),
    'min_samples_split': stats.randint(low=2, high=5),
    'learning_rate': [0.9, 0.5, 0.25, 0.1, 0.05, 0.01]
}

finalPredictions = np.zeros((nrows,2)) # for every impact, (real coord, predict coord)

starts_train_idx = np.where(trainVal == 0)
starts_train_idx = np.array(starts_train_idx)
starts_trainX = dataM[starts_train_idx,4:28] # the rest of the array is just 0 filler
starts_trainY = dataM[starts_train_idx,2]
starts_trainX = starts_trainX[0,:,:]
starts_trainY = starts_trainY.flatten()

n_fold = 5
# k_fold = StratifiedKFold(n_splits=n_fold,shuffle=True,random_state=13)

# first run the estimation on segment starts separately
k_fold = KFold(n_splits=n_fold, shuffle=True,random_state=13)
for k, (train, test) in enumerate(k_fold.split(starts_trainX, starts_trainY)):
    print("Running starts kfold")
    xTrain = starts_trainX[train]
    yTrain = starts_trainY[train]
    xTest = starts_trainX[test]
    yTest = starts_trainY[test]

    model_cv = RandomizedSearchCV(model_g, param_distributions=param_distributions_g, n_iter=25, verbose=0)
    model_cv.fit(xTrain, yTrain)
    
    yPredict = model_cv.predict(xTest)
  
    # saving values to finalPredictions
    test_idx = starts_train_idx[:,test]
    test_idx = np.array(test_idx)
    finalPredictions[test_idx,:] = np.transpose([yTest, yPredict])

# then run the estimation using previous impact data
# 20% test-train, so run 5 times:
for i in range(5):
    # find all training data
    print("Running other data with i: {:.1f}".format(i+1))

    test_idx = np.where(trainVal == i+1) # 0 indexing
    xTest = dataM[test_idx,4:]
    yTest = dataM[test_idx,2]
    xTest = xTest[0,:,:]
    yTest = yTest.flatten()
    trainData = np.delete(dataM,test_idx,axis=0) # all rows that are not train value i
    xTrain = trainData[:,4:]
    yTrain = trainData[:,2]
    model_cv = RandomizedSearchCV(model_g, param_distributions=param_distributions_g, n_iter=25, verbose=0)
    model_cv.fit(xTrain, yTrain)
    for n in range(np.size(test_idx)):
        yPredict = model_cv.predict(xTest[n].reshape(1,-1))
        test_idx = np.array(test_idx, dtype=object)
        # a = np.transpose([yTest[n], yPredict[0]])
        finalPredictions[test_idx[0,n],:] = np.transpose([yTest[n], yPredict[0]])
        if n != np.size(test_idx)-1: # not the last element
            if (test_idx[0,n+1]-test_idx[0,n]) == 1: # continuation of same walking cycle
                nextFVector = xTest[n+1]
                xTest[n+1] = np.append(nextFVector[0:-1],yPredict) # use current predicted value on next impact FVector

yTestAll = finalPredictions[:,0]
yPredictAll = finalPredictions[:,1]
mse_all = mean_squared_error(yTestAll, yPredictAll)
print("The RMSE on all tests: {:.4f}".format(np.sqrt(mse_all)))

DF = pd.DataFrame(finalPredictions)
DF.to_csv("jenny1_trackingloc_regular1_p1_results.csv") # CHANGE ***********************************

fig = plt.figure(figsize=(6,6))
plt.plot(yPredictAll,yTestAll,'o')
plt.xlabel('Predicted Values [m]')
plt.ylabel('Target Values [m]')
plt.title('X coordinate Performance')
plt.xlim([-4,4])
plt.ylim([-4,4])
plt.show()



# old code

# train_idx = np.where(dataM[:,0]==0)
# test_idx = np.where(dataM[:,0]==1)

# train_idx = np.ravel(train_idx)
# test_idx = np.ravel(test_idx)

# # define inputs/outputs
# train_output = dataM[train_idx,1]
# train_input = dataM[train_idx,2:]

# test_output = dataM[test_idx,1]
# test_input = dataM[test_idx,2:22]

# # GBR, 20 iter
# model_g = GradientBoostingRegressor(n_estimators=50)
# param_distributions_g = {
#     'n_estimators': stats.randint(low=10,high=1000),
#     'max_depth': stats.randint(low=2, high=6),
#     'min_samples_split': stats.randint(low=2, high=5),
#     'learning_rate': [0.9, 0.5, 0.25, 0.1, 0.05, 0.01]
# }

# # RF, 50 iter
# model_r = RandomForestRegressor(n_estimators=200)
# param_distributions_r = {
#     'n_estimators': stats.randint(low=50, high=300),
#     'min_samples_split': stats.randint(low=2, high=6),
#     'min_samples_leaf': stats.randint(low=1, high=5),
#     'min_weight_fraction_leaf': [0.5, 0.25, 0.1, 0.05, 0.01, 0.0]
# }

# model = model_g
# param_distributions = param_distributions_g
# n_it = 35
# name_str = 'GBR'

# yTestAll = [[1]]
# yPredictAll = [[1]]

# model_cv = RandomizedSearchCV(model, param_distributions=param_distributions, n_iter=n_it, verbose=0)
# model_cv.fit(train_input, train_output)


# prev_coord = dataM[test_idx[0],1]
# for i in range(np.size(test_idx)):
#     if i > 0:
#         # print('test')
#         # print(test_idx(i))
#         # print(test_idx(i-1))
#         if test_idx[i]-test_idx[i-1] == 1: # still in current intervention
#             yTest = dataM[test_idx[i],1]
#             xTest = dataM[test_idx[i],2:22]
#             xTest = np.append(xTest,prev_coord)

#             yPredict = model_cv.predict(xTest.reshape(1,-1))
#             prev_coord = yPredict

#              # saving values to plot
#             yTestAll = np.append(yTestAll,yTest)
#             yPredictAll = np.append(yPredictAll,yPredict)
#         else:
#             yTest = dataM[test_idx[i],1]
#             xTest = dataM[test_idx[i],2:22]
#             prev_coord = dataM[test_idx[i]-1,1]
#             xTest = np.append(xTest,prev_coord)

#             yPredict = model_cv.predict(xTest.reshape(1,-1))
#             prev_coord = yPredict

#              # saving values to plot
#             yTestAll = np.append(yTestAll,yTest)
#             yPredictAll = np.append(yPredictAll,yPredict)


# yTestAll = np.delete(yTestAll,0,0)
# yPredictAll = np.delete(yPredictAll,0,0)

# print('#####')
# print(name_str + ' results:')
# mse_all = mean_squared_error(yTestAll, yPredictAll)
# print("The RMSE on all tests: {:.4f}".format(np.sqrt(mse_all)))

# forprint =  np.column_stack((test_idx[1:],yPredictAll))
# # forprint = np.column_stack((forprint,yTestAll))
# DF = pd.DataFrame(forprint)
# DF.to_csv("trackingloc_1stquart_results.csv")

# fig = plt.figure(figsize=(6,6))
# plt.plot(yPredictAll,yTestAll,'o')
# plt.xlabel('Predicted Values [m]')
# plt.ylabel('Target Values [m]')
# plt.title('X coordinate Performance')
# plt.xlim([-4,4])
# plt.ylim([-4,4])
# plt.show()