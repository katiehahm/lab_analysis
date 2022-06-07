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
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler

# used to run localization prediction from experiment 3
# takes in csv with localization features
# outputs a csv with predicted locations for every impact

##################################################
# copied from exp4_notrecurisve_localization.py
# for using localization to redo footfall labeling

# read csv file

# data = pd.read_csv('C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/Jenny 1/ProcessedData/ExcelData/both_regular2_localization_p2_withta.csv',header=None)
# savefile = "both_regular2_localization_p2_results.csv"

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
yData = dataM[:,2]
yData = yData.flatten()

testdata = testdata.apply(pd.to_numeric, errors='coerce')
testdataM = pd.DataFrame(testdata).to_numpy()
testsh = testdataM.shape
print(testsh)
nTestrows = testsh[0]
nTestcols = testsh[1]
xTestData = testdataM[:,1:24]
print(xTestData.shape)
yTestData = testdataM[:,0]
yTestData = yTestData.flatten()

# GBR, 20 iter
model_g = GradientBoostingRegressor(n_estimators=50)
param_distributions_g = {
    'n_estimators': stats.randint(low=10,high=1000),
    'max_depth': stats.randint(low=2, high=6),
    'min_samples_split': stats.randint(low=2, high=5),
    'learning_rate': [0.9, 0.5, 0.25, 0.1, 0.05, 0.01]
}


model_cv = RandomizedSearchCV(model_g, param_distributions=param_distributions_g, n_iter=25, verbose=0)
model_cv.fit(xData, yData)
print("predicting 1....")
finalPredictions1 = model_cv.predict(xTestData)

data2 = data2.apply(pd.to_numeric, errors='coerce')
dataM2 = pd.DataFrame(data2).to_numpy()
sh2 = dataM2.shape
print(sh2)
nrows2 = sh2[0]
ncols2 = sh2[1]
xData2 = dataM2[:,4:27]
print(xData2.shape)
yData2 = dataM2[:,2]
yData2 = yData2.flatten()

model_cv2 = RandomizedSearchCV(model_g, param_distributions=param_distributions_g, n_iter=25, verbose=0)
model_cv2.fit(xData2, yData2)
print("predicting 2....")
finalPredictions2 = model_cv2.predict(xTestData)

finalPredictionsAll = np.transpose([finalPredictions1, finalPredictions2])
DF = pd.DataFrame(finalPredictionsAll)
DF.to_csv(savefile)