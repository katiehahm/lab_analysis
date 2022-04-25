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

# read csv file
data = pd.read_csv('C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/Praneeth experiment 3_11_22/ProcessedData/both_towards1_localization_person1.csv',header=None)
# data
# data_labels = data.iloc[0]
# take out the first row
# data = data.drop(data.index[0])
data = data.apply(pd.to_numeric, errors='coerce')
# convert from pandas dataframe to numpy array
dataM = pd.DataFrame(data).to_numpy()
sh = dataM.shape
nrows = sh[0]
ncols = sh[1]

xData = dataM[:,2:26]
yData = dataM[:,1]
# xData = xData[0,:,:]
yData = yData.flatten()

# normalization:
# min_max_scaler = MinMaxScaler().fit(xData)
# xData = min_max_scaler.transform(xData)

# standardization:
scaler = StandardScaler().fit(xData)
xData = scaler.transform(xData)

# GBR, 20 iter
model_g = GradientBoostingRegressor(n_estimators=50)
param_distributions_g = {
    'n_estimators': stats.randint(low=10,high=1000),
    'max_depth': stats.randint(low=2, high=6),
    'min_samples_split': stats.randint(low=2, high=5),
    'learning_rate': [0.9, 0.5, 0.25, 0.1, 0.05, 0.01]
}

finalPredictions = np.zeros((nrows,2)) # for every impact, (real coord, predict coord)

n_fold = 5
# k_fold = StratifiedKFold(n_splits=n_fold,shuffle=True,random_state=13)

k_fold = KFold(n_splits=n_fold, shuffle=True,random_state=13)
for k, (train, test) in enumerate(k_fold.split(xData, yData)):
    print("Running starts kfold")
    xTrain = xData[train]
    yTrain = yData[train]
    xTest = xData[test]
    yTest = yData[test]

    model_cv = RandomizedSearchCV(model_g, param_distributions=param_distributions_g, n_iter=25, verbose=0)
    model_cv.fit(xTrain, yTrain)
    
    yPredict = model_cv.predict(xTest)
  
    # saving values to finalPredictions
    test_idx = np.array(test)
    finalPredictions[test_idx,:] = np.transpose([yTest, yPredict])

yTestAll = finalPredictions[:,0]
yPredictAll = finalPredictions[:,1]
mse_all = mean_squared_error(yTestAll, yPredictAll)
print("The RMSE on all tests: {:.4f}".format(np.sqrt(mse_all)))

DF = pd.DataFrame(finalPredictions)
DF.to_csv("both_towards1_person1_results_norecursion.csv") # CHANGE ***********************************

fig = plt.figure(figsize=(6,6))
plt.plot(yPredictAll,yTestAll,'o')
plt.xlabel('Predicted Values [m]')
plt.ylabel('Target Values [m]')
plt.title('X coordinate Performance')
plt.xlim([-4,4])
plt.ylim([-4,4])
plt.show()
