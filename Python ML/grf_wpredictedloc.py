import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn import datasets, svm, metrics, preprocessing
# from sklearn.linear_model import LinearRegression
from scipy import stats 
from sklearn import datasets, ensemble
# from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_squared_error
# from sklearn.model_selection import train_test_split
# from sklearn.datasets import load_linnerud
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.ensemble import RandomForestRegressor
from tryGBR import GBR
from GBRtest import GBRtest
from RFtest import RFtest
from LinearRegtest import LinearRegtest
# from sklearn import linear_model
from sklearn.model_selection import train_test_split
import math
import statistics
import matplotlib.pyplot as plt
from sklearn.inspection import permutation_importance
from sklearn.model_selection import StratifiedKFold

from sklearn.linear_model import LinearRegression
from numpy.polynomial import Chebyshev
from numpy.polynomial.chebyshev import chebfit
from sklearn.preprocessing import PolynomialFeatures
from sklearn.model_selection import RandomizedSearchCV
from sklearn.model_selection import KFold


# read csv file
data = pd.read_csv('C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment3/ProcessedData/ExcelData/grf_features_subj4.csv',header=None)
data = data.apply(pd.to_numeric, errors='coerce')
# convert from pandas dataframe to numpy array
dataM = pd.DataFrame(data).to_numpy()
sh = dataM.shape
nrows = sh[0]
ncols = sh[1]

allX = dataM[:,2:]
ally = dataM[:,1]
takes = dataM[:,0]

regular_idx = np.where(takes == 1)
# regular_idx = np.append(regular_idx, np.where(takes== 5))

weight_idx = np.where(takes == 4)
# weight_idx = np.append(weight_idx, np.where(takes == 3))

X = dataM[np.append(regular_idx,weight_idx),2:]
y = dataM[np.append(regular_idx,weight_idx),1]

finalPredictions = np.zeros((nrows,3)) # for every impact, (take, real grf, predict grf)

n_fold = 5

# GBR, 20 iter
model_g = GradientBoostingRegressor(n_estimators=50)
param_distributions_g = {
    'n_estimators': stats.randint(low=10,high=1000),
    'max_depth': stats.randint(low=2, high=6),
    'min_samples_split': stats.randint(low=2, high=5),
    'learning_rate': [0.9, 0.5, 0.25, 0.1, 0.05, 0.01]
}

# RF, 50 iter
model_r = RandomForestRegressor(n_estimators=200)
param_distributions_r = {
    'n_estimators': stats.randint(low=50, high=300),
    'min_samples_split': stats.randint(low=2, high=6),
    'min_samples_leaf': stats.randint(low=1, high=5),
    'min_weight_fraction_leaf': [0.5, 0.25, 0.1, 0.05, 0.01, 0.0]
}

k_fold = KFold(n_splits=n_fold, shuffle=True,random_state=13)
for k, (train, test) in enumerate(k_fold.split(X, y)):
	print("Running starts kfold")
	xTrain = X[train]
	yTrain = y[train]
	xTest = X[test]
	yTest = y[test]

	# Polynomial regression
	# poly_reg = PolynomialFeatures(degree=3)
	# X_poly = poly_reg.fit_transform(xTrain)
	# lin_reg2 = LinearRegression()
	# lin_reg2.fit(X_poly,yTrain)

	# yPredict = lin_reg2.predict(poly_reg.fit_transform(xTest))

	model_cv = RandomizedSearchCV(model_g, param_distributions=param_distributions_g, n_iter=25, verbose=0)
	model_cv.fit(xTrain, yTrain)

	yPredict = model_cv.predict(xTest)

	# saving values to finalPredictions
	currTest = takes[test]
	finalPredictions[test,:] = np.transpose([currTest, yTest, yPredict])

DF = pd.DataFrame(finalPredictions)
DF.to_csv("grf_results_GBR_regularweightonly_subj4.csv") # CHANGE ***********************************

yPredictAll = finalPredictions[:,2]
yAll = finalPredictions[:,1]
mse = mean_squared_error(yPredictAll, yAll)
print("The root mean squared error (RMSE): {:.4f}".format(np.sqrt(mse)))
plt.plot(yPredictAll, yAll,'o')
plt.xlim([1,7])
plt.ylim([1,7])
plt.xlabel('Predicted tibial acceleration')
plt.ylabel('Measured tibial acceleration')
plt.show()
