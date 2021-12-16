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
from tryClassifierRandomizedCV import tryClassifierCV
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


# read csv file
data = pd.read_csv('12_7_21_grf_dataset_predictedloc.csv',header=None)
data = data.apply(pd.to_numeric, errors='coerce')
# convert from pandas dataframe to numpy array
dataM = pd.DataFrame(data).to_numpy()
sh = dataM.shape
nrows = sh[0]
ncols = sh[1]

# dataM[:,2:6] = np.divide(dataM[:,2:6],1000)

X = dataM[:,2:]
y = dataM[:,1]

# inputsTrain = X[0:172,:]
# outputTrain = y[0:172]
# inputsTest = X[172:,:]
# outputTest = y[172:]

# n_fold = 10

# kfold = StratifiedKFold(n_splits=n_fold,shuffle=True,random_state=13)
# k_fold = KFold(n_splits=n_fold, shuffle=True,random_state=13)

# for k, (train, test) in enumerate(k_fold.split(inp, out)):
#     print("Running kfold")
#     xTrain = X[train]
#     yTrain = y[train]
#     xTest = X[test]
#     yTest = y[test]

#     poly_reg = PolynomialFeatures(degree=2)
# 	X_poly = poly_reg.fit_transform(xTrain)
# 	lin_reg2 = LinearRegression()
# 	lin_reg2.fit(X_poly,yTrain)
# # z = np.polyfit(inputs,outputs,2)

poly_reg = PolynomialFeatures(degree=2)
X_poly = poly_reg.fit_transform(X)
lin_reg2 = LinearRegression()
lin_reg2.fit(X_poly,y)

# insoleWeight1_X = X[623:753,:]
# insoleWeight1_y = y[623:753]
# insoleWeight2_X = X[754:861,:]
# insoleWeight2_y = y[754:861]

# regular1_X = X[0:103,:]
# regular1_y = y[0:103]
# regular2_X = X[104:228,:]
# regular2_y = y[104:228]

# insoleWeight1_X = dataM[814:840,2:]
# insoleWeight1_y = dataM[814:840,1]
# insoleWeight2_X = dataM[840:862,2:]
# insoleWeight2_y = dataM[840:862,1]

# regular1_X = dataM[689:710,2:]
# regular1_y = dataM[689:710,1]
# regular2_X = dataM[710:735,2:]
# regular2_y = dataM[710:735,1]


predictions = lin_reg2.predict(poly_reg.fit_transform(X))
mse = mean_squared_error(predictions, y)
print("The root mean squared error (RMSE): {:.4f}".format(np.sqrt(mse)))
plt.plot(predictions, y,'o')
plt.xlim([1,7])
plt.ylim([1,7])
plt.xlabel('Predicted tibial acceleration')
plt.ylabel('Measured tibial acceleration')
plt.show()

# plt.subplot(2,2,1)
# predictions1 = lin_reg2.predict(poly_reg.fit_transform(insoleWeight1_X))
# mse = mean_squared_error(predictions1, insoleWeight1_y)
# print("The root mean squared error (RMSE) 1: {:.4f}".format(np.sqrt(mse)))
# plt.plot(predictions1, insoleWeight1_y, 'o')
# plt.title('InsoleWeight 1')
# plt.xlabel('Predicted Values')
# plt.ylabel('Target Values')
# plt.xlim([1,7])
# plt.ylim([1,7])

# plt.subplot(2,2,2)
# predictions2 = lin_reg2.predict(poly_reg.fit_transform(insoleWeight2_X))
# mse = mean_squared_error(predictions2, insoleWeight2_y)
# print("The root mean squared error (RMSE) 2: {:.4f}".format(np.sqrt(mse)))
# plt.plot(predictions2, insoleWeight2_y, 'o')
# plt.title('InsoleWeight 2')
# plt.xlabel('Predicted Values')
# plt.ylabel('Target Values')
# plt.xlim([1,7])
# plt.ylim([1,7])

# plt.subplot(2,2,3)
# predictions = lin_reg2.predict(poly_reg.fit_transform(regular1_X))
# mse = mean_squared_error(predictions, regular1_y)
# print("The root mean squared error (RMSE) 2: {:.4f}".format(np.sqrt(mse)))
# plt.plot(predictions, regular1_y, 'o')
# plt.title('Regular 1')
# plt.xlabel('Predicted Values')
# plt.ylabel('Target Values')
# plt.xlim([1,7])
# plt.ylim([1,7])

# plt.subplot(2,2,4)
# predictions = lin_reg2.predict(poly_reg.fit_transform(regular2_X))
# mse = mean_squared_error(predictions, regular2_y)
# print("The root mean squared error (RMSE) 2: {:.4f}".format(np.sqrt(mse)))
# plt.plot(predictions, regular2_y, 'o')
# plt.title('Regular 2')
# plt.xlabel('Predicted Values')
# plt.ylabel('Target Values')
# plt.xlim([1,7])
# plt.ylim([1,7])

# from sklearn.cluster import KMeans
# kmeans = KMeans(n_clusters=2)
# kmeans_limp = kmeans.fit_predict(predictions1.reshape(-1,1))

# kmeans2 = KMeans(n_clusters=2)
# kmeans_regular = kmeans2.fit_predict(predictions.reshape(-1,1))


# print(kmeans.cluster_centers_)
# print(kmeans2.cluster_centers_)

# plt.show()

# print(np.shape(outputs))
# print(np.shape(inputs))
# print(np.mean(outputs))
# inputs = preprocessing.normalize(inputs, norm='l1')

# # cross validation to increase amount of data
# n_fold = 10

# kfold = StratifiedKFold(n_splits=n_fold,shuffle=True,random_state=13)
# k_fold = KFold(n_splits=n_fold, shuffle=True,random_state=13)

# # GBR, 20 iter
# model_g = GradientBoostingRegressor(n_estimators=50)

# param_distributions_g = {
#     'n_estimators': stats.randint(low=100,high=800),
#     'max_depth': [2],
#     'min_samples_split': stats.randint(low=2, high=5),
#     'learning_rate': [0.01]
# }

# # RF, 50 iter
# model_r = RandomForestRegressor(n_estimators=200)
# param_distributions_r = {
#     'n_estimators': stats.randint(low=50, high=300),
#     'min_samples_split': stats.randint(low=2, high=6),
#     'min_samples_leaf': stats.randint(low=1, high=5),
#     'min_weight_fraction_leaf': [0.5, 0.25, 0.1, 0.05, 0.01, 0.0]
# }

# model_l = LinearRegression()
# param_distributions_l = {}
# graphit=True
# i = 1
# tryClassifierCV(inputs, outputs, model_r, param_distributions_r, 15, 'RF', k_fold,i,graphit)


