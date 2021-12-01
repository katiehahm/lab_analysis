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

import matplotlib.pyplot as plt
from sklearn.model_selection import RandomizedSearchCV
from sklearn.linear_model import LinearRegression

# read csv file
data = pd.read_csv('grf_features_yAccel_trainingset.csv',header=None)
data = data.apply(pd.to_numeric, errors='coerce')
# convert from pandas dataframe to numpy array
dataM = pd.DataFrame(data).to_numpy()
sh = dataM.shape
nrows = sh[0]
ncols = sh[1]

inputsall = preprocessing.normalize(dataM[:,2:],norm='l1')
inputTrain = inputsall[0:172,:]

inputTest = inputsall[173:,:]

# inputTrain = dataM[0:172,2:]
outputTrain = dataM[0:172,1]

print(np.shape(inputTrain))
print(np.shape(outputTrain))

# inputTest = dataM[173:,2:]
outputTest = dataM[173:,1]

# print(np.shape(outputs))
# print(np.shape(inputs))
# inputsNorm = preprocessing.normalize(inputs, norm='l1')

# GBR, 20 iter
model_g = GradientBoostingRegressor(n_estimators=50)

param_distributions_g = {
    'n_estimators': stats.randint(low=100,high=800),
    'max_depth': [2],
    'min_samples_split': stats.randint(low=2, high=5),
    'learning_rate': [0.01]
}

# RF, 50 iter
model_r = RandomForestRegressor(n_estimators=200)
param_distributions_r = {
    'n_estimators': stats.randint(low=50, high=300),
    'min_samples_split': stats.randint(low=2, high=6),
    'min_samples_leaf': stats.randint(low=1, high=5),
    'min_weight_fraction_leaf': [0.5, 0.25, 0.1, 0.05, 0.01, 0.0]
}

n_it = 20
model_cv = RandomizedSearchCV(model_r, param_distributions=param_distributions_r, n_iter=n_it, verbose=0)
model_cv.fit(inputTrain, outputTrain)

yPredict = model_cv.predict(inputTest)
mse = mean_squared_error(outputTest, yPredict)
print("The root mean squared error (RMSE): {:.4f}".format(np.sqrt(mse)))

# print(model_cv.get_params())
from pprint import pprint
# pprint(model_cv.get_params())
# pprint(model_cv.best_estimator_.get_params())

# fig = plt.figure(figsize=(6,6))
# plt.plot(yPredict,yTest,'o')
# plt.xlabel('Predicted Values [V]')
# plt.ylabel('Target Values [V]')
# subjnum = i + 1
# plt.title('Estimator performance')
# axes = plt.gca()
# y_min,y_max = axes.get_ylim()
# plt.xlim([y_min,y_max])
# # plt.xlim([-3,4])
# # plt.ylim([-3,4])

# plt.show()

separation_arr = np.array([173,256,356,459,566,671,776,861])
separation_arr = np.subtract(separation_arr, separation_arr[0])

fig = plt.figure()

y_min = 0
y_max = 7

plt.subplot(2,4,1)
x = yPredict[separation_arr[0]:separation_arr[1]]
y = outputTest[separation_arr[0]:separation_arr[1]]
plt.plot(x, y, 'o')
plt.title('Regular 1')
plt.xlabel('Predicted Values')
plt.ylabel('Target Values')
plt.xlim([1,4])
plt.ylim([1,4])

plt.subplot(2,4,2)
x = yPredict[separation_arr[1]:separation_arr[2]]
y = outputTest[separation_arr[1]:separation_arr[2]]
plt.plot(x, y, 'o')
plt.title('Regular 2')
plt.xlabel('Predicted Values')
plt.ylabel('Target Values')
plt.xlim([1,4])
plt.ylim([1,4])

plt.subplot(2,4,3)
x = yPredict[separation_arr[2]:separation_arr[3]]
y = outputTest[separation_arr[2]:separation_arr[3]]
plt.plot(x, y, 'o')
plt.title('Slow')
plt.xlabel('Predicted Values')
plt.ylabel('Target Values')
plt.xlim([1,3])
plt.ylim([1,3])

plt.subplot(2,4,4)
x = yPredict[separation_arr[3]:separation_arr[4]]
y = outputTest[separation_arr[3]:separation_arr[4]]
plt.plot(x, y, 'o')
plt.title('Stiff R knee')
plt.xlabel('Predicted Values')
plt.ylabel('Target Values')
plt.xlim([1,5])
plt.ylim([1,5])

plt.subplot(2,4,5)
x = yPredict[separation_arr[4]:separation_arr[5]]
y = outputTest[separation_arr[4]:separation_arr[5]]
plt.plot(x, y, 'o')
plt.title('Stiff L knee')
plt.xlabel('Predicted Values')
plt.ylabel('Target Values')
plt.xlim([1,4])
plt.ylim([1,4])

plt.subplot(2,4,6)
x = yPredict[separation_arr[5]:separation_arr[6]]
y = outputTest[separation_arr[5]:separation_arr[6]]
plt.plot(x, y, 'o')
plt.title('Insole R weight L 1')
plt.xlabel('Predicted Values')
plt.ylabel('Target Values')
plt.xlim([1,7])
plt.ylim([1,7])

plt.subplot(2,4,7)
x = yPredict[separation_arr[6]:separation_arr[7]]
y = outputTest[separation_arr[6]:separation_arr[7]]
plt.plot(x, y, 'o')
plt.title('Insole R weight L 2')
plt.xlabel('Predicted Values')
plt.ylabel('Target Values')
plt.xlim([1,7])
plt.ylim([1,7])

plt.show()