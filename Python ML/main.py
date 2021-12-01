import numpy as np  

import pandas as pd
# from IPython.display import Image
# import seaborn as sns
# from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import KFold
from sklearn import datasets, svm, metrics, preprocessing
# from sklearn.linear_model import LinearRegression
from scipy import stats 
# from sklearn.model_selection import RandomizedSearchCV
# from sklearn.neighbors import KNeighborsClassifier
# from sklearn.ensemble import RandomForestClassifier
# from sklearn.model_selection import StratifiedKFold
# from sklearn.metrics import confusion_matrix
# from xgboost import XGBClassifier
# from xgboost import plot_importance
# from matplotlib import pyplot
# from sklearn.model_selection import train_test_split
# from numpy import loadtxt
# from numpy import sort
# from sklearn.metrics import accuracy_score
# from sklearn.feature_selection import SelectFromModel
# from sklearn import datasets, ensemble
# from sklearn.inspection import permutation_importance
# from sklearn.metrics import mean_squared_error
# from sklearn.model_selection import train_test_split
# from sklearn.datasets import load_linnerud

from sklearn.ensemble import GradientBoostingRegressor
from sklearn.multioutput import MultiOutputRegressor
from sklearn.ensemble import RandomForestRegressor


from tryGBR import GBR
from tryClassifierCV import tryClassifierCV
from tryClassifierSingleCV import tryClassifierSingleCV2
from GBRtest import GBRtest
from RFtest import RFtest
from LinearRegtest import LinearRegtest
# from sklearn import linear_model

# read csv file
# data = pd.read_csv('06-04_06-10_finaldata.csv',header=None)
data = pd.read_csv('12-1-21_localization_dataset.csv',header=None)
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

# define inputs/outputs
outputs = dataM[:,0:1]
inputs = dataM[:,1:]
inputsNorm = preprocessing.normalize(inputs, norm='l1')


# cross validation to increase amount of data
n_fold = 5
k_fold = KFold(n_splits=n_fold, shuffle=True,random_state=30)

print(np.shape(outputs))
print(np.shape(inputsNorm))

param_dist = {
	'n_estimators' : 50,
	'max_depth' : 3,
	'min_samples_split' : 2,
	'learning_rate' : 0.1
}

# GBR(inputs, outputsM, param_dist, k_fold)

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

# print(model_g.get_params().keys())

tryClassifierSingleCV(inputs, outputs, model_g, param_distributions_g, 5, 'GBR', k_fold)
# GBRtest(inputs, outputsM, k_fold)
# RFtest(inputs, outputsM, k_fold)
# LinearRegtest(inputs, outputsM, k_fold)

# using previous location to estimate next one 8/23/21
# print(np.shape(outputs[1:,:]))
# print(np.shape(inputs[:-1,:]))
# inputs_wprev = np.concatenate((inputs[1:,:], outputsM[:-1,:]), axis=1)
# print(np.shape(outputsM[:-1,:]))
# print(np.shape(inputs_wprev))
# RFtest(inputs_wprev, outputsM[:-1,:], k_fold)