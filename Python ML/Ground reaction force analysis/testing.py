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

# read csv file
data = pd.read_csv('grf_features_yAccel_energyenvelope.csv',header=None)
data = data.apply(pd.to_numeric, errors='coerce')
# convert from pandas dataframe to numpy array
dataM = pd.DataFrame(data).to_numpy()
sh = dataM.shape
nrows = sh[0]
ncols = sh[1]

outputs = dataM[230:,1] 
inputs = dataM[230:,2:]

print(np.shape(outputs))
print(np.shape(inputs))
print(np.mean(outputs))
inputs = preprocessing.normalize(inputs, norm='l1')

# cross validation to increase amount of data
n_fold = 10

kfold = StratifiedKFold(n_splits=n_fold,shuffle=True,random_state=13)
k_fold = KFold(n_splits=n_fold, shuffle=True,random_state=13)

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

model_l = LinearRegression()
param_distributions_l = {}
graphit=True
i = 1
tryClassifierCV(inputs, outputs, model_r, param_distributions_r, 15, 'RF', k_fold,i,graphit)


