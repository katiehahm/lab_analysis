import numpy as np
import pandas as pd
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

# read csv file
data = pd.read_csv('distance_to_sensor_grf_features.csv',header=None)
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
# 0th column is the subject number
subject_split = np.array([0,826,1546,2501,3491,4268,5287,6276,7482,8667])

for i in range(3):
    i += 7
    subj_start = subject_split[i]
    if i == 9:
        subj_end = nrows
    else:
        subj_end = subject_split[i+1]
    outputs = dataM[subj_start:subj_end,1] # 826 starts subj 2
    inputs = dataM[subj_start:subj_end,2:]
    # inputsNorm = preprocessing.normalize(inputs, norm='l1')
    # outputsCM = np.divide(outputs,10)
    # outputsM = np.divide(outputsCM,100)

    # cross validation to increase amount of data
    n_fold = 5

    kfold = StratifiedKFold(n_splits=n_fold,shuffle=True,random_state=13)

    k_fold = KFold(n_splits=n_fold, shuffle=True,random_state=13)

    print(np.shape(outputs))
    print(np.shape(inputs))

    param_dist = {
    	'n_estimators' : 100,
    	'max_depth' : 3,
    	'min_samples_split' : 3,
    	'learning_rate' : 0.1
    }

    # GBR(inputs, outputsM, param_dist, k_fold)

    # GBR, 20 iter
    model_g = GradientBoostingRegressor(n_estimators=50)
    # all distributions to run randomcv
    # param_distributions_g = {
    #     'n_estimators': stats.randint(low=10,high=1000),
    #     'max_depth': stats.randint(low=2, high=6),
    #     'min_samples_split': stats.randint(low=2, high=5),
    #     'learning_rate': [0.9, 0.5, 0.25, 0.1, 0.05, 0.01]
    # }
    # modified distributions due to prior run
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

    # model_k = KNeighborsClassifier(n_neighbors=3)
    # param_distributions_k = {
    #     # 'sample_weights': stats.randint(Low=)
    # }

    print('################')
    print('Subject number: ')
    print(i+1)
    tryClassifierCV(inputs, outputs, model_g, param_distributions_g, 15, 'GBR', k_fold)

    # X_train, X_test, y_train, y_test = train_test_split(inputs, outputs, test_size=0.1, random_state=13)
    # params = {
    #     "n_estimators": 500,
    #     "max_depth": 4,
    #     "min_samples_split": 5,
    #     "learning_rate": 0.01,
    #     # "loss": "deviance",
    # }

# ################### this is just trying GBR quickly:
# gbr_model = GradientBoostingRegressor()
# gbr_model.set_params(**params)
# # reg = ensemble.GradientBoostingRegressor(**params)
# gbr_model.fit(X_train, y_train)

# y_predict = gbr_model.predict(X_test)

# # error_arr = y_test - y_predict
# # mse = np.sum(np.square(error_arr))
# mse = mean_squared_error(y_test, y_predict)
# rmse = math.sqrt(mse)
# print("The root mean squared error (MSE) on test set: {:.4f}".format(rmse))
# avg_fsr = statistics.mean(output)
# print("The average value of the fsr inputs: {:.4f}".format(avg_fsr))


# feature_importance = gbr_model.feature_importances_
# feature_names = np.array(['x coord','y coord','pkM1','pkM2','pkM3','pkM4','energy1','energy2','energy3','energy4'])
# sorted_idx = np.argsort(feature_importance)
# pos = np.arange(sorted_idx.shape[0]) + 0.5
# fig = plt.figure(figsize=(12, 6))
# plt.subplot(1, 2, 1)
# plt.barh(pos, feature_importance[sorted_idx], align="center")
# plt.yticks(pos, feature_names[sorted_idx])
# plt.title("Feature Importance (MDI)")

# result = permutation_importance(
#     gbr_model, X_test, y_test, n_repeats=10, random_state=42, n_jobs=2
# )
# sorted_idx = result.importances_mean.argsort()
# plt.subplot(1, 2, 2)
# plt.boxplot(
#     result.importances[sorted_idx].T,
#     vert=False,
#     labels=feature_names[sorted_idx],
# )
# plt.title("Permutation Importance (test set)")
# fig.tight_layout()
# plt.show()




# print(model_g.get_params().keys())

# tryClassifierCV(inputs, outputs, model_r, param_distributions_r, 20, 'RF', k_fold)
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

