# to run GBR with specific parameters (instead of using tryClassifierCV)
from sklearn.multioutput import MultiOutputRegressor
import numpy as np
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
from sklearn.model_selection import RandomizedSearchCV
from sklearn.ensemble import GradientBoostingRegressor

def GBRtest(inp, out, k_fold):
    yTestAll = [[1,1]]
    yPredictAll = [[1,1]]

    plt.rcParams.update({'font.size': 18})

    for k, (train, test) in enumerate(k_fold.split(inp, out)):
        print("Running kfold")
        xTrain = inp[train]
        yTrain = out[train]
        xTest = inp[test]
        yTest = out[test]

        param_dist = {'n_estimators' : 50, 'max_depth' : 3, 'min_samples_split' : 2, 'learning_rate' : 0.1}

        model_gbr = GradientBoostingRegressor()
        model_gbr.set_params(**param_dist)
        regr_multi = MultiOutputRegressor(model_gbr)
        regr_multi.fit(xTrain, yTrain)

        yPredict = regr_multi.predict(xTest)
      
        # saving values to plot
        yTestAll = np.append(yTestAll,yTest,axis=0)
        yPredictAll = np.append(yPredictAll,yPredict,axis=0)

        # print(regr_multi.get_params())
        # print("The best estimator across ALL searched params:", regr_multi.best_estimator_)
        # print("The best score across ALL searched params:", regr_multi.best_score_)
        # print("The best parameters across ALL searched params:", regr_multi.best_params_)

    yTestAll = np.delete(yTestAll,0,0)
    yPredictAll = np.delete(yPredictAll,0,0)

    print('#####')
    print('GBR Cross-validation results:')
    mse_all = mean_squared_error(yTestAll, yPredictAll)
    print("The RMSE on all tests: {:.4f}".format(np.sqrt(mse_all)))

    mse_x = mean_squared_error(yTestAll[:,0:1], yPredictAll[:,0:1])
    mse_z = mean_squared_error(yTestAll[:,1:], yPredictAll[:,1:])
    print("The RMSE on all tests X-dir: {:.4f}".format(np.sqrt(mse_x)))
    print("The RMSE on all tests Z-dir: {:.4f}".format(np.sqrt(mse_z)))

    fig = plt.figure(figsize=(6,6))
    plt.plot(yPredictAll[:,0:1],yTestAll[:,0:1],'o')
    plt.xlabel('Predicted Values [m]')
    plt.ylabel('Target Values [m]')
    plt.title('X coordinate Performance')
    plt.xlim([-4,4])
    plt.ylim([-4,4])
    plt.savefig('gbr1.png')

    fig2 = plt.figure(figsize=(6,6))
    plt.plot(yPredictAll[:,1:],yTestAll[:,1:],'o')
    plt.xlabel('Predicted Values [m]')
    plt.ylabel('Target Values [m]')
    plt.title('Z coordinate Performance')
    plt.xlim([-4,4])
    plt.ylim([-4,4])

    plt.savefig('gbr2.png')
    plt.show()
