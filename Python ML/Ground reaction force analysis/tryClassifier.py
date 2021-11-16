# THIS IS THE 1D VERSION
from sklearn.multioutput import MultiOutputRegressor
import numpy as np
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
import math

def tryClassifier(inp, out, model, param_distributions, name_str, k_fold):
    yTestAll = [1]
    yPredictAll = [1]
    
    for k, (train, test) in enumerate(k_fold.split(inp, out)):
        print("Running kfold")
        xTrain = inp[train]
        yTrain = out[train]
        xTest = inp[test]
        yTest = out[test]

        model_cv.set_params(param_distributions=param_distributions)
        model_cv.fit(xTrain, yTrain)
        
        yPredict = model_cv.predict(xTest)
      
        # saving values to plot
        yTestAll = np.append(yTestAll,yTest)
        yPredictAll = np.append(yPredictAll,yPredict)

        # finding accuracy for each location
        mse = mean_squared_error(yTest, yPredict)
        print("The root mean squared error (RMSE) on test set: {:.4f}".format(np.sqrt(mse)))

        # print(model_cv.get_params())

        # print("The best estimator across ALL searched params:", regr_multi.best_estimator_)
        # print("The best score across ALL searched params:", regr_multi.best_score_)
        # print("The best parameters across ALL searched params:", regr_multi.best_params_)

    yTestAll = np.delete(yTestAll,0)
    yPredictAll = np.delete(yPredictAll,0)

    print('#####')
    print(name_str + ' Cross-validation results:')
    mse_all = mean_squared_error(yTestAll, yPredictAll)
    print("The RMSE on all tests: {:.4f}".format(np.sqrt(mse_all)))

    fig = plt.figure(figsize=(6,6))
    plt.plot(yPredictAll,yTestAll,'o')
    plt.xlabel('Predicted Values [V]')
    plt.ylabel('Target Values [V]')
    plt.title('Ground reaction force estimator performance')
    axes = plt.gca()
    y_min,y_max = axes.get_ylim()
    plt.xlim([y_min,y_max])
    # plt.xlim([-3,4])
    # plt.ylim([-3,4])

    plt.show()

if __name__ == '__main__':
    tryClassifier()
