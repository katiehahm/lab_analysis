from sklearn.multioutput import MultiOutputRegressor
import numpy as np
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt
from sklearn.model_selection import RandomizedSearchCV

def tryClassifierSingleCV(inp, out, model, param_distributions, n_it, name_str, k_fold):
    yTestAll = [[1,1]]
    yPredictAll = [[1,1]]
    
    for k, (train, test) in enumerate(k_fold.split(inp, out)):
        print("Running kfold")
        xTrain = inp[train]
        yTrain = out[train]
        xTest = inp[test]
        yTest = out[test]

        
        model_cv = RandomizedSearchCV(model, param_distributions=param_distributions, n_iter=n_it, verbose=0)
        model_cv.fit(xTrain, yTrain)
        
        yPredict = model_cv.predict(xTest)
      
        # saving values to plot
        yTestAll = np.append(yTestAll,yTest,axis=0)
        yPredictAll = np.append(yPredictAll,yPredict,axis=0)

        # finding accuracy for each location
        mse = mean_squared_error(yTest, yPredict)
        print("The mean squared error (MSE) on test set: {:.4f}".format(np.sqrt(mse)))

        # print(regr_multi.get_params())
        # print("The best estimator across ALL searched params:", regr_multi.best_estimator_)
        # print("The best score across ALL searched params:", regr_multi.best_score_)
        # print("The best parameters across ALL searched params:", regr_multi.best_params_)

    yTestAll = np.delete(yTestAll,0,0)
    yPredictAll = np.delete(yPredictAll,0,0)

    print('#####')
    print(name_str + ' Cross-validation results:')
    mse_all = mean_squared_error(yTestAll, yPredictAll)
    print("The RMSE on all tests: {:.4f}".format(np.sqrt(mse_all)))

    fig = plt.figure(figsize=(6,6))
    plt.plot(yPredictAll,yTestAll,'o')
    plt.xlabel('Predicted Values [m]')
    plt.ylabel('Target Values [m]')
    plt.title('X coordinate Performance')
    plt.xlim([-3,4])
    plt.ylim([-3,4])
    plt.show()

if __name__ == '__main__':
    tryClassifier()
