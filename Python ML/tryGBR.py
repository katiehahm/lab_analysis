from sklearn.ensemble import GradientBoostingRegressor
from sklearn.multioutput import MultiOutputRegressor
import numpy as np
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt

def GBR(inp, out, param_dist, k_fold):
    yTestAll = [[1,1]]
    yPredictAll = [[1,1]]
    
    for k, (train, test) in enumerate(k_fold.split(inp, out)):
        xTrain = inp[train]
        yTrain = out[train]
        xTest = inp[test]
        yTest = out[test]

        model_gbr = GradientBoostingRegressor()
        model_gbr.set_params(**param_dist)
        regr_multi = MultiOutputRegressor(model_gbr)
        regr_multi.fit(xTrain, yTrain)

        print(regr_multi.get_params())
        
        yPredict = regr_multi.predict(xTest)
        
        # saving values to plot
        yTestAll = np.append(yTestAll,yTest,axis=0)
        yPredictAll = np.append(yPredictAll,yPredict,axis=0)

        # finding accuracy for each location
        mse = mean_squared_error(yTest, yPredict)
        print("The mean squared error (MSE) on test set: {:.4f}".format(mse))

    yTestAll = np.delete(yTestAll,0,0)
    yPredictAll = np.delete(yPredictAll,0,0)

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

    fig2 = plt.figure(figsize=(6,6))
    plt.plot(yPredictAll[:,1:],yTestAll[:,1:],'o')
    plt.xlabel('Predicted Values [m]')
    plt.ylabel('Target Values [m]')
    plt.title('Z coordinate Performance')

    plt.show()

if __name__ == '__main__':
    GBR()
