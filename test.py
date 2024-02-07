import LinReg
import numpy as np
import pandas as pd

data = pd.read_csv('dataset.txt', header=None)

regressor = LinReg.LinReg()

X = regressor.get_columns(data.values, np.ones(data.shape[1]))
RMSE_list = []
for _ in range(10):
    RMSE_list.append(regressor.get_fitness(X[:,:-1], X[:,-1]))
print("Max RMSE: ", max(RMSE_list))
print("Min RMSE: ", min(RMSE_list))
print("Mean RMSE: ", np.mean(RMSE_list))
print("Median RMSE: ", np.median(RMSE_list))

