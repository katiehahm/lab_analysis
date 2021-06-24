import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
# from sklearn.datasets import make_blobs
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
import math
from sklearn.metrics import confusion_matrix
import seaborn as sns

data = pd.read_csv('06-04_06-10_edited_finaldata_nonan_gridlabel.csv',header=None)
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
labels = dataM[:,ncols-1:]
features = dataM[:,0:ncols-1]

scaler = StandardScaler()
scaled_feat = scaler.fit_transform(features)

kmeans = KMeans(init="random", n_clusters = 18, n_init=1000, max_iter = 3000)
kmeans.fit(scaled_feat)

predict_labels = kmeans.labels_
predict_labels += 1

# convert back to coordinate for RMSE
labels_coord = [[1,1]]
predict_coord = [[1,1]]
for i in range(len(labels)):
	lbl = labels[i]
	x_lbl = 1 + math.floor( (lbl-1)/3 )
	y_lbl = 1 + (lbl-1)%3
	new_lbl = np.array([[x_lbl, y_lbl]], dtype="object")
	labels_coord = np.append(labels_coord,new_lbl,axis=0)
	p_lbl = predict_labels[i]
	px_lbl = 1 + math.floor( (p_lbl-1)/3 )
	py_lbl = 1 + (p_lbl-1)%3
	p_new_lbl = np.array([[px_lbl,py_lbl]], dtype="object")
	predict_coord = np.append(predict_coord,p_new_lbl, axis=0)
labels_coord = np.delete(labels_coord,0,0)
predict_coord = np.delete(predict_coord,0,0)

rmse = mean_squared_error(labels_coord, predict_coord, squared=False)
print("RMSE: {:.4f}".format(rmse*1.21)) # converting back to m from coordinates

# fig = plt.figure(figsize=(6,6))
# plt.plot(predict_labels,labels,'o')
# plt.xlabel('Predicted Values [m]')
# plt.ylabel('Target Values [m]')
# plt.title('X coordinate Performance')


confusion= confusion_matrix(labels,predict_labels)
    
lbls = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18']

df_cm = pd.DataFrame(confusion, index=lbls,
              columns=lbls)
plt.figure(figsize = (10,7))
ax1 = plt.axes()
ax1.set_title("KNN Confusion Matrix", fontsize = 18)
sns.heatmap(df_cm, annot=True)

fig2 = plt.figure(figsize=(6,6))
plt.hist(labels)

fig2 = plt.figure(figsize=(6,6))
plt.hist(predict_labels)

plt.show()