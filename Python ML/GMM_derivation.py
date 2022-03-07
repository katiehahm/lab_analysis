import imageio
import matplotlib.animation as ani
import matplotlib.cm as cmx
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from matplotlib.patches import Ellipse
from PIL import Image
from sklearn import datasets
from sklearn.cluster import KMeans

from sklearn.mixture import GaussianMixture

# 2/28/22 after 2nd com mtg, try derivation of GMM 

data = pd.read_csv('C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment3/ProcessedData/ExcelData/all_step_times1.csv',header=None)
data = data.apply(pd.to_numeric, errors='coerce')
# convert from pandas dataframe to numpy array
dataM = pd.DataFrame(data).to_numpy()
sh = dataM.shape
nrows = sh[0]
ncols = sh[1]

steps = dataM[:,1]
footlabel = dataM[:,2]
takes = dataM[:,0]

def gaussian(X, mu, cov):
    n = X.shape[1]
    diff = (X - mu).T
    return np.diagonal(1 / ((2 * np.pi) ** (n / 2) * np.linalg.det(cov) ** 0.5) * np.exp(-0.5 * np.dot(np.dot(diff.T, np.linalg.inv(cov)), diff))).reshape(-1, 1)

x0 = np.array([[0.05, 1.413, 0.212], [0.85, -0.3, 1.11], [11.1, 0.4, 1.5], [0.27, 0.12, 1.44], [88, 12.33, 1.44]])
mu = np.mean(x0, axis=0)
cov = np.dot((x0 - mu).T, x0 - mu) / (x0.shape[0] - 1)

y = gaussian(x0, mu=mu, cov=cov)
y

# step 1
def initialize_clusters(X, n_clusters):
    clusters = []
    idx = np.arange(X.shape[0])
    
    # We use the KMeans centroids to initialise the GMM
    # kmeans = KMeans(n_clusters).fit(X)
    # mu_k = kmeans.cluster_centers_

    # try using original GMM for cluster centers
    gmm = GaussianMixture(n_components=n_clusters, max_iter=50).fit(data)
    mu_k = gmm.means_
    mu_k = mu_k.flatten()
    
    for i in range(n_clusters):
        clusters.append({
            'pi_k': 1.0 / n_clusters,
            'mu_k': mu_k[i],
            'cov_k': np.identity(X.shape[1], dtype=np.float64)
        })
        
    return clusters

# step 2
def expectation_step(X, clusters):
    totals = np.zeros((X.shape[0], 1), dtype=np.float64)
    
    for cluster in clusters:
        pi_k = cluster['pi_k']
        mu_k = cluster['mu_k']
        cov_k = cluster['cov_k']
        
        gamma_nk = (pi_k * gaussian(X, mu_k, cov_k)).astype(np.float64)
        
        for i in range(X.shape[0]):
            totals[i] += gamma_nk[i]
        
        cluster['gamma_nk'] = gamma_nk
        cluster['totals'] = totals
        
    
    for cluster in clusters:
        cluster['gamma_nk'] /= cluster['totals']

# step 3
def maximization_step(X, clusters):
    N = float(X.shape[0])
  
    for cluster in clusters:
        gamma_nk = cluster['gamma_nk']
        cov_k = np.zeros((X.shape[1], X.shape[1]))
        
        N_k = np.sum(gamma_nk, axis=0)
        
        pi_k = 0.5 # CHANGED
        # pi_k = N_k / N # original
        mu_k = np.sum(gamma_nk * X, axis=0) / N_k
        
        for j in range(X.shape[0]):
            diff = (X[j] - mu_k).reshape(-1, 1)
            cov_k += gamma_nk[j] * np.dot(diff, diff.T)
            
        cov_k /= N_k
        
        cluster['pi_k'] = pi_k
        cluster['mu_k'] = mu_k
        cluster['cov_k'] = cov_k

def get_likelihood(X, clusters):
    sample_likelihoods = np.log(np.array([cluster['totals'] for cluster in clusters]))
    return np.sum(sample_likelihoods), sample_likelihoods

def train_gmm(X, n_clusters, n_epochs):
    clusters = initialize_clusters(X, n_clusters)
    likelihoods = np.zeros((n_epochs, ))
    scores = np.zeros((X.shape[0], n_clusters))
    history = []

    for i in range(n_epochs):
        clusters_snapshot = []
        
        # This is just for our later use in the graphs
        for cluster in clusters:
            clusters_snapshot.append({
                'mu_k': cluster['mu_k'].copy(),
                'cov_k': cluster['cov_k'].copy()
            })
            
        history.append(clusters_snapshot)
      
        expectation_step(X, clusters)
        maximization_step(X, clusters)

        likelihood, sample_likelihoods = get_likelihood(X, clusters)
        likelihoods[i] = likelihood

        # if np.max(np.array([cluster['pi_k'] for cluster in clusters])) > 0.55:
        # 	break

        # print('Epoch: ', i + 1, 'Likelihood: ', likelihood)
        
    for i, cluster in enumerate(clusters):
        scores[:, i] = np.log(cluster['gamma_nk']).reshape(-1)
        
    return clusters, likelihoods, scores, sample_likelihoods, history

# train model:
n_clusters = 2
n_epochs = 200 # sort of like iterations
for t in range(int(np.amax(takes))):
	data_idx = [i for i, value in enumerate(takes) if value == (t+1)]
	data = steps[data_idx]
	data = np.reshape(data, (-1,1))
	clusters, likelihoods, scores, sample_likelihoods, history = train_gmm(data, n_clusters, n_epochs)

	plt.figure(figsize=(10, 10))
	plt.title('Log-Likelihood')
	plt.plot(np.arange(1, n_epochs + 1), likelihoods)
	plt.show()

	gmm = GaussianMixture(n_components=n_clusters, max_iter=50).fit(data)
	gmm_scores = gmm.score_samples(data)

	data_footlabel = footlabel[data_idx]
	left_idx  = [i for i, value in enumerate(data_footlabel) if value == 0]
	right_idx  = [i for i, value in enumerate(data_footlabel) if value == 1]
	real1 = np.mean(data[left_idx,:])
	real2 = np.mean(data[right_idx,:])

	# print('Means by sklearn:\n', gmm.means_)
	# print('Means by our implementation:\n', np.array([cluster['mu_k'].tolist() for cluster in clusters]))
	# print('Weights by sklearn:\n', gmm.weights_)
	# # print('Weights by our implementation:\n', np.array([cluster['pi_k'].tolist() for cluster in clusters]))
	# print('Real step times of leg 1:\n', real1)
	# print('Real step times of leg 2:\n', real2)
	mean_results = np.array([cluster['mu_k'].tolist() for cluster in clusters])
	pi_results = np.array([cluster['pi_k'] for cluster in clusters])
	# print(pi_results)
	# for original: 
	# output = [gmm.means_[0,0], gmm.means_[1,0], gmm.weights_[0], gmm.weights_[1], mean_results[0,0], mean_results[1,0], pi_results[0,0], pi_results[1,0], real1, real2] # gmm1, gmm2, scal1, scal2, real1, real2
	# for edited version:
	output = [gmm.means_[0,0], gmm.means_[1,0], gmm.weights_[0], gmm.weights_[1], mean_results[0,0], mean_results[1,0], pi_results[0], pi_results[1], real1, real2] # gmm1, gmm2, scal1, scal2, real1, real2
	print(output)




# def create_cluster_animation(X, history, scores):
#     fig, ax = plt.subplots(1, 1, figsize=(10, 10))
#     colorset = ['blue', 'red', 'black']
#     images = []
    
#     for j, clusters in enumerate(history):
      
#         idx = 0
      
#         if j % 3 != 0:
#             continue
        
#         plt.cla()
        
#         for cluster in clusters:
#             mu = cluster['mu_k']
#             cov = cluster['cov_k']

#             eigenvalues, eigenvectors = np.linalg.eigh(cov)
#             order = eigenvalues.argsort()[::-1]
#             eigenvalues, eigenvectors = eigenvalues[order], eigenvectors[:, order]
#             vx, vy = eigenvectors[:,0][0], eigenvectors[:,0][1]
#             theta = np.arctan2(vy, vx)

#             color = colors.to_rgba(colorset[idx])

#             for cov_factor in range(1, 4):
#                 ell = Ellipse(xy=mu, width=np.sqrt(eigenvalues[0]) * cov_factor * 2, height=np.sqrt(eigenvalues[1]) * cov_factor * 2, angle=np.degrees(theta), linewidth=2)
#                 ell.set_facecolor((color[0], color[1], color[2], 1.0 / (cov_factor * 4.5)))
#                 ax.add_artist(ell)

#             ax.scatter(cluster['mu_k'][0], cluster['mu_k'][1], c=colorset[idx], s=1000, marker='+')
#             idx += 1

#         for i in range(X.shape[0]):
#             ax.scatter(X[i, 0], X[i, 1], c=colorset[np.argmax(scores[i])], marker='o')
        
#         fig.canvas.draw()
        
#         image = np.frombuffer(fig.canvas.tostring_rgb(), dtype='uint8')
#         image = image.reshape(fig.canvas.get_width_height()[::-1] + (3,))

#         images.append(image)
    
#     kwargs_write = {'fps':1.0, 'quantizer':'nq'}
#     imageio.mimsave('./gmm.gif', images, fps=1)
#     plt.show(Image.open('gmm.gif').convert('RGB'))
    
    
# create_cluster_animation(X, history, scores)



