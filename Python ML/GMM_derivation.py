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

data = pd.read_csv('C:/Users/Katie/Dropbox (MIT)/Lab/Analysis/Experiment4/Praneeth 5/ProcessedData/both_weight2_p2_GMM.csv',header=None)
data = data.apply(pd.to_numeric, errors='coerce')
# convert from pandas dataframe to numpy array
dataM = pd.DataFrame(data).to_numpy()
steps = dataM[:,0]

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
        
        pi_k = 0.5 # CHANGED #############################
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
data = steps
data = np.reshape(data, (-1,1))
clusters, likelihoods, scores, sample_likelihoods, history = train_gmm(data, n_clusters, n_epochs)

# plt.figure(figsize=(10, 10))
# plt.title('Log-Likelihood')
# plt.plot(np.arange(1, n_epochs + 1), likelihoods)
# plt.show()

gmm = GaussianMixture(n_components=n_clusters, max_iter=50).fit(data)
gmm_scores = gmm.score_samples(data)

mean_results = np.array([cluster['mu_k'].tolist() for cluster in clusters])
pi_results = np.array([cluster['pi_k'] for cluster in clusters])
# output = [gmm.means_[0,0], gmm.means_[1,0], gmm.weights_[0], gmm.weights_[1], mean_results[0,0], mean_results[1,0], pi_results[0], pi_results[1]] # gmm1, gmm2, scal1, scal2, real1, real2
output = [mean_results[0,0], mean_results[1,0]] # gmm1, gmm2, scal1, scal2, real1, real2
print(output)





