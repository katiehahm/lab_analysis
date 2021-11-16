%% Example - Simulated data (mixtures of multivariate Gaussian distributions)
clear;

rng(1);

%% Generate data
K = 4;          % number of classes
D = 2;          % dimension of each data point
n = ones(1,K)*50;   % sample size of each class

muMin = -5;    % Each mu \in [muMin, muMax]
muMax = +5;
muRange = muMax - muMin;

data = [];
Mu = cell(K,1);
Sigma = cell(K,1);
for i = 1:K
    Sigma{i} = randcorr(D);
    Mu{i} = rand(D,1)*muRange + muMin;
    
    data  = [data; mvnrnd(Mu{i}, Sigma{i}, n(i))];
end

%% Mixture model
mm = snob(data, {'mvg',1:D}, 'k', 1);
mm_Summary(mm);

%% Print matrix of KL divergences
mm_KLstats(mm, data);
