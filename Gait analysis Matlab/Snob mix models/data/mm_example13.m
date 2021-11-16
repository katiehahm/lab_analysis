%% Example - Simulated data (mixtures of single factor analysis models)
clear;

rng(1);

%% Generate data
K = 5;          % number of classes
D = 5;          % dimension of each data point
n = ones(1,K)*100;   % sample size of each class

muMin = -15;    % Each mu \in [muMin, muMax]
muMax = +15;
muRange = muMax - muMin;

data = [];
A = cell(K,1);
Mu = cell(K,1);
Sigma = cell(K,1);
for i = 1:K
    
    A     = randn(D,1);
    A     = (A / norm(A)) * (1+rand(1)*5);
    sigma1 = randi(2,D,1);
    Sigma1 = A*A' + diag(sigma1.^2);
    Mu{i} = rand(D,1)*muRange + muMin;
    data     = [data; mvnrnd(Mu{i}, Sigma1, n(i))];
    
end

%% Mixture model
mm = snob(data, {'sfa',1:D}, 'k', 3);
mm_Summary(mm);
