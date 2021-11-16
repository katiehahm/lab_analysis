%% Example - Acidity data (real data)
clear;

% Load the data. 
% The data consists of a single numerical continuous variable.
%
% Reference: 
% Richardson, S. and Green, P. J. (1997). 
% On Bayesian analysis of mixtures with unknown number of components (with discussion). Journal of the Royal Statistical Society, Series B, 59, 731--792
load data/acidity;  

% Run Snob with the following options: 
%     (1) the data is modelled using a univariate Gaussian distribution: {'norm',1}
%     (2) Snob will automatically attempt to discover the optimal number of
%     mixtures (subpopulations)
mm = snob(acidity, {'norm',1}, 'k', 5, 'varnames', {'Acidity'});

% Print a summary of all the components (parameters and structure) of the
% mixture model we have discovered. Snob discovered two classes; one of the
% classes [N(mu = 4.3, sigma = 0.4)] has n~92 data points, the other
% class [N(mu = 6.3, sigma = 0.5)] has n~63. 
%
% The total message length of this model is ~209 nits.
mm_Summary(mm);

% The mixing proportions of the mixture model are stored in mm.a
% The [mean, variance] parameters for class i are stored in mm.class{i}.model{1}.theta

% Plot the mixture distribution for the acidity data
mm_PlotModel1d(mm, acidity, 1);

% We now force Snob to use 3 classes and turn off subpopulation discovery. 
mm2 = snob(acidity, {'norm',1}, 'k', 3, 'fixedstructure', true, 'varnames', {'Acidity'});

% Print a summary of the new model. 
% The total message length of the 3-class model is ~217 nits; ~5.4 nits longer than
% the two class model. 
mm_Summary(mm2);

% The two class model is preferred as it has a smalled message length.
% Specifically, the two class model is approximately
exp( -(mm.msglen - mm2.msglen) )
% times more likely a posteriori than the three class model.

% Plot the CDF of the fitted mixture model with 2 classes
minX = min(acidity) - 1;    % range of values for plotting
maxX = max(acidity) + 1;
nPts = 1e3;
x = linspace(minX, maxX, nPts)';
y = zeros(nPts, 1);

for k = 1:mm.nClasses
    prop = mm.a(k); % mixing proportion
    theta = mm.class{k}.model{1}.theta; % parameters for class k
    mu = theta(1);  % mean and std. deviation
    sigma = sqrt(theta(2));
    y = y + (prop * normcdf(x,mu,sigma));   % CDF
end

figure;
plot(x, y, 'k-');
grid;
xlabel('X', 'fontsize', 16);
ylabel('F(X)', 'fontsize', 16);


% Now let's try fitting the same data using a mixutre model of
% inverse Gaussian distributions...
mm_ig = snob(acidity, {'igauss',1}, 'k', 5, 'varnames', {'Acidity'});

% Print a matrix of KL divergences for both the Gaussian and 
% inverse Gaussian models
mm_KLstats(mm, acidity);
mm_KLstats(mm_ig, acidity);
