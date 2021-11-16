%% Example - Simulated data (mixtures of beta distributions)
clear;

% Seed the random number generator
rng(1);

% Generate n=100 data points from a mixture of two beta distributions 
% x ~ 0.5 beta(30,10) + 0.5 beta(10,30)
x=[betarnd(30,10,50,1); betarnd(10,30,50,1)];

% Run Snob with the following options: 
%     (1) the data is modelled using a univariate beta distribution: {'beta',1}
%     (2) Snob will automatically attempt to discover the optimal number of
%     mixtures (subpopulations)
%     (3) We start the search with 5 mixtures {'k',5}
mm = snob(x, {'beta',1},'k',5);

% Print a summary of the discovered mixture model.
% The model discovered is x ~ 0.5 beta(32,11) + 0.5 beta(12,38)
mm_Summary(mm);

% Plot the PDF of the mixture model
mm_PlotModel1d(mm, x, 1);

% Plot the CDF of the fitted mixture model 
minX = 0;    % range of values for plotting
maxX = 1;
nPts = 1e3;
x = linspace(minX, maxX, nPts)';
y = zeros(nPts, 1);

for k = 1:mm.nClasses
    prop = mm.a(k); % mixing proportion
    theta = mm.class{k}.model{1}.theta; % parameters for class k
    a = theta(1);  % mean and std. deviation
    b = theta(2);
    y = y + (prop * betacdf(x,a,b));   % CDF
end

figure;
plot(x, y, 'k-');
grid;
xlabel('X', 'fontsize', 16);
ylabel('F(X)', 'fontsize', 16);