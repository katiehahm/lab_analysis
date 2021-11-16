clear;

%% Load the old failthful data set
% The data consists of two numerical continuous variables.
%
% Reference: 
% Azzalini, A. and Bowman, A. W. (1990). A look at some data on the Old Faithful geyser. Applied Statistics 39, 357–365.
load data/oldfaithful;

%% Cluster waiting times using snob
waiting = oldfaithful(:,2);
mm = snob(waiting, {'norm', 1},'k',2,'varnames',{'waiting'});

%% Summary
% Print a summary of all the components (parameters and structure) of the
% mixture model we have discovered. Snob discovered two classes
%
mm_Summary(mm);

% Plot the mixture distribution
mm_PlotModel1d(mm, waiting, 1);

%% Plot the CDF of the fitted mixture model with 2 classes
minX = min(waiting) - 10;    % range of values for plotting
maxX = max(waiting) + 10;
nPts = 1e3;
x = linspace(minX, maxX, nPts)';
y = zeros(nPts, 1);

for k = 1:mm.nClasses
    prop = mm.a(k); % mixing proportion
    theta = mm.class{k}.model{1}.theta;
    mu = theta(1);  % mean and std. deviation
    sigma = sqrt(theta(2));
    y = y + (prop * normcdf(x,mu,sigma));
end

figure;
plot(x, y, 'k-');
grid;
xlabel('X', 'fontsize', 16);
ylabel('F(X)', 'fontsize', 16);