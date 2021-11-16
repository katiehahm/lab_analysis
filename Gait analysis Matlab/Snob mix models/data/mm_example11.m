%% Example - Simulated data (mixtures of exponential distributions with fixed type I censoring)
clear;

% Seed the random number generator
rng(1);

%% We first generate some data with censoring.
% There are two classes:
C = 10;  % fixed censoring point
% (1) T ~ Exp(2), Y = min(T,C), Delta = I(T<=C)
n = 1e2;
T = exprnd(2, n, 1);
y1 = min(T, C);
delta1 = (T <= C)*1;    

% (2) T ~ Exp(1/5), Y = min(T,C), Delta = I(T<=C)
n = 50;
T = exprnd(1/5, n, 1);
y2 = min(T, C);
delta2 = (T <= C)*1;  

data = [[y1;y2], [delta1;delta2] ];

%% Run snob 

% We specify that we are dealing with censored exponential data.
% and start the search with k=1 classes.
mm = snob(data, {'cfixexp', [1,2]}, 'k', 1);

%% Lets look at the model snob discovered
mm_Summary(mm);

%% Plot the survival function
minX = 0;    % range of values for plotting
maxX = max(data(:,1)) + 10;
nPts = 1e3;
x = linspace(minX, maxX, nPts)';
y = zeros(nPts, 1);

for k = 1:mm.nClasses
    prop = mm.a(k); % mixing proportion
    theta = mm.class{k}.model{1}.theta;
    y = y + (prop * expcdf(x,theta));
end
Survival = 1 - y;

figure;
plot(x, Survival, 'k-');
grid;
xlabel('X', 'fontsize', 16);
ylabel('S(X)', 'fontsize', 16);
title('Mixture of exponential distributions', 'fontsize', 18);

%% Example - Simulated data (mixtures of Weibull distributions with fixed type I censoring)
% There are two classes:
C = 0.7;  % fixed censoring point
n = 100;
T = wblrnd(1, 1/2, n, 1);
y1 = min(T, C);
delta1 = (T <= C)*1;    

n = 150;
T = wblrnd(1, 5, n, 1);
y2 = min(T, C);
delta2 = (T <= C)*1;  

data = [[y1;y2], [delta1;delta2] ];

%% Run snob 
% We specify that we are dealing with censored exponential data.
% and start the search with k=1 classes.
mm = snob(data, {'cfixweibull', [1,2]}, 'k', 3);

%% Lets look at the model snob discovered
mm_Summary(mm);

%% Plot the survival function
minX = 0;    % range of values for plotting
maxX = max(data(:,1)) + 10;
nPts = 1e3;
x = linspace(minX, maxX, nPts)';
y = zeros(nPts, 1);

for k = 1:mm.nClasses
    prop = mm.a(k); % mixing proportion
    theta = mm.class{k}.model{1}.theta;
    scale = theta(1); shape = theta(2);
    y = y + (prop * wblcdf(x,scale,shape));
end
Survival = 1 - y;

figure;
plot(x, Survival, 'k-');
grid;
xlabel('X', 'fontsize', 16);
ylabel('S(X)', 'fontsize', 16);
title('Mixture of Weibull distributions', 'fontsize', 18);
