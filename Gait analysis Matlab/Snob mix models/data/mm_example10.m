%% Example - Simulated data (mixtures of exponential distributions with random type I censoring)
clear;

% Seed the random number generator
rng(1);

%% We first generate some data with censoring.
% There are two classes:
% (1) T ~ Exp(1), C ~ Exp(5), Y = min(T,C), Delta = I(T<=C)
n = 50;
a = 1; b = 5;       % censored proportion ~ 1/6
T = exprnd(b, n, 1);
C = exprnd(a, n, 1);
y1 = min(T, C);
delta1 = (T <= C)*1;    

% (2) T ~ Exp(5), C ~ Exp(10), Y = min(T,C), Delta = I(T<=C)
n = 75;
a = 5; b = 10;      % censored proportion ~ 1/3
T = exprnd(b, n, 1);
C = exprnd(a, n, 1);
y2 = min(T, C);
delta2 = (T <= C)*1;  

data = [[y1;y2], [delta1;delta2] ];

%% Run snob 
% We specify that we are dealing with censored exponential data.
% and start the search with k=1 classes.
mm = snob(data, {'crndexp', [1,2]}, 'k', 1);

%% Lets look at the model snob discovered
mm_Summary(mm);