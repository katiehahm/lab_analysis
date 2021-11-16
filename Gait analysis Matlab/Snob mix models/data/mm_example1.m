%% Example - Simulated data (mixtures of exponential distributions)
clear;

% Seed the random number generator
rng(1);

% Generate n=100 data points from a mixture of two exponential distributions 
% The first 60 points are from class #1 and the last 40 points from class #2
% x ~ 0.6 Exp(5) + 0.4 Exp(1)
x = [exprnd(5, 60, 1); exprnd(1, 40, 1)];

% Run Snob with the following options: 
%     (1) the data is modelled using a univariate exponential distribution: {'exp',1}
%     (2) Snob will automatically attempt to discover the optimal number of
%     mixtures (subpopulations)
%     (3) 'k', 5: is the starting number of classes
%     (4) 'display', false: disables printing during optimisation
mm = snob(x, {'exp',1},'k',5,'display',false);

% Print a summary of the discovered mixture model.
% The total message length of this model is 241.84 nits. 
% The model discovered is x ~ 0.42 Exp(0.9) + 0.58 Exp(6.4)
mm_Summary(mm);

% The mixing proportions of the mixture model are stored in mm.a
% The mean parameter for class i is stored in mm.class{i}.model{1}.theta

% Print Rand index and Adjusted Rand Index (ARI)
TrueInd = [ones(60,1); 2*ones(40,1)];   % True class assignments
[~,EstInd] = max(mm.r, [], 2);          % Estimated class assignments
[r, adjr] = RandIndex(TrueInd, EstInd); % Compute Rand index and ARI

fprintf('\n')
fprintf('*** Rand index = %5.3f, Adjusted rand index (ARI) = %5.3f\n', r, adjr);
fprintf('\n')

% What if we fit three sub-populations to this data?
% Run Snob with the following options: 
%     (1) the data is modelled using a univariate exponential distribution: {'exp',1}
%     (2) initial number of classes is 3: 'k',3
%     (2) Snob will NOT search for the best model structure: 'fixedstructure',true
mm2 = snob(x, {'exp',1},'k',3,'fixedstructure',true,'display',false);

% The total message length of the three class model is 243.19 nits.
mm_Summary(mm2);

% The two class model is preferred as it has a smalled message length.
fprintf('Two class exp model is %.2f more likely a posteriori than the three class exp model.\n', exp( -(mm.msglen - mm2.msglen) ));

% Next, we use a mixture of gamma distributions instead of exponentials.
mm3 = snob(x, {'gamma',1},'fixedstructure',true,'display',false);

% The two class exponential model is preferred as it has a smalled message length.
fprintf('Two class exp model is %.2f more likely a posteriori than the one-class gamma model.\n', exp( -(mm.msglen - mm3.msglen) ));
