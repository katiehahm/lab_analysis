%% Example - Diabetes data
clear;

% Load the data. 
% The data consists of 10 predictors (X = SEX, AGE, BMI, BP and six blood serum measurements S1, S2, S3, S4, S5, S6) of diabetes (Y). 
load data/diabetes;
X = [SEX, AGE, BMI, BP, S1, S2, S3, S4, S5, S6];
[~,p] = size(X);
data = [X, Y];

% Run Snob with the following options: 
%     (1) the data is modelled as follows: {'multi',1,'mvg',2:4,'mvg',5:10,'linreg',[11,1:p]}
%               'multi',1           column 1 (SEX) follows a multinomial distribution
%               'mvg',2:4           columns 2,3,4 (AGE,BMI,BP) follow a multivariate Gaussian distribution
%               'mvg',5:10          columns 5:10 (S1-S6) follow a multivaraite Gaussian distribution
%               linreg',[11,1:p]    column 11 is a Gaussian linear regression; covariates are columns 1:10
%     (2) initial number of classes is 4 ('k',4)
%     (2) Snob will automatically attempt to discover the optimal number of
%     mixtures (subpopulations)
VarNames = {'SEX', 'AGE', 'BMI', 'BP', 'S1', 'S2', 'S3', 'S4', 'S5', 'S6', 'DIABETES'};
mm = snob(data, {'multi',1,'mvg',2:4,'mvg',5:10,'linreg',[11,1:p]},'k',4,'varnames', VarNames);

% Print a summary of all the model parameters.
% There are two classes. The total message length is ~14,876 nits.
mm_Summary(mm);

% Now, suppose we model AGE,BMI and BP as uncorrelated [spherical Gaussian]
mm2 = snob(data, {'multi',1,'norm',2:4,'mvg',5:10,'linreg',[11,1:p]},'k',4,'varnames', VarNames);
mm_Summary(mm2);

% The codelengths are...
[mm.msglen, mm2.msglen]

% The uncorrelated model is 
exp( -(mm2.msglen - mm.msglen) )
% times more likely a posteriori than the MVG model.
