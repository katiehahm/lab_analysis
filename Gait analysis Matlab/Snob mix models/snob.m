%SNOB    Minimum Message Length mixture modelling
%   A finite mixture model is a statistical model that attempts to discover latent
%   (unobserved) subpopulations within an overall population. In finite
%   mixture models, we do not know the number of subpopulations nor do we
%   know which observations belong to which subpopulation. The aim is 
%   to discover the number of subpopulations and estimate parameters
%   associated with each subpopulation.
%
%   SNOB is a MATLAB implementation of finite mixture models that uses
%   the minimum message length criterion to determine the structure of the
%   mixture model (i.e., the number of subpopulations) and estimate all
%   corresponding parameters (see e.g., [1-5]). SNOB allows the user to
%   specify the desired number of subpopulations; if this is not specified, SNOB will
%   automatically try to discover this information. Currently, SNOB
%   supports mixtures of the following distributions: 
%
%   (.)    Univariate Gaussian distribution ('norm')
%               p(X|mu, sigma^2) = (1/sqrt(2*pi*sigma^2)) * exp(-(x-mu)/2/sigma^2),
%                                   X \in R, mu \in R, sigma > 0
%   (.)    Multivariate Gaussian distribution ('mvg')
%               p(X|mu, Sigma)   = (2*pi)^(-d/2) * det(Sigma)^(-1/2) * exp(-0.5*(x-mu)' inv(Sigma) (x-mu)),
%                                   X \in R^d, mu \in R^d, Sigma positive-definite
%   (.)    Weibull distribution ('weibull')
%               p(X|lambda, k)   = (k/lambda) * (x/lambda)^(k-1) * exp(-(x/lambda)^k),
%                                   X > 0, lambda > 0, k > 0
%   (.)    Weibull distribution with fixed type I censoring ('cfixweibull')
%               Data = [Y, delta], where
%               Y ~ min(T,C); delta = I(T <= C)
%               p(T|lambda, k)   = (k/lambda) * (t/lambda)^(k-1) * exp(-(t/lambda)^k),
%                                   T > 0, lambda > 0, k > 0
%   (.)    Exponential distribution ('exp')
%               p(X|lambda)      = (1/lambda) * exp(-x/lambda),
%                                   X > 0, lambda > 0
%   (.)    Exponential distribution with random type I censoring ('crndexp')
%               Data = [Y, delta], where
%               Y ~ min(T,C); delta = I(T <= C)
%               p(T|beta)        = (1/beta) * exp(-t/beta),
%               p(C|alpha)       = (1/alpha) * exp(-c/alpha),
%                                   T > 0, C > 0, alpha,beta > 0
%   (.)    Exponential distribution with fixed type I censoring ('cfixexp')
%               Data = [Y, delta], where
%               Y ~ min(T,C); delta = I(T <= C)
%               p(T|theta)        = (1/theta) * exp(-t/theta),
%                                   T > 0, C > 0, theta > 0
%   (.)    von Mises-Fisher distribution ('vmf')
%               p(X|kappa,mu)    = kappa^((d/2)-1)/(2*pi)^(d/2)/I_{d/2-1}(kappa) exp(kappa * mu' x)
%                                   ||X|| = ||mu|| = 1, kappa > 0, X,mu \in R^d
%   (.)    Inverse Gaussian distribution ('igauss')
%               p(X|mu,lambda)   = sqrt(1/2/pi/lambda/x^3) * exp(-(x-mu)^2/2/mu^2/x/lambda),
%                                   X > 0, mu > 0, lambda > 0
%   (.)    Poisson distribution ('poisson')
%               p(X|lambda)      = lambda^x * exp(-lambda) / x!,
%                                   X \in {0,1,2,3,...}, lambda > 0
%   (.)    Negative binomial ('negb')
%               p(X|r,p)         = nchoosek(X+r-1, X) (1-p)^r p^X
%                                   X \in {0,1,2,3,...}, r > 0, 0 < p < 1
%   (.)    Geometric distribution ('geometric')
%               p(X|theta )      = (1 - theta)^x * theta
%                                   X \in {0,1,2,3,...}, 0 < theta < 1
%   (.)    Multinomial distribution ('multi')
%               p(X|p)           = n! / (x1! * ... xk!) * p.^x, 
%                                   n = sum(x), sum(p) = 1
%   (.)    Gamma distribution ('gamma')
%               p(X|mu,phi)      = 1/gamma(phi)/(mu/phi)^phi * x^(phi-1) * exp(-x/(mu/phi)),
%                                  X > 0, mu > 0, phi > 0
%               Parametarization (alpha/shape,beta/rate):
%                                  alpha = phi, beta = phi / mu
%               Parametarization (k/shape,theta/scale):
%                                  k = phi, theta = mu / phi
%   (.)    Laplace distribution ('laplace')
%               p(X|mu,b)        = 1/2/b * exp(-abs(x - mu)/b),
%                                   X \in R, mu \in R, b > 0
%   (.)    Single factor analysis model ('sfa') 
%               x_nk             = mu_k + v_n a_k + sigma_k*r_nk,
%                                   {v_n, {r_nk, k=1,...,K},n=1,...,N} ~ N(0,1)
%   (.)    Beta distribution ('beta')
%               p(X|a,b) = x^(a-1) (1-x)^(b-1) / Beta(a,b)
% 
%   (.)    Gaussian linear regression ('linreg')
%               p(Y|X,theta)     = Gaussian(b0 + x'*b, sigma^2)
%                                   Y \in R, b0 \in R, b \ in R^d, sigma>0
%   (.)    Logistic regression ('logreg')
%               p(Y==1|X,theta)  = logit(b0 + x'*b)
%                                   Y \in {0,1}, b0 \in R, b \ in R^d              
%
%  The input arguments to SNOB() are:
%   
%   x           - [n x p] data set (n samples; p variables)
%   model_list  - [2k x 1] cell vector specifying k models for the data columns
%                 The format is: {'model',cols, 'model',cols, etc}
%       
%                 'model' is one of:
%                       'beta'      -> Beta distribution
%                       'exp'       -> Exponential distribution
%                       'crndexp'   -> Exponential distribution with random type I censoring
%                       'cfixexp'   -> Exponential distribution with fixed type I censoring
%                       'gamma'     -> Univariate gamma distribution
%                       'geometric' -> Geometric distribution
%                       'igauss'    -> Inverse Gaussian distribution
%                       'laplace'   -> Univariate Laplace distribution
%                       'linreg'    -> Gaussian linear regression
%                       'logreg'    -> Logistic regression
%                       'multi'     -> Multinomial distribution
%                       'mvg'       -> Multivariate normal distribution 
%                       'negb'      -> Negative binomial distribution
%                       'norm'      -> Univariate normal distribution
%                       'poisson'   -> Poisson distribution
%                       'sfa'       -> Multivariate normal distribution (single factor analysis)
%                       'vmf'       -> von Mises-Fisher distribution
%                       'weibull'   -> Weibull distribution
%
%                Except in the case of censored data, 'linreg' or 'logreg' models, the vector cols denotes which columns the
%                model applies to. That is, 
%                       {'norm', 1} -> the first column is data from a normal distribution
%                       {'exp',1,'poisson',[2,3],mvg,[4,5,6]} -> 
%                       column 1 is data from the exponential distribution;
%                       columns 2,3 are data from two Poisson distributions
%                       columns 4,5,6 are data from a single multivariate Gaussian distribution
%
%                For distributions with censored data
%                       {'cfixexp', [1,2], 'cfixexp', [3,4]} -> 
%                       column 1 is the data Y, 
%                       column 2 are the indicator variables with delta = 1
%                       (observed; T<=C) and delta = 0 (censored; T>=C)
%                       ... and similarly for columns 3 and 4.
%
%                In the case of the linear or logistic regression models, cols = [target, covariate(s)]. That is, 
%                       {'mvg',1:3,'linreg',[4,1,2,3]} -> the dependent variable is in column 4, 
%                                               the independent variables are in columns 1, 2 and 3
%                                               the independent variables are modelled by a single multivariate 
%                                               Gaussian distribution 
%                                               the dependent variable is follows a linear regression model 
%
%   args        - optional arguments:
%       'k', integer            - how many classes to start with [default=1]
%       'maxk', integer         - maximum number of classes [1 < k <= maxk]
%       'startmodel', struct    - start search from this mixture model [default=[]]
%       'fixedstructure', bool  - allow SNOB to attempt splitting/merging of classes? [default=false]
%       'display', bool         - show search progress [default=true]
%       'emmaxiter', integer    - maximum number of EM iterations [default=100]
%       'maxiter',integer       - maximum number of search iterations [default=100]
%       'maxtrycombine',integer - maximum number of classes to attempt to
%                                 split/combine during each search iteration [default=10]
%       'varnames', cell array  - cell array of strings, one for each variable
%                                 [default={'','',...}]
%       
%  Returns:
%   mm    - structure respresenting the complete mixture model
%
%   Key variables:
%   nClasses     - number of classes
%   nModelTypes  - number of attributes (eg of attributes: gender, height, etc)
%   a            - estimated mixing proportions (sum up to 1)
%   L            - negative log-likelihood of the data under this mixture
%   msglen       - MML codelength of the estimated mixture model
%   class{i}     - each class consists of M attributes which can be accesed
%                  with the model substructure. For example:
%                  mm.class{1}.model{1} is the statistical model for attribute 1 in class 1. 
%                  mm.class{2}.model{1} is the statistical model for attribute 1 in class 2. 
%                  etc.
%                  Parameters of the model are found in the field theta (see below). 
%
%  How to interpret the parameter field "theta":
%  'beta'        - [a, b] 
%  'exp'         - lambda 
%  'gamma'       - [mu, phi]
%                  Parametarization (alpha/shape,beta/rate):
%                      alpha = phi, beta = phi / mu
%                  Parametarization (k/shape,theta/scale):
%                      k = phi, theta = mu / phi
%  'geometric'   - theta
%  'igauss'      - [mu, lambda]
%  'norm'        - [mu, sigma^2]
%  'laplace'     - [mu, b]
%  'linreg'      - [sigma^2, b0, b']
%  'logreg'      - [b0, b']
%  'multi'       - [p0,p1,...] s.t. p0+p1+...=1
%  'mvg'         - [mu', Sigma(:)']
%  'negb'        - [mu, phi] where r = phi, p = 1-mu/(mu+phi)
%  'poisson'     - lambda
%  'sfa'         - [mu', sigma', a']
%  'weibull'     - [lambda,k]
%  'vmf'         - [kappa, mu']
%
%
%  Examples:
%  For examples of SNOB usage, please see data/mm_example1.m, data/mm_example2.m, etc.
%
% References:
% [1] Wallace, C. S. & Dowe, D. L.
%     MML clustering of multi-state, Poisson, von Mises circular and Gaussian distributions 
%     Statistics and Computing, 2000, 10, 73-83 
%
% [2] Wallace, C. S.
%     Intrinsic Classification of Spatially Correlated Data 
%     The Computer Journal, 1998, 41, 602-611 
%
% [3] Wallace, C. S.
%     Statistical and Inductive Inference by Minimum Message Length 
%     Springer, 2005
%
% [4] Schmidt, D. F. & Makalic, E.
%     Minimum Message Length Inference and Mixture Modelling of Inverse Gaussian Distributions 
%     AI 2012: Advances in Artificial Intelligence, Springer Berlin Heidelberg, 2012 , 7691 , 672-682 
%
% [5] Edwards, R. T. & Dowe, D. L.
%     Single factor analysis in MML mixture modelling 
%     Research and Development in Knowledge Discovery and Data Mining, Second Pacific-Asia Conference (PAKDD-98), 1998 
%
% [6] Parthan Kasarapu & Lloyd Allison
%     Minimum message length estimation of mixtures of multivariate Gaussian and von Mises-Fisher distributions
%     Mach Learn (2015) 100:333–378
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
function mm = snob(data, model_list, varargin)

%% Version number
VERSION = '0.50';

%% Parse options
inParser = inputParser;  

% Default parameter values
defaultMaxiter   = 100; 
defaultEMMaxiter = 100;
defaultKinit     = 1;
defaultToDisplay = true;
defaultInit      = 'kmeans++';
defaultGreedy    = false;
defaultFixedStructure = false;
defaultMaxTryCombine  = 10;
defaultStartModel     = [];
defaultNoCost         = []; % TODO: add option to avoid costing specified columns
defaultMaxK           = inf;
defaultVarNames = {};

expectedInit = {'random', 'kmeans++'};

% Required arguments
addRequired(inParser, 'model_list', @(x)(iscell(x) && ~isempty(x)));

% Optional arguments
addParameter(inParser, 'maxiter', defaultMaxiter, @(x) isnumeric(x) && isscalar(x) && (x > 1));
addParameter(inParser, 'emmaxiter', defaultEMMaxiter, @(x) isnumeric(x) && isscalar(x) && (x > 1));
addParameter(inParser, 'k',  defaultKinit, @(x) isnumeric(x) && isscalar(x) && (x > 0));
addParameter(inParser, 'display', defaultToDisplay, @islogical);
addParameter(inParser, 'init', defaultInit, @(x) any(validatestring(x, expectedInit)));
addParameter(inParser, 'fixedstructure', defaultFixedStructure, @islogical);
addParameter(inParser, 'greedy', defaultGreedy, @islogical);
addParameter(inParser, 'maxtrycombine',  defaultMaxTryCombine, @(x) isnumeric(x) && isscalar(x) && (x > 0));
addParameter(inParser, 'startmodel', defaultStartModel, @isstruct);
addParameter(inParser, 'nocost', defaultNoCost, @(x) isnumeric(x));
addParameter(inParser, 'maxk', defaultMaxK, @(x) isnumeric(x) && isscalar(x) && (x > 0));
addParameter(inParser, 'varnames', defaultVarNames, @iscell);

% Parse the input now
parse(inParser, model_list, varargin{:});  

% Get results of the parsing
models    = inParser.Results.model_list;  % list of models
mm        = inParser.Results.startmodel;  % did we get passed a model?

opts = struct;                                      % hold all options here
opts.nClasses       = inParser.Results.k;           % starting number of classes [default k=1]
opts.maxk           = inParser.Results.maxk;        % maximum number of classes [default kmax=inf]
opts.Initialisation = inParser.Results.init;        % initialisation strategy [default 'kmeans++']
opts.maxiter        = inParser.Results.maxiter;     % maximum number of search iterations [default maxiter=100]
opts.emmaxiter      = inParser.Results.emmaxiter;   % maximum number of EM iterations for a fixed model structure [default emmaxiter=100]
opts.display        = inParser.Results.display;     % verbose output [default display=true]
opts.fixedstructure = inParser.Results.fixedstructure;     % do we attempt merger/split combos? [default fixedstructure=false]
opts.greedy         = inParser.Results.greedy;      % if true, always pick the model with smallest message length; 
opts.ModelList      = model_list;
                                                    % if false, pick model stochastically [default]
opts.MaxTryCombines = inParser.Results.maxtrycombine;
opts.SearchOptions  = optimoptions('fminunc','display','off');
opts.NoCost         = inParser.Results.nocost;
opts.VarNames       = upper(inParser.Results.varnames);


%% Check Options
if(any(mod(opts.nClasses,1)))
    error('k must be an integer greater than 0 and less than maxk');
end
if(opts.nClasses > opts.maxk)
    error('k must be an integer greater than 0 and less than maxk');
end
if(any(mod(opts.maxk,1)))
    error('maxk must be an integer greater than 0');
end

%% Process model types
if(any(isinf(data(:))))
    error('Inf values detected in data')
end

[~,p] = size(data);
VarsUsed = false(p,1);
ModelTypes = {};
if (~isempty(models))
    i=1; Ix=1;
    while (i <= length(model_list))
        if(all(VarsUsed))
            error('Too many models specified');
        end
        [i, Ix, ModelTypes, VarsUsed] = mm_ProcessModelTypes(model_list, i, Ix, ModelTypes, data, VarsUsed);
    end
    
    % If insufficient models were added to cover all data columns, return an error
    if ~all(VarsUsed)
        error('Too few models specified -- some data columns are not modelled');
    end
end
opts.nModels = length(ModelTypes);

%% Check list of attribute names
if(~isempty(opts.VarNames))    
    if(length(opts.VarNames) ~= p)
        error('List of attribute names is incorrect size');
    end

    for i = 1:p
        if(~ischar(opts.VarNames{i}))
            error('All attribute names must be of type string');
        end
    end
    
    % Shorten strings to max 9 char length
    for i = 1:p
        chrLen = length(opts.VarNames{i});
        if(chrLen > 8)
            opts.VarNames{i} = [opts.VarNames{i}(1:8),'.'];
        end
    end
    
else % Create an empty list as no names were given
    for i = 1:p
        opts.VarNames{i} =  ['Attr.', num2str(i)];
    end
end


%% Create an initial model
if(isempty(mm))
    mm = mm_Create(data, ModelTypes, opts);
end

if(opts.display)
    maxlen = 0;
    fprintf('%s\n', repchar('=', maxlen + 85));
    fprintf('|%s|\n', centrestr(sprintf('MML mixture modelling ver. %s', VERSION), maxlen + 83));
    fprintf('|%s|\n', centrestr(sprintf('(c) Enes Makalic, Daniel F Schmidt. 2019-'), maxlen+83));
    fprintf('%s\n', repchar('=', maxlen + 85));    
end

%% Search for the best model
mm = mm_Search(mm, data);
    
end