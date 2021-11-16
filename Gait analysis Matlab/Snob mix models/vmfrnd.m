%VMFRND    Generate random numbers from the von Mises-Fisher distribution
%   VMFRND(kappa,mu,n) generates n random numbers from a p-dimensional
%   von Mises-Fisher (VMF) distribution with mean parameter mu 
%   and concentration parameter kappa > 0.
%
%  TODO: vectorize code
%
%  The input arguments are:
%   kappa    - [1 x 1] {kappa > 0} concentration parameter
%   mu       - [p x 1] {||mu|| = 1} mean parameter
%   n        - [1 x 1] {n > 0} number of parameters
%
%  Returns:
%   f    - [n x p] samples from the VMF distribution
%
% References:
%
% Wood, A. T. A., 
% Simulation of the von Mises Fisher distribution, 
% Communications in Statistics - Simulation and Computation, 23 , 157-164 (1994).
%
% Hornik, Kurt and Grün, Bettina, 
% movMF: An R Package for Fitting Mixtures of von Mises-Fisher Distributions, 
% Journal of Statistical Software, 58 (10), pp. 1-31 (2014).
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
function f = vmfrnd(kappa, mu, n)

%% Setup
mu = mu(:);
d = length(mu); % dimension

%% constants
b = (d - 1) / (2*kappa + sqrt(4*kappa*kappa + (d-1)*(d-1)));
x0 = (1 - b) / (1 + b);
c = kappa*x0 + (d-1)*log(1 - x0*x0);

%% Samples
f = zeros(n, d);
m = (d - 1)/2;

for i = 1:n
    t = -inf; 
    u = 1;
    while(t < log(u))
        z = betarnd(m, m);
        u = rand(1);
        w = (1 - (1 + b)*z) / (1 - (1 - b)*z);
        t = kappa*w + (d-1)*log(1 - x0*w) - c;
    end
    
    v = randn(d-1,1);
    v = v / norm(v);
    
    f(i,:) = [sqrt(1 - w*w)*v', w];    
end

Q = null(mu');
f = ([Q, mu]*f')';

end