%% MML_CONST    Computes the WF87 quantisation constant.
% 
% Computes (k/2) log(kappa_k + 1), where kappa_k is the mean squared
% quantisation error of an optimal quantising lattice in k
% dimensions.
% 
%  The input arguments are: 
%   k    - number of free parameters (k > 0)
% 
%  Returns:
%   f    - WF87 quantisation constant for dimension k (>0)
%
% References:
% [1] Agrell, Erik and Thomas Eriksson. 
% “Optimization of Lattices for Quantization.” IEEE Trans. Inf. Theory 44 (1998): 1814-1828.
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
function f = mml_const(k)

% MSE
kappa_k = [1 / 12               % k = 1
     5 / (36 * sqrt(3))         % k = 2
     19 / (192 * 2^(1/3))       % k = 3
     0.076603                   % k = 4
     0.075625                   % k = 5
     0.074244                   % k = 6
     0.073116                   % k = 7
     0.071682                   % k = 8
     0.071626                   % k = 9
     0.070814                   % k = 10
    ];

% WF87 constant
if(k < 10)
    f = (k/2)*(log(kappa_k(k)) + 1);
else
    f = -k*log(2*pi)/2 + log(k*pi)/2 + psi(1);
end

end