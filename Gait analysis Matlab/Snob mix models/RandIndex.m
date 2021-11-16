%RANDINDEX    Computes the Rand index and the Adjusted Rand Index (ARI)
%  RANDINDEX(.) computes the Rand index and the Adjusted Rand Index (ARI)
%  Code uses 
%  https://stackoverflow.com/questions/15793172/efficiently-generating-unique-pairs-of-integers
%  for generating all unique pairs.
%  
%  The input arguments are:
%   X    - partition, where each element of X is an integer in [1,numClassesX]
%   Y    - partition, where each element of Y is an integer in [1,numClassesY]
%
% (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
function [r, adjr] = RandIndex(X, Y)

%% Error checking
if(~isvector(X) || ~isvector(Y))
    error('The two inputs must be column or row vectors');
end
if(length(X) ~= length(Y))
    error('Input vectors must be the same length')
end

% Convert to column vectors
X = X(:); Y = Y(:);

%% Count pairs
N = length(X);
a = 0; b = 0; c = 0; d = 0;

if(N > 1e3) % slower, but uses little memory
    for i = 1:N-1
        for j = i+1:N
            a = a + (X(i) == X(j) && Y(i) == Y(j));
            b = b + (X(i) ~= X(j) && Y(i) ~= Y(j));
            c = c + (X(i) == X(j) && Y(i) ~= Y(j));
            d = d + (X(i) ~= X(j) && Y(i) == Y(j));
        end
    end
else
    % An alternative way to generate all pairs [p, q]
    % Copied from https://stackoverflow.com/questions/15793172/efficiently-generating-unique-pairs-of-integers
    % k = 1:21;
    % q = floor(sqrt(8*(k-1) + 1)/2 + 3/2);
    % p = k - (q-1).*(q-2)/2;    

    % Faster, but uses a lot of memory
    k = 1:(N/2*(N-1));
    q = floor(sqrt(8*(k-1) + 1)/2 + 3/2);
    p = k - (q-1).*(q-2)/2;

    a = sum(X(p) == X(q) & Y(p) == Y(q));
    b = sum(X(p) ~= X(q) & Y(p) ~= Y(q));
    c = sum(X(p) == X(q) & Y(p) ~= Y(q));
    d = sum(X(p) ~= X(q) & Y(p) == Y(q));

end

M = N/2*(N-1);  % Total number of unique pairs

%% Rand index
% Formula: https://en.wikipedia.org/wiki/Rand_index
r = (a+b) / M;

%% Adjusted rand index (ARI)
% Formula: https://www.sciencedirect.com/topics/computer-science/adjusted-rand-index
num = (a - (a + c)*(a + d) / M);
den = (c + d + 2*a) / 2 - (a + c)*(a + d)/M;
adjr = num / den;

end