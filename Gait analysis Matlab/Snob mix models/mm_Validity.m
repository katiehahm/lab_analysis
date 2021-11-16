%MM_VALIDITY    Computes cluster validity indices for evaluating clustering 
%   MM_VALIDITY(mm,data,distance) computes two internal cluster
%   validity indices (Dunn Index; DI and an index based on the KS test [1]; KS_sym)
%   for evaluating the quality of clustering.
%  
%   Roughly speaking, the higher the DI/KS_Sym, the better the clustering.
%
%   The input arguments are:
%   mm       - [struct] mixture model
%   data     - [n x p] data set the model was trained on
%   distance - [string] type of distance; same as the distance argument to
%              Matlab function pdist().
%              e.g. 'euclidean','cityblock','correlation',etc.
%
%   Returns:
%   DI        - [1 x 2] Dunn Index based on the mean and min distance
%   KS_Sym    - [1 x 1] KS similarity test [1]
%
%   References:
%   [1] An Internal Cluster Validity Index Using a Distance-based Separability Measure
%       Shuyue Guan; Murray Loew
%       2020 IEEE 32nd International Conference on Tools with Artificial Intelligence (ICTAI)
%       DOI: 10.1109/ICTAI50040.2020.00131
%
%  (c) Copyright Enes Makalic and Daniel F. Schmidt, 2019-
function [DI,KS_sim] = mm_Validity(mm, data, distance)   

[~,cid] = max(mm.r,[],2);                           % determine which class each point belongs to

%% Compute validity indices
% Numerator:
% 1. Minimum mean distance between all pairs of points belonging to different
% clusters
% 2. Minimum distance between all pairs of points belonging to different
% clusters
%
% Denominator:
% 1. Maximum mean distance between all pairs of points in the same cluster
% 2. Maximum distance between all pairs of points in the same cluster
%
numerator = inf(1,2);
denominator = -inf(1,2);
s = zeros(mm.nClasses,1);

for k = 1:mm.nClasses
    ix1 = cid == k; % all points in cluster k
    ix2 = (~ix1);   % all other points NOT in cluster k
    
    Tbw = pdist2(data(ix1,:),data(ix2,:), distance);  % between cluster distance
    numerator(1) = min([nanmean(Tbw(:)), numerator(1)]);
    numerator(2) = min([nanmin(Tbw(:)), numerator(2)]);
    
    Twn = pdist(data(ix1,:), distance);               % within cluster distance
    denominator(1) = max([nanmean(Twn(:)), denominator(1)]);
    denominator(2) = max([nanmax(Twn(:)), denominator(2)]);  
    
    % KS Test [1]
    [~,~,s(k)] = kstest2(Tbw(:), Twn(:));
    
end

DI = numerator ./ denominator;

% Note: unlike in [1], we average the KS test statistics
KS_sim = nanmean(s);   

end