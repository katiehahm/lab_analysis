function [r, ep] = mm_EstimateR(mm, data)

% Normalize to get the posterior probability of assignment to classes
% ep = exp(-mm_Likelihood(mm, data, 1:mm.nModelTypes));
% ep(ep==0) = realmin; 
% r = bsxfun(@rdivide, ep, sum(ep,2));

% Normalize to get the posterior probability of assignment to classes
p = -mm_Likelihood(mm, data, 1:mm.nModelTypes);
p = max(min(p, 700), -700); % log(p) \in [-700, 700]
ep = exp(p);
%ep(ep == 0) = realmin;
if(mm.nClasses > 1)
    r = bsxfun(@rdivide, ep, exp(logsumexp(p)));
else
    r = ones(mm.N, 1);
end

end