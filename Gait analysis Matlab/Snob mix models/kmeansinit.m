% Based on Matlabs kmeans() function
function [idx, D] = kmeansinit(X, K)

[n,d] = size(X);
centres = zeros(K, d);

%% Choose initial centre
centres(1,:) = datasample(X,1,1);
minDist = inf(1,n); % Distance to the nearest centroid

%% Choose remaining centres according to kmeans++
for j = 2:K
    minDist = min(minDist,distance(centres((j-1),:)', X'));
    minDist(minDist<eps) = eps;
    denominator = sum(minDist);
    
    if denominator==0 || isinf(denominator) || isnan(denominator)
        centres(j:K,:) = datasample(X, K-j+1, 1,'Replace',false);
        break;
    end
    sampleProbability = minDist/denominator;    
    centres(j,:) = datasample(X, 1, 1, 'Replace', false, 'Weights', sampleProbability);            
end

%% Compute distance of every point to every centre and perform initial assignment
D = distance(X', centres');
[~, idx] = min(D, [], 2);

end

function f = distance(A, B)

f = bsxfun(@plus,dot(A,A,1)',dot(B,B,1))-2*(A'*B);

end

