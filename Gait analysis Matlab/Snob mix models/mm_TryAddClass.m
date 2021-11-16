function [mm, msglen] = mm_TryAddClass(mm, data)

if(mm.opts.display)
    fprintf('             Attempting to add and re-assign...\n');
end  

K = mm.nClasses + 1;

%% Add a new class
mm.class{K} = mm_CreateClass(mm.ModelTypes);
mm.nClasses = K;

% run kmeans 
[~, Dist] = kmeansinit(data, K);

R  = bsxfun(@rdivide, Dist, sum(Dist,2));
% Heuristic: stop points being assigned 100% to a single class
R = bsxfun(@rdivide, 0.1+R, sum(0.1+R,2)); 
ix = sum(isnan(Dist),2) > 0;
if(any(ix))                
    t       = rand(sum(ix), K);
    R(ix,:) = bsxfun(@rdivide, t, sum(t,2));
end

mm.r  = R;
mm.Nk = sum(mm.r,1)';
mm.a = (mm.Nk+1/2) ./ (mm.N+K/2); %mm.a  = mm.Nk / sum(mm.Nk);

%% Re-estimate the model parameters for the two classes using the new memberships
if(min(mm.Nk) >= mm.MinMembers)
    mm = mm_EstimateTheta(mm, data, 1:mm.nClasses);
end

%% Run the EM procedure
mm = mm_EM(mm, data);
msglen = mm.msglen;

if(mm.opts.display)
    fprintf('             Best      [ADD]: msglen = %10.2f nits\n', msglen);
end  

end