function mm = mm_Create(data, ModelTypes, opts)

%% Structure to holde all mixture model information
mm.name = 'Mixture Model';
mm.nClasses = opts.nClasses;
mm.ModelTypes = ModelTypes;
mm.nModelTypes = length(ModelTypes);
mm.N = size(data,1);
mm.opts = opts;

%% Get minimum members required for each class
minmembers = zeros(mm.nModelTypes, 1);
for i = 1:mm.nModelTypes
    minmembers(i) = mm.ModelTypes{i}.MinMembers;
end
mm.MinMembers = max(minmembers);    % minimum items for each class

K = opts.nClasses;

%% Initialise all parameters using specified algorithm
switch opts.Initialisation
    
    %% Random assignment of data
    case 'random'
        
        % Mixing proportions
        mm.a = rand(K,1); 
        mm.a = mm.a ./ sum(mm.a);

        % Create parameters for each model type and class
        for k = 1:K
            mm.class{k} = mm_CreateClass(ModelTypes);
        end

        % Randomly assign data equally to all K classes
        mm.r  = rand(size(data,1), K);
        mm.r  = bsxfun(@rdivide, mm.r, sum(mm.r,2));
        mm.Nk = sum(mm.r,1)';
    
    %% Assignment based on the k-means++ algorithm
    case 'kmeans++'
        
        if(K > 1)
            nRep = 25; BestDist = []; target = inf;
            for t = 1:nRep
                [~, Dist] = kmeansinit(data, K);
                S = sum(Dist(:));
                if(S < target)
                    BestDist = Dist; target = S;
                end
            end
            
            R = bsxfun(@rdivide, BestDist, sum(BestDist,2));
            % Heuristic: stop points being assigned 100% to a single class
            R = bsxfun(@rdivide, 0.1+R, sum(0.1+R,2)); 
            ix = sum(isnan(BestDist),2) > 0;
            if(any(ix))                
                t       = rand(sum(ix), K);
                R(ix,:) = bsxfun(@rdivide, t, sum(t,2));
            end            
        else
            R = ones(mm.N, 1);  % all data belongs to one class fully
        end
                
        mm.r  = R;
        mm.Nk = sum(mm.r, 1)';
        mm.a  = (mm.Nk+1/2) ./ (mm.N+K/2); % MML87 estimate vs ML: mm.Nk / sum(mm.Nk);
        
        for k = 1:K
            mm.class{k} = mm_CreateClass(ModelTypes);
        end            
end

%% negative log-likeliood and message length of mixture model
mm.L      = inf;
mm.msglen = inf;        

%% Run estimation functions for all the classes to seed models with parameters (using random assignments)
if(min(mm.Nk) >= mm.MinMembers) % do this only if we have enough items in each class
    mm = mm_EstimateTheta(mm, data, 1:K);
end

end