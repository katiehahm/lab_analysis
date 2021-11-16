% trying snob data folder to try mixture models other than gaussian
% 11/9/21

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj'; 
subj = '8';
takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'normal2'};
Mmodels = zeros(length(takes),9);

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, subj, '_', intervention];
    load(string(filename))
    Fs = 12800;

    % for original data
    differences = zeros(1,length(arrival_idx)-1);
    for i = 1:(length(arrival_idx)-1)
        curr = min(arrival_idx(i+1,:));
        prev = min(arrival_idx(i,:));
        differences(i) = (curr - prev)/Fs;
        if curr - prev < 0
            differences(i) = [];
        else
            differences(i) = (curr - prev)/Fs;
        end
    end
    
    % extracting outliers (step lengths that are too large)
    % that are most likely bc of extraight_straight_paths
    mu=mean(differences);
    sig = std(differences);
    outliers = find(differences>(mu*1.5) );
    differences(outliers) = [];
    differences(find(differences == 0)) = [];
    
    mm = snob(transpose(differences), {'igauss',1},'k',2,'display',false);
    mm_Summary(mm);
    
    proportion = mm.a;
    nll = mm.L;
%     
% %     disp(intervention)
% %     GMModel = fitgmdist(transpose(differences),2, 'RegularizationValue',0.1)
% %     S = struct('mu',[mu; mu+0.1],'Sigma',[sig; sig+0.01],'ComponentProportion',[1/2,1/2]);
% %     GM = fitgmdist(transpose(differences),2,'SharedCovariance',true,'CovarianceType','diagonal','Start',S);
%     GM = fitgmdist(transpose(differences),2);
%     proportion = GM.ComponentProportion;
%     mu = GM.mu;
%     sig = GM.Sigma;
%     GMmodels(take,:) = [mu(1), mu(2), sig(1), sig(2), proportion(1), proportion(2), GM.AIC, GM.BIC, GM.NegativeLogLikelihood];
%     
end
