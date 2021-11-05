%% 10/25/21 plotting to histograms and extracting distribution parameters

clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
subj = '2'; % number of subject
% intervention = 'normal1';
% takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'normal2'};
takes = {'normal1', 'insole', 'weight', 'count', 'normal2'}; % for subj2

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, subj, '_', intervention, '_extract_straight_paths'];
    % filename = [data_root_katie, subj, '_', intervention];

    load(string(filename))
    Fs = 12800;

    % for extract_straight_paths:
    differences = [];
    for i = 1:(length(data)-1)
        if data(i+1,19) ~= -1 % edge doesn't say it's the start of seg
            curr = min(data(i+1,1:4));
            prev = min(data(i,1:4));
            differences(end+1) = (curr - prev)/Fs;
        end
    end

    % for original data
    % differences = zeros(1,length(arrival_idx)-1);
    % for i = 1:(length(arrival_idx)-1)
    %     curr = min(arrival_idx(i+1,:));
    %     prev = min(arrival_idx(i,:));
    %     differences(i) = (curr - prev)/Fs;
    % end

    mu=mean(differences);
    outliers = find(differences>(mu*1.5) );
    differences(outliers) = [];

    bin = round(1+3.22*log(numel(differences)));
    figure
    hf=histfit(differences,bin,'kernel');
    figure
    histfit(differences,bin)
    hold on
    x=get(hf(2),'XData'); 
    y=get(hf(2),'YData');
    plot(x,y,'Color','b','LineWidth', 2)
    titleroot = 'Step time distribution ';
    title([titleroot, intervention])
    xlabel('Step Time')
    ylabel('Occurances')
    mu=mean(differences);
    sigma=std(differences);
    hold on
    line([mu, mu], ylim, 'Color', 'c', 'LineWidth', 1); 
    line([mu + sigma, mu + sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
    line([mu - sigma, mu - sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
    legend('Histogram','Normal','Fitted','1 std')
    k = kurtosis(differences);
    s = skewness(differences);
    bimodality = (s^2 + 1)/k;
    iqrange = iqr(differences);

    disp([intervention, ' ','mu ','sigma ','bimodality ','kurtosis ','skewness ', 'iqr'])
    disp([mu, sigma, bimodality, k, s, iqrange])

    % figure;
    % plot(ones(1,length(differences)), differences, 'kx')
end

%% gaussian mixture model 10/25/21

clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
subj = '1'; % number of subject
% intervention = 'normal1';
takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'normal2'};
% takes = {'normal1', 'insole', 'weight', 'count', 'normal2'}; % subject 2
GMmodels = zeros(length(takes),9);

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, subj, '_', intervention, '_extract_straight_paths'];
    % filename = [data_root_katie, subj, '_', intervention];

    load(string(filename))
    Fs = 12800;

    % for extract_straight_paths:
    differences = [];
    for i = 1:(length(data)-1)
        if data(i+1,19) ~= -1 % edge doesn't say it's the start of seg
            curr = min(data(i+1,1:4));
            prev = min(data(i,1:4));
            differences(end+1) = (curr - prev)/Fs;
        end
    end

    % for original data
    % differences = zeros(1,length(arrival_idx)-1);
    % for i = 1:(length(arrival_idx)-1)
    %     curr = min(arrival_idx(i+1,:));
    %     prev = min(arrival_idx(i,:));
    %     differences(i) = (curr - prev)/Fs;
    % end
    
    % extracting outliers (step lengths that are too large)
    % that are most likely bc of extraight_straight_paths
    mu=mean(differences);
    outliers = find(differences>(mu*1.5) );
    differences(outliers) = [];
    
    
%     disp(intervention)
%     GMModel = fitgmdist(transpose(differences),2, 'RegularizationValue',0.1)
    GM = fitgmdist(transpose(differences),2);
    proportion = GM.ComponentProportion;
    mu = GM.mu;
    sig = GM.Sigma;
    GMmodels(take,:) = [mu(1), mu(2), sig(1), sig(2), proportion(1), proportion(2), GM.AIC, GM.BIC, GM.NegativeLogLikelihood];
    
end

GMmodels