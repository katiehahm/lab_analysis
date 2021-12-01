%% 10/25/21 plotting to histograms and extracting distribution parameters

clear all
close all
% data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
% subj = '1'; % number of subject
% intervention = 'normal1';
% takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'normal2'};
takes = {'regular1', 'regular2', 'slow', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};
% takes = {'normal1', 'insole', 'weight', 'count', 'normal2'}; % for subj2

for take = 1:length(takes)
    intervention = char(takes(take));
%     filename = [data_root_katie, subj, '_', intervention, '_extract_straight_paths'];
%     filename = [data_root_katie, subj, '_', intervention];
    filename = [data_root_katie, intervention];

    load(string(filename))
    Fs = 12800;

%     for extract_straight_paths:
    differences = [];
    for i = 2:length(whichfoot)
        if walk_segments(i) ~= -1 % not the start of episode
            curr = min(arrival_idx(i,1:4));
            prev = min(arrival_idx(i-1,1:4));
            differences(end+1) = (curr - prev)/Fs;
        end
    end

%     % for original data
%     differences = zeros(1,length(arrival_idx)-1);
%     for i = 1:(length(arrival_idx)-1)
%         curr = min(arrival_idx(i+1,:));
%         prev = min(arrival_idx(i,:));
%         differences(i) = (curr - prev)/Fs;
%     end

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

%% 11/16/21 single intervention plot

clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
subj = '1'; % number of subject
intervention = 'weight';

filename = [data_root_katie, subj, '_', intervention, '_extract_straight_paths'];
load(string(filename))
Fs = 12800;

differences = [];
for i = 2:length(whichfoot)
    if walk_episodes(i) ~= -1 % not the start of episode
        curr = min(arrival_idx(i,1:4));
        prev = min(arrival_idx(i-1,1:4));
        differences(end+1) = (curr - prev)/Fs;
    end
end

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

figure;
plot(ones(1,length(differences)), differences, 'kx')

figure;
plot(differences,'b.')

%% gaussian mixture model 10/25/21

clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\'; 
% subj = '1'; % number of subject
% intervention = 'normal1';
% takes = {'normal1', 'slow', 'insole', 'weight', 'count','normal2'};
takes = {'regular1', 'regular2', 'slow', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};
GMmodels = zeros(length(takes),9);

for take = 1:length(takes)
    intervention = char(takes(take));
%     filename = [data_root_katie, subj, '_', intervention, '_extract_straight_paths'];
%     filename = [data_root_katie, subj, '_', intervention];
    filename = [data_root_katie,intervention];

    load(string(filename))
    Fs = 12800;

    % for extract_straight_paths:
    differences = [];
    for i = 2:length(whichfoot)
        if walk_segments(i) ~= -1 % not the start of episode
            curr = min(arrival_idx(i,1:4));
            prev = min(arrival_idx(i-1,1:4));
            differences(end+1) = (curr - prev)/Fs;
        end
    end

    % for original data
%     differences = zeros(1,length(arrival_idx)-1);
%     for i = 1:(length(arrival_idx)-1)
%         curr = min(arrival_idx(i+1,:));
%         prev = min(arrival_idx(i,:));
%         differences(i) = (curr - prev)/Fs;
%         if curr - prev < 0
%             differences(i) = [];
%         else
%             differences(i) = (curr - prev)/Fs;
%         end
%     end
    
    % extracting outliers (step lengths that are too large)
    % that are most likely bc of extraight_straight_paths
    mu=mean(differences);
    sig = std(differences);
    outliers = find(differences>(mu*1.5) );
    differences(outliers) = [];
    differences(find(differences == 0)) = [];
    
%     disp(intervention)
%     GMModel = fitgmdist(transpose(differences),2, 'RegularizationValue',0.1)
    mu_s = [mu; mu+0.01];
    sigma_s = zeros(1,1,2);
    sigma_s(1,1,:) = [sig; sig];
%     sigma_s(:,:,2) = [sig; sig];
    pcomponents = [1/2,1/2];
%     S = struct('mu',mu_s,'Sigma',sigma_s,'ComponentProportion',pcomponents);
%     GM = fitgmdist(transpose(differences),2,'Start',S);
    GM = fitgmdist(transpose(differences),2,'RegularizationValue',0.0001);
    proportion = GM.ComponentProportion;
    mu = GM.mu;
    sig = GM.Sigma;
    GMmodels(take,:) = [mu(1), mu(2), sig(1), sig(2), proportion(1), proportion(2), GM.AIC, GM.BIC, GM.NegativeLogLikelihood];
    
end

GMmodels

%% 11/10/21 scaling computation

[r,c] = size(GMmodels);
scaled_means = zeros(r,3);

for i = 1:r
    [prop1,~] = max(GMmodels(i,5:6));
    if (1-prop1) < abs(0.5-prop1)
        scaled_means(i,1) = GMmodels(i,1) + (GMmodels(i,2)-GMmodels(i,1))*(1- (GMmodels(i,5))^2);
        scaled_means(i,2) = GMmodels(i,2) + (GMmodels(i,1)-GMmodels(i,2))*(1- (GMmodels(i,6))^2);
        scaled_means(i,3) = abs(scaled_means(i,1)-scaled_means(i,2));
    else
        scaled_means(i,1) = GMmodels(i,1) + (GMmodels(i,2)-GMmodels(i,1))*abs(0.5- GMmodels(i,5));
        scaled_means(i,2) = GMmodels(i,2) + (GMmodels(i,1)-GMmodels(i,2))*abs(0.5- GMmodels(i,6));
        scaled_means(i,3) = abs(scaled_means(i,1)-scaled_means(i,2));
    end
end
scaled_means
% comparison_matrix = [real_times(:,1),scaled_means(:,1),scaled_means(:,2),real_times(:,2),real_times(:,3)];

%% calculation of how accurate GMM scaled guess is for each step time 11/23/21

real_time = all_params(:,1:2);
gmm_time = scaled_means(:,1:2);

max_real = max(real_time,[],2);
min_real = min(real_time,[],2);
max_gmm = max(gmm_time,[],2);
min_gmm = min(real_time,[],2);

performance = sum(abs(max_real - max_gmm)) + sum(abs(min_real-min_gmm));
performance = performance/(2*length(real_time))

%% comparison_matrix calculations 11/10/21

small = 0;
large = 0;
count = 0;
for i = 1:length(comparison_matrix)
    if comparison_matrix(i,1) ~= 0
        large = large + abs( max(comparison_matrix(i,2:3)) - max(comparison_matrix(i,4:5)) );
        small = small + abs( min(comparison_matrix(i,2:3)) - min(comparison_matrix(i,4:5)) );
        count = count + 1;
    end
end
large/count
small/count

figure; plot(max(comparison_matrix(:,2),comparison_matrix(:,3)),max(comparison_matrix(:,4),comparison_matrix(:,5)))
figure; plot(max(comparison_matrix(:,2),comparison_matrix(:,3)),max(comparison_matrix(:,4),comparison_matrix(:,5)),'.')
hold on
plot(min(comparison_matrix(:,2),comparison_matrix(:,3)),min(comparison_matrix(:,4),comparison_matrix(:,5)),'r.')
hold on
plot(linspace(0,2.5,100),linspace(0,2.5,100))
title('Performance of L/R step time estimation')
xlabel('Estimated step time')
ylabel('Real step time')
legend(['Larger step','Smaller step','y=x line'])
legend('Larger step','Smaller step','y=x line')
ylim([0.4 2.5])
xlim([0.4 2.5])

%%
idx = kmeans(transpose(differences),2);
group1 = find(idx==1);
group2 = find(idx == 2);
mean(differences(group1))
mean(differences(group2))

%% gaussian mixture model write to csv 11/9/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'box', 'normal2'};
GMmodels = zeros(70,10);
for j = 1:10
    subj = int2str(j);

    for take = 1:length(takes)
        if j < 6 && take == 6
            continue
        elseif j == 2 && take == 2
            continue
        else
            intervention = char(takes(take));
            filename = [data_root_katie, subj, '_', intervention];

            load(string(filename))
            Fs = 12800;

            % for extract_straight_paths:
        %     differences = [];
        %     for i = 1:(length(data)-1)
        %         if data(i+1,19) ~= -1 % edge doesn't say it's the start of seg
        %             curr = min(data(i+1,1:4));
        %             prev = min(data(i,1:4));
        %             differences(end+1) = (curr - prev)/Fs;
        %         end
        %     end

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
            mu =mean(differences);
            sig = std(differences);
            outliers = find(differences>(mu*1.5) );
            differences(outliers) = [];
            differences(find(differences == 0)) = [];

            mu_s = [mu; mu+0.01];
            sigma_s = zeros(1,1,2);
            sigma_s(1,1,:) = [sig; sig];
            pcomponents = [1/2,1/2];
            S = struct('mu',mu_s,'Sigma',sigma_s,'ComponentProportion',pcomponents);
            GM = fitgmdist(transpose(differences),2,'Start',S,'RegularizationValue',0.001);
            proportion = GM.ComponentProportion;
            mu = GM.mu;
            sig = GM.Sigma;
            GMmodels((j-1)*7 + take,:) = [j,mu(1), mu(2), sig(1), sig(2), proportion(1), proportion(2), GM.AIC, GM.BIC, GM.NegativeLogLikelihood];
        end
        
    end
end
GMmodels_reg = GMmodels;