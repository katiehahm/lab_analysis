%% 10/19/21 to analyze step time of experiment
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
subj = '2'; % number of subject
intervention = 'count'; % normal1 slow insole weight count normal2
filename = [data_root_katie, subj, '_', intervention];

load(filename)
Fs = 12800;
differences = zeros(1,length(arrival_idx)-1);
for i = 1:(length(arrival_idx)-1)
    curr = min(arrival_idx(i+1,:));
    prev = min(arrival_idx(i,:));
    differences(i) = (curr - prev)/Fs;
end

bin = round(1+3.22*log(numel(differences)));
% figure
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

mu
sigma
k = kurtosis(differences);
s = skewness(differences);
bimodality = (s^2 + 1)/k

%% 10/20/21 for running data files that have straight paths extracted

clear all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
subj = '1'; % number of subject
intervention = 'weight'; % normal1 slow insole weight count normal2
filename = [data_root_katie, subj, '_', intervention, '_extract_straight_paths'];

load(filename)
Fs = 12800;

differences = [];
for i = 1:(length(data)-1)
    if data(i+1,19) ~= -1 % edge doesn't say it's the start of seg
        curr = min(data(i+1,1:4));
        prev = min(data(i,1:4));
        differences(end+1) = (curr - prev)/Fs;
    end
end

bin = round(1+3.22*log(numel(differences)));
% figure
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

mu
sigma
k = kurtosis(differences);
s = skewness(differences);
bimodality = (s^2 + 1)/k
