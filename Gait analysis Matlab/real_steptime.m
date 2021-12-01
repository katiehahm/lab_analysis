% 10/21/21
% to get real left vs right step times to compare against
% prediction from hist_steptime

% clear all
close all
% data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\';
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_21_21\ProcessedData\';
% subj = '1'; % number of subject
% takes = {'normal1', 'slow', 'insole', 'weight', 'count','normal2'};
takes = {'regular1', 'regular2', 'slow', 'stiffRknee', 'stiffLknee', 'insoleRweightL1', 'insoleRweightL2'};

% Fs = 518.5;
Fs = 296.3;

all_params = zeros(length(takes),4); % stores left_right mean, std, right_left mean, std

for take = 1:length(takes)
    intervention = char(takes(take));
%     filename = [data_root_katie, 'Subj', subj, '_', intervention, '_extract_straight_paths'];
%     filename = [data_root_katie, 'Subj',subj, '_', intervention];
    filename = [data_root_katie,intervention];

    load(string(filename))
    
%     if take == 1
%         intervention = 'Regular 1';
%     end
%     if take == length(takes)
%         intervention = 'Regular 2';
%     end
    
    left_right_diff = [];
    right_left_diff = [];

    for i = 2:length(whichfoot)
        % this if statement is for extract_straight_paths
        if walk_segments(i) ~= -1 % not the start of episode
            if whichfoot(i) == 1 % right foot
                if whichfoot(i-1) == 0 % left foot
                    left_right_diff(end+1) = impacts(i,1)-impacts(i-1,1);
                end
            elseif whichfoot(i) == 0 % left foot
                if whichfoot(i-1) == 1 % right foot
                    right_left_diff(end+1) = impacts(i,1)-impacts(i-1,1);
                end
            end
        end
    end
    
    left_right_diff = left_right_diff./Fs;
    right_left_diff = right_left_diff./Fs;
    
    left_right_mean = mean(left_right_diff);
    left_right_std = std(left_right_diff);
    right_left_mean = mean(right_left_diff);
    right_left_std = std(right_left_diff);
%     all_params(take,:) = [left_right_mean, left_right_std, right_left_mean, right_left_std];
    all_params(take,:) = [left_right_mean, right_left_mean, left_right_std, right_left_std];
    
%     color1 = [0.7 0.7 1];
    color1 = [1 0.5 0.5];
    color2 = [88/255 162/255 115/255];
    
    
%     fig = figure;
%     h1 = histogram(left_right_diff);
%     h1.Normalization = 'count';
%     h1.BinWidth = 0.05;
%     h1.FaceColor = color1;
%     alpha(0.7)
%     hold on
%     line([left_right_mean, left_right_mean], [0, max(h1.Values)], 'Color', 'r', 'LineWidth', 1.5);
%     
%     x1 = left_right_mean - 3*left_right_std:0.01:left_right_mean + 3*left_right_std;
%     y1 = gaussmf(x1,[left_right_std left_right_mean]);
%     plot(x1,y1.*max(h1.Values), 'r')
%     
%     h2 = histogram(right_left_diff);
%     h2.Normalization = 'count';
%     h2.BinWidth = 0.05;
%     h2.FaceColor = color2;
%     alpha(0.7)
%     line([right_left_mean, right_left_mean], [0, max(h2.Values)], 'Color', 'g', 'LineWidth', 1.5);
%     x2 = right_left_mean - 3*right_left_std:0.01:right_left_mean + 3*right_left_std;
%     y2 = gaussmf(x2,[right_left_std right_left_mean]);
%     plot(x2,y2.*max(h2.Values), 'g')
%     
%     titleroot = ['Subject ',subj,' step time distribution ', intervention];
%     title(titleroot)
%     xlabel('Step Time')
%     ylabel('Occurances')
%     legend('Left to right','L2R mean','L2R Normal','Right to left','R2L mean','R2L Normal')
%     
%     savefig(fig,[data_root_katie,'Figs\',titleroot,'.fig'])
%     saveas(fig,[data_root_katie,'Figs\',titleroot,'.png'])
end

all_params

%% 11/16/21
% single intervention plot

clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\';
subj = '10'; % number of subject
intervention = 'weight';

Fs = 518.5;

filename = [data_root_katie, 'Subj', subj, '_', intervention, '_extract_straight_paths'];
load(string(filename))

left_right_diff = [];
right_left_diff = [];

for i = 2:length(whichfoot)
    % this if statement is for extract_straight_paths
    if walk_episodes(i) ~= -1 % not the start of episode
        if whichfoot(i) == 1 % right foot
            if whichfoot(i-1) == 0 % left foot
                left_right_diff(end+1) = impacts(i,1)-impacts(i-1,1);
            end
        elseif whichfoot(i) == 0 % left foot
            if whichfoot(i-1) == 1 % right foot
                right_left_diff(end+1) = impacts(i,1)-impacts(i-1,1);
            end
        end
    end
end

left_right_diff = left_right_diff./Fs;
right_left_diff = right_left_diff./Fs;

left_right_mean = mean(left_right_diff);
left_right_std = std(left_right_diff);
right_left_mean = mean(right_left_diff);
right_left_std = std(right_left_diff);

splitN = 5;
progression = zeros(splitN,2);
lrN = round(length(left_right_diff)/splitN);
rlN = round(length(right_left_diff)/splitN);
progression(1,:) = [mean(left_right_diff(1:lrN)),mean(right_left_diff(1:rlN))];

for i = 1:splitN-1 % split into splitN sections
    lastidx = min(length(left_right_diff),lrN*(i+1));
    progression(i+1,1) = mean(left_right_diff(lrN*i+1:lastidx));
    lastidx = min(length(right_left_diff),rlN*(i+1));
    progression(i+1,2) = mean(right_left_diff(rlN*i+1:lastidx));
end

progression


%% 11/10/21 trying gaussian mix on original fsr data

total_diff = [left_right_diff,right_left_diff];
total_diff = total_diff(randperm(numel(total_diff)));
for i = 1:length(total_diff)
    
    mu_s = [left_right_mean; right_left_mean];
    sigma_s = zeros(1,1,2);
    sigma_s(1,1,:) = [left_right_std; right_left_std];
%     sigma_s(:,:,2) = [sig; sig];
    pcomponents = [1/2,1/2];
    S = struct('mu',mu_s,'Sigma',sigma_s,'ComponentProportion',pcomponents);
    GM = fitgmdist(transpose(total_diff),2,'Start',S);
%     GM = fitgmdist(transpose(total_diff),2);
%     GM = fitgmdist(transpose(differences),2);
    proportion = GM.ComponentProportion;
    mu = GM.mu;
    sig = GM.Sigma;
    GMmodels(take,:) = [mu(1), mu(2), sig(1), sig(2), proportion(1), proportion(2), GM.AIC, GM.BIC, GM.NegativeLogLikelihood];
    
end

GMmodels

%% storing real step time 11/10

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\';
takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'box', 'normal2'};
Fs = 518.5;
real_times = zeros(70,5);

for j = 1:10
    subj = int2str(j);

    for take = 1:length(takes)
        if j < 6 && take == 6
            continue
        elseif j == 2 && take == 2
            continue
        else
            intervention = char(takes(take));
            filename = [data_root_katie, 'Subj',subj, '_', intervention];
            load(string(filename))
            new_impacts = impacts; % for non straight paths data

            left_right_diff = [];
            right_left_diff = [];

            for i = 2:length(new_impacts)
                if new_impacts(i,7) == 1 % right foot
                    if new_impacts(i-1,7) == 0 % left foot
                        left_right_diff(end+1) = new_impacts(i,1)-new_impacts(i-1,1);
                    end
                elseif new_impacts(i,7) == 0 % left foot
                    if new_impacts(i-1,7) == 1 % right foot
                        right_left_diff(end+1) = new_impacts(i,1)-new_impacts(i-1,1);
                    end
                end
            end

            left_right_diff = left_right_diff./Fs;
            right_left_diff = right_left_diff./Fs;

            left_right_mean = mean(left_right_diff);
            left_right_std = std(left_right_diff);
            right_left_mean = mean(right_left_diff);
            right_left_std = std(right_left_diff);
            real_times((j-1)*7 + take,:) = [j,left_right_mean,right_left_mean,left_right_std,right_left_std];
        end
    end
end
save('C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\real_steptimes.mat','real_times')