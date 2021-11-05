% 10/21/21
% to get real left vs right step times to compare against
% prediction from hist_steptime

clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\ProcessedData\Subj';
subj = '1'; % number of subject
takes = {'normal1', 'slow', 'insole', 'weight', 'count', 'normal2'};
% takes = {'normal1', 'insole', 'weight', 'count', 'normal2'};

Fs = 519;

all_params = zeros(length(takes),4); % stores left_right mean, std, right_left mean, std

for take = 1:length(takes)
    intervention = char(takes(take));
%     filename = [data_root_katie, subj, '_', intervention, '_extract_straight_paths'];
    filename = [data_root_katie, subj, '_', intervention];

    load(string(filename))
    new_impacts = impacts; % for non straight paths data
    
    if take == 1
        intervention = 'Regular 1';
    end
    if take == length(takes)
        intervention = 'Regular 2';
    end
    
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
    all_params(take,:) = [left_right_mean, left_right_std, right_left_mean, right_left_std];
    
%     color1 = [0.7 0.7 1];
    color1 = [1 0.5 0.5];
    color2 = [88/255 162/255 115/255];
    
    
    figure;
    h1 = histogram(left_right_diff);
    h1.Normalization = 'count';
    h1.BinWidth = 0.05;
    h1.FaceColor = color1;
    alpha(0.7)
    hold on
    line([left_right_mean, left_right_mean], [0, max(h1.Values)], 'Color', 'r', 'LineWidth', 1.5);
    
    x1 = left_right_mean - 3*left_right_std:0.01:left_right_mean + 3*left_right_std;
    y1 = gaussmf(x1,[left_right_std left_right_mean]);
    plot(x1,y1.*max(h1.Values), 'r')
    
    h2 = histogram(right_left_diff);
    h2.Normalization = 'count';
    h2.BinWidth = 0.05;
    h2.FaceColor = color2;
    alpha(0.7)
    line([right_left_mean, right_left_mean], [0, max(h2.Values)], 'Color', 'g', 'LineWidth', 1.5);
    x2 = right_left_mean - 3*right_left_std:0.01:right_left_mean + 3*right_left_std;
    y2 = gaussmf(x2,[right_left_std right_left_mean]);
    plot(x2,y2.*max(h2.Values), 'g')
    
    titleroot = 'Step time distribution ';
    title([titleroot, intervention])
    xlabel('Step Time')
    ylabel('Occurances')
    legend('Left to right','L2R mean','L2R Normal','Right to left','R2L mean','R2L Normal')
    
end

all_params