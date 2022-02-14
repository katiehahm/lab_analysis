%% ordering impacts 12/19/21
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '2'; % number of subject
takes = {'regular1', 'brace1','weight1','weight2','regular2'};

featureV = zeros(1,24);
for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie,'subj',subj,'_',intervention];
    load(string(filename))
    
    [~,sortI] = sort(arrival_idx(:,1));
    
    acc_pk_idx = acc_pk_idx(sortI,:);
    acc_pks = acc_pks(sortI,:);
    arrival_idx = arrival_idx(sortI,:);
    coordinates = coordinates(sortI,:);
    energy_envelope = energy_envelope(sortI,:);
    energy_squared = energy_squared(sortI,:);
    % EXTRACTED_PTS ARE NOT SORTED
    impacts = impacts(sortI,:);
    peak_idx = peak_idx(sortI,:);
    peak_mag = peak_mag(sortI,:);
    walk_edges = walk_edges(sortI);
    whichfoot = whichfoot(sortI);
    
    save(filename,'pcbTime','filt_pcbD','arrival_idx','peak_idx','peak_mag', ...
    'fsrTime','fsrData','impacts', ...
    'acc_pks','acc_pk_idx',...
    'mocapT','mocapR','mocapL','extracted_pts_R','extracted_pts_L','coordinates','whichfoot')
%     save(filename,'pcbTime','filt_pcbD','arrival_idx','peak_idx','peak_mag', ...
%     'fsrTime','fsrData','impacts', ...
%     'acc_pks','acc_pk_idx',...
%     'mocapT','mocapR','mocapL')
    disp(append("Saved as ", filename))
end

%% add start/end marker for all takes 12/13/21
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '2'; % number of subject
takes = {'regular1', 'brace1', 'brace2','weight2','regular2'};

Fs = 296.3;

all_params = zeros(length(takes),4); % stores left_right mean, std, right_left mean, std

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie,'subj',subj,'_',intervention];
    load(string(filename))
    
    impactN = length(impacts(:,1));
    walk_edges = zeros(1,impactN);
    walk_edges(1) = -1; 
    
    for i = 2:impactN
        curr_i = fsrTime(impacts(i,1));
        prev_i = fsrTime(impacts(i-1,1));
        if curr_i - prev_i > 2 % segments are separated by at least 1.5 secs
            walk_edges(i) = -1;
            walk_edges(i-1) = 1;
        end
    end
    
    walk_edges(impactN) = 1; % last one is always edge
    
    save(filename,'walk_edges','-append')
    
    % visualize
    figure;
    plot(pcbTime, filt_pcbD(:,1))
    hold on
    starts = find(walk_edges == -1);
    plot(fsrTime(impacts(starts,1)),0,'r.','MarkerSize',12)
    lasts = find(walk_edges == 1);
    plot(fsrTime(impacts(lasts,1)),0,'b.','MarkerSize',12)
end

%% energy extraction 12/13/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '2'; % number of subject
takes = {'regular1', 'brace1', 'weight1','weight2','regular2'};
Fs = 12800;
noise_window = Fs*0.2; % take the 0.2s before the arrival of impact to get noise

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie,'subj',subj,'_',intervention];
    load(string(filename))
    
    energy = zeros(length(arrival_idx),4);
    energy_envelope = zeros(length(arrival_idx),4);
    window_end_idx = zeros(length(arrival_idx),4);
    ups = zeros(length(arrival_idx),4);
    for j = 1:4
        for i = 1:length(arrival_idx)
            if i == length(arrival_idx)
                window_end = length(filt_pcbD);
            else
                window_end = arrival_idx(i+1,j);
            end
            curr_i = arrival_idx(i,j);
            window = filt_pcbD(curr_i:window_end,j);
            noise_thresh = max(abs(filt_pcbD(curr_i - noise_window:curr_i-noise_window*0.1)))+0.0000005;
            if ~isempty(window)
%                 [up,lo] = envelope(abs(window),300,'rms');
                [up,lo] = envelope(window,300,'peak');
                indeces = find(up < noise_thresh);
                while isempty(indeces) % keep raising threshold until find the signal end
                    noise_thresh = noise_thresh + 0.000001;
                    indeces = find(up < noise_thresh);
                end
                window = window(1:indeces(1));
                window_end_idx(i,j) = indeces(1) + curr_i;
                window_energy = sum(abs(window).^2);
                energy_envelope(i,j) = sum(up(1:indeces(1)).^2);
                energy(i,j) = window_energy;
            end
        end

    end
    energy_squared = energy;
    save(filename,'energy_squared','energy_envelope','-append')
    
    findzeros = find(energy_squared == 0);
    if ~isempty(findzeros)
        disp("There is a zero in energy matrix")
        take
        findzeros
    end
end



%% real step time 12/13/21
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '2'; % number of subject
takes = {'regular1', 'brace1', 'brace2', 'weight1','weight2','regular2'};

Fs = 296.3;

all_params = zeros(length(takes),4); % stores left_right mean, std, right_left mean, std

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie,'subj',subj,'_',intervention];
    load(string(filename))
    
    left_right_diff = [];
    right_left_diff = [];

    for i = 2:length(walk_edges)
        if walk_edges(i) ~= -1 % not the start of episode
            if impacts(i,4) == 1 % right foot
                if impacts(i-1,4) == 0 % left foot
                    left_right_diff(end+1) = impacts(i,1)-impacts(i-1,1);
                end
            elseif impacts(i,4) == 0 % left foot
                if impacts(i-1,4) == 1 % right foot
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
    
    % deleting outliers: over 2 std
    left_min = left_right_mean - 2*left_right_std;
    left_max = left_right_mean + 2*left_right_std;
    left_idx = find(left_right_diff > left_max | left_right_diff < left_min);
    left_right_diff(left_idx) = [];
    right_min = right_left_mean - 2*right_left_std;
    right_max = right_left_mean + 2*right_left_std;
    right_idx = find(right_left_diff > right_max | right_left_diff < right_min);
    right_left_diff(right_idx) = [];
    
    all_params(take,:) = [left_right_mean, right_left_mean, left_right_std, right_left_std];
    
    color1 = [0.7 0.7 1];
    color1 = [1 0.5 0.5];
    color2 = [88/255 162/255 115/255];
    
    
    fig = figure;
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
    
    titleroot = ['Subject ',subj,' step time distribution ', intervention];
    title(titleroot)
    xlabel('Step Time')
    ylabel('Occurances')
    legend('Left to right','L2R mean','L2R Normal','Right to left','R2L mean','R2L Normal')
    
%     savefig(fig,[data_root_katie,'Figs\',titleroot,'.fig'])
%     saveas(fig,[data_root_katie,'Figs\',titleroot,'.png'])
end

all_params

%% GMM analysis 12/13/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '2'; % number of subject
takes = {'regular1', 'brace1', 'brace2','weight1', 'weight2','regular2'};
GMmodels = zeros(length(takes),9);

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie,'subj',subj,'_',intervention];
    load(string(filename))
    Fs = 12800;

    differences = [];
    for i = 2:length(walk_edges)
        if walk_edges(i) ~= -1 % not the start of episode
            curr = min(arrival_idx(i,1:4));
            prev = min(arrival_idx(i-1,1:4));
            differences(end+1) = (curr - prev)/Fs;
        end
    end
    mu=mean(differences);
    sig = std(differences);
    
    % extracting outliers (step lengths that are too large)
    % that are most likely bc of extraight_straight_paths
    outliers = find(differences>(mu+1.5*sig) | differences<(mu-1.5*sig));
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
    GM = fitgmdist(transpose(differences),2,'RegularizationValue',0.000001);
    proportion = GM.ComponentProportion;
    mu = GM.mu;
    sig = GM.Sigma;
    GMmodels(take,:) = [mu(1), mu(2), sig(1), sig(2), proportion(1), proportion(2), GM.AIC, GM.BIC, GM.NegativeLogLikelihood];
    
end

GMmodels

% scaling computation

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

% calculation of how accurate GMM scaled guess is for each step time 

real_time = all_params(:,1:2);
gmm_time = scaled_means(:,1:2);

max_real = max(real_time,[],2);
min_real = min(real_time,[],2);
max_gmm = max(gmm_time,[],2);
min_gmm = min(gmm_time,[],2);

performance = sum(abs(max_real - max_gmm));
performance = performance + sum(abs(min_real - min_gmm));
performance = performance/(2*length(real_time))

%% GMM and real step time figures 12/13/21

figure;
y = [max_real,max_gmm,min_real,min_gmm];
b = bar(y);
set(gca,'xticklabel',takes)
legend('Larger real step','Larger predicted step','Shorter real step','Shorter predicted step')
title('Subj 2 predicted and real step times for shorter and longer step')

% for point figure
% figure;
% xticks = [0,10,20,30,40,50];
% plot(xticks, max_real,'b.')
% hold on
% plot(xticks + 1, max_gmm, 'r*')
% plot(xticks + 5, min_real,'y*')
% plot(xticks + 6, min_gmm, 'c*')
% legend('Larger real step','Larger predicted step','Shorter real step','Shorter predicted step')
% 
% xmidticks = [5,15,25,35,45,55];
% allticks = sort([xticks,xmidticks]);
% xlabelWords = {'regular1 longer step','regular1 shorter step','brace1 longer step','brace1 shorter step','brace2 longer step','brace2 shorter step',...
%     'weight1 longer step','weight1 shorter step','weight2 longer step','weight2 shorter step','regular2 longer step','regular2 shorter step'};
% set(gca,'xtick',allticks,'xticklabel',xlabelWords)

%% real GRF distribution visualization 12/13/21

clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '4'; % number of subject
takes = {'regular1', 'brace1', 'brace2', 'weight2','regular2'};

all_params = zeros(length(takes),4);
accel_bar_data = zeros(length(takes),2);

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie,'subj',subj,'_',intervention];
    load(string(filename))
    
    left_acc = [];
    right_acc = [];
    for i = 1:length(walk_edges)
        if impacts(i,4) == 1 % right foot
            right_acc(end+1) = acc_pks(i,2);
        elseif impacts(i,4) == 0 % left foot
            left_acc(end+1) = acc_pks(i,2);
        end
    end
    
    accel_bar_data(take,:) = [mean(left_acc),mean(right_acc)];

end
accel_bar_data
figure;
bar(accel_bar_data)
set(gca,'XTickLabel',takes(1:end));
% title(['Subj ', subj, ' differences in GRF from 1st regular walk to interventions'])
title('Average tibial acceleration of L and R legs')
legend('Left foot','Right foot')



%% localization csv 12/14/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '4'; % number of subject
takes = {'regular1', 'brace1','brace2','weight2','regular2'};

featureV = zeros(1,24);
for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie,'subj',subj,'_',intervention];
    load(string(filename))
    
    Nsegments = find(walk_edges == -1); % find how many segments there are
    train = round(length(Nsegments)*0.2); % # segments in one test set
    segN = 0;
    
    for i = 1:length(impacts)
        xcoord = coordinates(i,1)/1000;
        curr_mag = peak_mag(i,:);
        curr_energy = energy_squared(i,:);
        min_arr = min(arrival_idx(i,:));
        arrivals = arrival_idx(i,:) - min_arr;
        
        if walk_edges(i) == -1 % start of segment
            feature = [take, 0, xcoord, curr_mag, arrivals, curr_energy, zeros(1,9)];
            featureV(end+1,:) = feature;
            segN = segN + 1;
        else % don't count prev step if it's a new walking segment
            if mod(segN,train) == 0
                trainValue = floor(segN/train);
            else
                trainValue = floor(segN/train)+1;
            end
            
            if trainValue > 5
                trainValue = 5;
            end
            
            prev_x = coordinates(i-1,1)/1000;
            prev_mag = peak_mag(i-1,:);
            prev_energy = energy_squared(i-1,:);
            
            % xy coord, arrival idx, peak mag, energy, xy prev coord
            feature = [take, trainValue, xcoord,curr_mag, arrivals, curr_energy, curr_mag ./ prev_mag, curr_energy ./ prev_energy, prev_x];
            featureV(end+1,:) = feature;
        end
    end
end
featureV(1,:) = [];
filename = ['C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\ExcelData\subj',subj,'_localization.csv'];
writematrix(featureV,filename) 
% then run recursive_localization.py

%% using recursive localization results to make grf excel 12/20/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '4'; % number of subject
takes = {'regular1', 'brace1','brace2','weight2','regular2'};

T = readtable([data_root_katie, 'ExcelData\trackingloc_subj',subj,'_results']);
A = table2array(T);
real_loc = A(:,1);
predict_loc = A(:,2);

count = 1;

s1 = [-3.590,-3.343];
s2 = [-3.580,2.61];
s3 = [3.639,2.11];
s4 = [3.650,-3.412];

featureV = zeros(1,14);
for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie,'subj',subj,'_',intervention];
    load(string(filename))
    
    for i = 1:length(impacts)
        if round(coordinates(i,1)/1000,2) == round(real_loc(count),2)
            xcoord = predict_loc(count);
            ycoord = coordinates(i,2);

            dist1 = sqrt( (xcoord-s1(1)).^2 + (ycoord-s1(2)).^2 );
            dist2 = sqrt( (xcoord-s2(1)).^2 + (ycoord-s2(2)).^2 );
            dist3 = sqrt( (xcoord-s3(1)).^2 + (ycoord-s3(2)).^2 );
            dist4 = sqrt( (xcoord-s4(1)).^2 + (ycoord-s4(2)).^2 );
            
            accY = acc_pks(i,2);

            % take #, acc, dist 1-4, pk mag 1-4, energy 1-4
            feature = [take,abs(accY),dist1,dist2,dist3,dist4,peak_mag(i,:),energy_squared(i,:)]; %last two are just to check
            featureV(end+1,:) = feature;
            count = count + 1;
        else
            disp("wrong here")
            take
            i
        end
    end
end
featureV(1,:) = [];
filename = [data_root_katie, 'ExcelData\grf_features_subj',subj,'.csv'];
writematrix(featureV,filename) 

% then run grf_wpredictedloc.py

%% use grf predicted results to show difference with kmeans

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '1';
T = readtable([data_root_katie, 'ExcelData\grf_results_GBR_regularweightonly_subj',subj]);
A = table2array(T);
real_grf = A(:,2);
predict_grf = A(:,3);
take = A(:,1);

count = 1;

% these are indeces
N = length(nonzeros(take));
N = round(N/4);

regular1 = predict_grf(1:N);
regular2 = predict_grf(N+1:N*2);
weight1 = predict_grf(N*2+1:N*3);
weight2 = predict_grf(N*3+1:N*4);

real_reg1 = real_grf(1:N);
real_reg2 = real_grf(N+1:N*2);
real_w1 = real_grf(N*2+1:N*3);
real_w2 = real_grf(N*3+1:N*4);

[~,regular1_C] = kmeans(regular1,2)
[~,regular2_C] = kmeans(regular2,2)
[~,weight1_C] = kmeans(weight1,2)
[~,weight2_C] = kmeans(weight2,2)

xmax = 5;

figure;
plot(regular1,real_reg1,'.')
title('Regular 1 GRF prediction')
xlim([1 xmax])
ylim([1 xmax])
xlabel('Predicted GRF')
ylabel('Real GRF')

figure;
plot(regular2,real_reg2,'.')
title('Regular 2 GRF prediction')
xlim([1 xmax])
ylim([1 xmax])
xlabel('Predicted GRF')
ylabel('Real GRF')

figure;
plot(weight1,real_w1,'.')
title('Weight 1 GRF prediction')
xlim([1 xmax])
ylim([1 xmax])
xlabel('Predicted GRF')
ylabel('Real GRF')

figure;
plot(weight2,real_w2,'.')
title('Weight2 GRF prediction')
xlim([1 xmax])
ylim([1 xmax])
xlabel('Predicted GRF')
ylabel('Real GRF')

