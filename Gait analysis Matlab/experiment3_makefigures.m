%% example of raw signal for paper figure 12/21/21

close all
clear all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\subj2_regular1';
data_root_katie2 = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\Subj 2\subj2_regular1';
data_root_katie3 = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\Subj 2\subj2_fsr_regular1';
load(string(data_root_katie2))
load(string(data_root_katie3))
load(string(data_root_katie))

[accData, accTime, fsrData, fsrTime] = clip_accelfsr_fromMocap(Data, Time, Fs);
[pcbData, pcbTime] = clip_pcb_fromMocap(datas, times,4);

t = 5; % 30 secs in
t2 = 9;
fsri = findTindex(t,fsrTime);
pcbi = findTindex(t,pcbTime);
mci = findTindex(t,mocapT);
acci = findTindex(t,accTime);
fsri2 = findTindex(t2,fsrTime);
pcbi2 = findTindex(t2,pcbTime);
mci2 = findTindex(t2,mocapT);
acci2 = findTindex(t2,accTime);

figure;
subplot(4,1,1)
plot(pcbTime(pcbi:pcbi2),pcbData(pcbi:pcbi2,:))
xlabel('Time (s)')
ylabel('Volts (V)')
title('Accelerometer data')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 13)

subplot(4,1,2)
plot(fsrTime(fsri:fsri2),fsrData(fsri:fsri2,1),'b')
hold on
plot(fsrTime(fsri:fsri2),fsrData(fsri:fsri2,2),'r')
legend('Left','Right')
xlabel('Time (s)')
ylabel('Volts (V)')
title('FSR data')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 13)

subplot(4,1,3)
plot(accTime(acci:acci2),accData(acci:acci2,2),'b')
hold on
plot(accTime(acci:acci2),accData(acci:acci2,5),'r')
legend('Left','Right')
xlabel('Time (s)')
ylabel('Acceleration (g)')
title('Tibial acceleration data')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 13)

subplot(4,1,4)
% plot(mocapT(mci:mci2),mocapL(mci:mci2,1),'b')
% hold on
% plot(mocapT(mci:mci2),mocapR(mci:mci2,1),'r')
% legend('Left','Right')
% xlabel('Time (s)')
% ylabel('Meters (m)')
% title('X coordinate mocap data')
% subplot(5,1,5)
plot(mocapT(mci:mci2),mocapL(mci:mci2,1),'b')
hold on
plot(mocapT(mci:mci2),mocapR(mci:mci2,1),'r')
legend('Left','Right')
xlabel('Time (s)')
ylabel('X-Coordinate (mm)')
title('Foot location from mocap data')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 13)

%% avg step time and cadence 12/22/21

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '4'; % number of subject
takes = {'regular1', 'brace1', 'brace2', 'weight2','regular2'};

% real_time = 0;
% real_count = 0;
% 
% calc_time = 0;
% calc_count = 0;

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie,'subj',subj,'_',intervention];
    load(string(filename))
    
    for i = 1:length(impacts)
        if walk_edges(i) ~= -1 % is not the start of the walk 
            time_diff = fsrTime(impacts(i,1)) - fsrTime(impacts(i-1,1));
            real_time = real_time + time_diff;
            real_count = real_count + 1;
            
            pcb_diff = pcbTime(arrival_idx(i,1)) - pcbTime(arrival_idx(i-1,1));
            calc_time = calc_time + pcb_diff;
            calc_count = calc_count + 1;
        end
    end
end

real_time/real_count
calc_time/calc_count

%% start/stop frequently figure 12/22/21
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\subj1_regular1';
load(string(data_root_katie))
first = 0;
sec = 30;

pcbi = findTindex(first, pcbTime);
pcbi2 = findTindex(sec, pcbTime);

figure; 
plot(pcbTime(pcbi:pcbi2),filt_pcbD(pcbi:pcbi2,:))

%% GMM summary figure 12/22/21

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\GMM_results.xlsx';
T = readtable(data_root_katie);
A = table2array(T);

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
takes = {'regular1', 'brace1', 'brace2', 'weight1','weight2','regular2'};

all_params = zeros(length(takes),4);
accel_bar_data = zeros(length(takes),2);

subj = '1';

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

values = [0.0188 0.0134 0.2246 0.2177;...
    0.0177 0.0337 0.3235 0.3007;...
    0.0145 0.0185 0.1629 0.1922;...
    0.0892 0.0434 0.1602 0.1085;...
    0.0431 0.0239 0.3393 0.3291];

figure;
bar(values.*100)
legend('SI regular','Predicted SI regular','SI brace','Predicted SI brace')
xlabel('Subject number')
ylabel('Symmetry Index (%)')
title('Comparison of measured and predicted SI')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 15)

%% show variation in step time data (for exponent slides) 5/3/22

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\subj1_brace2';
load(string(data_root_katie))


left_step = [];
right_step = [];
for i = 2:length(impacts)
    diff = (impacts(i,1) - impacts(i-1,1))/296.2963;
    if diff < 1
        if impacts(i,4) == 0
            left_step(end+1) = diff;
        else
            right_step(end+1) = diff;
        end
    end
end
        
figure;
plot(left_step)
hold on
plot(right_step,'r')
title('Observed braced walking step time of each foot')
xlabel('Time (s)')
ylabel('Step time (s)')
legend('Left step time','Right step time')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 15)
xlim([0 90])

%% signal dimming 12/23/21

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\subj1_regular1';
load(string(data_root_katie))

pcbi = 12800*3;
pcbi2 = 12800*8;

figure; 
plot(pcbTime(pcbi:pcbi2),filt_pcbD(pcbi:pcbi2,1))
xlim([3 8])
ylim([-0.01 0.01])
xlabel('Time (s)')
ylabel('Volts (V)')
title('Example of signal attenuation due to distance')
set(gca, 'FontName', 'Times New Roman')
set(gca, 'FontSize', 15)

%% find average 

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '1';
T = readtable([data_root_katie, 'ExcelData\grf_results_GBR_regularweightonly_subj',subj]);
A = table2array(T);
real_grf = A(:,2);
take = A(:,1);
N = length(nonzeros(take));
avg_TA1 = sum(nonzeros(real_grf))/N;

subj = '2';
T = readtable([data_root_katie, 'ExcelData\grf_results_GBR_regularweightonly_subj',subj]);
A = table2array(T);
take = A(:,1);
real_grf = A(:,2);
N = N + length(nonzeros(take));
avg_TA2 = mean(nonzeros(real_grf));

subj = '3';
T = readtable([data_root_katie, 'ExcelData\grf_results_GBR_regularweightonly_subj',subj]);
A = table2array(T);
take = A(:,1);
real_grf = A(:,2);
N = N + length(nonzeros(take));
avg_TA3 = mean(nonzeros(real_grf));

subj = '4';
T = readtable([data_root_katie, 'ExcelData\grf_results_GBR_regularweightonly_subj',subj]);
A = table2array(T);
take = A(:,1);
real_grf = A(:,2);
N = N + length(nonzeros(take));
avg_TA4 = mean(nonzeros(real_grf));

(avg_TA1+avg_TA2+avg_TA3+avg_TA4)/4

%% TA average / range plot

right_arr = ones(length(TA_reg1_R),1);
left_arr = zeros(length(TA_reg1_R),1);
data_root = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\subj';
takes = {'regular1', 'brace1','brace2','weight1','weight2','regular2'};
M = [];
for i = 1:6
    intervention = char(takes(i));
    load([data_root,'2_',intervention])
    
    idx = find(acc_pks(:,4) == 1);
    TA_R = acc_pks(idx,2);
    idx = find(acc_pks(:,4) == 0);
    TA_L = acc_pks(idx,2);
    
    right_arr = ones(length(TA_R),1);
    left_arr = zeros(length(TA_L),1);
    
    name_arr_R = right_arr;
    for x = 1:length(name_arr_R)
        name_arr_R(x) = i;
    end
    
    name_arr_L = left_arr;
    for x = 1:length(name_arr_L)
        name_arr_L(x) = i;
    end
    
    M1 = [right_arr,name_arr_R,TA_R];
    M2 = [left_arr,name_arr_L,TA_L];
    
    M = [M;M1;M2];

end

writematrix(M,'TAplot.csv')

%% uncomment later
tbl = readtable('TAplot.csv');
takeOrder = {'1','2','3','4','5','6'};
% tbl.Take = categorical(tbl.Take,takeOrder);
figure;
boxchart(tbl.Take,tbl.TA,'GroupByColor',tbl.Foot,'JitterOutliers','off')
ylabel('Acceleration (g)')
legend('Right Foot','Left Root')
xticks([1,2,3,4,5,6]);
xticklabels({'regular1','brace1','brace2','weight1','weight2','regular2'})
title('Observed TA of each foot for one subject')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14)

%% load subj3_brace1_sorted

diff = [];
for i = 2:length(arrival_idx(:,1))
    difference = min(arrival_idx(i,:)) - min(arrival_idx(i-1,:));
    if difference < 12800 % less than 1 sec bc idk where walk_edges is
        diff(end+1) = difference;
    end
end

diff = diff ./ 12800;

fig = figure;
h1 = histogram(diff);
h1.Normalization = 'count';
h1.BinWidth = 0.05;
% h1.FaceColor = color1;
alpha(0.7)
hold on
line([0.4658, 0.4658], [0, max(h1.Values)], 'Color', 'r', 'LineWidth', 1.5);
line([0.5710, 0.5710], [0, max(h1.Values)], 'Color', 'k', 'LineWidth', 1.5);
% line([0.4638, 0.4638], [0, max(h1.Values)], 'Color', 'k', 'LineWidth', 1.5);
% line([0.6028, 0.6028], [0, max(h1.Values)], 'Color', 'k', 'LineWidth', 1.5);

% legend('Step time differences','Predicted avg step time (leg 1)','Predicted avg step time (leg 2)','Real avg step time (leg 1)','Real avg step time (leg 2)')
legend('Step time differences','Predicted avg step time (leg 1)','Predicted avg step time (leg 2)')
xlabel('Number of occurances')
ylabel('Step time (s)')
title('Step time difference distribution and predicted step time for each leg')

%% localization

data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\ExcelData\trackingloc_subj1_results.csv';
T = readtable(data_root_katie);
A = table2array(T);

predicted = A(:,2);
target = A(:,1);

figure;
plot(predicted, target,'.','MarkerSize',12)
title('Tracking localization performance')
xlabel('Predicted Values [m]')
ylabel('Target Values [m]')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14)

%% GRF and kmeans

close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '2';
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

[reg_idx,regular1_C] = kmeans(regular1,2);
[~,regular2_C] = kmeans(regular2,2);
[weight_idx,weight1_C] = kmeans(weight1,2);
[~,weight2_C] = kmeans(weight2,2);

reg_idx1 = find(reg_idx == 1);
reg_idx2 = find(reg_idx == 2);
weight_idx1 = find(weight_idx == 1);
weight_idx2 = find(weight_idx == 2);

takes = {'regular1','weight1','weight2','regular2'};

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


xmax = 3.5;

% old plots
figure;
plot(regular1,real_reg1,'k.')
hold on
plot(min(regular1_C),min(accel_bar_data(1,:)),'rx','MarkerSize',12)
plot(max(regular1_C),max(accel_bar_data(1,:)),'b*','MarkerSize',12)
title('Regular walking TA prediction performance')
xlim([1 xmax])
ylim([1 xmax])
xlabel('Predicted GRF')
ylabel('Real GRF')
legend('TA','Leg 1 mean','Leg 2 mean')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',16)

figure;
plot(weight1,real_w1,'k.')
hold on
plot(min(weight1_C),min(accel_bar_data(2,:)),'rx','MarkerSize',12)
plot(max(weight1_C),max(accel_bar_data(2,:)),'b*','MarkerSize',12)
title('Walking with weight TA prediction performance')
xlim([1 xmax])
ylim([1 xmax])
xlabel('Predicted GRF')
ylabel('Real GRF')
legend('TA','Leg 1 mean','Leg 2 mean')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',16)

% new plots 2/14/22 (for com mtg)
figure;
plot(regular1(reg_idx1),real_reg1(reg_idx1),'r.')
hold on
plot(regular1(reg_idx2),real_reg1(reg_idx2),'b.')
plot(min(regular1_C),min(accel_bar_data(1,:)),'kx','MarkerSize',20,'LineWidth',3)
plot(max(regular1_C),max(accel_bar_data(1,:)),'k*','MarkerSize',20,'LineWidth',3)
plot(1:0.01:xmax,1:0.01:xmax,'k--')
title('Regular walking TA prediction performance')
xlim([1 xmax])
ylim([1 xmax])
xlabel('Predicted TA')
ylabel('Real TA')
legend('Cluster 1','Cluster 2','Cluster 1 centroid','Cluster 2 centroid')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',16)

figure;
plot(weight1(weight_idx1),real_w1(weight_idx1),'r.')
hold on
plot(weight1(weight_idx2),real_w1(weight_idx2),'b.')
plot(min(weight1_C),min(accel_bar_data(2,:)),'kx','MarkerSize',20,'LineWidth',3)
plot(max(weight1_C),max(accel_bar_data(2,:)),'k*','MarkerSize',20,'LineWidth',3)
plot(1:0.01:xmax,1:0.01:xmax,'k--')
title('Walking with weight TA prediction performance')
xlim([1 xmax])
ylim([1 xmax])
xlabel('Predicted TA')
ylabel('Real TA')
legend('Cluster 1','Cluster 2','Cluster 1 centroid','Cluster 2 centroid')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',16)



%% testing if step time data is normally distributed with qq plot 2/1/22

step_times = [];
for i = 2:length(impacts(:,1))
    if walk_edges(i) ~= -1
        curr_time = impacts(i,1) - impacts(i-1,1);
        if curr_time < 175
            step_times(end+1) = curr_time;
        end
    end
end

qqplot(step_times)

%% GMM graphic 2/13/22


close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '2'; % number of subject
takes = {'regular1', 'brace1'};
GMmodels = zeros(length(takes),9);
differences1 = [];
differences2 = [];

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
    
    if take == 1
        differences1 = differences;
    else
        differences2 = differences;
    end
    
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


%% got differences from prev section
std = [0.0005, 0.0002, 0.0006, 0.0008];
means = [0.5419, 0.5905, 0.5489, 0.6491];
pptions = [0.8763, 0.1237, 0.5441, 0.4559];
scaled_means = [0.5532, 0.5426, 0.5533, 0.6447];

sort_diff1 = sort(differences1);
sort_diff2 = sort(differences2);

% regular modified
data1 = sort_diff1(1:2:end);
data1 = data1(1:round(length(data1)*0.95));
data2 = sort_diff1(2:2:end);
data2 = data2(round(length(data2)*0.05):end);
figure;
h1 = histfit(data1,10);
h1(1).FaceColor = [.8 .8 1];
h1(1).BarWidth = 1;
% h1(1).FaceAlpha = 0.5;
h1(2).Color = 'blue';
hold on
h2 = histfit(data2,11);
h2(1).FaceColor = [1 .8 .8];
h2(1).FaceAlpha = 0.5;
h2(1).BarWidth = 1;
h2(2).Color = 'red';
title('Regular walking step time distribution')
xlabel('Step time (s)')
ylabel('Occurances')
h = [h1(2);h2(2)];
legend(h,'GMM Mixture 1','GMM Mixture 2')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14)

%% regular
pp_start1 = round(length(sort_diff1)*pptions(1));
data1 = sort_diff1(1:pp_start1);
data2 = sort_diff1(pp_start1+1:end);
figure;
h1 = histfit(data1);
h1(1).FaceColor = [.8 .8 1];
h1(2).Color = 'blue';
hold on
h2 = histfit(data2);
h2(1).FaceColor = [1 .8 .8];
h2(2).Color = 'red';
title('Regular walking step time distribution')
xlabel('Step time (s)')
ylabel('Occurances')
h = [h1(2);h2(2)];
legend(h,'GMM Mixture 1','GMM Mixture 2')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14)

% brace
pp_start1 = round(length(sort_diff2)*pptions(3));
data1 = sort_diff2(1:pp_start1);
data2 = sort_diff2(pp_start1+1:end);
figure;
h1 = histfit(data1);
h1(1).FaceColor = [.8 .8 1];
h1(2).Color = 'blue';
hold on
h2 = histfit(data2);
h2(1).FaceColor = [1 .8 .8];
h2(2).Color = 'red';
title('Brace walking step time distribution')
xlabel('Step time (s)')
ylabel('Occurances')
h = [h1(2);h2(2)];
legend(h,'GMM Mixture 1','GMM Mixture 2')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14)

% regular w/ scaling
pp_start1 = round(length(sort_diff1)*pptions(1));
data1 = sort_diff1(1:pp_start1);
data2 = sort_diff1(pp_start1+1:end);
figure;
h1 = histfit(data1);
h1(1).FaceColor = [.8 .8 1];
h1(2).Color = 'blue';
hold on
h2 = histfit(data2);
h2(1).FaceColor = [1 .8 .8];
h2(2).Color = 'red';
title('Regular walking step time distribution')
xlabel('Step time (s)')
ylabel('Occurances')
h = [h1(2);h2(2)];
legend(h,'GMM Mixture 1','GMM Mixture 2')
set(findobj(gcf,'type','axes'),'FontName','Times New Roman','FontSize',14)

