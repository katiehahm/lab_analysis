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
[pcbData, pcbTime] = clip_pcb_fromMocap(datas, times);

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
plot(mocapT(mci:mci2),mocapL(mci:mci2,2),'b')
hold on
plot(mocapT(mci:mci2),mocapR(mci:mci2,2),'r')
legend('Left','Right')
xlabel('Time (s)')
ylabel('Meters (m)')
title('Foot height from mocap data')
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

% close all
% data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\GMM_results.csv';
% T = readtable(data_root_katie);
% A = table2array(T);
% 
% all_params = zeros(length(takes),4);
% accel_bar_data = zeros(length(takes),2);
% 
% for take = 1:length(takes)
%     intervention = char(takes(take));
%     filename = [data_root_katie,'subj',subj,'_',intervention];
%     load(string(filename))
%     
%     left_acc = [];
%     right_acc = [];
%     for i = 1:length(walk_edges)
%         if impacts(i,4) == 1 % right foot
%             right_acc(end+1) = acc_pks(i,2);
%         elseif impacts(i,4) == 0 % left foot
%             left_acc(end+1) = acc_pks(i,2);
%         end
%     end
%     
%     accel_bar_data(take,:) = [mean(left_acc),mean(right_acc)];
% 
% end
% accel_bar_data
% figure;
% bar(accel_bar_data)
% set(gca,'XTickLabel',takes(1:end));
% % title(['Subj ', subj, ' differences in GRF from 1st regular walk to interventions'])
% title('Average tibial acceleration of L and R legs')
% legend('Left foot','Right foot')

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



