clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\Hammering\';

load([data_root_katie, 'take1'])
load([data_root_katie, 'take1_fsr'])

% storing variables
fsr1_data = Data(1,:);
fsr1_trigger = Data(5,:);
fsr1_time = Time(1,:);
fsr1_trigger_time = Time(5,:);
data1 = datas(:,1:4);
data1_trigger = datas(:,8);
time1 = times;

load([data_root_katie, 'take2'])
load([data_root_katie, 'take2_fsr'])

% storing variables
fsr2_data = Data(1,:);
fsr2_trigger = Data(5,:);
fsr2_time = Time(1,:);
fsr2_trigger_time = Time(5,:);
data2 = datas(:,1:4);
data2_trigger = datas(:,8);
time2 = times;

%% clipping
% figure; plot(fsr1_trigger)
% figure; plot(fsr2_trigger)
% figure; plot(data1_trigger)
% figure; plot(data2_trigger)

% got these values from above 4 plots
fsr1_start = 11080;
fsr1_end = 181363;
fsr2_start = 7752;
fsr2_end = 363511;
pcb1_start = 132033;
pcb1_end = 1112798;
pcb2_start = 193379;
pcb2_end = 2242480;

% clip fsr
t_start = fsr1_trigger_time(fsr1_start);
t_end = fsr1_trigger_time(fsr1_end);
i_start = findTindex(t_start,fsr1_time);
i_end = findTindex(t_end,fsr1_time);
fsr1_data = fsr1_data(i_start:i_end);
fsr1_time = linspace(0,t_end-t_start,length(fsr1_data));

t_start = fsr2_trigger_time(fsr2_start);
t_end = fsr2_trigger_time(fsr2_end);
i_start = findTindex(t_start,fsr2_time);
i_end = findTindex(t_end,fsr2_time);
fsr2_data = fsr2_data(i_start:i_end);
fsr2_time = linspace(0,t_end-t_start,length(fsr2_data));

% clip pcb
data1 = data1(pcb1_start:pcb1_end,:);
time1 = linspace(0,time1(pcb1_end)-time1(pcb1_start),length(data1(:,1)));
data2 = data2(pcb2_start:pcb2_end,:);
time2 = linspace(0,time2(pcb2_end)-time2(pcb2_start),length(data2(:,1)));

% check this is correct. These numbers should be nearly the same:
fsr1_time(end)
time1(end)

fsr2_time(end)
time2(end)

%% overall plot for visual check
figure(1)
plotN = 5;
subplot(plotN,1,1)
hold on
plot(fsr1_time, fsr1_data)
for i = 2:5
    subplot(plotN,1,i)
    hold on
    plot(time1, data1(:,i-1))
end

figure(2)
plotN = 5;
subplot(plotN,1,1)
hold on
plot(fsr2_time, fsr2_data)
for i = 2:5
    subplot(plotN,1,i)
    hold on
    plot(time2, data2(:,i-1))
end

% 1st dataset is weird, so ignore from here ########### 10/15/21

%% extract arrival times
% clean and filter pcb data
% data2_ds = downsample(data2, 10); % downsampled by 10
% time2_ds = downsample(time2, 10);
data_ds = data1;
time_ds = time1;
data_h = hilbert(data_ds);
data_h = imag(data_h);
% data2_h = data2_ds;
data_filt = lpf_data(data_h);
% data2_filt = data2_h;
fsr_data = fsr1_data;
fsr_time = fsr1_time;
% finding impacts based on fsr
figure;
plot(fsr_data)
findpeaks(fsr_data,'MinPeakProminence',1.7,'Annotate','extents','MinPeakDistance',1000)
[pks,locs,~,~] = findpeaks(fsr_data,'MinPeakProminence',1.7,'Annotate','extents','MinPeakDistance',1000);

pk_time = fsr_time(locs);
arrival_idx = zeros(length(pk_time),4);
peak_mg = zeros(length(pk_time),4);
for i = 1:length(pk_time)
    for j = 1:4
        idx = findTindex(pk_time(i),time2_ds);
        window = data_filt(idx-5000:idx+5000,j);
        arrival_idx(i,j) = aic_pick(window, 'to_peak') + idx;
        peak_mg(i,j) = max(abs(window));
%         figure; plot(window)
%         hold on
%         plot(arrival_idx(i,j)-idx,0,'rx')
%         plot(0,max(abs(window)),'gx')
    end
end

%% plot to visualize if onset times are consistent

figure;
plot(pks(2:end))
title('Peak magnitudes of fsr')

figure;
hold on
diff31 = arrival_idx(2:end,3)-arrival_idx(2:end,1);
diff32 = arrival_idx(2:end,3)-arrival_idx(2:end,2);
diff34 = arrival_idx(2:end,3)-arrival_idx(2:end,4);
plot(diff31,'x-')
plot(diff32,'x-')
plot(diff34,'x-')
title('pcb arrival differences')
legend('3-1','3-2','3-4')

figure;
hold on
plot(peak_mg(2:end,3)/peak_mg(2:end,1),'x-')
plot(peak_mg(2:end,3)/peak_mg(2:end,2),'x-')
plot(peak_mg(2:end,3)/peak_mg(2:end,4),'x-')
title('pcb peak magnitude ratios')
legend('3/1','3/2','3/4')

figure;
hold on
plot(diff31./diff32)
plot(diff31./diff34)
plot(diff32./diff34)
title('pcb arrival differences ratios')
legend('3-1/3-2','3-1/3-4','3-2/3-4')

figure;
hold on
plot(diff31-diff32)
plot(diff31-diff34)
plot(diff32-diff34)
title('pcb arrival differences differences')
legend('(3-1)-(3-2)','(3-1)-(3-4)','(3-2)-(3-4)')

% figure;
% % this needs to be proper start of impact not the location of the 
% % peak! ########################## 10/15/21
% hold on
% plot(fsr2_time(locs2)-time2_ds(arrival_idx(:,1)))
% plot(fsr2_time(locs2)-time2_ds(arrival_idx(:,2)))
% plot(fsr2_time(locs2)-time2_ds(arrival_idx(:,3)))
% plot(fsr2_time(locs2)-time2_ds(arrival_idx(:,4)))
% title('fsr to pcb differences')
% legend('1','2','3','4')
