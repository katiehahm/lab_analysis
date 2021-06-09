% load data 5/26/21
clear all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\';

% datestr = '05-09-2021_15-08-29'; % Rn 110
datestr = '06-04-2021_13-47-54'; % Rn 116
load([data_root_katie, 'DelsysFSRs\Run_number_116_', datestr])
load([data_root_katie, 'ImmersionFloorSensors\data_', datestr])
T = readtable([data_root_katie, 'MotionCapture\', datestr]);

num_sensors = 4;
Fs = 12800;

% map of fsr sensor location to array number
keySet = {'Rballi','Rtoe','Rheel','Rballo','Ltoe','Lheel','Lballi','Lballo'};
% valueSet = [1 2 3 4 5 6 7 8]; % 05-09-2021_15-08-29
valueSet = [2 3 1 2 7 5 6 6]; % 06-04-2021_13-47-54
Mfsr = containers.Map(keySet, valueSet);

% map of mocap sensor marker to location
% b = ball, t = toe, h = heel
% y is up/down
keySet = {'Lbx','Lby','Lbz','Lhx','Lhy','Lhz','Ltx','Lty','Ltz','Rbx','Rby','Rbz','Rhx','Rhy','Rhz','Rtx','Rty','Rtz'};
% 06-04-2021_13-47-54
% valueSet = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
% 05-09-2021_15-08-29
valueSet = [6 7 8 3 4 5 9 10 11 15 16 17 12 13 14 18 19 20];
Mmocap = containers.Map(keySet, valueSet);

% clip accel data from square triggers in fsr
Time = Time(1,:);
[clipped_time, clipped_data] = clip_to_fsr(times, datas);

% convert and clip stomp mocap csv
[mocapT, mocapL, mocapR] = convertMocap(T, Mmocap);

% clip stomp fsr and pcb 
[clip_stomp_data,clip_stomp_time,clip_fsrD,clip_fsrT] = clip_stomp(clipped_time,clipped_data,Time,Data);


% clip end of data
[pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL] = clip_end(clip_stomp_data,clip_stomp_time,clip_fsrT,clip_fsrD,mocapT,mocapR,mocapL);

% visualize how much to clip
figure; plot(pcbT,pcbD(:,2))

%% clip 1 sec from start of data to filter out 
timeseg = 8.2;
[pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL] = clip_tosettle(timeseg,pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL);


% overall plot for visual check
plot_3data(pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL,Mfsr)

% clean and filter pcb data
filt_pcbD = lpf_data(pcbD);
filt_pcbD = noise_to_zero(filt_pcbD,Fs);

% Final workable data arrays:
% fsrD, fsrT, filt_pcbD, pcbT, mocapT, mocapR, mocapL

% finding footfalls based on fsr heel data 6/7/21
[impactT, impacts] = findimpacts_fsr(fsrT,fsrD,Mfsr);

% extract pcb segments and mocap locations from fsr heel impacts

% pcb extract
[arrival_idx, peak_idx, peak_mag, impactT, impacts] = findimpacts_pcb(impactT,impacts,pcbT,filt_pcbD,Fs,num_sensors,true);

% mocap extract
[coordinates, whichfoot] = findimpacts_mocap(impactT,mocapT,mocapR,mocapL,true);

% final imapct results:
% coordinates, arrival_idx, peak_idx, peak_mag
DataT = [arrival_idx, peak_idx, peak_mag, coordinates];
% eliminate = find(isnan(DataT(:,13))); % eliminate any row where coordinates = NaN
eliminate = find(DataT(:,13)==0); % eliminate any row where coordinates = 0
DataT(eliminate,:) = [];
arrival_idx(eliminate,:) = [];

% final visualization
figure;
plot(pcbT,filt_pcbD)
hold on
arrival_idx_nonan = arrival_idx(~isnan(arrival_idx(:,1)),1);
plot(pcbT(arrival_idx_nonan),0,'k.','MarkerSize',15)
title('PCB final extraction arrival times')

%% 6/9/21 for saving data
filename = [data_root_katie, 'ProcessedData\', datestr];
save(filename, 'DataT', 'pcbT', 'filt_pcbD', 'fsrT', 'fsrD', 'mocapT', 'mocapR', 'mocapL')
disp(append("Saved as ", filename))












%% probably trash everything under this: was used to make a small plot for mini group
%% walking segments extract
% 
% wWidth = 1000; % number of datapoints in a typical impact window width
% [~,seg1] = max(accel(:,2));
% [~,seg2] = max(accel(:,3));
% [~,seg3] = max(accel(seg1+wWidth:end,2));
% seg3 = seg3 + seg1 + wWidth;
% 
% window1 = accel(seg1-wWidth:seg1+wWidth,2);
% window2 = accel(seg2-wWidth:seg2+wWidth,3);
% window3 = accel(seg3-wWidth:seg3+wWidth,2);
% 
% arrival1 = aic_pick(window1, 'to_peak') + seg1 - wWidth;
% arrival2 = aic_pick(window2, 'to_peak') + seg2 - wWidth;
% arrival3 = aic_pick(window3, 'to_peak') + seg3 - wWidth;
% 
% start1 = arrival1+wWidth;
% end1 = arrival2-wWidth;
% start2 = arrival2+wWidth;
% end2 = arrival3-wWidth;
% 
% % final accelerometer segments
% accelW1 = accel(start1:end1,:);
% accelT1 = accelT(start1:end1);
% accelW2 = accel(start2:end2,:);
% accelT2 = accelT(start2:end2);
% 
% startT1 = accelT1(1);
% endT1 = accelT1(end);
% startT2 = accelT2(1);
% endT2 = accelT2(end);
% 
% fsr_startI1 = (find(fsrT>startT1,1));
% fsr_endI1 = (find(fsrT>endT1,1));
% fsr_startI2 = (find(fsrT>startT2,1));
% fsr_endI2 = (find(fsrT>endT2,1));
% 
% % final fsr segments
% fsrT1 = fsrT(fsr_startI1:fsr_endI1);
% fsrT2 = fsrT(fsr_startI2:fsr_endI2);
% fsrW1 = fsr(:,fsr_startI1:fsr_endI1);
% fsrW2 = fsr(:,fsr_startI2:fsr_endI2);
% 
% %% process pcb data
% 
% impactN = 1;
% impactWidth = 10000; % number of indeces an impact lasts in data
% mass = m; % mass of object if dropped
% 
% arrival_idx = zeros(n,sn);
% peak_idx = zeros(n,sn);
% peak_mag = zeros(n,sn);
% % prominence = zeros(n,sn);
% % echo = zeros(n,sn);
% % area_under_curve = zeros(n,sn);
% 
% datas = filt_data;
% 
% for j = 1:sn
%     for i = 1:impactN
%         [~, maxidx] = max(datas(:,j));
%         peak_idx(i,j) = maxidx;
%         min_idx = maxidx - impactWidth*0.2;
%         max_idx = min_idx + impactWidth;
%         datas(min_idx:max_idx,j) = 0;
%     end
% end
% 
% %% making an example plot 5/12/21
% 
% % figure
% % subplot(8,1,1)
% % plot(fsrT1,fsrW1(1,:))
% % subplot(8,1,2)
% % plot(fsrT1,fsrW1(4,:))
% % % subplot(8,1,3)
% % % plot(accelT1,accelW1(:,1))
% % % subplot(8,1,4)
% % % plot(accelT1,accelW1(:,2))
% % subplot(8,1,5)
% % plot(accelT1,accelW1(:,3))
% % subplot(8,1,6)
% % plot(accelT1,accelW1(:,4))
% % subplot(8,1,7)
% % plot(fsrT1,fsrW1(5,:))
% % subplot(8,1,8)
% % plot(fsrT1,fsrW1(8,:))
% 
% % fsr = [Rheel;Rball1;Rball2;Rtoe;Lheel;Lball1;Lball2;Ltoe];
% 
% figure;
% subplot(6,1,1)
% plot(fsrT1(4234:4857),fsrW1(1,4234:4857))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Heel strike FSR')
% subplot(6,1,2)
% plot(fsrT1(4234:4857),fsrW1(4,4234:4857))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Toe strike FSR')
% subplot(6,1,3)
% plot(fsrT1(4234:4857),fsrW2(1,4234:4857))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Heel strike FSR')
% subplot(6,1,4)
% plot(fsrT1(4234:4857),fsrW2(4,4234:4857))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Toe strike FSR')
% subplot(6,1,5)
% plot(accelT1(104535:119895),accelW1(104535:119895,3))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Sensor 3 accelerometer')
% subplot(6,1,6)
% plot(accelT1(104535:119895),accelW1(104535:119895,4))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Sensor 4 accelerometer')
% 
% %% example plot 5/13/21
% 
% figure;
% subplot(8,1,1)
% plot(fsrT1(4234:4857),fsrW1(1,4234:4857))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Right Heel')
% subplot(8,1,2)
% plot(fsrT1(4234:4857),fsrW1(4,4234:4857))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Right Toe')
% subplot(8,1,3)
% plot(fsrT1(4234:4857),fsrW2(1,4234:4857))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Left Heel')
% subplot(8,1,4)
% plot(fsrT1(4234:4857),fsrW2(4,4234:4857))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Right Toe')
% 
% subplot(8,1,5)
% plot(accelT1(104535:119895),accelW1(104535:119895,1))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Sensor 1 accelerometer')
% subplot(8,1,6)
% plot(accelT1(104535:119895),accelW1(104535:119895,2))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Sensor 2 accelerometer')
% subplot(8,1,7)
% plot(accelT1(104535:119895),accelW1(104535:119895,3))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Sensor 3 accelerometer')
% subplot(8,1,8)
% plot(accelT1(104535:119895),accelW1(104535:119895,4))
% hold on
% xline(24.5237, 'b')
% xline(25.056, 'r')
% xline(25.2392, 'g')
% xlabel('Time (s)')
% title('Sensor 4 accelerometer')
% %%
% peak_idx = sort(peak_idx, 'ascend');
% 
% for j = 1:sn
%     for i = 1:n
%         pk_i = peak_idx(i,j);
%         min_idx = pk_i - impactWidth*0.2;
%         max_idx = min_idx + impactWidth;
%         window = filt_data(min_idx:max_idx,j);
%         arrival_idx(i,j) = aic_pick(window, 'to_peak')+min_idx;
%         peak_mag(i,j) = max(window);
%         % prominence calc
%         % echo calc
%         % area under curve calc
%     end
% end
% 
% figure;
% for i = 1:sn
%     hold on
%     subplot(sn,1,i)
%     title(append('Filtered data at ', loc_names(i)))
%     xlabel('Impact number')
%     ylabel('Volts (V)')
%     hold on
%     plot(filt_data(:,i))
%     plot(arrival_idx(:,i),0,'bx')
%     plot(peak_idx(:,i),peak_mag(:,i),'ro')
% end
% 
% % save([data_root_katie, filename], 'raw_data','loc_names','Fs','filt_data','impactN','arrival_idx','peak_idx','peak_mag','mass', '-append')
%    
% %%
% fsr = Data(1,:);
% Time = Time(1,:);
% 
% for i = 1:length(fsr)
%     if fsr(i) < 0.02
%         fsr(i) = 0;
%     end
% end
% 
% 
% 
% 
