% load data 5/26/21
clear all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\';

% ################# change ########################
% datestr = '05-09-2021_15-08-29'; % Rn 110
% datestr = '06-04-2021_13-47-54'; % Rn 116
datestr = '06-10-2021_14-17-29'; % Rn 131-134
load([data_root_katie, 'DelsysFSRs\Run_number_134_', datestr]) % ### change
load([data_root_katie, 'ImmersionFloorSensors\data_', datestr])
T = readtable([data_root_katie, 'MotionCapture\', datestr]);

num_sensors = 4;
Fs = 12800;

% ############### change #######################
% map of fsr sensor location to array number
keySet = {'Rballi','Rtoe','Rheel','Rballo','Ltoe','Lheel','Lballi','Lballo'};
% valueSet = [1 2 3 4 5 6 7 8]; % 05-09-2021_15-08-29
% valueSet = [2 3 1 2 7 5 6 6]; % 06-04-2021_13-47-54
valueSet = [5 6 7 8 1 2 3 4]; % 06-10-2021_
Mfsr = containers.Map(keySet, valueSet);

% ################# change ####################
% map of mocap sensor marker to location
% b = ball, t = toe, h = heel
% y is up/down
keySet = {'Lbx','Lby','Lbz','Lhx','Lhy','Lhz','Ltx','Lty','Ltz','Rbx','Rby','Rbz','Rhx','Rhy','Rhz','Rtx','Rty','Rtz'};
% 06-04-2021_13-47-54
% valueSet = [3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
% 05-09-2021_15-08-29
%valueSet = [6 7 8 3 4 5 9 10 11 15 16 17 12 13 14 18 19 20];
% 06-10-2021_13-23-05
valueSet = [6 7 8 3 4 5 6 7 8 12 13 14 9 10 11 12 13 14];
Mmocap = containers.Map(keySet, valueSet);

% clip accel data from square triggers in fsr
Time = Time(1,:);
[clipped_time, clipped_data] = clip_to_fsr(times, datas);

% convert and clip stomp mocap csv
[mocapT, mocapL, mocapR] = convertMocap(T, Mmocap);

% clip stomp fsr and pcb 
[clip_stomp_data,clip_stomp_time,clip_fsrD,clip_fsrT] = clip_stomp(clipped_time,clipped_data,Time,Data);


% clip end of data
[pre_pcbD,pre_pcbT,pre_fsrT,pre_fsrD,pre_mocapT,pre_mocapR,pre_mocapL] = clip_end(clip_stomp_data,clip_stomp_time,clip_fsrT,clip_fsrD,mocapT,mocapR,mocapL);

% visualize how much to clip
figure; plot(pre_pcbT,pre_pcbD(:,2))

%% clip 1 sec from start of data to filter out 
% ############# change ###################
timeseg = 14; % find time right before first real impact. -1.
[pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL] = clip_tosettle(timeseg,pre_pcbD,pre_pcbT,pre_fsrT,pre_fsrD,pre_mocapT,pre_mocapR,pre_mocapL);


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
[extracted_pts_R, extracted_pts_L, coordinates, whichfoot] = findimpacts_mocap(impactT,mocapT,mocapR,mocapL,true);

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

