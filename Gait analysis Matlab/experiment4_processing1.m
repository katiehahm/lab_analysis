%% copied from experiment4_allprocessingcompiled.m, for processing1-6

clear all
close all
filepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment4\Praneeth 5\';
% ################# change ########################
intervention = 'both_weight2'; % regular1 limp1 limp2 weight1 weight2 regular2
% #################################################

% load 3 datasets
load([filepath, intervention])
load([filepath, 'fsr_', intervention])
T = readtable([filepath, 'mocap_', intervention]);

% constants
sensorN = 6;
Fs_pcb = 12800;
MCkeySet = {'Lx1','Ly1','Lz1','Rx1','Ry1','Rz1','Lx2','Ly2','Lz2','Rx2','Ry2','Rz2'};
MCvalueSet = [3 4 5 6 7 8 9 10 11 12 13 14]; % order of Mocap left/right based on keyset
Mmocap = containers.Map(MCkeySet, MCvalueSet);

Fs_mocap = 240;
[mocapT, mocapL1, mocapR1, mocapL2, mocapR2] = convertMocap2pp(T, Mmocap);
allmocap = cat(3, mocapL1,mocapR1,mocapL2,mocapR2);

[accData, accTime, fsrData, fsrTime] = clip_accelfsr_fromMocap2pp_ver2(fsrData, accData, trigger_signal, Fs_fsr, Fs_acc, Fs_trigger);
datas(:,10) = []; % just for Praneeth 5 bc trigger is in row 9
[pcbData, pcbTime] = clip_pcb_fromMocap(datas, times, sensorN);

% check this is correct. These numbers should be nearly the same:
mocapT(end)
fsrTime(end)
pcbTime(end)
accTime(end)

%% delete start and end of data when walking on/off
figure; plot(mocapT, mocapL2(:,1),'b.')
figure; plot(pcbTime, pcbData)

%% cont.

start_time = 11; % change ##
startidx_pcb = findTindex(start_time,pcbTime);
startidx_fsr = findTindex(start_time,fsrTime);
startidx_acc = findTindex(start_time,accTime);
startidx_mocap = findTindex(start_time,mocapT);

fsrData = fsrData(startidx_fsr:end,:);
fsrTime = fsrTime - fsrTime(startidx_fsr);
fsrTime = fsrTime(startidx_fsr:end);
pcbTime = pcbTime - pcbTime(startidx_pcb);
pcbTime = pcbTime(startidx_pcb:end,:);
pcbData = pcbData(startidx_pcb:end,:);
mocapT = mocapT - mocapT(startidx_mocap);
mocapT = mocapT(startidx_mocap:end);
allmocap = allmocap(startidx_mocap:end,:,:);
accTime = accTime - accTime(startidx_acc);
accTime = accTime(startidx_acc:end);
accData = accData(startidx_acc:end,:);


last_time = 136 - start_time; % change ##
lastidx_pcb = findTindex(last_time,pcbTime);
lastidx_fsr = findTindex(last_time,fsrTime);
lastidx_acc = findTindex(last_time,accTime);
lastidx_mocap = findTindex(last_time,mocapT);

fsrData = fsrData(1:lastidx_fsr,:);
fsrTime = fsrTime(1:lastidx_fsr);
pcbTime = pcbTime(1:lastidx_pcb);
pcbData = pcbData(1:lastidx_pcb,:);
mocapT = mocapT(1:lastidx_mocap);
allmocap = allmocap(1:lastidx_mocap,:,:);
accTime = accTime(1:lastidx_acc);
accData = accData(1:lastidx_acc,:);

mocapT(end)
fsrTime(end)
pcbTime(end)
accTime(end)

figure; plot(mocapT, allmocap(:,1,3),'b.')
figure; plot(pcbTime, pcbData)

%% clean up foot drops error in sensor jenny second half 5/21/22 (limp2)
% first run next section. Then run "clean up foot drops", then run this
% section, then run next section again.

% old_fsrData = fsrData;
% % set any drops under 5 = prev value
% for i = 1:length(fsrData(:,4))
%     curr_val = fsrData(i,4);
%     if curr_val < 6.5
%         fsrData(i,4) = fsrData(i-1,4);
%     end
% end

%% finding footfalls based on fsr heel data
L_dist = 180; % min distance between heel peaks
R_dist = 180;
min_threshL = 25; % value between peaks threshold, if not lower than this then omit peak
min_threshR = 25;

% min distance between heel peaks
dist = [180, 180, 180, 180]; % L1, R1, L2, R2
min_thresh = [25, 25, 25, 25]; % L1, R1, L2, R2

% nan_idx = find(fsrData(:,4) < 5);
% fsrData(nan_idx,4) = NaN;

impacts = findimpacts_fsr_accel2pp(fsrTime, mocapT, allmocap, fsrData,dist,min_thresh);

%% to clean up foot drops in sensor 5/21/22 

% % set all value after a peak to be flat zero
% wrong_impacts = find(impacts(:,4) == 4); % foot 4
% fsrL = length(fsrData(:,4));
% for i = 1:length(wrong_impacts)
%     curr_peak_idx = impacts(wrong_impacts(i),2);
%     fsrData(min(curr_peak_idx + 5,fsrL):min(curr_peak_idx + 140,fsrL),4) = 8.5;
% end
% figure; plot(fsrData(:,4))
% then run above section again

%% fix small errors in impacts 11/5/21
heel_start_wrong = [20932,18375]; % these need to be same length
heel_start_right = [20992,18505];

heel_pk_wrong = []; % index, these need to be same length
heel_pk_right = [];

delete_pks = []; % index of peaks to delete

impacts = manual_fix_fsr2pp(impacts,fsrData,heel_start_wrong,heel_start_right,heel_pk_wrong,heel_pk_right,delete_pks);
impacts = sortrows(impacts,1);
%% mocap extract
[coordinates] = findimpacts_mocap2pp(impacts,fsrTime,mocapT,allmocap,true);

%% extract accel data
% gets the peak of abs value of accelerometer values for x,y,z directions
% 4th col in variable 'acc_pks' indicates whether it is LorR leg.
window = 0.3; % seconds
offset = 0.25;

[acc_pks,acc_pk_idx] = find_accel_impacts2pp(impacts,fsrTime,window,offset,accTime,accData,fsrData);

%% saving data to matlab

processedfilepath = [filepath,'ProcessedData\',intervention];
save(processedfilepath,'processedfilepath','pcbTime','datas','pcbData', ...
    'fsrTime','fsrData','impacts','Fs_pcb','Fs_fsr','Fs_acc','Fs_trigger', ...
    'acc_pks','acc_pk_idx','accTime','accData',...
    'mocapT','allmocap','coordinates','Fs_mocap','sensorN')
disp(append("Saved as ", processedfilepath))
close all

%% find walk edges and segments
% coordinates(:,1) = coordinates(:,1)./1000; % only for praneeth 5
% eliminate impacts that are off the floor
elim_idx = [];
for i = 1:length(impacts)
    curr_coord = coordinates(i,1); % x value
    if curr_coord > 3.5 | curr_coord < -3.5
        elim_idx(end+1) = i;
    end
end
impacts(elim_idx,:) = [];
coordinates(elim_idx,:) = [];
acc_pk_idx(elim_idx,:) = [];
acc_pks(elim_idx,:) = [];


% label walk edges and segments
impactN = length(impacts);
walk_edges = zeros(1,impactN);
walk_edges(1) = -1;
walk_edges(end) = 1;
ID_labels = zeros(1,impactN);

% ID labels
o_idx = find(impacts(:,4) == 1 | impacts(:,4) == 2);
x_idx = find(impacts(:,4) == 3 | impacts(:,4) == 4);
ID_labels(o_idx) = 1;
ID_labels(x_idx) = 2;

% find walk_edges
for i = 2:impactN
    if (impacts(i,1) - impacts(i-1,1))/Fs_fsr > 1.5 % start of new walking segment
        walk_edges(i-1) = 1;
        walk_edges(i) = -1;
    end
end

% find segments
segments = find(walk_edges == -1);
segments(end+1) = impactN+1;

% find average step time for person o
step_o1_diff = [];
step_o2_diff = [];
real_o_idx = find(impacts(:,4) == 1 | impacts(:,4) == 2);
real_o_times = impacts(real_o_idx,1)./Fs_fsr;
for x = 2:length(real_o_times)
    diff = real_o_times(x)-real_o_times(x-1);
    if diff < 1.5 % not a turning point
        if impacts(real_o_idx(x),4) == 1 % left foot
            step_o1_diff(end+1) = diff;
        else
            step_o2_diff(end+1) = diff;
        end
    end
end

% same for person x
step_x1_diff = [];
step_x2_diff = [];
real_x_idx = find(impacts(:,4) == 3 | impacts(:,4) == 4);
real_x_times = impacts(real_x_idx,1)./Fs_fsr;
for x = 2:length(real_x_times)
    diff = real_x_times(x)-real_x_times(x-1);
    if diff < 1.5 % not a turning point
        if impacts(real_x_idx(x),4) == 3 % left foot
            step_x1_diff(end+1) = diff;
        else
            step_x2_diff(end+1) = diff;
        end
    end
end

step_o1 = mean(step_o1_diff)
step_o2 = mean(step_o2_diff)
step_x1 = mean(step_x1_diff)
step_x2 = mean(step_x2_diff)

save(processedfilepath,'impactN','impacts','coordinates','acc_pk_idx','acc_pks',...
    'ID_labels','walk_edges','segments','step_o1','step_o2','step_x1','step_x2','-append')

% sanity check walk edges

for i = 1:4
    figure; plot(fsrData(:,i))
    hold on
    impact_idx = find(impacts(:,4) == i);
    plot(impacts(impact_idx,1),fsrData(impacts(impact_idx,1),i),'rx')
end
%% wiener filter noise 5/17/22

% make sure these noise segments are pure noise
for s = 1:6
    figure;
    plot(datas(1:(Fs_pcb*15),s))
end

%% wiener filter noise cont. 5/17/22

wien_pcbD = zeros(length(pcbData),sensorN);
for s = 1:6
    noise_clip = datas(Fs_pcb*12:Fs_pcb*15,s);
    data = pcbData(:,s);
    [wiener_pcb,~] = WienerNoiseReduction([noise_clip;data],Fs_pcb,length(noise_clip));
    wiener_pcb = wiener_pcb(length(noise_clip)+1:end);
    wien_pcbD(:,s) = wiener_pcb;
end

% sanity check
figure; plot(pcbData(:,1))
title('Original')
figure; plot(wien_pcbD(:,1))
title('Wiener')

save(processedfilepath,'wien_pcbD','-append')