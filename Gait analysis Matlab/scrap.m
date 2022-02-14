N = length(arrival_idx(:,1));
down_arrival_idx = round(arrival_idx(:,1)./10);
for i = 1:10
    % just sensor 1
    clip = filt_pcbD(arrival_idx(i):arrival_idx(i)+2560,1); % 0.5s window
    clip_d = iddata(clip,[],0.000078125); % sample time is 0.5
    na = 20;
    nc = 20;
%     sys = armax(clip_d,[na nc]);
    sys = ar(clip_d, na);
    figure; compare(clip_d, sys)
end

%% for experiments 12/13/21-12/16/21 initial file to run

clear all
close all
filepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\';
subj = '1'; 
% ################# change ########################
intervention = 'regular1'; % regular1 brace1 brace2 weight1 weight2 regular2
MCvalueSet = [3 4 5 6 7 8]; % order of Mocap left/right based on keyset
FSRvalueSet = [1 2]; % order of FSR inputs based on keyset
FSRkeySet = {'Lheel','Rheel'};
% #################################################

filepath = [filepath, 'Subj ',subj,'\subj',subj','_'];
% load 3 datasets
load([filepath, intervention])
load([filepath, 'fsr_', intervention])
T = readtable([filepath, 'mocap_', intervention]);

% constants
num_sensors = 4;
Fs_pcb = 12800;
Mfsr = containers.Map(FSRkeySet, FSRvalueSet);
MCkeySet = {'Lx','Ly','Lz','Rx','Ry','Rz'};
Mmocap = containers.Map(MCkeySet, MCvalueSet);

[mocapT, mocapL, mocapR] = convertMocap(T, Mmocap);

% for i = 239144:266129 % for error in fsr trigger
%     Data(9,i) = -3;
% end

%% data clipping to mocap length 8/30/21
[accData, accTime, fsrData, fsrTime] = clip_accelfsr_fromMocap(Data, Time, Fs);
[pcbData, pcbTime] = clip_pcb_fromMocap(datas, times);

% check this is correct. These numbers should be nearly the same:
mocapT(end)
fsrTime(end)
pcbTime(end)


%% extracting data 9/9/21
% overall plot for visual check
% plot_3data(pcbData,pcbTime,fsrData,fsrTime,mocapR,mocapL,mocapT,Mfsr)

% clean and filter pcb data
downData = downsample(pcbData, 10);
filt_pcbD = lpf_data(pcbData);
filt_down_pcbD = lpf_data(downData);

% finding footfalls based on fsr heel data
% [impacts, Rheel,- Lheel, Rtoe, Ltoe] = findimpacts_fsr(fsrTime,fsrData,Mfsr);
L_dist = 180; % min distance between heel peaks
R_dist = 180;
min_threshL = 25; % value between peaks threshold, if not lower than this then omit peak
min_threshR = 25;

% sometimes the fsr values are 0
for i = 1:length(fsrData)
    if fsrData(i,1) < 10
        fsrData(i,1) = fsrData(i-1,1);
    end
    if fsrData(i,2) < 1
        fsrData(i,2) = fsrData(i-1,2);
    end
end

impacts = findimpacts_fsr_accel(fsrTime,fsrData,Mfsr,L_dist,R_dist,min_threshL,min_threshR);

%% fix small errors in impacts 11/5/21
heel_start_wrong = [27739,37126]; % these need to be same length
heel_start_right = [27688,36986];

heel_pk_wrong = []; % index, these need to be same length
heel_pk_right = [];

impacts = manual_fix_fsr(impacts,fsrData,Mfsr,heel_start_wrong,heel_start_right,heel_pk_wrong,heel_pk_right);

%% to fix subj4_regular1 bc fsr data was bad
% old_impacts = impacts;
% keepimpacts = [];
% for i = 1:length(impacts)
%     if impacts(i,1) < 14085
%         keepimpacts(end+1) = i;
%     elseif impacts(i,1) > 24344
%         keepimpacts(end+1) = i;
%     end
% end
% impacts = old_impacts(keepimpacts,:);
% 
% figure;
% subplot(2,1,1)
% hold on
% plot(fsrData(:,Mfsr('Lheel')))
% lefts = find(impacts(:,4) == 0);
% plot(impacts(lefts,1),fsrData(impacts(lefts,1),Mfsr('Lheel')),'rx','MarkerSize',12)
% plot(impacts(lefts,2),fsrData(impacts(lefts,2),Mfsr('Lheel')),'bx','MarkerSize',12)
% title('Left heel')
% 
% subplot(2,1,2)
% hold on
% plot(fsrData(:,Mfsr('Rheel')))
% rights = find(impacts(:,4) == 1);
% plot(impacts(rights,1),fsrData(impacts(rights,1),Mfsr('Rheel')),'rx','MarkerSize',12)
% plot(impacts(rights,2),fsrData(impacts(rights,2),Mfsr('Rheel')),'bx','MarkerSize',12)
% title('Right heel')

%% extract turning points from looking at x cord 12/13/21

% eliminate any impacts off the floor
x_min = -3500;
x_max = 3600;
elim_i = [1]; % first step is always on the resonance

for i = 1:length(impacts)
    impact_time = fsrTime(impacts(i,1));
    mocap_i = findTindex(impact_time,mocapT);
    if impacts(i,4) == 1 % right foot
        coord = mocapR(mocap_i);
        if coord > x_max || coord < x_min
            elim_i(end+1) = i;
        end
    else % left foot
        coord = mocapL(mocap_i);
        if coord > x_max || coord < x_min
            elim_i(end+1) = i;
        end
    end
end

impacts(elim_i,:) = [];

% visually check
% fsr
RorL = length(impacts(1,:)); % always the last column
figure;
subplot(2,1,1)
hold on
plot(fsrData(:,Mfsr('Lheel')))
lefts = find(impacts(:,RorL) == 0);
plot(impacts(lefts,1),fsrData(impacts(lefts,1),Mfsr('Lheel')),'rx','MarkerSize',12)
plot(impacts(lefts,2),fsrData(impacts(lefts,2),Mfsr('Lheel')),'bx','MarkerSize',12)
title('Left heel')

subplot(2,1,2)
hold on
plot(fsrData(:,Mfsr('Rheel')))
rights = find(impacts(:,RorL) == 1);
plot(impacts(rights,1),fsrData(impacts(rights,1),Mfsr('Rheel')),'rx','MarkerSize',8)
plot(impacts(rights,2),fsrData(impacts(rights,2),Mfsr('Rheel')),'bx','MarkerSize',12)
title('Right heel')

% pcb
figure;
plot(pcbTime,filt_pcbD(:,1))
hold on
plot(fsrTime(impacts(:,1)),0,'r.','MarkerSize',12)

%% manual eliminate turns bc mocap doesn't work

% % eliminate any impacts off the floor
% starts_elim = [3.08578,10.646,17.9223,25.1478,33.1033,40.7816,48.5384,56.3171,63.9938,71.0226,79.1123,86.9313,95.4166,103.019,111.287,119.404,127.321,135.015,142.831,151.43];
% ends_elim = [7.43836,15.1642,22.3268,29.7812,37.9017,45.4147,53.2272,61.1672,68.7294,75.8789,84.0383,91.6504,99.5757,107.9,116.212,124.351,132.146,139.729,147.674,156.196];
% 
% fsr_idx = findTindex(starts_elim(1),fsrTime);
% elim_idx = find(impacts(:,1) < fsr_idx);
% impacts(elim_idx,:) = [];
% 
% for i = 2:length(starts_elim)
%     idx_start = findTindex(starts_elim(i),fsrTime);
%     idx_end = findTindex(ends_elim(i-1),fsrTime);
%     
%     elim_idx = find(impacts(:,1) > idx_end & impacts(:,1) < idx_start);
%     impacts(elim_idx,:) = [];
% end
% 
% fsr_idx = findTindex(ends_elim(end),fsrTime);
% elim_idx = find(impacts(:,1) > fsr_idx);
% impacts(elim_idx,:) = [];
% 
% % visually check
% % fsr
% RorL = length(impacts(1,:)); % always the last column
% figure;
% subplot(2,1,1)
% hold on
% plot(fsrData(:,Mfsr('Lheel')))
% lefts = find(impacts(:,RorL) == 0);
% plot(impacts(lefts,1),fsrData(impacts(lefts,1),Mfsr('Lheel')),'rx','MarkerSize',12)
% plot(impacts(lefts,2),fsrData(impacts(lefts,2),Mfsr('Lheel')),'bx','MarkerSize',12)
% title('Left heel')
% 
% subplot(2,1,2)
% hold on
% plot(fsrData(:,Mfsr('Rheel')))
% rights = find(impacts(:,RorL) == 1);
% plot(impacts(rights,1),fsrData(impacts(rights,1),Mfsr('Rheel')),'rx','MarkerSize',12)
% plot(impacts(rights,2),fsrData(impacts(rights,2),Mfsr('Rheel')),'bx','MarkerSize',12)
% title('Right heel')
% 
% % pcb
% figure;
% plot(pcbTime,filt_pcbD(:,1))
% hold on
% plot(fsrTime(impacts(:,1)),0,'r.','MarkerSize',12)
%% pcb adjust params to capture whole impact 11/5/21
% look at where the impacts are and adjust window_width and offset
% parameters in the next section

close all
Fs = 12800;
secs = 20; % plot 10 secs of data
figure;
for i = 1:4
    subplot(4,1,i)
    plot(pcbTime(1:secs*Fs), filt_pcbD(1:secs*Fs,i))
    hold on
    indeces = find(fsrTime(impacts(:,1)) < secs);
    plot(fsrTime(impacts(indeces,1)),0,'r.','MarkerSize',10)
end

%% pcb and mocap extract

window_width = Fs*0.3; % shaking lasts < N seconds
offset = 0.1; % start window N behind fsr start
[arrival_idx, peak_idx, peak_mag] = findimpacts_pcb(window_width,offset,impacts,fsrTime,pcbTime,filt_pcbD,num_sensors,true);

% mocap extract
[extracted_pts_R, extracted_pts_L, coordinates, whichfoot] = findimpacts_mocap(impacts,fsrTime,mocapT,mocapR,mocapL,true);

%% fix small errors in pcb 11/6/21
arrival_wrong = [13.7388,13.7409,13.7386,13.7367]; % these are matrixes, 4 cols for sensorN and each row is a different wrong impact
arrival_right = [13.6563,13.6555,13.6635,13.6618];
peak_idx_right = [13.6748,13.668,13.6674,13.697];

[arrival_idx,peak_idx,peak_mag] = manual_fix_pcb(arrival_wrong,arrival_right,peak_idx_right,pcbTime,arrival_idx,peak_idx,peak_mag,filt_pcbD);


%% fix small errors in mocap 11/5/21
r_wrong = [24475]; % these need to be same length
r_right = [24713];

l_wrong = []; % index, these need to be same length
l_right = [];

[extracted_pts_R,extracted_pts_L,coordinates] = manual_fix_mocap(r_wrong,r_right,l_wrong,l_right,mocapR,mocapL,extracted_pts_R,extracted_pts_L,coordinates,whichfoot);

%% extract accel data 11/22/21
% gets the peak of abs value of accelerometer values for x,y,z directions
% 4th col in variable 'acc_pks' indicates whether it is LorR leg.
window = 0.3; % seconds
offset = 0.25;

[acc_pks,acc_pk_idx] = find_accel_impacts(impacts,fsrTime,window,offset,accTime,accData,fsrData,Mfsr);

%% saving data to matlab and excel 9/9/21
processedfilepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
filename = [processedfilepath, 'subj',subj,'_',intervention];

save(filename,'pcbTime','filt_pcbD','arrival_idx','peak_idx','peak_mag', ...
    'fsrTime','fsrData','impacts', ...
    'acc_pks','acc_pk_idx','accTime','accData',...
    'mocapT','mocapR','mocapL','extracted_pts_R','extracted_pts_L','coordinates','whichfoot')
disp(append("Saved as ", filename))

