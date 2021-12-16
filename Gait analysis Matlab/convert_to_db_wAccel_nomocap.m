clear all
close all
filepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\12-2-21\';

% ################# change ########################
intervention = 'thickbrace_fast'; % regular metal doufurt doufurt_fast thickbrace thickbrace_fast
MCvalueSet = [3 4 5 6 7 8]; % order of Mocap left/right based on keyset
FSRvalueSet = [1 2]; % order of FSR inputs based on keyset
FSRkeySet = {'Lheel','Rheel'};
% #################################################

% load 3 datasets
load([filepath, intervention])
load([filepath, intervention, '_fsr'])

% constants
num_sensors = 4;
Fs_pcb = 12800;
Mfsr = containers.Map(FSRkeySet, FSRvalueSet);

% for i = 239144:266129 % for error in fsr trigger
%     Data(9,i) = -3;
% end

%% data clipping to mocap length 8/30/21
[accData, accTime, fsrData, fsrTime] = clip_accelfsr_fromMocap(Data, Time, Fs);
[pcbData, pcbTime] = clip_pcb_fromMocap(datas, times);

% check this is correct. These numbers should be nearly the same:
fsrTime(end)
pcbTime(end)


%% clip end when walking off the floor at the end 9/10/21
% don't need this for human subject testing 10/18/21
figure; plot(pcbTime, pcbData(:,2)) % use this to see when end is
%%
clip_end = 94.5; % delete data after this timestamp for pcb,fsr ##### change ####

fsr_end = findTindex(clip_end,fsrTime);
fsrTime = fsrTime(1:fsr_end);
fsrData = fsrData(1:fsr_end,:);

acc_end = findTindex(clip_end,accTime);
accTime = accTime(1:acc_end);
accData = accData(1:acc_end,:);

pcb_end = findTindex(clip_end,pcbTime);
pcbTime = pcbTime(1:pcb_end);
pcbData = pcbData(1:pcb_end,:);


%%
clip_start = 1.4; % delete data before this timestamp for pcb,fsr ##### change ####

fsr_start = findTindex(clip_start,fsrTime);
fsrTime = fsrTime(fsr_start:end);
fsrData = fsrData(fsr_start:end,:);

acc_start = findTindex(clip_start,accTime);
accTime = accTime(acc_start:end);
accData = accData(acc_start:end,:);

pcb_start = findTindex(clip_start,pcbTime);
pcbTime = pcbTime(pcb_start:end);
pcbData = pcbData(pcb_start:end,:);
%% extracting data 9/9/21
% overall plot for visual check
% plot_3data(pcbData,pcbTime,fsrData,fsrTime,mocapR,mocapL,mocapT,Mfsr)

% clean and filter pcb data
filt_pcbD = lpf_data(pcbData);

% finding footfalls based on fsr heel data
% [impacts, Rheel,- Lheel, Rtoe, Ltoe] = findimpacts_fsr(fsrTime,fsrData,Mfsr);
L_dist = 180; % min distance between heel peaks
R_dist = 180;
min_threshL = 5; % value between peaks threshold, if not lower than this then omit peak
min_threshR = 5;

impacts = findimpacts_fsr_accel(fsrTime,fsrData,Mfsr,L_dist,R_dist,min_threshL,min_threshR);

%% fix small errors in impacts 11/5/21
heel_start_wrong = [2561,13969,4670,6765,8992,18025,26150]; % these need to be same length
heel_start_right = [2786,13959,4604,7031,9188,18010,26348];

heel_pk_wrong = []; % index, these need to be same length
heel_pk_right = [];

impacts = manual_fix_fsr(impacts,fsrData,Mfsr,heel_start_wrong,heel_start_right,heel_pk_wrong,heel_pk_right);

%% extract turning points 11/19/21

pkdist = 3000;
pkprom = 0.001;
figure;
findpeaks(filt_pcbD(:,1),'MinPeakDistance',pkdist,'MinPeakProminence',pkprom)
[pks, locs] = findpeaks(filt_pcbD(:,1),'MinPeakDistance',pkdist,'MinPeakProminence',pkprom);

%% extract turning points pt2 11/19/21

add_pks = [];
omit_pks = [];

start_walk = [pcbTime(locs(1))];
end_walk = [];

turn_threshold = 1; % peaks are within this # of datapoints 
for i = 2:length(locs)
    currT = pcbTime(locs(i));
    prevT = pcbTime(locs(i-1));
    if currT - prevT > turn_threshold % then person got off floor
        start_walk(end+1) = currT;
        end_walk(end+1) = prevT;
    end
end
end_walk(end+1) = pcbTime(locs(end));

% give them buffer bc peak loc doesn't include some parts of the impulse
buffer = 0.2;
start_walk = start_walk - buffer;
end_walk = end_walk + buffer;

figure;
plot(pcbTime, filt_pcbD(:,1))
hold on
plot(start_walk, 0, 'r.','MarkerSize',12)
plot(end_walk, 0, 'b.', 'MarkerSize',12)

%% extract turning points manual fix

start_wrong = []; % these need to be same length
start_right = [];

end_wrong = []; % timestamp, these need to be same length
end_right = [];

if ~isempty(start_wrong)
    for i = 1:length(start_wrong)
        indeces = abs(start_walk - start_wrong(i));
        [~,change_idx] = min(indeces);
        start_walk(change_idx) = start_right(i);
    end
end

if ~isempty(end_wrong)
    for i = 1:length(end_wrong)
        indeces = abs(end_walk - end_wrong(i));
        [~,change_idx] = min(indeces);
        end_walk(change_idx) = end_right(i);
    end
end

% visually check
figure;
plot(pcbTime, filt_pcbD(:,1))
hold on
plot(start_walk, 0, 'r.','MarkerSize',12)
plot(end_walk, 0, 'b.', 'MarkerSize',12)

%% extract turning points pt 3 11/19/21
% delete parts of impacts var that are outside of the walking episode
% ranges

fsrI = findTindex(start_walk(1),fsrTime);
impacts(find(impacts(:,1)<fsrI),:) = [];
for i = 2:length(start_walk)
    start_fsrI = findTindex(start_walk(i),fsrTime);
    end_fsrI = findTindex(end_walk(i-1),fsrTime);
    delete_idx = find(impacts(:,1)<start_fsrI & impacts(:,1)>end_fsrI);
    impacts(delete_idx,:) = [];
end
fsrI = findTindex(end_walk(end),fsrTime);
impacts(find(impacts(:,1)>fsrI),:) = [];

% visually check
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

window_width = Fs*0.35; % shaking lasts < N seconds
offset = 0.1; % start window N behind fsr start
[arrival_idx, peak_idx, peak_mag] = findimpacts_pcb(window_width,offset,impacts,fsrTime,pcbTime,filt_pcbD,num_sensors,true);

%% fix small errors in pcb 11/6/21
arrival_wrong = [NaN,NaN,NaN,8.10383;NaN,NaN,NaN,86.8053]; % these are matrixes, 4 cols for sensorN and each row is a different wrong impact
arrival_right = [NaN,NaN,NaN,7.98969;NaN,NaN,NaN,87.3008];
peak_idx_right = [NaN,NaN,NaN,7.99328;NaN,NaN,NaN,87.3051];

[arrival_idx,peak_idx,peak_mag] = manual_fix_pcb(arrival_wrong,arrival_right,peak_idx_right,pcbTime,arrival_idx,peak_idx,peak_mag,filt_pcbD);

%% extract accel data 11/22/21
% gets the peak of abs value of accelerometer values for x,y,z directions
% 4th col in variable 'acc_pks' indicates whether it is LorR leg.
window = 0.3; % seconds
offset = 0.25;

[acc_pks,acc_pk_idx] = find_accel_impacts(impacts,fsrTime,window,offset,accTime,accData,fsrData,Mfsr);

%% saving data to matlab and excel 9/9/21
processedfilepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\12-2-21\ProcessedData\';
filename = [processedfilepath, intervention];
% for findimpacts_fsr:
% save(filename,'pcbTime','filt_pcbD','arrival_idx','peak_idx','peak_mag', ...
%     'fsrTime','fsrData','impacts','Rheel','Lheel','Rtoe','Ltoe', ...
%     'mocapT','mocapR','mocapL','extracted_pts_R','extracted_pts_L','coordinates','whichfoot')

% for findimpacts_fsr_compact:
whichfoot = impacts(:,4);
save(filename,'pcbTime','filt_pcbD','arrival_idx','peak_idx','peak_mag', ...
    'fsrTime','fsrData','impacts', ...
    'acc_pks','acc_pk_idx','whichfoot') % walk_segments
disp(append("Saved as ", filename))

% EXCEL
% excelfilename = [filename, '_excel.xlsx'];
% T = table(arrival_idx,peak_idx,peak_mag,impacts,coordinates);
% writetable(T,excelfilename)

