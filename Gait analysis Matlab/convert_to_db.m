clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\';

% ################# change ########################
subj = '10'; % number of subject
intervention = 'insole'; % normal1 slow insole weight count box normal2
% box (for subj 6+)
FSRvalueSet = [1 2 5 6]; % order of FSR inputs based on keyset
MCvalueSet = [3 4 5 6 7 8]; % order of Mocap left/right based on keyset
% #################################################

% load 3 datasets
new_root = [data_root_katie, 'Subject ', subj, '\subj', subj, '_'];
load([new_root, intervention])
load([new_root, 'fsr_', intervention])
T = readtable([new_root, 'mocap_', intervention]);

% constants
num_sensors = 4;
Fs = 12800;
FSRkeySet = {'Lheel','Ltoe','Rheel','Rtoe'};
Mfsr = containers.Map(FSRkeySet, FSRvalueSet);
MCkeySet = {'Lx','Ly','Lz','Rx','Ry','Rz'};
Mmocap = containers.Map(MCkeySet, MCvalueSet);

[mocapT, mocapL, mocapR] = convertMocap(T, Mmocap);

for i = 239144:266129 % for error in fsr trigger
    Data(9,i) = -3;
end



%% data clipping to mocap length 8/30/21
[fsrData, fsrTime] = clip_fsr_fromMocap(Data, Time);
[pcbData, pcbTime] = clip_pcb_fromMocap(datas, times);

% check this is correct. These numbers should be nearly the same:
mocapT(end)
fsrTime(end)
pcbTime(end)


%% clip end when walking off the floor at the end 9/10/21
% don't need this for human subject testing 10/18/21
figure; plot(pcbTime, pcbData(:,2)) % use this to see when end is
%%
clip_end = 116; % delete data after this timestamp for pcb,fsr ##### change ####

fsr_end = findTindex(clip_end,fsrTime);
fsrTime = fsrTime(1:fsr_end);
fsrData = fsrData(1:fsr_end,:);

pcb_end = findTindex(clip_end,pcbTime);
pcbTime = pcbTime(1:pcb_end);
pcbData = pcbData(1:pcb_end,:);

mocap_end = findTindex(clip_end,mocapT);
mocapT = mocapT(1:mocap_end);
mocapR = mocapR(1:mocap_end,:);
mocapL = mocapL(1:mocap_end,:);
%% extracting data 9/9/21
% overall plot for visual check
% plot_3data(pcbData,pcbTime,fsrData,fsrTime,mocapR,mocapL,mocapT,Mfsr)

% clean and filter pcb data
filt_pcbD = lpf_data(pcbData);

% finding footfalls based on fsr heel data
% [impacts, Rheel,- Lheel, Rtoe, Ltoe] = findimpacts_fsr(fsrTime,fsrData,Mfsr);
L_dist = 350; % min distance between heel peaks
R_dist = 400;
min_threshL = 15; % value between peaks threshold, if not lower than this then omit peak
min_threshR = 15;
toe_threshL = 5.5; % value after toe peak where spike ends (toe off)
toe_threshR = 3.1;

impacts = findimpacts_fsr_compact(fsrTime,fsrData,Mfsr,L_dist,R_dist,min_threshL,min_threshR,toe_threshL,toe_threshR);

%% fix small errors in impacts 11/5/21
heel_start_wrong = [37804,41249,49332,1981,4285,12082,12267,14875,16688,20394,23683,26009,30041,32875,35430,42394,44442,46985,48524,48545,49416]; % these need to be same length
heel_start_right = [38316,41240,49393,1959,4268,11917,12559,14860,16661,20763,23678,25993,30025,32860,35766,42734,44426,47368,48503,49098,49676];

heel_pk_wrong = []; % index, these need to be same length
heel_pk_right = [];

impacts = manual_fix_fsr(impacts,fsrData,Mfsr,heel_start_wrong,heel_start_right,heel_pk_wrong,heel_pk_right);

% % toe peak/off data 11/6/21
% 
% toe_threshL = 1.3; % value after toe peak where spike ends (toe off)
% toe_threshR = 1.3;

%% pcb adjust params to capture whole impact 11/5/21
% look at where the impacts are and adjust window_width and offset
% parameters in the next section

close all
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

window_width = Fs*0.45; % shaking lasts < N seconds
offset = 0.2; % start window N behind fsr start
[arrival_idx, peak_idx, peak_mag] = findimpacts_pcb(window_width,offset,impacts,fsrTime,pcbTime,filt_pcbD,num_sensors,true);

% mocap extract
[extracted_pts_R, extracted_pts_L, coordinates, whichfoot] = findimpacts_mocap(impacts,fsrTime,mocapT,mocapR,mocapL,true);

%% fix small errors in pcb 11/6/21
arrival_wrong = [NaN,NaN,NaN,8.10383;NaN,NaN,NaN,86.8053]; % these are matrixes, 4 cols for sensorN and each row is a different wrong impact
arrival_right = [NaN,NaN,NaN,7.98969;NaN,NaN,NaN,87.3008];
peak_idx_right = [NaN,NaN,NaN,7.99328;NaN,NaN,NaN,87.3051];

[arrival_idx,peak_idx,peak_mag] = manual_fix_pcb(arrival_wrong,arrival_right,peak_idx_right,pcbTime,arrival_idx,peak_idx,peak_mag,filt_pcbD);


%% fix small errors in mocap 11/5/21
r_wrong = []; % these need to be same length
r_right = [];

l_wrong = [10255,15125,21727]; % index, these need to be same length
l_right = [10365,15223,21825];

[extracted_pts_R,extracted_pts_L,coordinates] = manual_fix_mocap(r_wrong,r_right,l_wrong,l_right,mocapR,mocapL,extracted_pts_R,extracted_pts_L,coordinates,whichfoot);

%% saving data to matlab and excel 9/9/21
filename = [data_root_katie, 'ProcessedData\Subj', subj, '_', intervention];
% for findimpacts_fsr:
% save(filename,'pcbTime','filt_pcbD','arrival_idx','peak_idx','peak_mag', ...
%     'fsrTime','fsrData','impacts','Rheel','Lheel','Rtoe','Ltoe', ...
%     'mocapT','mocapR','mocapL','extracted_pts_R','extracted_pts_L','coordinates','whichfoot')

% for findimpacts_fsr_compact:
save(filename,'pcbTime','filt_pcbD','arrival_idx','peak_idx','peak_mag', ...
    'fsrTime','fsrData','impacts', ...
    'mocapT','mocapR','mocapL','extracted_pts_R','extracted_pts_L','coordinates','whichfoot')
disp(append("Saved as ", filename))

% EXCEL
excelfilename = [filename, '_excel.xlsx'];
T = table(arrival_idx,peak_idx,peak_mag,impacts,coordinates);
writetable(T,excelfilename)

%% clip end when walking off the floor at the end 9/10/21
% don't need this for human subject testing 10/18/21
figure; plot(pcbTime, pcbData(:,1)) % use this to see when end is
clip_end = 54; % (secs) ##### change ####

fsr_end = findTindex(clip_end,fsrTime);
fsrTime = fsrTime(1:fsr_end);
fsrData = fsrData(1:fsr_end,:);

pcb_end = findTindex(clip_end,pcbTime);
pcbTime = pcbTime(1:pcb_end);
pcbData = pcbData(1:pcb_end,:);

mocap_end = findTindex(clip_end,mocapT);
mocapT = mocapT(1:mocap_end);
mocapR = mocapR(1:mocap_end,:);
mocapL = mocapL(1:mocap_end,:);