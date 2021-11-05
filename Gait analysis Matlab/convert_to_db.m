clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\';

% ################# change ########################
subj = '5'; % number of subject
intervention = 'count'; % normal1 slow insole weight count normal2
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

% for i = 280252:283765 % for error in fsr trigger
%     Data(9,i) = -3;
% end
% for i = 103656:113140 % for error in fsr trigger
%     Data(9,i) = -3;
% end
%% data clipping to mocap length 8/30/21
[fsrData, fsrTime] = clip_fsr_fromMocap(Data, Time);
[pcbData, pcbTime] = clip_pcb_fromMocap(datas, times);

% check this is correct. These numbers should be nearly the same:
mocapT(end)
fsrTime(end)
pcbTime(end)

%% extracting data 9/9/21
% overall plot for visual check
% plot_3data(pcbData,pcbTime,fsrData,fsrTime,mocapR,mocapL,mocapT,Mfsr)

% clean and filter pcb data
filt_pcbD = lpf_data(pcbData);

% finding footfalls based on fsr heel data
% [impacts, Rheel, Lheel, Rtoe, Ltoe] = findimpacts_fsr(fsrTime,fsrData,Mfsr);
impacts = findimpacts_fsr_compact(fsrTime,fsrData,Mfsr);

%% fix small errors in impacts 11/5/21
heel_start_wrong = [15099]; % these need to be same length
heel_start_right = [15082];

heel_pk_wrong = [20452]; % index, these need to be same length
heel_pk_right = [20924];

impacts = manual_fix_fsr(impacts,fsrData,Mfsr,heel_start_wrong,heel_start_right,heel_pk_wrong,heel_pk_right);

%% pcb adjust params to capture whole impact 11/5/21
% look at where the impacts are and adjust window_width and offset
% parameters in the next section

secs = 10; % plot 10 secs of data
figure;
for i = 1:4
    subplot(4,1,i)
    plot(pcbTime(1:secs*Fs), filt_pcbD(1:secs*Fs,i))
    hold on
    indeces = find(fsrTime(impacts(:,1)) < secs);
    plot(fsrTime(impacts(indeces,1)),0,'r.','MarkerSize',10)
end

%% pcb and mocap extract

window_width = Fs*0.4; % shaking lasts < N seconds
offset = 0.05; % start window N behind fsr start
[arrival_idx, peak_idx, peak_mag] = findimpacts_pcb(window_width,offset,impacts,fsrTime,pcbTime,filt_pcbD,num_sensors,true);

% mocap extract
[extracted_pts_R, extracted_pts_L, coordinates, whichfoot] = findimpacts_mocap(impacts,fsrTime,mocapT,mocapR,mocapL,true);

%% fix small errors in mocap 11/5/21
r_wrong = [15099]; % these need to be same length
r_right = [15082];

l_wrong = [20452]; % index, these need to be same length
l_right = [20924];

impacts = manual_fix_mocap(mocap);% ######################## stopped here 11/5/21

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