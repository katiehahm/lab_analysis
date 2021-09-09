clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\';

% ################# change ########################
datestr = '08_30_1';
run_num = '_221';
FSRvalueSet = [1 2 3 4 5 6 7 8]; % order of FSR inputs based on keyset
MCvalueSet = [3 4 5 6 7 8]; % order of Mocap left/right based on keyset
% #################################################

% load 3 datasets
load([data_root_katie, datestr, run_num])
load([data_root_katie, datestr])
T = readtable([data_root_katie, 'Take ', datestr]);

% constants
num_sensors = 4;
Fs = 12800;
FSRkeySet = {'Rballi','Rtoe','Rheel','Rballo','Ltoe','Lheel','Lballi','Lballo'};
Mfsr = containers.Map(FSRkeySet, FSRvalueSet);
MCkeySet = {'Lx','Ly','Lz','Rx','Ry','Rz'};
Mmocap = containers.Map(MCkeySet, MCvalueSet);

[mocapT, mocapL, mocapR] = convertMocap(T, Mmocap);

%% data clipping to mocap length 8/30/21
[fsrData, fsrTime] = clip_fsr_fromMocap(Data, Time);
[pcbData, pcbTime] = clip_pcb_fromMocap(datas, times);

% check this is correct. These numbers should be nearly the same:
mocapT(end)
fsrTime(end)
pcbTime(end)

%% extracting data 9/9/21
% overall plot for visual check
plot_3data(pcbData,pcbTime,fsrData,fsrTime,mocapR,mocapL,mocapT,Mfsr)

% clean and filter pcb data
filt_pcbD = lpf_data(pcbData);

% finding footfalls based on fsr heel data
[impacts, Rheel, Lheel, Rtoe, Ltoe] = findimpacts_fsr(fsrTime,fsrData,Mfsr);

% pcb extract
[arrival_idx, peak_idx, peak_mag] = findimpacts_pcb(impacts,fsrTime,pcbTime,filt_pcbD,Fs,num_sensors,true);

% mocap extract
[extracted_pts_R, extracted_pts_L, coordinates, whichfoot] = findimpacts_mocap(impacts,fsrTime,mocapT,mocapR,mocapL,true);

%% saving data to matlab and excel 9/9/21
filename = [data_root_katie, 'ProcessedData\', datestr];
save(filename,'pcbTime','filt_pcbD','arrival_idx','peak_idx','peak_mag', ...
    'fsrTime','fsrData','impacts','Rheel','Lheel','Rtoe','Ltoe', ...
    'mocapT','mocapR','mocapL','extracted_pts_R','extracted_pts_L','coordinates','whichfoot')
disp(append("Saved as ", filename))

excelfilename = [data_root_katie, 'ProcessedData\', datestr, '_excel.xlsx'];
T = table(arrival_idx,peak_idx,peak_mag,impacts,coordinates);
writetable(T,excelfilename)

