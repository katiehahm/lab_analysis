clear all
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

%%
% overall plot for visual check
plot_3data(pcbD,pcbT,fsrT,fsrD,mocapT,mocapR,mocapL,Mfsr)

% clean and filter pcb data
filt_pcbD = lpf_data(pcbD);

% finding footfalls based on fsr heel data 6/7/21
[impactT, impactsM, RheelT, RtoeT, LheelT, LtoeT] = findimpacts_fsr(fsrT,fsrD,Mfsr);
[RheelT, RtoeT, LheelT, LtoeT] = heel_toe_fsr(impactsT, fsrT, fsrD, Mfsr);


% pcb extract
[arrival_idx, peak_idx, peak_mag, impactT, impacts] = findimpacts_pcb(impactT,pcbT,filt_pcbD,Fs,num_sensors,true);

% mocap extract
[extracted_pts_R, extracted_pts_L, coordinates, whichfoot] = findimpacts_mocap(impactT,mocapT,mocapR,mocapL,true);
