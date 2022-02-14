% 1/12/22
% used to estimate TA from armax model

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

[accData, accTime, fsrData, fsrTime] = clip_accelfsr_fromMocap(Data, Time, Fs);
[pcbData, pcbTime] = clip_pcb_fromMocap(datas, times);

% check this is correct. These numbers should be nearly the same:
mocapT(end)
fsrTime(end)
pcbTime(end)

% clean and filter pcb data
downData = downsample(pcbData, 10);
filt_pcbD = lpf_data(pcbData);
filt_down_pcbD = lpf_data(downData);

%% load processed data to 
filepath = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment3\ProcessedData\';
subj = '1'; 
intervention = 'regular1'; % regular1 brace1 brace2 weight1 weight2 regular2
filepath = [filepath, 'subj',subj,'_',intervention];
load(filepath)


N = length(arrival_idx(:,1));
down_arrival_idx = round(arrival_idx(:,1)./10);
for i = 1:10
    % just sensor 1
    clip = filt_down_pcbD(down_arrival_idx(i):down_arrival_idx(i)+256,1); % 0.5s window
    clip_d = iddata(clip,[],0.00078125); % sample time is 0.5
    na = 100;
    nc = 50;
%     sys = armax(clip_d,[na nc]);
    sys = ar(clip_d, na);
    figure; compare(clip_d, sys)
end