clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_18_21\';

% ################# change ########################
intervention = 'stompL_sameloc'; 
FSRvalueSet = [1 2 3 4 5 6 7 8]; % order of FSR inputs based on keyset
MCvalueSet = [3 4 5 6 7 8]; % order of Mocap left/right based on keyset
% #################################################

% load 3 datasets
new_root = [data_root_katie];
load([new_root, intervention])
load([new_root, intervention, '_fsr'])
T = readtable([new_root, intervention, '_mocap']);

% constants
num_sensors = 4;
Fs = 12800;
FSRkeySet = {'Ltoe','Lheel','Lheel2','Lheel3','Rtoe','Rheel','Rheel2','Rheel3'};
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


%% clip end when walking off the floor at the end 9/10/21
% don't need this for human subject testing 10/18/21
figure; plot(pcbTime, pcbData(:,2)) % use this to see when end is
%%
clip_end = 52; % delete data after this timestamp for pcb,fsr ##### change ####

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

%%
clip_start = 13.75; % delete data before this timestamp for pcb,fsr ##### change ####

fsr_start = findTindex(clip_start,fsrTime);
fsrTime = fsrTime(fsr_start:end);
fsrData = fsrData(fsr_start:end,:);

pcb_start = findTindex(clip_start,pcbTime);
pcbTime = pcbTime(pcb_start:end);
pcbData = pcbData(pcb_start:end,:);

mocap_start = findTindex(clip_start,mocapT);
mocapT = mocapT(mocap_start:end);
mocapR = mocapR(mocap_start:end,:);
mocapL = mocapL(mocap_start:end,:);
%% extracting data 9/9/21
% overall plot for visual check
% plot_3data(pcbData,pcbTime,fsrData,fsrTime,mocapR,mocapL,mocapT,Mfsr)

% clean and filter pcb data
filt_pcbD = lpf_data(pcbData);

% finding footfalls based on fsr heel data
% [impacts, Rheel,- Lheel, Rtoe, Ltoe] = findimpacts_fsr(fsrTime,fsrData,Mfsr);
L_dist = 350; % min distance between heel peaks
R_dist = 300;
min_threshL = 15; % value between peaks threshold, if not lower than this then omit peak
min_threshR = 15;
toe_threshL = 5.5; % value after toe peak where spike ends (toe off)
toe_threshR = 3.1;

impacts = findimpacts_fsr_compact(fsrTime,fsrData,Mfsr,L_dist,R_dist,min_threshL,min_threshR,toe_threshL,toe_threshR);


%% fix small errors in impacts 11/5/21
heel_start_wrong = [42042,42610]; % these need to be same length
heel_start_right = [42009,42585];

heel_pk_wrong = []; % index, these need to be same length
heel_pk_right = [];

impacts = manual_fix_fsr(impacts,fsrData,Mfsr,heel_start_wrong,heel_start_right,heel_pk_wrong,heel_pk_right);

% % toe peak/off data 11/6/21
% 
% toe_threshL = 1.3; % value after toe peak where spike ends (toe off)
% toe_threshR = 1.3;

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

turn_threshold = 1.3; % peaks are within this # of datapoints 
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

window_width = Fs*0.3; % shaking lasts < N seconds
offset = 0.06; % start window N behind fsr start
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
filename = [data_root_katie, 'ProcessedData\', intervention];
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
% excelfilename = [filename, '_excel.xlsx'];
% T = table(arrival_idx,peak_idx,peak_mag,impacts,coordinates);
% writetable(T,excelfilename)

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
%%
% 9/9/21, 11/16/21
% uses change in direction in x-dir to determine walking episode separation
% plots the original points and overlays with x with the deleted points
clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_18_21\ProcessedData\';
takes = {'regular','rkneestiff1','rkneestiff2','run','insoleRweightL','stompL'};
deleted_pts = cell(length(takes),1);

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, intervention];
    load(filename)

    % 10/20/21 
    % assumes that subj always walks in the long direction, so 
    % any change in x direction is a turn
    curr_x = sign( coordinates(2,1)-coordinates(1,1) );
    for i = 2:(length(coordinates)-1)
        new_x = sign( coordinates(i+1,1)-coordinates(i,1) );
        if curr_x ~= new_x % change in direction
            curr_x = new_x;
            deleted_pts{take}(end+1) = i;
        end
    end
    % compare these two figures with the two at the end of this section
    % the turning points should've been eliminated
    figure; 
    plot(deleted_pts{take},coordinates(deleted_pts{take},1),'bx')
    hold on
    plot(coordinates(:,1),'r.')
end

%% based on above plots, add extra points
% 11/16/21

deleted_pts{1} = sort([deleted_pts{1} ]);
deleted_pts{2} = sort([deleted_pts{2} ]);
deleted_pts{3} = sort([deleted_pts{3} ]);
deleted_pts{4} = sort([deleted_pts{4} ]);
deleted_pts{5} = sort([deleted_pts{5} ,101]);
deleted_pts{6} = sort([deleted_pts{6} ]);

% omit_pts{1} = [];
% omit_pts{2} = [];
% omit_pts{3} = [];
% omit_pts{4} = [];
% omit_pts{5} = [];
% % omit_pts{6} = [];


for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, intervention];
    load(filename)
    
    figure; 
    plot(deleted_pts{take},coordinates(deleted_pts{take},1),'bx')
    hold on
    plot(coordinates(:,1),'r.')
    title(intervention)
end

%% then save the new data
% 11/16/21

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, intervention];
    load(filename)
    
    delete_idx = deleted_pts{take};
    
    % start of segment is -1, end of segment is 1, 0 in between
    walk_episodes = zeros(length(arrival_idx),1);
    if delete_idx(1) ~= 1
        walk_episodes(1) = -1;
        walk_episodes(delete_idx(1)-1) = 1;
    end
    for i = 2:length(delete_idx)
        % assume if more than one is deleted, they are always consecutive
        if delete_idx(i)-delete_idx(i-1) > 1 % not consecutive
            walk_episodes(delete_idx(i)-1) = 1;
            walk_episodes(delete_idx(i-1)+1) = -1;
        end
    end
    if delete_idx(end) < length(walk_episodes)
        walk_episodes(end) = 1;
        walk_episodes(delete_idx(end)+1) = -1;
    end
    arrival_idx(delete_idx,:) = [];
    coordinates(delete_idx,:) = [];
    impacts(delete_idx,:) = [];
    peak_idx(delete_idx,:) = [];
    peak_mag(delete_idx,:) = [];
    whichfoot(delete_idx,:) = [];
    walk_episodes(delete_idx,:) = [];
    
    figure;
    plot(coordinates(:,1),'k.')
    title(intervention)
    hold on
    plot(find(walk_episodes==1),coordinates(find(walk_episodes==1),1),'bx')
    plot(find(walk_episodes==-1),coordinates(find(walk_episodes==-1),1),'rx')
    
    filename = [filename, '_extract_straight_paths'];
    save(filename,'arrival_idx','coordinates','impacts','peak_idx','peak_mag','whichfoot','walk_episodes','filt_pcbD','fsrData','fsrTime','mocapL','mocapR','mocapT','pcbTime')
end
%% real_steptime analysis
clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_18_21\ProcessedData\';
takes = {'regular','rkneestiff1','rkneestiff2','run','insoleRweightL','stompL'};

Fs = 518.5;

all_params = zeros(length(takes),4); % stores left_right mean, std, right_left mean, std

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, intervention, '_extract_straight_paths'];
    load(string(filename))
    
    left_right_diff = [];
    right_left_diff = [];

    for i = 2:length(whichfoot)
        % this if statement is for extract_straight_paths
        if walk_episodes(i) ~= -1 % not the start of episode
            if whichfoot(i) == 1 % right foot
                if whichfoot(i-1) == 0 % left foot
                    left_right_diff(end+1) = impacts(i,1)-impacts(i-1,1);
                end
            elseif whichfoot(i) == 0 % left foot
                if whichfoot(i-1) == 1 % right foot
                    right_left_diff(end+1) = impacts(i,1)-impacts(i-1,1);
                end
            end
        end
    end
    
    left_right_diff = left_right_diff./Fs;
    right_left_diff = right_left_diff./Fs;
    
    left_right_mean = mean(left_right_diff);
    left_right_std = std(left_right_diff);
    right_left_mean = mean(right_left_diff);
    right_left_std = std(right_left_diff);
    all_params(take,:) = [left_right_mean, right_left_mean, left_right_std, right_left_std];
    
end

all_params

%% 10/25/21 plotting to histograms and extracting distribution parameters

clear all
close all
data_root_katie = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_18_21\ProcessedData\';
takes = {'regular','rkneestiff1','rkneestiff2','run','insoleRweightL','stompL'};
GMmodels = zeros(length(takes),9);

for take = 1:length(takes)
    intervention = char(takes(take));
    filename = [data_root_katie, intervention, '_extract_straight_paths'];
    load(string(filename))
    Fs = 12800;

%     for extract_straight_paths:
    differences = [];
    for i = 2:length(whichfoot)
        if walk_episodes(i) ~= -1 % not the start of episode
            curr = min(arrival_idx(i,1:4));
            prev = min(arrival_idx(i-1,1:4));
            differences(end+1) = (curr - prev)/Fs;
        end
    end

    mu=mean(differences);
    outliers = find(differences>(mu*1.5) );
    differences(outliers) = [];

    bin = round(1+3.22*log(numel(differences)));
    figure
    hf=histfit(differences,bin,'kernel');
    figure
    histfit(differences,bin)
    hold on
    x=get(hf(2),'XData'); 
    y=get(hf(2),'YData');
    plot(x,y,'Color','b','LineWidth', 2)
    titleroot = 'Step time distribution ';
    title([titleroot, intervention])
    xlabel('Step Time')
    ylabel('Occurances')
    mu=mean(differences);
    sig=std(differences);
    sigma = sig;
    hold on
    line([mu, mu], ylim, 'Color', 'c', 'LineWidth', 1); 
    line([mu + sigma, mu + sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
    line([mu - sigma, mu - sigma], ylim, 'Color', 'c', 'LineWidth', 0.5); 
    legend('Histogram','Normal','Fitted','1 std')
    k = kurtosis(differences);
    s = skewness(differences);
    bimodality = (s^2 + 1)/k;
    iqrange = iqr(differences);

    disp([intervention, ' ','mu ','sigma ','bimodality ','kurtosis ','skewness ', 'iqr'])
    disp([mu, sigma, bimodality, k, s, iqrange])
    
    mu_s = [mu; mu+0.01];
    sigma_s = zeros(1,1,2);
    sigma_s(1,1,:) = [sig; sig];
    pcomponents = [1/2,1/2];
    S = struct('mu',mu_s,'Sigma',sigma_s,'ComponentProportion',pcomponents);
%     GM = fitgmdist(transpose(differences),2,'Start',S);
    GM = fitgmdist(transpose(differences),2,'RegularizationValue',0.0001);
    proportion = GM.ComponentProportion;
    mu = GM.mu;
    sig = GM.Sigma;
    GMmodels(take,:) = [mu(1), mu(2), sig(1), sig(2), proportion(1), proportion(2), GM.AIC, GM.BIC, GM.NegativeLogLikelihood];
end

GMmodels

%% accel imu data analysis 11/20/21

load('C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_18_21\hammer_imu');
load('C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_18_21\hammer_imu_fsr');

[fsrData, fsrTime] = clip_fsr_fromMocap(Data, Time, 148.1);
[pcbData, pcbTime] = clip_pcb_fromMocap(datas, times);

% check this is correct. These numbers should be nearly the same:
fsrTime(end)
pcbTime(end)

accX = fsrData(:,5);
accY = fsrData(:,6);
accZ = fsrData(:,7);
fsr = fsrData(:,1);

start_clip = findTindex(15,fsrTime);
end_clip = findTindex(67.4,fsrTime);

accX = accX(start_clip:end_clip);
accY = accY(start_clip:end_clip);
accZ = accZ(start_clip:end_clip);
fsrTime = linspace(1,52.4,length(accX));

start_clip = findTindex(15,pcbTime);
end_clip = findTindex(67.4,pcbTime); 

pcbData = pcbData(start_clip:end_clip,:);
pcbTime = linspace(1,52.4,length(pcbData(:,1)));

figure;
subplot(4,1,1)
plot(fsrTime,accX)
subplot(4,1,2)
plot(fsrTime,accY)
subplot(4,1,3)
plot(fsrTime,accZ)
subplot(4,1,4)
plot(pcbTime,pcbData)

%% accel imu analysis cont.
[pksY,locsY,~,~] = findpeaks(abs(accY),'MinPeakDistance',80,'MinPeakProminence',0.8);
N = length(locsY);
pksX = zeros(N,1);
locsX = zeros(N,1);
pksZ = zeros(N,1);
locxZ = zeros(N,1);
pksP = zeros(N,4);
locsP = zeros(N,4);

for i = 1:N
    loc = locsY(i);
    [pksX(i),locsX(i)] = max(abs(accX(loc-40:loc+40)));
    locsX(i) = locsX(i) + loc-40;
    [pksZ(i),locsZ(i)] = max(abs(accZ(loc-40:loc+40)));
    locsZ(i) = locsZ(i) + loc-40;
    pcb_i = findTindex(fsrTime(loc),pcbTime);
    for j = 1:4
        [pksP(i,j),locsP(i,j)] = max(pcbData(pcb_i-1000:pcb_i+1000,j));
        locsP(i,j) = locsP(i,j) + pcb_i-1000;
    end
end
    
figure;
subplot(4,1,1)
plot(abs(accX))
hold on
plot(locsX,pksX,'rx','MarkerSize',10)
subplot(4,1,2)
plot(abs(accY))
hold on
plot(locsY,pksY,'rx','MarkerSize',10)
subplot(4,1,3)
plot(abs(accZ))
hold on
plot(locsZ,pksZ,'rx','MarkerSize',10)
subplot(4,1,4)
plot(abs(pcbData))
hold on
plot(locsP(:,1),pksP(:,1),'rx','MarkerSize',10)

%% energy extraction 11/20/21
% get arrival idx
arrival_idx = zeros(N,4);
offset = 0.2;
window_width = 12800*0.4;
filt_pcbD = pcbData;
for i = 1:N
    starti = locsP(i,1) - round(offset*window_width);
    for j = 1:4
        endi = min(starti+window_width, length(filt_pcbD(:,1)));
        window = filt_pcbD(starti:endi,j);
        if any(window > 0)
            arrival_idx(i,j) = aic_pick(window, 'to_peak')+starti;
        else
            arrival_idx(i,j) = NaN;
        end
    end
end

% get energy
noise_thresh_arr = [0.002,0.00125,0.0012,0.0015];

energy = zeros(length(arrival_idx),4);
for j = 1:4
    for i = 1:length(arrival_idx)
        noise_thresh = noise_thresh_arr(j);
        if i == length(arrival_idx)
            window_end = length(filt_pcbD);
        else
            window_end = arrival_idx(i+1,j);
        end
        window = filt_pcbD(arrival_idx(i,j):window_end,j);
        if ~isempty(window)
            [up,lo] = envelope(abs(window),500,'rms');
            indeces = find(up < noise_thresh);
            while isempty(indeces) % keep raising threshold until find the signal end
                noise_thresh = noise_thresh + 0.0005;
                indeces = find(up < noise_thresh);
            end
            window = window(1:indeces(1));
            window_energy = sum(abs(window));
            energy(i,j) = window_energy; % energy might need to be squared sig but still ok?
        end
    end
end

filename = 'C:\Users\Katie\Dropbox (MIT)\Lab\Analysis\Experiment2\11_18_21\ProcessedData\hammer_imu';
save(filename,'pcbTime','filt_pcbD','arrival_idx','pksP','pksX','pksY','pksZ','locsX','locsY','locsZ', ...
    'fsrTime','fsrData','energy')

%% different types of excel saving
excelfilename = [filename, '_excel_sqrtZImu.csv'];
all_imu = pksZ;
T = table(all_imu,pksP,energy);
writetable(T,excelfilename)

